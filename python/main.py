# from osgeo import ogr, osr, gdal
import os.path

import gdal
import numpy as np
import heapq
from time import time
from osgeo import osr
import sys

# values according to:
# https://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/flow-direction.htm
direction_values = {
    (0, 0): 0,
    (0, 1): 1,
    (1, 1): 2,
    (1, 0): 4,
    (1, -1): 8,
    (0, -1): 16,
    (-1, -1): 32,
    (-1, 0): 64,
    (-1, 1): 128,
}

# reversed mapping to get the draining neighbor for the cell by its flow value
deltas = {v: k for k, v in direction_values.items()}

nodata = 0


class Raster:
    def __init__(self, nrows, ncols, origin, pixelsize):
        self.pixels = np.empty((nrows, ncols))
        self.ncols = ncols
        self.nrows = nrows
        self.total_pixels = ncols * nrows
        self.origin = origin
        self.pixelsize = pixelsize
        self.insertion_counter = 0
        self.visit_order = []

    def copy(self):
        return Raster(self.nrows, self.ncols, self.origin, self.pixelsize)

    def fill(self, value):
        self.pixels.fill(value)

    def set_pixels(self, array):
        self.pixels = array

    def new_order(self):
        """Insertion order counter for the points added to the queue
        Later used to resolve comparison of equal height points by which was added first"""
        self.insertion_counter += 1
        return self.insertion_counter

    def boundary_to_queue(self):
        """Returns a queue with boundary pixels"""
        queue = []

        for i in range(self.ncols):
            r1 = RasterCell(0, i, self.pixels[0][i], self.new_order())
            queue.append(r1)
            r2 = RasterCell(self.nrows - 1, i, self.pixels[-1][i], self.new_order())
            queue.append(r2)

        # skip first row and last row to eliminate duplicate pixels in corners
        for i, (v1, v2) in enumerate([(row[0], row[-1]) for row in self.pixels[1:-1]], start=1):
            r1 = RasterCell(i, 0, v1, self.new_order())
            queue.append(r1)
            r2 = RasterCell(i, self.ncols - 1, v2, self.new_order())
            queue.append(r2)

        heapq.heapify(queue)
        return queue

    def get_neighbors(self, cell):
        """Returns an iterator of the cells neighbors and the flow value from neighbor to cell"""
        for dy, dx in deltas.values():
            if dx == dy == 0:
                continue
            new_row = cell.row + dy
            new_col = cell.col + dx

            if 0 <= new_row < self.nrows and 0 <= new_col < self.ncols:
                yield (RasterCell(new_row, new_col, self.pixels[new_row][new_col]),
                       direction_values[(-dy, -dx)])

    def flow_direction(self):
        """To be called on dtm raster, returns flow_direction raster"""
        queue = self.boundary_to_queue()
        result_raster = self.copy()
        result_raster.fill(0)
        processed = np.ones(self.pixels.shape, dtype=bool)
        processed[self.pixels.ndim * (slice(1, -1),)] = False
        # boundary points are in queue so counted as processed

        while len(queue) > 0:
            cell = heapq.heappop(queue)
            result_raster.visit_order.append(cell)
            for neighbor, direction in self.get_neighbors(cell):
                if not processed[neighbor.row, neighbor.col]:
                    # set drainage direction for the neighbor to the searched cell
                    result_raster.pixels[neighbor.row, neighbor.col] = direction
                    # insertion order to resolve equal height cells - first in first out if heights are equal
                    neighbor.insertion_order = self.new_order()
                    heapq.heappush(queue, neighbor)
                    processed[neighbor.row, neighbor.col] = True

        return result_raster

    def flow_accumulation(self):
        """To be called on flow direction raster. Returns flow accumulation raster"""
        result_raster = self.copy()
        result_raster.fill(1)

        # opposite order than in flow direction
        for cell in reversed(self.visit_order):
            dy, dx = deltas[self.pixels[cell.row, cell.col]]
            result_raster.pixels[cell.row + dy, cell.col + dx] += result_raster.pixels[cell.row, cell.col]

        return result_raster

    def save(self, filename, gdaltype, epsg=4326):
        originX, originY = self.origin
        pixelWidth, pixelHeight = self.pixelsize

        driver = gdal.GetDriverByName('GTiff')
        outRaster = driver.Create(filename, self.ncols, self.nrows, 1, gdaltype)
        outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
        outband = outRaster.GetRasterBand(1)
        outband.SetNoDataValue(nodata)
        outband.WriteArray(self.pixels)
        outRasterSRS = osr.SpatialReference()
        outRasterSRS.ImportFromEPSG(epsg)
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()

    def set_edges_nodata(self):
        for i in range(self.ncols):
            self.pixels[0][i] = nodata
            self.pixels[self.nrows - 1][i] = nodata
        for i in range(1, self.nrows - 1):
            self.pixels[i][0] = nodata
            self.pixels[i][self.ncols - 1] = nodata


def load_raster(filename):
    ds = gdal.Open(filename)
    cols = ds.RasterXSize
    rows = ds.RasterYSize

    gtrf = ds.GetGeoTransform()
    origin_x, origin_y = gtrf[0], gtrf[3]
    pixel_width, pixel_height = gtrf[1], gtrf[5]

    dtm = Raster(rows, cols, (origin_x, origin_y), (pixel_width, pixel_height))
    arr = np.array(ds.GetRasterBand(1).ReadAsArray())
    dtm.set_pixels(arr)

    return dtm


class RasterCell:
    def __init__(self, row, col, height, insertion_order=0):
        self.col = col
        self.row = row
        self.height = height
        self.insertion_order = insertion_order

    def __eq__(self, other):
        return self.height == other.height

    def __lt__(self, other):
        if self.height == other.height:
            return self.insertion_order < other.insertion_order
        return self.height < other.height

    def __repr__(self):
        return "(h={}, x={}, y={}, o={})".format(self.height, self.col, self.row, self.insertion_order)


def main():
    filename = sys.argv[1]
    direction_filename = sys.argv[2]
    accumulation_filename = sys.argv[3]

    dtm = load_raster(filename)
    flow_direction = dtm.flow_direction()
    flow_accumulation = flow_direction.flow_accumulation()

    flow_accumulation.set_edges_nodata()

    flow_direction.save(direction_filename, gdal.GDT_Byte)
    flow_accumulation.save(accumulation_filename, gdal.GDT_UInt32)


start = time()
main()
stop = time()

print('Run time: {:.2f}s'.format(stop - start))
