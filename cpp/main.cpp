#include <iostream>
#include <queue>
#include <chrono>
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "cassert"
#include "stack"

using namespace std;

// direction values from neighbor to the cell
map<pair<int, int>, int> direction_values;

//the same but reversed, so that you can get a neighbor (ex. (1, 1) or (-1, 0) by flow value
map<int, pair<int, int>> deltas;

// default nodata value
int nodata = 0;

// cpp type to gdal datatype mapping
template <typename T>
struct GdalTypeMap { static const GDALDataType type; };

template <>
const GDALDataType GdalTypeMap<uint8_t>::type = GDT_Byte;
template <>
const GDALDataType GdalTypeMap<uint32_t>::type = GDT_UInt32;
template <>
const GDALDataType GdalTypeMap<int16_t>::type = GDT_Int16;
template <>
const GDALDataType GdalTypeMap<int32_t>::type = GDT_Int32;

// Raster cell with a defined sorting operator
struct RasterCell {
    int row, col, height; // row and column of the cell
    unsigned int insertion_order;

    // Defines a new link to a cell
    RasterCell(int row, int col, int height, unsigned int insertion_order=0) {
        this->row = row;
        this->col = col;
        this->insertion_order = insertion_order;
        this->height = height;
    }

    // reversed because priority queue sorts descending
    bool operator<(const RasterCell &other) const {
        if (this->height == other.height){
            // resolve equal heights by insertion order
            return this->insertion_order > other.insertion_order;
        }
        return this->height > other.height;
    }
};

// Write the values in a linked raster cell (useful for debugging)
std::ostream& operator<<(std::ostream& os, const RasterCell& c) {
    os << "{h=" << c.height << ", o=" << c.insertion_order << ", x=" << c.col << ", y=" << c.row << "}";
    return os;
}

// Data storage, access, calculations and saving to file
struct Raster {
  std::vector<vector<int>> pixels; // where everything is stored
  std::vector<RasterCell> visit_order; //to be used later for accumulation
  int nrows, ncols; // number of columns and rows
  unsigned int insertion_counter, total_pixels; // counter to keep track of order of cell insertion
  pair<float, float> origin, pixelsize; // raster properties
  vector<vector<bool>> processed; // mark the processed cells
  const char *projection; // projection of the raster

  // Initialise a raster with ncols and nrows
  Raster(int nrows, int ncols, pair<double, double> origin, pair<double, double> pixelsize, const char *projection) {
    this->nrows = nrows;
    this->ncols = ncols;
    total_pixels = ncols * nrows;
    pixels.reserve(nrows);
    visit_order.reserve(total_pixels);
    this->origin = origin;
    this->pixelsize = pixelsize;
    insertion_counter = 0;
    this->projection = projection;
  }

  Raster copy() const {
      return {this->nrows, this->ncols, this->origin, this->pixelsize, this->projection};
  }
  
  // Fill values of an entire row
  void add_scanline(const int *line) {
      vector<int> row;
      row.reserve(ncols);
    for (int i = 0; i < ncols; ++i) {
        row.push_back(line[i]);
    }
    pixels.push_back(row);
  }
  
  // Fill entire raster with chosen value
  void fill(int value) {
      for (int i = 0; i < nrows; i++) {
          vector<int> row;
          row.reserve(ncols);
          for (int j = 0; j < ncols; j++)
              row.push_back(value);
          pixels.push_back(row);
      }
  }
  
  // Access the value of a raster cell to read or write it
  int &operator()(int row, int col) {
    assert(col >= 0 && col < ncols);
    assert(row >= 0 && row < nrows);
    return pixels[row][col];
  }
  
  // Access the value of a raster cell to read it
  int operator()(int row, int col) const {
      assert(col >= 0 && col < ncols);
      assert(row >= 0 && row < nrows);
      return pixels[row][col];
  }

  // Manage insertion order. Returns current counter and increments
  unsigned int new_order(){
      return ++insertion_counter;
  }

  // Returns a priority queue of boundary cells and marks them as processed
  priority_queue<RasterCell, vector<RasterCell>> boundary_to_queue(){
      priority_queue<RasterCell, vector<RasterCell>> q;
      vector<vector<bool>> p(nrows, std::vector<bool>(ncols, false));
      processed = p;

      for (int i = 0; i < ncols; i++){
          auto r1 = RasterCell(0, i, pixels[0][i], new_order());
          auto r2 = RasterCell(nrows - 1, i, pixels.back()[i], new_order());
          q.push(r1);
          q.push(r2);

          processed[0][i] = true;
          processed[nrows - 1][i] = true;
      }

      // skip first and last row to avoid corner duplicates
      for (int i = 1; i < nrows - 1; i++){
          auto r1 = RasterCell(i, 0, pixels[i][0], new_order());
          auto r2 = RasterCell(i, ncols - 1, pixels[i][ncols - 1], new_order());
          q.push(r1);
          q.push(r2);

          processed[i][0] = true;
          processed[i][ncols - 1] = true;
      }

      return q;
  }

  // returns a vector of neighbors and the direction values to flow towards the given cell
  vector<pair<RasterCell, int>> get_neighbors(const RasterCell &cell){
      vector<pair<RasterCell, int>> neighbors;
      neighbors.reserve(8);
      for (const auto & delta: deltas){
          int dx = delta.second.second;
          int dy = delta.second.first;
          if (dy == 0 && dx == 0)
              continue;

          int new_row = cell.row + dy;
          int new_col = cell.col + dx;

          // check if valid neighbor
          if (new_row >= 0 && new_row < nrows && new_col >= 0 && new_col < ncols){
              RasterCell new_cell = RasterCell(new_row, new_col, pixels[new_row][new_col]);
              int direction = direction_values[make_pair(-dy, -dx)];
              neighbors.emplace_back(new_cell, direction);
          }
      }
      return neighbors;
  }

  // To be called on a dtm raster. Returns flow direction raster.
  Raster flow_direction(){
      auto q = boundary_to_queue();
      auto result_raster = copy();
      result_raster.fill(0);

      while(!q.empty()){
          RasterCell cell = q.top();
          q.pop();

          result_raster.visit_order.push_back(cell);
          for (auto it: get_neighbors(cell)){
              RasterCell neighbor = it.first;
              int direction = it.second;

              // if neighbor hasn't been processed add it to the queue
              if (not processed[neighbor.row][neighbor.col]){
                  result_raster(neighbor.row, neighbor.col) = direction;
                  neighbor.insertion_order = new_order();
                  q.push(neighbor);
                  processed[neighbor.row][neighbor.col] = true;
              }
          }
      }
    return result_raster;
  }

  // to be called on flow direction raster. Returns raster of flow accumulation.
  Raster flow_accumulation(){
      auto result_raster = copy();
      result_raster.fill(1);

      // reversed visit order
      for (auto cell = visit_order.rbegin();
       cell != visit_order.rend(); ++cell ) {

          auto delta = deltas[pixels[cell->row][cell->col]];
          int dy = delta.first;
          int dx = delta.second;

          result_raster.pixels[cell->row + dy][cell->col + dx] += result_raster.pixels[cell->row][cell->col];
      }

      return result_raster;
  }

  // write to a file
  template <typename datatype>
  void save(const string &filename){
      GDALDataType gdt = GdalTypeMap<datatype>::type;

      float origin_x = origin.first;
      float origin_y = origin.second;
      float pixel_width = pixelsize.first;
      float pixel_height = pixelsize.second;

      double transform[6]{origin_x, pixel_width, 0, origin_y, 0, pixel_height};

      auto driverGeotiff = GetGDALDriverManager()->GetDriverByName("GTiff");
      auto geotiffDataset = driverGeotiff->Create(filename.c_str(),ncols,nrows,1,gdt,nullptr);
      geotiffDataset->SetGeoTransform(transform);
      geotiffDataset->SetProjection(projection);
      geotiffDataset->GetRasterBand(1)->SetNoDataValue(nodata);

      auto *rowBuff = (datatype*) CPLMalloc(sizeof(datatype)*ncols);
      for(int row=0; row<nrows; row++) {
          for(int col=0; col<ncols; col++) {
              rowBuff[col] = (datatype)pixels[row][col];
          }
          CPLErr e = geotiffDataset->GetRasterBand(1)->RasterIO(GF_Write,0,row,ncols,1,rowBuff,ncols,1,gdt,0,0);
      }
      GDALClose(geotiffDataset) ;
      // for some reason this causes the program to stop
//      GDALDestroyDriverManager();
  }

  // Sets the edge pixels to nodata, because they are unreliable.
  void set_edges_nodata(){
      for (int i = 0; i < ncols; i++){
          pixels[0][i] = nodata;
          pixels[nrows - 1][i] = nodata;
      }
      for (int i = 1; i < nrows - 1; i++){
          pixels[i][0] = nodata;
          pixels[i][ncols - 1] = nodata;
      }
  }
};

Raster load_raster(GDALDataset *input_dataset){
    // Print dataset info
    double geo_transform[6];
    std::cout << "Driver: " << input_dataset->GetDriver()->GetDescription() << "/" << input_dataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME) << std::endl;
    std::cout << "Size is " << input_dataset->GetRasterXSize() << "x" << input_dataset->GetRasterYSize() << "x" << input_dataset->GetRasterCount() << std::endl;
    if (input_dataset->GetProjectionRef() != NULL) std::cout << "Projection is '" << input_dataset->GetProjectionRef() << "'" << std::endl;
    if (input_dataset->GetGeoTransform(geo_transform) == CE_None) {
        std::cout << "Origin = (" << geo_transform[0] << ", " << geo_transform[3] << ")" << std::endl;
        std::cout << "Pixel Size = (" << geo_transform[1] << ", " << geo_transform[5] << ")" << std::endl;
    }

    // Print Band 1 info
    GDALRasterBand *input_band;
    int nBlockXSize, nBlockYSize;
    int bGotMin, bGotMax;
    double adfMinMax[2];
    input_band = input_dataset->GetRasterBand(1);
    input_band->GetBlockSize(&nBlockXSize, &nBlockYSize);
    std::cout << "Band 1 Block=" << nBlockXSize << "x" << nBlockYSize << " Type=" << GDALGetDataTypeName(input_band->GetRasterDataType()) << " ColorInterp=" << GDALGetColorInterpretationName(input_band->GetColorInterpretation()) << std::endl;
    adfMinMax[0] = input_band->GetMinimum(&bGotMin);
    adfMinMax[1] = input_band->GetMaximum(&bGotMax);
    if (!(bGotMin && bGotMax)) GDALComputeRasterMinMax((GDALRasterBandH)input_band, TRUE, adfMinMax);
    std::cout << "Min=" << adfMinMax[0] << " Max=" << adfMinMax[1] << std::endl;

    // Read Band 1 line by line
    int nXSize = input_band->GetXSize();
    int nYSize = input_band->GetYSize();
    pair<double, double> origin = make_pair(geo_transform[0], geo_transform[3]);
    pair<double, double> pixelsize = make_pair(geo_transform[1], geo_transform[5]);
    auto projection = input_dataset->GetProjectionRef();
    Raster input_raster(nYSize, nXSize, origin, pixelsize, projection);
    for (int current_scanline = 0; current_scanline < nYSize; ++current_scanline) {
        int *scanline = (int *)CPLMalloc(sizeof(float)*nXSize);
        if (input_band->RasterIO(GF_Read, 0, current_scanline, nXSize, 1,
                                 scanline, nXSize, 1, GDT_Int32,
                                 0, 0) != CPLE_None) {
            std::cerr << "Couldn't read scanline " << current_scanline << std::endl;
            exit(1);
        } input_raster.add_scanline(scanline);
        CPLFree(scanline);
    } std::cout << "Created raster: " << input_raster.ncols << "x" << input_raster.pixels.size()/input_raster.nrows << " = " << input_raster.pixels.size() << std::endl;

    return input_raster;
}

int main(int argc, const char * argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    string filename = argv[1];
    string direction_filename = argv[2];
    string accumulation_filename = argv[3];

    // globals
    direction_values = {
            {make_pair(0, 0), 0},
            {make_pair(0, 1), 1},
            {make_pair(1, 1), 2},
            {make_pair(1, 0), 4},
            {make_pair(1, -1), 8},
            {make_pair(0, -1), 16},
            {make_pair(-1, -1), 32},
            {make_pair(-1, 0), 64},
            {make_pair(-1, 1), 128}
    };

    for (auto & dir_val : direction_values)
        deltas[dir_val.second] = dir_val.first;

  // Open dataset
    GDALDataset  *input_dataset;
    GDALAllRegister();
    input_dataset = (GDALDataset *)GDALOpen(filename.c_str(), GA_ReadOnly);
    if (input_dataset == NULL) {
    std::cerr << "Couldn't open file" << std::endl;
    return 1;
    }

    Raster input_raster = load_raster(input_dataset);

  // Flow direction
    auto flow_direction = input_raster.flow_direction();
    auto flow_accumulation = flow_direction.flow_accumulation();
    flow_accumulation.set_edges_nodata();

    // flow direction needs less storage because we know the value range to be 0-128
    flow_direction.save<uint8_t>(direction_filename);
    flow_accumulation.save<uint32_t>(accumulation_filename);
  
  // Close input dataset
    GDALClose(input_dataset);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout.precision(3);
    std::cout << "Run time: " << duration.count() / 1000000.0<< " seconds" << std::endl;
    return 0;
}