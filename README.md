# runoff-modelling

<p align="center">
  <img src="https://user-images.githubusercontent.com/38667147/154864849-9987f979-176a-4c88-8952-e31b32f39aed.png" />
</p>


Implementation of single flow (no dispersion) run off modelling using a Least Cost Paths algorithm to solve flats and sinks.
The code was written in Python first and then rewritten in CPP.

Command line parameters:
- input elevation raster
- output flow direction raster file name
- output flow accumulation raster file name

Example values:
_/home/S33W070.hgt myflow_dir.tif myflow_acc.tif_

This implementation uses the knowledge described in the 3D terrain book: https://github.com/tudelft3d/terrainbook

More about the LCP algorithm for single flow modelling can be found here: https://www.hydrol-earth-syst-sci.net/15/667/2011/

