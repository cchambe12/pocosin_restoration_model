// Integrate Lidar instead of NED - notes from Github Copilot

// Load DEM
var dem = ee.Image('USGS/3DEP/10m');

// Define AOI
var aoi = ee.Geometry.Rectangle([-77.9, 34.4, -76.0, 36.8]);

// Export as RASTER (GeoTIFF)
Export.image.toDrive({
  image: dem,
  description: 'NC_DEM_10m',
  scale: 10,
  region: aoi,
  fileFormat: 'GeoTIFF',
  maxPixels: 1e13
});

