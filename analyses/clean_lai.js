// Earth Engine: export annual MODIS LAI (MCD15A3H) for 2007-2025 to Drive
// Robust to years with no images in the geometry/time window.

var geometry = ee.Geometry.Polygon(
  [[[-78.783740, 33.726053],
    [-78.783740, 36.8],
    [-75.322985, 36.8],
    [-75.322985, 33.726053]]], null, false);

// Parameters
var startYear = 2007;
var endYear = 2025;
var modisId = 'MODIS/006/MCD15A3H';
var laiBand = 'Lai';
var scaleFactor = 1.0;   // adjust if product uses integer encoding (e.g., 0.001 or 0.1)
var exportScale = 500;   // native resolution for MCD15A3H
var exportFolder = 'MODIS_LAI'; // Drive folder name
var exportPrefix = 'modis_lai_'; // will result in filenames like modis_lai_2007

// Load and narrow collection once (spatial + band)
// Don't filter date here, we filter per-year below
var modis = ee.ImageCollection(modisId)
  .select([laiBand])
  .filterBounds(geometry);

// Create a template image (one band named 'LAI') to use when a yearly collection is empty.
// We derive it from modis.first() so it has a valid projection/CRS/resolution/mask.
// If modis.first() is null (collection empty overall) this will fail; in that case check the dataset ID.
var firstImg = ee.Image(modis.first());
var template = firstImg
  .select([0])        // select the first band (should exist if collection not empty)
  .multiply(0)        // gives zeros
  .rename('LAI')
  .updateMask(firstImg.select([0]).mask().multiply(0)); // fully masked (all NA)

// Build list of years and map to images (annual mean)
var years = ee.List.sequence(startYear, endYear);

var annualCollection = ee.ImageCollection(
  years.map(function(y) {
    y = ee.Number(y).int();
    var start = ee.Date.fromYMD(y, 1, 1);
    var end = start.advance(1, 'year');
    var yearly = modis.filterDate(start, end);
    var count = yearly.size();
    // If count > 0 use yearly.mean(), otherwise use masked template
    var img = ee.Image(ee.Algorithms.If(
      count.gt(0),
      yearly.mean().multiply(scaleFactor).rename('LAI'),
      template
    ));
    return img.set({
      'year': y,
      'system:time_start': start.millis(),
      'description': ee.String('MODIS MCD15A3H annual mean LAI for ').cat(ee.Number(y).format())
    });
  })
);

// Quick checks
print('Annual LAI collection (size):', annualCollection.size());
print('Sample image properties:', annualCollection.first().propertyNames());
Map.centerObject(geometry, 7);
Map.addLayer(annualCollection.first(), {min:0, max:6, palette:['white','yellow','green']}, 'LAI sample 1st year');

// Export using geetools batch exporter (optional). Uncomment if you have geetools available.
// var batch = require('users/fitoprincipe/geetools:batch');
// var options = { name: 'modis_lai_2007_2025', scale: exportScale, maxPixels: 1e13, region: geometry };
// batch.Download.ImageCollection.toDrive(annualCollection, exportFolder, options);

// Alternate: Export each year via Export.image.toDrive (no external dependency)
years.getInfo().forEach(function(y){
  var img = ee.Image(annualCollection.filter(ee.Filter.eq('year', y)).first());
  Export.image.toDrive({
    image: img,
    description: exportPrefix + y,
    folder: exportFolder,
    fileNamePrefix: exportPrefix + y,
    scale: exportScale,
    region: geometry,
    maxPixels: 1e13
  });
});

// Notes:
// - If modis.first() is null (the collection ID is wrong or MCD15A3H has no images for this area), the script will error on creating the template. Check modis.size() and modis.first() in the Console.
// - Adjust scaleFactor if the catalog indicates a scaling (some products require multiplying by e.g. 0.001).
// - This script gives one masked image for years with no data; masked images will be exported but contain no valid pixels.