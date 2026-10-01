//// Calculate NDVI in Earth Engine
// https://code.earthengine.google.com/
// 1) This is the region for our area of interest in NC
// Update the Geometry
var geometry = 
  ee.Geometry.Polygon(
    [[[-78.783740, 33.726053],
      [-78.783740, 36.608362],
      [-75.322985, 36.608362],
      [-75.322985, 33.726053]]], null, false);


// 2) Below is the Script for GEE

// Specify your start date as the start of a season
// Winter (December 1), Spring (March 1), Summer (June 1) or Autumn (September 1)
var start = '2025-04-01'
// Also, specify the number of three-month seasons to collect data on.
var numSeasons = 1

// Calculate the end date and the index of the last month in the sequence
var end = ee.Date(start).advance(7 * numSeasons, 'month')
var lastStartingMonth = (numSeasons - 1) * 7

// Create a list where each item of the list is the start date for a season
var seasonStartMonths = ee.List.sequence(0, lastStartingMonth, 7).map(function(n){
  return ee.Date(start).advance(n, 'month')
})

// Add NDVI and remove all other bands
function calculateNDVI(image) {
  var ndvi = image.normalizedDifference(['SR_B5','SR_B4']).rename('NDVI')
  return image.select().addBands(ndvi);
}

var mycollection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
.filterBounds(geometry)
.filterDate(start, end)
.map(calculateNDVI)

// Map over the list of seasons
var ndvimax_spring = ee.ImageCollection(seasonStartMonths.map(function(d) {
  var startDate = ee.Date(d)
  return mycollection.filterDate(startDate, startDate.advance(7, 'month'))
  .max()
  .set({
    'system:id': startDate.format("YYYYMMdd"),
    'system:time_start': startDate.millis()
  })
}))

print(ndvimax_spring)

// Import geetools
var batch = require('users/fitoprincipe/geetools:batch');

// Export the ImageCollection
var options = {
  name: 'ndvi2025',
  scale: 30,
  maxPixels: 1e13,
  region: geometry
}


// Export the image, specifying scale and region.
batch.Download.ImageCollection.toDrive(ndvimax_spring,
                                       'NC NDVI', options);