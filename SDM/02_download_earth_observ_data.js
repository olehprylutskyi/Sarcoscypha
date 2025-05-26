// Author: Oleh Prylutskyi, GNU General Public License - 3.0

// This script acquires Earth observation data from various sources and produce 
// georeferenced rasters (*.tiff) to be downloaded to the Google Drive.
// This script is ment to be run within Google Earth Engine Code editor (JavaScript).

// Load the national boundary shapefile of Ukraine from the public Asset
var Ukraine = ee.FeatureCollection('users/olegpril12/Ukraine');

// Define Area Of Interest (AOI) as the entire country of Ukraine
var AOI = Ukraine;

// Set output resolution in metres (MODIS native resolution is ~500 m)
var scale = 500;

// Centre the map display over the AOI at zoom level 5
Map.centerObject(AOI, 5);


/////////// MODIS Summer and Winter NDVIs ///////////////
// MOD13A1.061 Terra Vegetation Indices 16-Day Global 500m
// https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13A1

// --------- Summer NDVI ---------
// Filter for years 2013–2022 and summer months (May to September)
var image = ee.ImageCollection('MODIS/061/MOD13A1')
            .filter(ee.Filter.calendarRange(2013,2022,'year'))
            .filter(ee.Filter.calendarRange(5,9,'month'))
            .select('NDVI')
            .median();

var dataset_name = 'MODIS_NDVI_Summer';

// Define visualisation parameters
var ndviVis = {
  min: 0.0,
  max: 9000.0,
  palette: [
    'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
    '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
    '012E01', '011D01', '011301'
  ],
};

// Optional: Add to the interactive map (layer is off by default)
Map.addLayer(image, ndviVis, 'NDVIs', false);

// Export as GeoTIFF to Google Drive
Export.image.toDrive({
  image: image,
  description: dataset_name,
  folder: 'GEE_data',
  scale: scale,
  region: AOI,
  crs: 'EPSG:4326',
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF'
});


// --------- Winter NDVI ---------
// Note: calendarRange(11, 3, 'month') wraps around year-end (Nov–Mar)
var image = ee.ImageCollection('MODIS/061/MOD13A1')
            .filter(ee.Filter.calendarRange(2013,2022,'year'))
            .filter(ee.Filter.calendarRange(11,3,'month'))     // winter
            .select('NDVI')
            .median();

var dataset_name = 'MODIS_NDVI_Winter';

var ndviVis = {
  min: 0.0,
  max: 9000.0,
  palette: [
    'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
    '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
    '012E01', '011D01', '011301'
  ],
};

Map.addLayer(image, ndviVis, 'NDVIw', false);


Export.image.toDrive({
  image: image,
  description: dataset_name,
  folder: 'GEE_data',
  scale: scale,
  region: AOI,
  crs: 'EPSG:4326',
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF'
});


/////////// MODIS EVI: Summer and Winter ///////////////
// EVI: Enhanced Vegetation Index, more sensitive in high-biomass areas
// MOD13A1.061 Terra Vegetation Indices 16-Day Global 500m
// https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13A1

// --------- Summer EVI ---------
var image = ee.ImageCollection('MODIS/061/MOD13A1')
            .filter(ee.Filter.calendarRange(2013,2022,'year'))
            .filter(ee.Filter.calendarRange(5,9,'month'))     // summer
            .select('EVI')
            .median();

var dataset_name = 'MODIS_EVI_Summer';

var ndviVis = {
  min: 0.0,
  max: 9000.0,
  palette: [
    'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
    '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
    '012E01', '011D01', '011301'
  ],
};

Map.addLayer(image, ndviVis, 'EVIs', false);


Export.image.toDrive({
  image: image,
  description: dataset_name,
  folder: 'GEE_data',
  scale: scale,
  region: AOI,
  crs: 'EPSG:4326',
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF'
});



// --------- Winter EVI ---------
var image = ee.ImageCollection('MODIS/061/MOD13A1')
            .filter(ee.Filter.calendarRange(2013,2022,'year'))
            .filter(ee.Filter.calendarRange(11,3,'month'))     // winter
            .select('EVI')
            .median();

var dataset_name = 'MODIS_EVI_Winter';

var ndviVis = {
  min: 0.0,
  max: 9000.0,
  palette: [
    'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
    '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
    '012E01', '011D01', '011301'
  ],
};

Map.addLayer(image, ndviVis, 'EVIw', false);


Export.image.toDrive({
  image: image,
  description: dataset_name,
  folder: 'GEE_data',
  scale: scale,
  region: AOI,
  crs: 'EPSG:4326',
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF'
});



/////////// MODIS NDWI: Summer and Winter //////////////////////
// NDWI from MODIS MCD43A4 product (surface reflectance-based index)
// https://developers.google.com/earth-engine/datasets/catalog/MODIS_MCD43A4_006_NDWI

// --------- Summer NDWI ---------
var image = ee.ImageCollection('MODIS/MCD43A4_006_NDWI')
            .filter(ee.Filter.calendarRange(2013,2022,'year'))
            .filter(ee.Filter.calendarRange(5,9,'month'))     // summer
            .select('NDWI')
            .median();

var dataset_name = 'MODIS_NDWI_Summer';

var ndwiVis = {
  min: 0.0,
  max: 1.0,
  palette: ['0000ff', '00ffff', 'ffff00', 'ff0000', 'ffffff'],
};

Map.addLayer(image, ndwiVis, 'NDWIs', false);

Export.image.toDrive({
  image: image,
  description: dataset_name,
  folder: 'GEE_data',
  scale: scale,
  region: AOI,
  crs: 'EPSG:4326',
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF'
});



// --------- Winter NDWI ---------
var image = ee.ImageCollection('MODIS/MCD43A4_006_NDWI')
            .filter(ee.Filter.calendarRange(2013,2022,'year'))
            .filter(ee.Filter.calendarRange(11,3,'month'))     // winter
            .select('NDWI')
            .median();

var dataset_name = 'MODIS_NDWI_Winter';

var ndwiVis = {
  min: 0.0,
  max: 1.0,
  palette: ['0000ff', '00ffff', 'ffff00', 'ff0000', 'ffffff'],
};

Map.addLayer(image, ndwiVis, 'NDWIw', false);

Export.image.toDrive({
  image: image,
  description: dataset_name,
  folder: 'GEE_data',
  scale: scale,
  region: AOI,
  crs: 'EPSG:4326',
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF'
});


// End of the script
// Note: Make sure to run this script in the Google Earth Engine Code Editor.
//       You can adjust the date ranges, visualisation parameters, and export settings as needed.
//       The exported files will be saved in the 'GEE_data' folder in your Google Drive.
//       You can access the exported files from your Google Drive account.