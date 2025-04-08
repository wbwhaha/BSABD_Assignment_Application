Map.setCenter(86.9250, 27.9881, 15);
Map.setOptions("satellite");

// Define the study area: a buffered point around Mount Everest.
var everestPoint = ee.Geometry.Point([86.9250, 27.9881]);
var everestBuffer = everestPoint.buffer(15000); // 15 km buffer

Map.centerObject(everestBuffer, 10);
Map.addLayer(everestBuffer, {color: 'blue'}, 'Everest Buffer', 1, 0.2);
Map.addLayer(everestPoint, {color: 'red'}, 'The Everest');

// Define the time period
var startDate = '2023-01-01';
var endDate = '2023-03-31';

// Load Sentinel-2 Surface Reflectance imagery.
var s2 = ee.ImageCollection("COPERNICUS/S2_SR")
    .filterBounds(everestBuffer)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

// Cloud masking using the SCL band.
function maskS2clouds(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(3).and(scl.neq(8)).and(scl.neq(9));
  return image.updateMask(mask);
}
var s2Masked = s2.map(maskS2clouds);

// Create a median composite image and clip to the study area.
var composite = s2Masked.median().clip(everestBuffer);

// Compute the NDSI using Sentinel-2 bands (B3: Green, B11: SWIR) and clip.
var ndsi = composite.normalizedDifference(['B3', 'B11']).rename('NDSI').clip(everestBuffer);

// Create a binary snow mask (adjust the threshold as needed).
var ndsiThreshold = 0.45;
var snowMask = ndsi.gt(ndsiThreshold).clip(everestBuffer);

// Compute a local snow cover percentage via a neighborhood mean.
var kernel = ee.Kernel.square({radius: 10, units: 'pixels'});
var snowFraction = snowMask.reduceNeighborhood({
  reducer: ee.Reducer.mean(),
  kernel: kernel
}).rename('snowFraction');
var snowPercentage = snowFraction.multiply(100).rename('snowPercentage');

// Classify the snow cover into four classes: Class 1: 0–25%, Class 2: 25–50%, Class 3: 50–75%, Class 4: 75–100%.
var classImage = snowPercentage.expression(
  "(b('snowPercentage') < 25) ? 1" +
  " : (b('snowPercentage') < 50) ? 2" +
  " : (b('snowPercentage') < 75) ? 3" +
  " : 4"
).rename('snowClass').clip(everestBuffer);

// For reference, add layers to the map.
Map.addLayer(composite, {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000}, 'Composite');
Map.addLayer(classImage, {min: 1, max: 4, palette: ['#D3D3D3', '#CCCCFF', '#4169E1', '#E0FFFF']}, 'Snow Classification');

// Merge route segments by unique route name.
var uniqueNames = routes.aggregate_array('NAME').distinct();
var mergedRoutes = ee.FeatureCollection(uniqueNames.map(function(routeName) {
  var parts = routes.filter(ee.Filter.eq('NAME', routeName));
  // Merge geometries by taking the union of all parts.
  var mergedGeometry = parts.geometry();
  return ee.Feature(mergedGeometry, {'NAME': routeName});
}));
print("Merged Routes:", mergedRoutes);

// Define a function to calculate the dangerous index for each merged route.
// Dangerous index is defined as a weighted average of the snow classes: dangerous_index = (1*f1 + 2*f2 + 3*f3 + 4*f4) / (f1+f2+f3+f4)
function calculateDangerousIndex(route) {
  // Buffer the merged route to capture adjacent conditions.
  var routeBuffer = route.geometry().buffer(50);
  
  // Compute the frequency histogram of snow classes within the route.
  var histDict = classImage.reduceRegion({
    reducer: ee.Reducer.frequencyHistogram(),
    geometry: routeBuffer,
    scale: 10,
    maxPixels: 1e9
  }).get('snowClass');
  
  histDict = ee.Dictionary(histDict);

  var f1 = ee.Number(histDict.get('1', 0));
  var f2 = ee.Number(histDict.get('2', 0));
  var f3 = ee.Number(histDict.get('3', 0));
  var f4 = ee.Number(histDict.get('4', 0));
  
  var total = f1.add(f2).add(f3).add(f4);

  var dangerousIndex = ee.Algorithms.If(total.eq(0),
    999,
    f1.multiply(1).add(f2.multiply(2))
      .add(f3.multiply(3)).add(f4.multiply(4)).divide(total)
  );
  
  return route.set({ dangerous_index: dangerousIndex });
}

// Calculate the dangerous index for each merged route.
var routesWithDanger = mergedRoutes.map(calculateDangerousIndex);

print("Routes with Dangerous Index:", routesWithDanger);

// Identify the safest route (lowest dangerous index).
var safestRoute = routesWithDanger.sort('dangerous_index').first();
var minDangerIndex = ee.Number(safestRoute.get('dangerous_index'));
print('Safest Route Dangerous Index:', minDangerIndex);

// Color the safest route in green and the others in red.
var tolerance = 1e-6;

var styledRoutes = routesWithDanger.map(function(feature) {
  var index = ee.Number(feature.get('dangerous_index'));
  var isSafest = index.subtract(minDangerIndex).abs().lt(tolerance);
  var color = ee.Algorithms.If(isSafest, 'green', 'red');
  return feature.set({styleColor: { color: color, width: 3, fillColor: '00000000' } });
});

// Use the styleColor property for stroke color when visualizing.
var routeVis = styledRoutes.style({
  styleProperty: 'styleColor'
  });
  
Map.addLayer(routeVis, {}, 'Routes (Safest in Green)');