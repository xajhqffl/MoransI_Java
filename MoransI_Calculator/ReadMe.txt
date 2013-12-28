1. Calculation of Moran’s I index
There are two components in calculating Moran’s I. The first one is the spatial weight matrix of the input spatial features. There are many different ways computing the spatial weights. This program provides the calculation of inverse spatial weight based on point features. The second component is the calcualation of Moran’s I using spatial weights, points locations and one attribute data for those locations. The calculation formula is available online, such as Wikipedia.

So first, with input points, you can use SpatialWeight.calc_inverse_spatial_weigths(points,threshold,neighbors,weights,ifNormalize)
to calculate spatial weights. It returns the spatial neighbors and weights into the arrays.
Then you can use Moran.I(attributes,weights,neighbors) to calculate Moran’s I. It returns Moran’s I index and p-value.

2. Spatial weights class
The SpatialWeight class also provides functionalities such as writing and reading spatial weights files. 
So after you calculate spatial weights once, you can store the spatial weights into a gwt file (spatial weight file). It can be read later into this program and reused to calculate Moran’s I for the same features, but different attributes.

3. Utility function for SOM cells
This program also has a utility function for converting SOM cell ID to spatial x,y locations. It uses the same coordinates from the SOMAnalyst toolbox.
