import java.lang.Math;
import org.apache.commons.math3.distribution.NormalDistribution;

void setup(){
  
String[] lines = loadStrings("/Users/fangmingdu/Documents/mallet_training/Ph3_100k.cod");
double[][] points = new double[10000][2];
SpatialWeight sp = new SpatialWeight();

double[] attr = new double[10000];
for(int i=1;i<lines.length;i++)
{
  String[] attributes = lines[i].split(" ");
  attr[i-1]=Double.parseDouble(attributes[100]);

  double[] xy = sp.calcSOMCellLoc(i-1,100,100);
  points[i-1]=xy;
}
//double[][] points=new double[][]{{10, 10}, {20, 10}, {40, 10}, {15, 20}, {30, 20}, {30, 30}};

/*----------------- COMPUTE spatial weights and SAVE to file ----------------------------*/
//int[][] neighbors = new int[10000][2];
//double[][] weights = new double[10000][2];
//sp.calc_inverse_spatial_weigths(points,1.5,neighbors,weights,true);
//sp.writeWeightsFile(neighbors,weights,"test");

/*----------------- READ Spaital weights From files ----------------------------*/
Object[] obj = sp.readWeightsFile("/Users/fangmingdu/MoransI_Java/MoransI_Calculator/test.gwt");
int[][] neighbors = (int[][])obj[0];
double[][] weights = (double[][])obj[1];

/*----------------- COMPUTE Moran's I and p-value ----------------------------*/
println(Moran.I(attr,weights,neighbors));

}
