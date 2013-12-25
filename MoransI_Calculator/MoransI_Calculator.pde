import java.lang.Math;
import org.apache.commons.math3.distribution.NormalDistribution;

void setup(){
  
String[] lines = loadStrings("/Users/fangmingdu/Documents/mallet_training/Ph3_100k.cod");
double[][] points = new double[10000][2];

double[] attr = new double[10000];
//String[] output = new String[10000];
for(int i=1;i<lines.length;i++)
{
  String[] attributes = lines[i].split(" ");
  attr[i-1]=Double.parseDouble(attributes[0]);

  double[] xy = calcSOMCellLoc(i-1,100,100);
  points[i-1]=xy;
}
int[][] neighbors = new int[10000][];
double[][] weights = new double[10000][];
calc_inverse_spatial_weigths(points,1.5,neighbors,weights,true);
println(Moran.I(attr,weights,neighbors));
//NormalDistribution n = new NormalDistribution();
}
