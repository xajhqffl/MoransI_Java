//Step 1. location generator based on SOM neuron index

//Step 2. create inverse distance matrix of all the neurons

import org.jblas.*;

import java.util.Random;
import java.lang.Math;

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
println(neighbors[0]);
println(weights[0]);
println(Moran.I(attr,weights,neighbors));
}

double calcMoransI(double[] attribtues,double[][] weights,int[][] neighbors)
{
  //1.calculate mean value of the attributes
  double mean = calcMean(attribtues);
  
  float sum1=(float)0,sum2=(float)0,sum3=(float)0;
  
  //2.calculate moran's I based on wiki formula http://en.wikipedia.org/wiki/Moran%27s_I    http://resources.arcgis.com/en/help/main/10.1/index.html#/How_Spatial_Autocorrelation_Global_Moran_s_I_works/005p0000000t000000/
  for(int i=0;i<weights.length;i++)
  {
    sum3+=Math.pow((attribtues[i]-mean),2);
    
    for(int j=0;j<weights.length;j++)
    {
      double zi=attribtues[i]-mean;
      double zj=attribtues[j]-mean;
      double wij = get_spatialweights_fromNeighbors(weights,neighbors,i,j);
      
      //self weights equal zero
      if(i==j)
        wij=0;
      
      sum1+=zi*zj*wij;
      sum2+=wij;
      
    }
  }
  println(attribtues.length);
  println(sum1);
  println(sum2);
  println(sum3);
  float morani = attribtues.length/sum2 * sum1/sum3;
  return morani;
  
}

double get_spatialweights_fromNeighbors(double[][] weights,int[][] neighbors,int host,int neighbor)
  {
    for(int i=0;i<neighbors[host].length;i++)
    {
      if(neighbors[host][i] == neighbor)
        return weights[host][i];
    }
    
    return 0;
  }

