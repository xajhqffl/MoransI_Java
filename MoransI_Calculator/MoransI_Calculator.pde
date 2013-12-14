//Step 1. location generator based on SOM neuron index

//Step 2. create inverse distance matrix of all the neurons

import org.jblas.*;

import java.util.Random;
import java.lang.Math;

void setup(){
  
//read SOM training file
println(calcSOMCellLoc(1910,100,100));
String[] lines = loadStrings("C:\\Users\\Fangming\\Documents\\thesis\\new_map\\data\\Ph3_100k.cod");
double[] attr = new double[10000];
float[] x = new float[10000];
float[] y = new float[10000];
for(int i=1;i<lines.length;i++)
{
  println(i);
  String[] attributes = lines[i].split(" ");
  attr[i-1]=Double.parseDouble(attributes[0]);
  println(attr[i-1]);
  
  float[] xy = calcSOMCellLoc(i,100,100);
  x[i-1]=xy[0];
  y[i-1]=xy[1];
  println(x[i-1]);
}
FloatMatrix dm = getInverseDistMatrix(x,y);

//testing distance calculators
//println(calcDist((double)0,(double)0,(double)2,(double)2));

//testing inverse distance matrix calculator
//double[] x = new double[]{34.13583,34.17611,33.82361,34.19944,34.06694,33.92917,34.01500,34.06722,34.08333,34.38750};
//double[] y = new double[]{-117.9236,-118.3153,-118.1875,-118.5347,-117.7514,-118.2097,-118.0597,-118.2264,-118.1069,-118.5347};
//DoubleMatrix dm = getInverseDistMatrix(x,y);
////dm.print();
//
//double[] attr = new double[]{7.225806,5.899194,4.052885,7.181452,6.076613,3.157258,5.201613,4.717742,6.532258,7.540323};
////println(calcMean(attr));
//
//println(calcMoransI(attr,dm));
}

double calcMoransI(double[] attribtues,DoubleMatrix weights)
{
  //1.calculate mean value of the attributes
  double mean = calcMean(attribtues);
  
  double sum1=(double)0,sum2=(double)0,sum3=(double)0;
  
  //2.calculate moran's I based on wiki formula http://en.wikipedia.org/wiki/Moran%27s_I    http://resources.arcgis.com/en/help/main/10.1/index.html#/How_Spatial_Autocorrelation_Global_Moran_s_I_works/005p0000000t000000/
  for(int i=0;i<weights.getRows();i++)
  {
    for(int j=0;j<weights.getColumns();j++)
    {
      double zi=attribtues[i]-mean;
      double zj=attribtues[j]-mean;
      double wij = weights.get(i,j);
      
      sum1+=zi*zj*wij;
      sum2+=wij;
      sum3+=Math.pow(zi,2);
      
    }
  }
  
  double morani = attribtues.length/sum2 * sum1/sum3;
  return morani;
  
}

FloatMatrix getInverseDistMatrix(float[] X,float[] Y)
{
  FloatMatrix inverseDistMx = FloatMatrix.zeros(X.length,X.length);
  println("initialized");
   
  if(X.length != Y.length){
    println("ERROR!Points do not have equal number of x, y coordinates!");
    return null;
  }
  else{
    for(int i=0;i<X.length;i++){
      for(int j=i+1;j<X.length;j++){
        float dist = 1/calcDist(X[i],Y[i],X[j],Y[j]);
        inverseDistMx.put(i,j,dist);
        inverseDistMx.put(j,i,dist);
        println(dist);
       // dist[j] = calcDist(X[i],Y[i]);
      }
    }
  }
   
  return inverseDistMx;
}

float calcDist(float x1,float y1,float x2,float y2){
  return (float)(Math.sqrt(Math.pow((x1-x2),2)+Math.pow((y1-y2),2)));
}

double calcMean(double[] values)
{
  double sum = (double)0;
  for(int i=0;i<values.length;i++)
    sum += values[i];
    
  return sum/values.length;
  
}

/**
* The function calcuates the SOM cell X,Y locations in the spatial space by its index.
* This calculation also follows the same rule how SOM cells are created in the SOM Analyst toolbox in ArcGIS.
*
* @author:Fangming Du
* @param index the index of SOM cells,here index is defined to start from 1, which simplified the calculation
* @param xdim the dimensions of x direction, horizontal direction
* @param ydim the dimensions of y direction, vertial direction
* @return array of x y locations, array[0] x location, array[1] y location
*
*/
float[] calcSOMCellLoc(int index,int xdim,int ydim)
{
  float deltaY=(float)(1.0*Math.pow(0.75,0.5));
  float deltaX=1.0;
    
  int rowindex = index/xdim;//row index starts from 0, so the first row would be 0
  int columnindex = index-rowindex*xdim;//the column index starts from 1
  
  float xloc = (columnindex-1)*deltaX+0+rowindex%2*0.5;
  float yloc = rowindex*deltaY+0;
  
  return new float[]{xloc,yloc};
  
}
