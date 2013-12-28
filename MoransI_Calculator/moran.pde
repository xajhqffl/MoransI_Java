public static class Moran{
  
  /**
  *
  * calculate spatial autocorelation value Moran's I 
  * It uses the formula from Arcgis and wikipedia webstes as references.
  * http://en.wikipedia.org/wiki/Moran's_I
  * http://resources.arcgis.com/en/help/main/10.1/index.html#/How_Spatial_Autocorrelation_Global_Moran_s_I_works/005p0000000t000000/
  * 
  * It gets the same result as ArcGIS Moran's I toolbox (use Inverse Distance as spatial weight measure)
  * 
  * @param attributes one array of the attributes to be computed for Moran's I
  * @param weights  results from inverse spatial weights
  * @param neighbors results from inverse spatial weights
  * @return double[2]  the first parameter is the Moran's I, the second the parameter is the p-value
  *
  */
  public static double[] I(double[] attribtues,double[][] weights,int[][] neighbors){
    
      //1.calculate mean value of the attributes
    double mean = calcMean(attribtues);
    int n = attribtues.length;
    
    
    double sum1=(double)0,sum2=(double)0,sum3=(double)0,s1=(double)0.0,s2=(double)0.0;
    
    //2.calculate moran's I based on wiki formula http://en.wikipedia.org/wiki/Moran%27s_I    http://resources.arcgis.com/en/help/main/10.1/index.html#/How_Spatial_Autocorrelation_Global_Moran_s_I_works/005p0000000t000000/
    for(int i=0;i<weights.length;i++)
    {
      sum3+=Math.pow((attribtues[i]-mean),2);
      
      double tempwij_all=0.0;
      double tempwji_all=0.0;
      
      for(int j=0;j<weights.length;j++)
      {
        double zi=attribtues[i]-mean;
        double zj=attribtues[j]-mean;
        double wij = get_spatialweights_fromNeighbors(weights,neighbors,i,j);
        double wji = get_spatialweights_fromNeighbors(weights,neighbors,j,i);
        
        //self weights equal zero
        if(i==j)
          wij=0;
        
        sum1+=zi*zj*wij;
        sum2+=wij;
        
        s1+=(wij+wji)*(wij+wji);
        
        tempwij_all+=wij;
        tempwji_all+=wji;
        
      }
      
      s2+=(tempwij_all+tempwji_all)*(tempwij_all+tempwji_all);
    }
//    println(attribtues.length);
//    println(sum1);
//    println(sum2);
//    println(sum3);
    double morani = n/sum2 * sum1/sum3;
    double moranEI = -1. / (n - 1);
    
    //--------calculate p value------------ 
    //-------refers to the souce code in pysal library-----
    //--------http://pysal.org/library/esda/moran.html---
    //--------pysal.esda.moran.Moran--------
    s1/=2;
//    println("s0 " + sum2);
//    println("s1 " + s1);
//    println("s2 " + s2);
    double s0=sum2;
    double v_num = n * n * s1 - n * s2 + 3 * s0 * s0;
    double v_den = (n - 1) * (n + 1) * sum2 * sum2;
    double VI_norm = v_num / v_den - (1.0 / (n - 1)) * (1.0 / (n - 1));
    double seI_norm = Math.pow(VI_norm,0.5);
    double z_norm = (morani - moranEI) / seI_norm ;
    
    NormalDistribution n_distr = new NormalDistribution(0,1);
    double p_norm = 2.0 * (1 - n_distr.cumulativeProbability(Math.abs(z_norm)));
    
    return new double[]{morani,p_norm};
    
  }

  /**
  *
  * get spatial weights from i's j neighbor
  *
  */
  private static double get_spatialweights_fromNeighbors(double[][] weights,int[][] neighbors,int host,int neighbor)
  {
    //println(host);
    for(int i=0;i<neighbors[host].length;i++)
    {
      if(neighbors[host][i] == neighbor)
        return weights[host][i];
    }
    
    return 0;
  }
  

}





public class SpatialWeight{
  
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
  public double[] calcSOMCellLoc(int index,int xdim,int ydim)
  {
    double deltaY=(1.0*Math.pow(0.75,0.5));
    float deltaX=1.0;
    double xloc;
      
    int rowindex = index/xdim;//row index starts from 0, so the first row would be 0
    int columnindex = index-rowindex*xdim;//the column index starts from 0
    
    if(rowindex%2 == 0)
      xloc = columnindex*deltaX+0;
    else
      xloc = columnindex*deltaX+0.5;
    double yloc = rowindex*deltaY+0;
    
    return new double[]{xloc,yloc};
    
  }
  
  
  /**
  * Caculate Inverse spatial weights within a given distance threshold
  * refer to pysal python library calculations of spatial weigths
  * http://www.pysal.org/library/weights/weights.html#weights-spatial-weights
  * The result of the spatial weights have two arrays, one describe the neighborhood configuration for each point, and the other array stores the corresponding weights for each neighbor
  *
  * @author:Fangming Du 
  * @param points array, each dimension has two values, x and y
  * @param threshold filter out irrelevant neighbors
  * @param neighbors empty int arrays for storing neighbors
  * @param weights empty double array for storing weigths
  * @param rowStandardization whether apply row standardization for weights
  * @return 
  */
  public void  calc_inverse_spatial_weigths(double[][] points,float threshold,int[][] neighbors,double[][] weights,boolean rowStandardization)
  {
    //iterate through each points
    int len = points.length;
    //neighbors = new int[len][];
    //weights
    for(int i=0;i<len;i++)
    {
      ArrayList tempweights= new ArrayList();
      ArrayList tempneibhors=new ArrayList();
      for(int j=0;j<len;j++)
      {
        if(i==j)
          continue;
        double dist = calcDist(points[i][0],points[i][1],points[j][0],points[j][1]);
        if(dist<threshold)
        {
          tempweights.add(1/dist);
          tempneibhors.add(j);
        }
      }
      
      neighbors[i] = convertInts(tempneibhors);
      weights[i] = convertDoubles(tempweights);
      
      if(rowStandardization)
      {
        double sum_weight = calcSum(weights[i]);
        for(int iter=0;iter<weights[i].length;iter++)
        {
          weights[i][iter] = weights[i][iter]/sum_weight;
        }
      }
      
    }
  }
  
  
  /**
  *
  * write spatial neighbors and weights to .gwt file format, which is sparse matrix file. defined at http://www.biomedware.com/files/documentation/spacestat/data/export/Spatial_Weights_Files.htm
  * It's one of the most common used file format for storing spatial weights
  *
  * @param neighbors the spatial configuration for neighors
  * @param weights the spatial weights for each neighbor
  * @param fielPath the path for the file to be saved
  */
  public void writeWeightsFile(int[][] neighbors,double[][] weights,String filePath)
  {
    ArrayList results = new ArrayList();
    results.add("1\t"+neighbors.length);
    
    for(int i=0;i<neighbors.length;i++)
    {
      //has neighbors
      if(neighbors[i].length>0)
      {
        for(int j=0;j<neighbors[i].length;j++)
        {
          results.add(i+"\t"+neighbors[i][j]+"\t"+weights[i][j]);
        }
      }
    }
    
    saveStrings(filePath+".gwt",convertStrings(results));
  }
  
  
  /**
  *
  * read spatial neighbors and weights to arrays, which is sparse matrix file. defined at http://www.biomedware.com/files/documentation/spacestat/data/export/Spatial_Weights_Files.htm
  *
  * @param neighbors the spatial configuration for neighors
  * @param weights the spatial weights for each neighbor
  * @param fielPath the path for the file to be read
  */
  public Object[] readWeightsFile(String filePath)
  {
    String[] lines = loadStrings(filePath);
    int num = Integer.parseInt(lines[0].split("\t")[1]);
    //println(num);
    int[][] neighbors = new int[num][];
    double[][] weights = new double[num][];
    int oldhost = Integer.parseInt(lines[1].split("\t")[0]);
    
    //println(neighbors[neighbors.length-1]);
    
    ArrayList temp_neighbor = new ArrayList();
    ArrayList temp_weight = new ArrayList();
    
    for(int i=1;i<lines.length;i++)
    {
      String[] data = lines[i].split("\t");
      int host = Integer.parseInt(data[0]);
      
      
      if(oldhost == host)
      {
        temp_neighbor.add(Integer.parseInt(data[1]));
        temp_weight.add(Double.parseDouble(data[2]));
      }
      else
      {
        neighbors[oldhost] = convertInts(temp_neighbor);
        weights[oldhost] = convertDoubles(temp_weight);
        
        temp_neighbor.clear();
        temp_weight.clear();
        
        temp_neighbor.add(Integer.parseInt(data[1]));
        temp_weight.add(Double.parseDouble(data[2]));
        
        oldhost=host;
      }
    }
    
    neighbors[oldhost] = convertInts(temp_neighbor);
    weights[oldhost] = convertDoubles(temp_weight);
    
    return new Object[]{neighbors, weights};
   
  }
  
}


/**
*
* convert arraylist of double type and int type to arrays
*
*/
static double[] convertDoubles(ArrayList doubles)
{
    double[] ret = new double[doubles.size()];
    for (int i = 0; i < ret.length; i++)
    {
        ret[i] = (Double)(doubles.get(i));
    }
    return ret;
}
static int[] convertInts(ArrayList ints)
{
    int[] ret = new int[ints.size()];
    for (int i = 0; i < ret.length; i++)
    {
        ret[i] = (Integer)(ints.get(i));
    }
    return ret;
}
static String[] convertStrings(ArrayList strings)
{
    String[] ret = new String[strings.size()];
    for (int i = 0; i < ret.length; i++)
    {
        ret[i] = (String)(strings.get(i));
    }
    return ret;
}


//=============== UTILITY FUNCTIONS ================

/**
*
* calculate distance between two points
*
*/
static double calcDist(double x1,double y1,double x2,double y2){
  double x_pow2 = Math.pow((x1-x2),2); //<>//
  double y_pow2 = Math.pow((y1-y2),2);
  //println(x1+" " + y1 + " " + x2 + " " +y2);
  //println(x_pow2 + " "+y_pow2); //<>//
  //println(Math.sqrt(x_pow2+y_pow2));
  return Math.sqrt(x_pow2+y_pow2);
  //return (float)(Math.sqrt(+Math));
}

/**
*
* calculate mean values of double arrays
*
*/
static double calcMean(double[] values)
{
  double sum = (double)0;
  for(int i=0;i<values.length;i++)
    sum += values[i];
    
  return sum/values.length;
  
}
/**
*
* calculate the sum values of double arrays
*
*/
static double calcSum(double[] values)
{
  double sum = (double)0;
  for(int i=0;i<values.length;i++)
  {
    sum+= values[i];
  }
  
  return sum;
}
