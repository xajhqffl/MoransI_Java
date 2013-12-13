//Step 1. location generator based on SOM neuron index

//Step 2. create inverse distance matrix of all the neurons

import org.jblas.*;

import java.util.Random;
import java.lang.Math;

void setup(){
  
//  DoubleMatrix A = DoubleMatrix.randn(100,100);
//  A.transpose();
//test it in windows sourcetree

Random generator = new Random();
Double[] A = new Double[100];
Double[] B = new Double[100];
for(int i=0;i<100;i++){
  A[i] = generator.nextDouble();
  B[i] = generator.nextDouble();
}
  
}

Float[][] getInverseDistMatrix(Double[] X,Double[] Y)
{
  if(X.length != Y.length){
    println("ERROR!Points do not have equal number of x, y coordinates!");
    return null;
  }
  else{
    Double[] dist = new Double[X.length];
    for(int i=0;i<X.length;i++){
      for(int j=i+1;j<X.length;j++){
        dist[i] = calcDist(X[i],Y[i])
      }
    }
  }
   
  return new Float[][]{{4321.0232},{0.32}};
}

Double calcDist(Double x1,Double y1,Double x2,Double y2){
  return Math.sqrt(Math.pow((x1-x2),2)+Math.pow((y1-y2),2));
}
