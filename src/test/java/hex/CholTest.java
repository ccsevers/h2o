package hex;

import static org.junit.Assert.*;
import hex.gram.Gram.InPlaceCholesky;

import java.io.File;
import java.util.concurrent.ExecutionException;
import java.util.Random;
import java.util.Arrays;
import java.lang.Math;
import org.junit.Test;

import Jama.CholeskyDecomposition;
import Jama.Matrix;

import jsr166y.*;
import water.H2O;
import water.util.Utils;
import water.util.Log;
import water.TestUtil;

public class CholTest extends TestUtil{

  @Test public void testPar1Blk10 () {
    System.out.println("Running testPar1Blk10");
    for (int sz = 2000; sz < 4000; sz+=1000) {
      DataSetup data = new DataSetup(sz, 12345);
      new ForkJoinPool(2).invoke(new TestSetup(data,1,100));
    }
  } 
  
  private final static class DataSetup {
    public double xx[][];
    public DataSetup(int N, int rseed) {
      xx = new double[N][];
      Random r = new Random(rseed);
      double x[][] = new double[N][];
      for (int i = 0; i < N; i++) {
        x[i] = new double[N]; 
        for (int j = 0; j < N; j++)
          x[i][j] = r.nextGaussian();
      }
      for (int i = 0; i < N; i++) {
        xx[i] = new double[N];
        for (int j = 0; j < N; j++)
          xx[i][j] = 0;
      }
      for (int k = 0; k < N; k++)
        for (int i = 0; i < N; i++)
          for (int j = 0; j <= i; j++) {
            xx[i][j] += x[k][i]*x[k][j];
            xx[j][i] = xx[i][j];
          }
    }
  }

  private static  void print_matrix(double matrix[][]) {
    for (int i = 0; i < matrix.length; i++) {
      for (int j = 0; j <= i; j++)
        System.out.print(String.format("%.2f ",matrix[i][j]));
      System.out.println();
    }
  }

  private final static class TestSetup extends RecursiveAction {
    DataSetup data;
    int par; 
    int blk;
    public TestSetup(DataSetup data, int par, int blk) {
      this.data = data; this.par = par; this.blk = blk;
    }
    public void compute() {
      System.out.println("CREATING MATRIX FOR JAMA. ");
      double[][] jamaxx = new double[data.xx.length][];
      for (int i = 0; i < jamaxx.length; i++) 
        jamaxx[i] = Arrays.copyOfRange(data.xx[i], 0, jamaxx.length);
      long start = System.currentTimeMillis();
      CholeskyDecomposition jamaChol = new Matrix(jamaxx).chol();
      Log.err("JAMA CHOLESKY TAKES " + (System.currentTimeMillis() - start) + " MILLISECONDS.");
      if (!jamaChol.isSPD())  return;
      double[][] jama = jamaChol.getL().getArray();
      System.out.println("JAMA CHOLESKY -- "); // print_matrix(jama);
      start = System.currentTimeMillis();
      double[][] chol = new InPlaceCholesky(data.xx, par, blk).getL();
      Log.err("H2O CHOLESKY TAKES " + (System.currentTimeMillis() - start) + " MILLISECONDS.");
      System.out.println("H2O CHOLESKY -- "); // print_matrix(chol);
      assertEquals(jama.length, chol.length);
      for (int i = 0; i < chol.length; i++)
        for( int j = 0; j <= i; j++)
          assertEquals(jama[i][j], chol[i][j], 0.0001);
    }
  }
}

