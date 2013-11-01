package hex;

import static org.junit.Assert.*;
import hex.gram.Gram.InPlaceCholesky;

import java.io.File;
import java.util.concurrent.ExecutionException;
import java.util.Random;
import java.lang.Math;
import org.junit.Test;

import Jama.CholeskyDecomposition;
import Jama.Matrix;

import water.util.Utils;
import water.TestUtil;



public class CholTest extends TestUtil{

  @Test public void testPar1Blk10 () {
    System.out.println("Running testPar1Blk10");
    DataSetup data = new DataSetup(1000, 12345);
    new TestSetup(data,1,10).test();
  } 
  
  private final static class DataSetup {
    public double xx[][];
    public DataSetup(int N, int rseed) {
      xx = new double[N][];
      Random r = new Random(rseed);
      for (int i = 0; i < N; i++) {
        xx[i] = new double[i+1];
        for (int j = 0; j < i; j++)
        xx[i][j] = r.nextGaussian();
        double x = r.nextGaussian();
        xx[i][i] = Math.max(x, -x);
      }
    }
  }

  private final static class TestSetup {
    DataSetup data;
    int par; 
    int blk;
    public TestSetup(DataSetup data, int par, int blk) {
      this.data = data; this.par = par; this.blk = blk;
    }
    public void test() {
      double[][] jamaChol = new Matrix(data.xx).chol().getL().getArray();
      double[][] chol = new InPlaceCholesky(data.xx, par, blk).getL();
      assertEquals(jamaChol.length, chol.length);
      for (int i = 0; i < chol.length; i++) {
        assertEquals(jamaChol[i].length, chol[i].length);
        for( int j = 0; j < chol[i].length; j++)
          assertEquals(jamaChol[i][j], chol[i][j], 0.0001);
      }
    }
  }
}

