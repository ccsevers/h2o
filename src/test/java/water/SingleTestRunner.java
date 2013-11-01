package water;

import java.io.File;
import java.io.IOException;
import java.lang.annotation.*;
import java.lang.reflect.Method;
import java.net.ServerSocket;
import java.util.*;

import org.apache.commons.lang.ArrayUtils;
import org.junit.Ignore;
import org.junit.Test;
import org.junit.runner.Result;
import org.junit.runner.notification.Failure;

import water.H2O;
import water.deploy.Node;
import water.deploy.NodeVM;
import water.parser.ParseFolderTestBig;
import water.util.Log;
import water.util.Utils;

public class SingleTestRunner {
  // TODO
  @Retention(RetentionPolicy.RUNTIME)
  public @interface Nightly {
  }

  public static void main(String[] args) throws Exception {
    // Can be necessary to run in parallel to other clouds, so find open ports
    String flat = "";
    // Force all IPs to local so that users can run with a firewall
    String[] a = new String[] { "-ip", "127.0.0.1", "-flatfile", Utils.writeFile(flat).getAbsolutePath() };
    H2O.OPT_ARGS.ip = "127.0.0.1";
    args = (String[]) ArrayUtils.addAll(a, args);

    ArrayList<Node> nodes = new ArrayList<Node>();
//     nodes.add(new NodeVM(args));

    args = Utils.append(new String[] { "-mainClass", Master.class.getName() }, args);
    Node master = new NodeVM(args);
    nodes.add(master);

    File out = null, err = null, sandbox = new File("sandbox");
    sandbox.mkdirs();
    Utils.clearFolder(sandbox);
    for( int i = 0; i < nodes.size(); i++ ) {
      out = File.createTempFile("junit-" + i + "-out-", null, sandbox);
      err = File.createTempFile("junit-" + i + "-err-", null, sandbox);
      nodes.get(i).persistIO(out.getAbsolutePath(), err.getAbsolutePath());
      nodes.get(i).start();
    }

    int exit = master.waitFor();
    if( exit != 0 ) {
      Log.log(out, System.out);
      Thread.sleep(100); // Or mixed (?)
      Log.log(err, System.err);
    }
    for( Node node : nodes )
      node.kill();
    if( exit == 0 )
      System.out.println("OK");
    System.exit(exit);
  }

  private static boolean isOpen(int port) throws Exception {
    ServerSocket s = null;
    try {
      s = new ServerSocket(port);
      return true;
    } catch( IOException ex ) {
      return false;
    } finally {
      if( s != null )
        s.close();
    }
  }

  public static class Master {
    public static void main(String[] args) {
      try {
        List<Class> tests = new ArrayList<Class>();
        boolean classArg = false;
        String[] realArgs = args;
        for (int i = 0; i < args.length; i++) {
          String arg = args[i];
          if (arg.equals("--")) { 
            classArg = true; 
            realArgs = new String[i]; System.arraycopy(args, 0, realArgs, 0, i); continue; 
          }
          if (classArg) tests.add(Class.forName(arg)); 
        }
        
        H2O.main(realArgs);
        TestUtil.stall_till_cloudsize(1);
        
        Result r = org.junit.runner.JUnitCore.runClasses(tests.toArray(new Class[0]));
        if( r.getFailureCount() == 0 ) {
          System.out.println("Successfully ran the following tests in " + (r.getRunTime() / 1000) + "s");
          for( Class c : tests )
            System.out.println(c.getName());
        } else {
          for( Failure f : r.getFailures() ) {
            System.err.println(f.getDescription());
            if( f.getException() != null )
              f.getException().printStackTrace();
          }
        }
        System.exit(r.getFailureCount());
      } catch( Throwable t ) {
        t.printStackTrace();
        System.exit(1);
      }
    }
  }
}