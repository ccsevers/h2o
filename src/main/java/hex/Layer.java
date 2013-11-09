package hex;

import hex.rng.MersenneTwisterRNG;

import java.util.Random;

import water.*;
import water.api.DocGen;
import water.api.Request.API;
import water.fvec.Chunk;
import water.fvec.Vec;

/**
 * Neural network layer.
 *
 * @author cypof
 */
public abstract class Layer extends Iced {
  static final int API_WEAVER = 1;
  public static DocGen.FieldDoc[] DOC_FIELDS;

  @API(help = "Number of neurons")
  public int units;

  @API(help = "Learning rate")
  public float rate;

  @API(help = "Learning rate annealing")
  public float rate_annealing;

  @API(help = "L2 regularisation")
  public float l2;

  // TODO disabled for now, not enough testing
  @ParamsSearch.Info(origin = 1)
  float _momentum;
  float _momentumAnnealing;

  // TODO
  float _perWeight;
  float _perWeightAnnealing;

  // Current rate and momentum
  transient float _r, _m;

  // Weights, biases, activity, error
  transient float[] _w, _b, _a, _e;

  // Last weights & per-weight rate data
  transient float[] _wPrev, _wInit, _wMult;
  transient float[] _bPrev, _bInit, _bMult;

  transient float[] _wSpeed, _bSpeed;

  // Previous layer
  transient Layer _in;

  // Optional visible units bias, e.g. for pre-training
  transient float[] _v, _gv;

  public final void init(Layer[] ls, int index) {
    init(ls, index, true, 0);
  }

  public void init(Layer[] ls, int index, boolean weights, long step) {
    _a = new float[units];
    _e = new float[units];
    _in = ls[index - 1];

    if( weights ) {
      _w = new float[units * _in.units];
      _b = new float[units];

      // deeplearning.net tutorial (TODO special ones for rectifier & softmax?)
      // TODO only subset of inputs?
      Random rand = new MersenneTwisterRNG(MersenneTwisterRNG.SEEDS);
      float min = (float) -Math.sqrt(6. / (_in.units + units));
      float max = (float) +Math.sqrt(6. / (_in.units + units));
      for( int i = 0; i < _w.length; i++ )
        _w[i] = rand(rand, min, max);
    }

    if( _momentum != 0 ) {
      _wPrev = new float[_w.length];
      _bPrev = new float[_b.length];
      for( int i = 0; i < _w.length; i++ )
        _wPrev[i] = _w[i];
      _wSpeed = new float[_w.length];
      _bSpeed = new float[_b.length];
    }

    if( _perWeight != 0 ) {
      _wInit = new float[_w.length];
      _wMult = new float[_w.length];
      for( int i = 0; i < _w.length; i++ ) {
        _wInit[i] = _w[i];
        _wMult[i] = 1;
      }
      _bInit = new float[_b.length];
      _bMult = new float[_b.length];
      for( int i = 0; i < _b.length; i++ ) {
        _bInit[i] = _b[i];
        _bMult[i] = 1;
      }
    }

    anneal(step);
  }

  abstract void fprop();

  void dropout_fprop() {
      this.fprop();
  }

  abstract void bprop();

  public final void anneal(long n) {
    _r = rate(n);
    _m = _momentum * (n + 1) / ((n + 1) + _momentumAnnealing);
  }

  public float rate(long n) {
    return rate / (1 + rate_annealing * n);
  }

  public final void momentum(long n) {
    for( int i = 0; i < _w.length; i++ )
      adjust(i, _w, _wPrev, _wInit, _wMult);

    for( int i = 0; i < _b.length; i++ )
      adjust(i, _b, _bPrev, _bInit, _bMult);
  }

  private final void adjust(int i, float[] w, float[] prev, float[] init, float[] mult) {
    float coef = 1;

    if( init != null ) {
      float g = w[i] - init[i];
      boolean sign = g > 0;
      boolean last = mult[i] > 0;
      coef = Math.abs(mult[i]);
      // If the gradient kept its sign, increase
      if( sign == last ) {
        if( coef < 4 )
          coef += _perWeight;
      } else
        coef *= 1 - _perWeight;
      mult[i] = sign ? coef : -coef;
      w[i] = init[i] + coef * g;
    }

    if( prev != null ) {
      // Nesterov's Accelerated Gradient
      float v = (w[i] - prev[i]) * _m;
      prev[i] = w[i];
      w[i] += coef * v;
      if( w == _w )
        _wSpeed[i] = v;
      else
        _bSpeed[i] = v;
    }

    if( init != null )
      init[i] = w[i];
  }

  public static abstract class Input extends Layer {
    long _pos, _len;

    @Override public void init(Layer[] ls, int index, boolean weights, long step) {
      _a = new float[units];
    }

    @Override void bprop() {
      throw new UnsupportedOperationException();
    }

    public final long move() {
      return _pos = _pos == _len - 1 ? 0 : _pos + 1;
    }
  }

  public static class VecsInput extends Input {
    static final int API_WEAVER = 1;
    public static DocGen.FieldDoc[] DOC_FIELDS;

    public Vec[] vecs;

    @API(help = "Categoricals identified on the training set")
    int[] categoricals;

    @API(help = "Normalisation stats used during training")
    float[] subs, muls;

    transient Chunk[] _chunks;

    private float prob;
    private Random rand;
    private float scale;


    VecsInput() {
    }

    public VecsInput(Vec[] vecs, VecsInput stats, boolean dropout, float prob) {
      assert prob <= 1.0f;
      assert prob > 0 ;
      this.prob = prob;
      this.scale = 1/prob;
      this.rand = new MersenneTwisterRNG(MersenneTwisterRNG.SEEDS);
      units = stats != null ? stats.subs.length : expand(vecs);
      this.vecs = vecs;
      _len = vecs[0].length();

      if( stats != null ) {
        assert stats.categoricals.length == vecs.length;
        categoricals = stats.categoricals;
        assert stats.subs.length == units;
        subs = stats.subs;
        muls = stats.muls;
      } else {
        categoricals = new int[vecs.length];
        for( int i = 0; i < vecs.length; i++ )
          categoricals[i] = categories(vecs[i]);
        subs = new float[units];
        muls = new float[units];
        stats(vecs);
      }
    }

    public VecsInput(Vec[] vecs, VecsInput stats) {
      this(vecs, stats, false, 1.0f);
    }

    static int categories(Vec vec) {
      if( vec.domain() == null )
        return 1;
      return (int) (vec.max() - vec.min());
    }

    static int expand(Vec[] vecs) {
      int n = 0;
      for( int i = 0; i < vecs.length; i++ )
        n += categories(vecs[i]);
      return n;
    }

    private void stats(Vec[] vecs) {
      Stats stats = new Stats();
      stats._units = units;
      stats._categoricals = categoricals;
      stats.doAll(vecs);
      for( int i = 0; i < vecs.length; i++ ) {
        subs[i] = (float) stats._means[i];
        double sigma = Math.sqrt(stats._sigms[i] / (stats._rows - 1));
        muls[i] = (float) (sigma > 1e-6 ? 1 / sigma : 1);
      }
    }

    @Override void fprop() {
      if( _chunks == null )
        _chunks = new Chunk[vecs.length];
      for( int i = 0; i < vecs.length; i++ ) {
        Chunk c = _chunks[i];
        if( c == null || c._vec != vecs[i] || _pos < c._start || _pos >= c._start + c._len )
          _chunks[i] = vecs[i].chunk(_pos);
      }
      ChunksInput.set(_chunks, _a, (int) (_pos - _chunks[0]._start), subs, muls, categoricals);
    }

    @Override void dropout_fprop() {
      if( _chunks == null )
          _chunks = new Chunk[vecs.length];
      for( int i = 0; i < vecs.length; i++ ) {
        Chunk c = _chunks[i];
        if( c == null || c._vec != vecs[i] || _pos < c._start || _pos >= c._start + c._len )
            _chunks[i] = vecs[i].chunk(_pos);
      }
      ChunksInput.set(_chunks, _a, (int) (_pos - _chunks[0]._start), subs, muls, categoricals);
      for(int o = 0; o < _a.length; o++)
        _a[o] *= (rand.nextFloat() < prob ? 1 : 0) * scale;
      }
  }

  /**
   * Stats with expanded categoricals.
   */
  static class Stats extends MRTask2<Stats> {
    int _units;
    int[] _categoricals;
    double[] _means, _sigms;
    long _rows;
    transient float[] _subs, _muls;

    @Override protected void setupLocal() {
      _subs = new float[_units];
      _muls = new float[_units];
      for( int i = 0; i < _muls.length; i++ )
        _muls[i] = 1;
    }

    @Override public void map(Chunk[] cs) {
      _means = new double[_units];
      _sigms = new double[_units];
      float[] a = new float[_means.length];
      for( int r = 0; r < cs[0]._len; r++ ) {
        ChunksInput.set(cs, a, r, _subs, _muls, _categoricals);
        for( int c = 0; c < a.length; c++ )
          _means[c] += a[c];
      }
      for( int c = 0; c < a.length; c++ )
        _means[c] /= cs[0]._len;
      for( int r = 0; r < cs[0]._len; r++ ) {
        ChunksInput.set(cs, a, r, _subs, _muls, _categoricals);
        for( int c = 0; c < a.length; c++ )
          _sigms[c] += (a[c] - _means[c]) * (a[c] - _means[c]);
      }
      _rows += cs[0]._len;
    }

    @Override public void reduce(Stats rs) {
      reduce(_means, _sigms, _rows, rs._means, rs._sigms, rs._rows);
      _rows += rs._rows;
    }

    static void reduce(double[] ma, double[] sa, long ra, double[] mb, double[] sb, long rb) {
      for( int c = 0; c < ma.length; c++ ) {
        double delta = ma[c] - mb[c];
        ma[c] = (ma[c] * ra + mb[c] * rb) / (ra + rb);
        sa[c] = sa[c] + sb[c] + delta * delta * ra * rb / (ra + rb);
      }
    }

    @Override public boolean logVerbose() {
      return !H2O.DEBUG;
    }
  }

  static class ChunksInput extends Input {
    transient Chunk[] _chunks;
    float[] _subs, _muls;
    int[] _categoricals;

    public ChunksInput(Chunk[] chunks, VecsInput stats) {
      units = stats.subs.length;
      _chunks = chunks;
      _subs = stats.subs;
      _muls = stats.muls;
      _categoricals = stats.categoricals;
    }

    @Override void fprop() {
      set(_chunks, _a, (int) _pos, _subs, _muls, _categoricals);
    }

    static void set(Chunk[] chunks, float[] a, int row, float[] subs, float[] muls, int[] categoricals) {
      int n = 0;
      for( int i = 0; i < categoricals.length; i++ ) {
        double d = chunks[i].at0(row);
        d = Double.isNaN(d) ? 0 : d;
        if( categoricals[i] == 1 ) {
          d -= subs[n];
          d *= muls[n];
          a[n++] = (float) d;
        } else {
          int cat = categoricals[i];
          for( int c = 0; c < cat; c++ )
            a[n + c] = -subs[n + c];
          int c = (int) d - (int) chunks[i]._vec.min() - 1;
          if( c >= 0 )
            a[n + c] = (1 - subs[n + c]) * muls[n + c];
          n += cat;
        }
      }
      assert n == a.length;
    }
  }

  public static abstract class Output extends Layer {
    Input _input;

    @Override public void init(Layer[] ls, int index, boolean weights, long step) {
      super.init(ls, index, weights, step);
      _input = (Input) ls[0];
    }
  }

  public static abstract class Softmax extends Output {
    abstract int label();

    @Override void fprop() {
      float max = Float.NEGATIVE_INFINITY;
      for( int o = 0; o < _a.length; o++ ) {
        _a[o] = 0;
        for( int i = 0; i < _in._a.length; i++ )
          _a[o] += _w[o * _in._a.length + i] * _in._a[i];
        _a[o] += _b[o];
        if( max < _a[o] )
          max = _a[o];
      }
      float scale = 0;
      for( int o = 0; o < _a.length; o++ ) {
        _a[o] = (float) Math.exp(_a[o] - max);
        scale += _a[o];
      }
      for( int o = 0; o < _a.length; o++ )
        _a[o] /= scale;
    }

    @Override void bprop() {
      int label = label();
      for( int o = 0; o < _a.length; o++ ) {
        float t = o == label ? 1 : 0;
        float e = t - _a[o];
        // Gradient is error * derivative of Softmax: (1 - x) * x
        float g = e * (1 - _a[o]) * _a[o];
        for( int i = 0; i < _in._a.length; i++ ) {
          int w = o * _in._a.length + i;
          _in._e[i] += g * _w[w];
          _w[w] += _r * (g * _in._a[i] - _w[w] * l2);
        }
        _b[o] += _r * g;
      }
    }
  }

  public static class VecSoftmax extends Softmax {
    static final int API_WEAVER = 1;
    public static DocGen.FieldDoc[] DOC_FIELDS;

    public Vec vec;

    @API(help = "Min response value on the training set")
    int min;

    VecSoftmax() {
    }

    public VecSoftmax(Vec vec, VecSoftmax stats) {
      this.units = stats != null ? stats.units : (int) (vec.max() - vec.min() + 1);
      this.vec = vec;
      this.min = stats != null ? stats.min : (int) vec.min();
    }

    @Override int label() {
      return (int) vec.at8(_input._pos) - min;
    }
  }

  static class ChunkSoftmax extends Softmax {
    transient Chunk _chunk;
    int _min;

    public ChunkSoftmax(Chunk chunk, VecSoftmax stats) {
      units = stats.units;
      _chunk = chunk;
      _min = stats.min;

      // TODO extract layer info in separate Ice
      rate = stats.rate;
      rate_annealing = stats.rate_annealing;
      _momentum = stats._momentum;
      _momentumAnnealing = stats._momentumAnnealing;
      _perWeight = stats._perWeight;
      _perWeightAnnealing = stats._perWeightAnnealing;
      l2 = stats.l2;
    }

    @Override int label() {
      return (int) _chunk.at80((int) _input._pos) - _min;
    }
  }

  public static abstract class Linear extends Output {
    abstract float value();

    @Override void fprop() {
      for( int o = 0; o < _a.length; o++ ) {
        _a[o] = 0;
        for( int i = 0; i < _in._a.length; i++ )
          _a[o] += _w[o * _in._a.length + i] * _in._a[i];
        _a[o] += _b[o];
      }
    }

    @Override void bprop() {
      float v = value();
      for( int o = 0; o < _a.length; o++ ) {
        float g = v - _a[o];
        for( int i = 0; i < _in._a.length; i++ ) {
          int w = o * _in._a.length + i;
          _in._e[i] += g * _w[w];
          _w[w] += _r * (g * _in._a[i] - _w[w] * l2);
        }
        _b[o] += _r * g;
      }
    }
  }

  public static class VecLinear extends Linear {
    Vec _vec;

    VecLinear() {
    }

    public VecLinear(Vec vec, VecLinear stats) {
      this.units = stats != null ? stats.units : 1;
      _vec = vec;
    }

    @Override float value() {
      double d = _vec.at(_input._pos);
      return Double.isNaN(d) ? 0 : (float) d;
    }
  }

  static class ChunkLinear extends Linear {
    transient Chunk _chunk;

    public ChunkLinear(Chunk chunk, VecLinear stats) {
      units = stats.units;
      _chunk = chunk;

      // TODO extract layer info in separate Ice
      rate = stats.rate;
      rate_annealing = stats.rate_annealing;
      _momentum = stats._momentum;
      _momentumAnnealing = stats._momentumAnnealing;
      _perWeight = stats._perWeight;
      _perWeightAnnealing = stats._perWeightAnnealing;
      l2 = stats.l2;
    }

    @Override float value() {
      double d = _chunk.at0((int) _input._pos);
      return Double.isNaN(d) ? 0 : (float) d;
    }
  }

  public static class Tanh extends Layer {
        Tanh() {
        }

        public Tanh(int units) {
            this.units = units;
        }

        @Override void fprop() {
            for( int o = 0; o < _a.length; o++ ) {
                _a[o] = 0;
                for( int i = 0; i < _in._a.length; i++ )
                    _a[o] += _w[o * _in._a.length + i] * _in._a[i];
                _a[o] += _b[o];

                // tanh approx, slightly faster, untested
                // float a = Math.abs(_a[o]);
                // float b = 12 + a * (6 + a * (3 + a));
                // _a[o] = (_a[o] * b) / (a * b + 24);

                _a[o] = (float) Math.tanh(_a[o]);
            }
        }

        @Override void bprop() {
            for( int o = 0; o < _a.length; o++ ) {
                // Gradient is error * derivative of hyperbolic tangent: (1 - x^2)
                float g = _e[o] * (1 - _a[o] * _a[o]);
                for( int i = 0; i < _in._a.length; i++ ) {
                    int w = o * _in._a.length + i;
                    if( _in._e != null )
                        _in._e[i] += g * _w[w];
                    _w[w] += _r * (g * _in._a[i] - _w[w] * l2);
                }
                _b[o] += _r * g;
            }
        }
    }

    public static class Maxout extends Layer {
        Maxout() {
        }

        private float prob;
        private Random rand;
        private float scale;

        public Maxout(int units, float prob) {
            this.units = units;
            assert prob <= 1;
            assert prob > 0 ;
            this.prob = prob;
            this.scale = 1/prob;
            this.rand = new MersenneTwisterRNG(MersenneTwisterRNG.SEEDS);
        }

        @Override void fprop() {
            for( int o = 0; o < _a.length; o++ ) {
                _a[o] = Float.MIN_VALUE;
                float cur = 0.0f;
                for( int i = 0; i < _in._a.length; i++ )
                    cur = _w[o * _in._a.length + i] * _in._a[i];
                    _a[o] = java.lang.Math.max(cur, _a[o]);
                _a[o] += _b[o];
            }
        }

        @Override void dropout_fprop() {
            for( int o = 0; o < _a.length; o++ ) {
                _a[o] = Float.MIN_VALUE;
                float cur = 0.0f;
                for( int i = 0; i < _in._a.length; i++ )
                    // is it better to pre-fill the mask with rands?
                    cur = _w[o * _in._a.length + i] * _in._a[i];
                    _a[o] = java.lang.Math.max(cur, _a[o]);
                _a[o] += _b[o];
                _a[o] *= (rand.nextFloat() < prob ? 1 : 0) * scale;
            }
        }

        @Override void bprop() {
            for( int o = 0; o < _a.length; o++ ) {
                float g = _e[o];
//                if( _a[o] < 0 )   Not sure if we should be using maxout with a hard zero bottom
//                    g = 0;
                for( int i = 0; i < _in._a.length; i++ ) {
                    int w = o * _in._a.length + i;
                    if( _in._e != null )
                        _in._e[i] += g * _w[w];
                    _w[w] += _r * (g * _in._a[i] - _w[w] * l2);
                }
                _b[o] += _r * g;
            }
        }
    }


  /**
   * Apply tanh to the weights' transpose. Used for auto-encoders.
   */
  public static class TanhPrime extends Layer {
    TanhPrime() {
    }

    public TanhPrime(int units) {
      this.units = units;
    }

    @Override void fprop() {
      for( int o = 0; o < _a.length; o++ ) {
        _a[o] = 0;
        for( int i = 0; i < _in._a.length; i++ )
          _a[o] += _w[i * _in._a.length + o] * _in._a[i];
        _a[o] += _b[o];
        _a[o] = (float) Math.tanh(_a[o]);
      }
    }

    @Override void bprop() {
      for( int o = 0; o < _a.length; o++ ) {
        assert _in._in.units == units;
        float g = _in._in._a[o] - _a[o];
        // TODO derivative?
        for( int i = 0; i < _in._a.length; i++ ) {
          int w = i * _in._a.length + o;
          if( _in._e != null )
            _in._e[i] += g * _w[w];
          _w[w] += _r * (g * _in._a[i] - _w[w] * l2);
        }
        _b[o] += _r * g;
      }
    }
  }

  public static class Rectifier extends Layer {
    Rectifier() {
    }

    public Rectifier(int units) {
      this.units = units;
    }

    @Override public void init(Layer[] ls, int index, boolean weights, long step) {
      super.init(ls, index, weights, step);

      if( weights ) {
        Random rand = new MersenneTwisterRNG(MersenneTwisterRNG.SEEDS);
        int count = Math.min(15, _in.units);
        float min = -.1f, max = +.1f;
        for( int o = 0; o < units; o++ ) {
          for( int n = 0; n < count; n++ ) {
            int i = rand.nextInt(_in.units);
            int w = o * _in.units + i;
            _w[w] = rand(rand, min, max);
          }
        }
        for( int i = 0; i < _b.length; i++ )
          _b[i] = 1;
        for( int i = 0; _v != null && i < _v.length; i++ )
          _v[i] = 1;
      }
    }

    @Override void fprop() {
      for( int o = 0; o < _a.length; o++ ) {
        _a[o] = 0;
        for( int i = 0; i < _in._a.length; i++ )
          _a[o] += _w[o * _in._a.length + i] * _in._a[i];
        _a[o] += _b[o];
        if( _a[o] < 0 )
          _a[o] = 0;
      }
    }

    @Override void bprop() {
      for( int o = 0; o < _a.length; o++ ) {
        float g = _e[o];
        if( _a[o] < 0 )
          g = 0;
        for( int i = 0; i < _in._a.length; i++ ) {
          int w = o * _in._a.length + i;
          if( _in._e != null )
            _in._e[i] += g * _w[w];
          _w[w] += _r * (g * _in._a[i] - _w[w] * l2);
        }
        _b[o] += _r * g;
      }
    }
  }

  //

  @Override public Layer clone() {
    return (Layer) super.clone();
  }

  public static void shareWeights(Layer src, Layer dst) {
    dst._w = src._w;
    dst._b = src._b;
  }

  public static void shareWeights(Layer[] src, Layer[] dst) {
    for( int y = 1; y < src.length; y++ )
      shareWeights(src[y], dst[y]);
  }

  // If layer is a RBM

  /**
   * TODO inject noise in units <br>
   * mean 0 and variance 1 / ( 1 + e-x )
   */
  void contrastiveDivergence(float[] in) {
//    float[] v1 = in;
//    float[] h1 = new float[_b.length];
//    fprop(v1, h1);
//    float[] v2 = generate(h1);
//    float[] h2 = new float[_b.length];
//    fprop(v2, h2);

//    for( int o = 0; o < _b.length; o++ )
//      for( int i = 0; i < _v.length; i++ )
//        _gw[o * _v.length + i] += _rate * ((h1[o] * v1[i]) - (h2[o] * v2[i]));
//
//    for( int o = 0; o < _gb.length; o++ )
//      _gb[o] += _rate * (h1[o] - h2[o]);

//    for( int i = 0; i < _gv.length; i++ )
//      _gv[i] += _rate * (v1[i] - v2[i]);
  }

  final void adjustVisible() {
    if( _gv != null ) {
      for( int v = 0; v < _gv.length; v++ ) {
        _v[v] += _gv[v];
        _gv[v] *= 1 - _momentum;
      }
    }
  }

  float[] generate(float[] hidden) {
    assert hidden.length == _b.length;
    float[] visible = new float[_v.length];
    for( int o = 0; o < hidden.length; o++ )
      for( int i = 0; i < _in._a.length; i++ )
        visible[i] += _w[o * _in._a.length + i] * hidden[o];
    for( int i = 0; i < visible.length; i++ ) {
      visible[i] += _v[i];
      if( visible[i] < 0 )
        visible[i] = 0;
    }
    return visible;
  }

  float error(float[] in1) {
    float[] out1 = new float[_b.length];
//    fprop(in1, out1);
    float[] in2 = generate(out1);
    float error = 0;
    for( int i = 0; i < in1.length; i++ ) {
      float d = in2[i] - in1[i];
      error += d * d;
    }
    return error;
  }

  float freeEnergy(float[] in) {
    float energy = 0.0f;
    for( int i = 0; i < in.length; i++ )
      energy -= in[i] * _v[i];
    for( int o = 0; o < _b.length; o++ ) {
      float out = 0;
      for( int i = 0; i < in.length; i++ )
        out += _w[o * in.length + i] * in[i];
      out += _b[o];
      energy -= Math.log(1 + Math.exp(out));
    }
    return energy;
  }

  private static float rand(Random rand, float min, float max) {
    return min + rand.nextFloat() * (max - min);
  }

  @Override public AutoBuffer writeJSON(AutoBuffer bb) {
    bb.put1('{');
    bb.putJSONStr("type").put1(':').putJSONStr(getClass().getName());
    bb.put1(',');
    writeJSONFields(bb);
    bb.put1('}');
    return bb;
  }
}
