using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;
using ArcFrame.Core.Results;

namespace ArcFrame.Core.Experimental
{
    /// <summary>
    /// Embeds a d dimensional curve into R^N with a basis that can vary with arc length s.
    /// Not sure if this has a purpose or not. 
    /// </summary>
    public sealed class FrameFieldPromotedCurve : IArcLengthCurve
    {
        private readonly IArcLengthCurve _inner;   // dim d
        private readonly Func<double, double[,]> _E;   // N x d ON columns at s
        private readonly Func<double, double[]> _p0;  // N offset at s
        private readonly int _N, _d;

        // NOTE: This simple version assumes the inner parameter s is "close enough"
        // to the output arc-length parameter. For true unit-speed output, you'd
        // precompute a lookup of s_out(s) by integrating ||E(s) T_in(s) + E'(s) P_in(s) + p0'(s)||,
        // then invert it (monotone interpolation) inside Evaluate.
        public FrameFieldPromotedCurve(
            IArcLengthCurve inner,
            Func<double, double[,]> E_of_s,
            Func<double, double[]> p0_of_s)
        {
            _inner = inner ?? throw new ArgumentNullException(nameof(inner));
            _E = E_of_s ?? throw new ArgumentNullException(nameof(E_of_s));
            _p0 = p0_of_s ?? (_ => new double[(_E(0.0)).GetLength(0)]);

            _d = inner.Dimension;
            var E0 = _E(0.0);
            _N = E0.GetLength(0);
            if (E0.GetLength(1) != _d) throw new ArgumentException("E(s) must be N x d.");
        }

        public int Dimension => _N;
        public double Length => _inner.Length; // NOT exact if E, p0 vary; see note above.

        public Sample Evaluate(double s)
        {
            var q = _inner.Evaluate(s);           // q in R^d
            var E = _E(s);                        // N x d
            var p = _p0(s);                       // N

            var P = new double[_N];
            //var T = new double[_N];
            var R = Helpers.Multiply(E, q.R);
            for (int i = 0; i < _N; i++)
            {
                double px = 0, tx = 0;
                for (int j = 0; j < _d; j++) 
                { 
                    px += E[i, j] * q.P[j]; 
                    tx += E[i, j] * q.T[j]; 
                }
                P[i] = px + p[i];
                //T[i] = tx; // if you need exact unit-speed, renormalize T after reparam
            }

            // Curvatures: preserve inner generalized curvatures; pad zeros
            var K = new double[System.Math.Max(0, _N - 1)];
            if (q.k != null) 
            {
                for (int j = 0; j < System.Math.Min(q.k.Length, K.Length); j++)
                {
                    K[j] = q.k[j];
                }
            }

            return new Sample(P, R, q.s, K);
        }

        public double[] Position(double s) => Evaluate(s).P;

        public double[] Tangent(double s) => Evaluate(s).T;
    }
}
