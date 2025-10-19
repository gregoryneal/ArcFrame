using ArcFrame.Core.Results;

namespace ArcFrame.Core.Geometry
{
    /// <summary>
    /// Embeds the curve into a higher dimension N using axisMap (default: map i→i). No pose change.
    /// Normal use case this simply embeds the curve in a higher dimension without extra transformations.
    /// You could use axisMap to remap say the inner curves X axis to the PromotedCurves Z axis.
    /// </summary>
    public sealed class PromotedCurve : IArcLengthCurve
    {
        private readonly IArcLengthCurve _inner;
        private readonly int _N;
        private readonly int[] _map; // length = inner.Dimension, indices in [0,_N)

        public PromotedCurve(IArcLengthCurve inner, int targetN, int[]? axisMap = null)
        {
            _inner = inner; _N = targetN;
            if (targetN < inner.Dimension) throw new ArgumentException();
            _map = axisMap ?? DefaultMap(inner.Dimension, targetN);
        }
        public int Dimension => _N;
        public double Length => _inner.Length;
        public double[] Position(double s) => Evaluate(s).P;
        public double[] Tangent(double s) => Evaluate(s).T;

        public Sample Evaluate(double s)
        {
            var q = _inner.Evaluate(s);
            var P = new double[_N];
            var T = new double[_N];
            var R = new double[_N, _N];
            // place inner dims into selected axes, rest = 0
            for (int i = 0; i < _map.Length; i++)
            {
                P[_map[i]] = q.P[i];
                for (int j = 0; j < _map.Length; j++)
                {
                    R[_map[i], j] = q.R[i, j];
                }
            }

            // Place a 1 at each i,i position on the new rotation matrix, so it is SO(_N)
            for (int i = _map.Length; i < _N; i++)
            {
                R[i, i] = 1;
            }

            var K = new double[System.Math.Max(0, _N - 1)];
            // carry kappa_1; higher generalized curvatures stay 0 (planar inner stays planar)
            if (q.k != null && q.k.Length > 0 && K.Length > 0) K[0] = q.k[0];
            return new Sample(P, R, q.s, K);
        }

        /// <summary>
        /// Get a number of evenly spaced samples along the arc length
        /// </summary>
        /// <param name="count"></param>
        /// <returns></returns>
        public Sample[] GetSamples(int count)
        {
            Sample[] samples = new Sample[count];
            if (count <= 0) count = 1;
            double s;
            for (int i = 0; i < count; i++)
            {
                s = (i == count - 1) ? Length : i * (Length / System.Math.Max(1, count - 1));
                samples[i] = Evaluate(s);
            }
            return samples;
        }

        /// <summary>
        /// Default axis mapping for static embedding of axis.
        /// </summary>
        /// <param name="d"></param>
        /// <param name="N"></param>
        /// <returns>[0, 1, 2, ..., d-1]</returns>
        static int[] DefaultMap(int d, int N) { var m = new int[d]; for (int i = 0; i < d; i++) m[i] = i; return m; }
    }
}
