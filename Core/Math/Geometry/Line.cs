using ArcFrame.Core.Math;
using ArcFrame.Core.Results;

namespace ArcFrame.Core.Geometry
{
    /// <summary>
    /// A line class for N dimensional space.
    /// </summary>
    public sealed class Line : IArcLengthCurve
    {
        private readonly double[] _a, _t;
        private readonly double[,] R; //the orthonormal basis vectors of the moving frame, doesn't change on the line.
        private readonly double _len;
        private readonly int _n;

        public Line(double[] a, double[] b)
        {
            _n = a.Length;
            _a = (double[])a.Clone();
            //distance between a and b

            var d = Helpers.Subtract(b, a);// new double[_n];
            double ss = 0;
            for (int i = 0; i < _n; i++)
            {
                d[i] = b[i] - a[i];
                ss += d[i] * d[i];
            }
            _len = System.Math.Sqrt(ss);
            _t = new double[_n];
            if (_len > 0)
            {
                for (int i = 0; i < _n; i++)
                {
                    _t[i] = d[i] / _len;
                }
            }
            else _t[0] = 1;

            double[] T = Helpers.Normalize(Helpers.Subtract(b, a));
            R = ONFrame.R0_FromT_Complete(T);
        }

        public static Line From2D(double x0, double z0, double x1, double z1)
        {
            return new Line(new[] { x0, 0, z0 }, new[] { x1, 0, z1 });
        }

        public static Line From3D(double x0, double y0, double z0, double x1, double y1, double z1)
        {
            return new Line(new[] { x0, y0, z0 }, new[] { x1, y1, z1 });
        }

        public int Dimension => _n;

        public double Length => _len;
        public double[] Position(double s) => Evaluate(s).P;
        public double[] Tangent(double s) => Evaluate(s).T;

        public Sample Evaluate(double s)
        {
            if (s < 0) s = 0;
            else if (s > _len) s = _len;
            var p = new double[_n];
            for (int i = 0; i < _n; i++)
            {
                p[i] = _a[i] + (_t[i] * s);
            }
            return new Sample(p, (double[,])R.Clone(), s, new double[System.Math.Max(0, _n - 1)]);
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
    }
}
