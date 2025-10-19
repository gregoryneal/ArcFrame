namespace ArcFrame.Core.Geometry.Splines
{
    /// <summary>
    /// Class for evaluating a cubic spline. This sits inside the curvature law for a spline curve.
    /// </summary>
    public sealed class CubicSpline1D
    {
        readonly double[] x, a, b, c, d; // intervals ([x[i], x[i+1]): S_i(t)=a + b*τ + c*τ^2 + d*τ^3
        /// <summary>
        /// Create the spline and solve for cubic parameters.
        /// </summary>
        /// <param name="xs">Knots (sub interval positions along the arc length of the curve)</param>
        /// <param name="ys">Curvature offset per knot</param>
        /// <exception cref="ArgumentException"></exception>
        public CubicSpline1D(double[] xs, double[] ys)
        {
            if (xs.Length != ys.Length || xs.Length < 2) throw new ArgumentException();
            int n = xs.Length;
            x = (double[])xs.Clone();
            a = (double[])ys.Clone();
            b = new double[n - 1];
            c = new double[n];
            d = new double[n - 1];

            // Natural cubic via tridiagonal solve for c (second derivs)            
            var h = new double[n - 1];      //knot lengths

            for (int i = 0; i < n - 1; i++)
            {
                h[i] = x[i + 1] - x[i];
            }

            var alpha = new double[n];

            for (int i = 1; i < n - 1; i++)
            {
                alpha[i] = 3 * ((a[i + 1] - a[i]) / h[i] - (a[i] - a[i - 1]) / h[i - 1]);
            }

            var l = new double[n];
            var mu = new double[n];
            var z = new double[n];
            l[0] = 1;
            mu[0] = z[0] = 0;
            for (int i = 1; i < n - 1; i++)
            {
                l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
                mu[i] = h[i] / l[i];
                z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
            }
            l[n - 1] = 1; z[n - 1] = c[n - 1] = 0;
            for (int j = n - 2; j >= 0; j--)
            {
                c[j] = z[j] - mu[j] * c[j + 1];
                b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (2 * c[j] + c[j + 1]) / 3.0;
                d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
            }
        }
        int Find(double s)
        {
            int i = Array.BinarySearch(x, s);
            return i >= 0 ? System.Math.Max(0, i - 1) : System.Math.Max(0, ~i - 1);
        }
        public double Eval(double s)
        {
            int i = Find(s);
            double t = s - x[i];
            return a[i] + b[i] * t + c[i] * t * t + d[i] * t * t * t;
        }

        public void EvalJet(double s, int order, out double v0, out double v1, out double v2, out double v3)
        {
            int i = Find(s);
            double t = s - x[i];
            double ai = a[i], bi = b[i], ci = c[i], di = d[i];
            v0 = ai + bi * t + ci * t * t + di * t * t * t;
            v1 = bi + 2 * ci * t + 3 * di * t * t;
            v2 = 2 * ci + 6 * di * t;
            v3 = 6 * di;
        }
    }

}
