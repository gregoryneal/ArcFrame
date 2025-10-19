using ArcFrame.Core.Math;
using ArcFrame.Core.Results;

namespace ArcFrame.Core.Geometry.Splines
{
    /// <summary>
    /// Disguise a catmullrom spline as an arc length parameterized curve.
    /// Fake it till you make it, this is probably not all that accurate compared to 
    /// most generators. Make sure you duplicate the end point.
    /// 
    /// TODO: Fix, this does not do what i want it to do at all. It is probably the arc length table
    /// </summary>
    public sealed class CatmullRomCurve : IArcLengthCurve
    {
        /// <summary>
        /// Our spline points
        /// </summary>
        readonly double[][] P;
        readonly double alpha;
        readonly ArcLengthTable table;
        public int Dimension => P[0].Length;
        public double Length => table.Length;

        public CatmullRomCurve(IReadOnlyList<double[]> points, double alpha = 0.5, int samples = 1024)
        {
            if (points.Count < 4) throw new ArgumentException("Need ≥4 points (with end duplicates if open).");
            this.alpha = alpha;
            P = new double[points.Count][];
            for (int i = 0; i < P.Length; i++)
            {
                P[i] = (double[])points[i].Clone();
            }
            table = new ArcLengthTable(t => Pos(t), samples);
        }

        // Evaluate basis over segment k with local u ∈ [0,1]
        (int k, double u) Locate(double t)
        {
            t = System.Math.Max(0, System.Math.Min(1, t));
            double seg = t * (P.Length - 3); // segments = n-3
            int k = System.Math.Min(P.Length - 4, (int)System.Math.Floor(seg));
            return (k, seg - k);
        }

        double[] Pos(double t)
        {
            var (k, u) = Locate(t);
            var p0 = P[k];
            var p1 = P[k + 1];
            var p2 = P[k + 2];
            var p3 = P[k + 3];
            // centripetal parameterization
            double[] ts = new double[4];
            ts[0] = 0;
            ts[1] = ts[0] + System.Math.Pow(Helpers.Len(Helpers.Subtract(p0, p1)), alpha);
            ts[2] = ts[1] + System.Math.Pow(Helpers.Len(Helpers.Subtract(p1, p2)), alpha);
            ts[3] = ts[2] + System.Math.Pow(Helpers.Len(Helpers.Subtract(p2, p3)), alpha);
            double tau = ts[1] + u * (ts[2] - ts[1]);

            var A1 = Helpers.Lerp(p0, p1, (tau - ts[0]) / (ts[1] - ts[0]));
            var A2 = Helpers.Lerp(p1, p2, (tau - ts[1]) / (ts[2] - ts[1]));
            var A3 = Helpers.Lerp(p2, p3, (tau - ts[2]) / (ts[3] - ts[2]));
            var B1 = Helpers.Lerp(A1, A2, (tau - ts[1]) / (ts[2] - ts[1]));
            var B2 = Helpers.Lerp(A2, A3, (tau - ts[2]) / (ts[3] - ts[2]));
            return Helpers.Lerp(B1, B2, (tau - ts[2]) / (ts[3] - ts[2]));
        }

        public Sample Evaluate(double s)
        {
            double t = table.MapStoT(s);
            // finite diff for derivatives (simple and robust)
            double dt = 1E-4;
            double t0 = System.Math.Max(0, t - dt);
            double t1 = System.Math.Min(1, t + dt);
            var p = Pos(t);
            var p0 = Pos(t0);
            var p1 = Pos(t1);
            var T = Helpers.Normalize(Helpers.Subtract(p1, p0));
            // curvature magnitude κ1 in ND: || (I - T T^T) p'' || / ||p'||^2  ~ via finite differences
            var v1 = Helpers.Multiply(Helpers.Subtract(p1, p), 1.0 / (t1 - t));
            var v0 = Helpers.Multiply(Helpers.Subtract(p, p0), 1.0 / (t - t0));
            var aApprox = Helpers.Multiply(Helpers.Subtract(v1, v0), 1.0 / ((t1 - t0) / 2.0));
            var K = new double[System.Math.Max(0, Dimension - 1)];
            K[0] = Helpers.Len(Helpers.Reject(aApprox, T)) / System.Math.Pow(Helpers.Len(Helpers.Multiply(Helpers.Subtract(p1, p0), 1.0 / (t1 - t0))), 2);

            //read this
            //https://people.engr.tamu.edu/sueda/courses/CSCE450/2023F/labs/L01/index.html
            throw new NotImplementedException();
            //return new Sample(p, T, s, K);
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