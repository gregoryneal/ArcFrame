using System;
using ArcFrame.Core.Math;

namespace ArcFrame.Core.Geometry.Splines
{
    /// <summary>
    /// For extrinsic splines defined by the position of their nodes in space, we map the arc length to the spline parameterization [0, L] -> [0, 1]
    /// </summary>
    public sealed class ArcLengthTable
    {
        readonly double[] t, s; // t in [0,1], s cumulative in [0,L]
        /// <summary>
        /// Total length of the cached polyline.
        /// </summary>
        public double Length => s[^1];

        /// <summary>
        /// Create the mapping from [0, L] -> [0, 1]
        /// </summary>
        /// <param name="pos">A function mapping the spline parameter in [0, 1] to the coordinate position along that spline.</param>
        /// <param name="samples">Number of samples to take (larger is slower but more accurate as usual)</param>
        public ArcLengthTable(Func<double, double[]> pos, int samples = 512)
        {
            t = new double[samples + 1];
            s = new double[samples + 1];
            double accum = 0;
            double[] prev = pos(0);
            for (int i = 0; i <= samples; i++)
            {
                double ti = (double)i / samples;
                t[i] = ti;
                var pi = pos(ti);
                if (i > 0)
                {
                    accum += Helpers.Len(Helpers.Subtract(prev, pi)); //distance between prev and pi 
                }
                s[i] = accum;
                prev = pi;
            }
            // normalize to arc length (optionally refine with Simpson/RK)
        }

        /// <summary>
        /// Maps s in [0, L] to [0, 1]
        /// </summary>
        /// <param name="sQuery"></param>
        /// <returns></returns>
        public double MapStoT(double sQuery)
        {
            if (sQuery > Length) sQuery = Length;
            if (sQuery < 0) sQuery = 0;
            return sQuery / Length;
        }

        /// <summary>
        /// Maps arc length s in [0, L] to the normalized spline parameter t in [0, 1].
        /// </summary>
        /// <param name="sQuery">The arc length (s) value to map.</param>
        /// <returns>The normalized spline parameter t corresponding to the given arc length.</returns>
        public double MapSToTUniform(double sQuery)
        {
            // Clamp sQuery to be within the valid arc length range [0, Length]
            if (sQuery < 0) sQuery = 0;
            if (sQuery > Length) sQuery = Length;

            // Binary search for the segment in the arc length table
            int low = 0, high = s.Length - 1;
            while (low < high)
            {
                int mid = (low + high) / 2;
                if (s[mid] < sQuery)
                    low = mid + 1;
                else
                    high = mid;
            }

            // Low will be the index of the smallest s[i] that is greater than or equal to sQuery
            // Interpolate between s[low-1] and s[low]
            if (low == 0)
                return t[0]; // Edge case for the first element
            else if (s[low] == sQuery)
                return t[low]; // Exact match, no interpolation needed
            else
            {
                // Interpolate between the two closest points
                double sPrev = s[low - 1];
                double sNext = s[low];
                double tPrev = t[low - 1];
                double tNext = t[low];

                // Linear interpolation of t for the given sQuery
                double tInterpolated = tPrev + (sQuery - sPrev) * (tNext - tPrev) / (sNext - sPrev);
                return tInterpolated;
            }
        }

        /// <summary>
        /// Export the currently generated arc length mapping from s to t
        /// </summary>
        /// <param name="arcS"></param>
        /// <param name="arcT"></param>
        public void Export(out double[] arcS, out double[] arcT)
        {
            arcS = s;
            arcT = t;
        }
    }
}
