﻿using ArcFrame.Core.Math;

namespace ArcFrame.Core.Geometry.Splines
{
    /// <summary>
    /// For extrinsic splines defined by the position of their nodes in space, we map the arc length to the spline parameterization [0, L] -> [0, 1]
    /// </summary>
    public sealed class ArcLengthTable
    {
        readonly double[] t, s; // t in [0,1], s cumulative in [0,L]
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

        public double MapStoT(double sQuery)
        {
            if (sQuery > Length) sQuery = Length;
            if (sQuery < 0) sQuery = 0;
            return sQuery / Length;
        }
    }

}
