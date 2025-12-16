using System;

namespace ArcFrame.Core.Math.Geometry.Splines
{
    /// <summary>
    /// Cubic hermite spline
    /// </summary>
    public class HermiteSpline : Spline
    {
        /// <summary>
        /// Packed layout: [Pi, Pi+1, Mi, Mi+1] per segment ⇒ 4 columns per segment
        /// </summary>
        public override int ComputeSegmentCount() => System.Math.Max(1, ControlPoints.Length / 4);

        /// <inheritdoc/>
        public override int GetControlPointStartForSegment(int segIndex)
        {
            // in precomputed _gbCache, use blocks 0, 4, 8, ..
            return 4 * segIndex;
        }

        /// <summary>
        /// Pack the points and tangents into a single array since tangents are used as control points in the matrix multiplication.
        /// </summary>
        /// <param name="points"></param>
        /// <param name="tangents"></param>
        /// <param name="frame"></param>
        /// <param name="fastMode"></param>
        /// <param name="cacheSamplesOverride"></param>
        /// <exception cref="ArgumentException"></exception>
        public HermiteSpline(double[][] points, 
            double[][] tangents, 
            FrameModel frame = FrameModel.Frenet,
            bool fastMode = false,
            int cacheSamplesOverride = 0)
        : base(Pack(points, tangents), ConstructBasis(), frame, fastMode, cacheSamplesOverride)
        {
            if (tangents.Length != points.Length)
                throw new ArgumentException("HermiteSpline: tangents length must match points length.");
        }

        /// <summary>
        /// Use control points that are already pre packed.
        /// </summary>
        /// <param name="ptsAndTangents"></param>
        /// <param name="frame"></param>
        /// <param name="fastMode"></param>
        /// <param name="cacheSamplesOverride"></param>
        public HermiteSpline(double[][] ptsAndTangents, 
            FrameModel frame = FrameModel.Frenet,
            bool fastMode = false,
            int cacheSamplesOverride = 0)
        : base(ptsAndTangents, ConstructBasis(), frame, fastMode, cacheSamplesOverride) { }

        /// <summary>
        /// Locate the start index of the control point array given the parameter t in [0, 1].
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        /*
        protected override (int k, double u) Locate(double t)
        {
            t = System.Math.Clamp(t, 0, 1);

            // Packed layout: [Pi, Pi+1, Mi, Mi+1] per segment ⇒ 4 columns per segment
            int segCount = ComputeSegmentCount();

            // Map to segment index in [0..segCount-1] and local u in [0,1)
            double seg = t * segCount;

            // Keep the last sample on the last segment with u=1 handled by evaluation, not by rolling over
            int i = System.Math.Min(segCount - 1, (int)System.Math.Floor(seg));
            double u = seg - i;

            // Start column index for this segment's 4-tuple
            int cpi = 4 * i;

            return (cpi, u);
        }*/

        /// <summary>
        /// Pack the points and tangents together.
        /// </summary>
        /// <param name="pts"></param>
        /// <param name="tans"></param>
        /// <returns></returns>
        private static double[][] Pack(double[][] pts, double[][] tans)
        {
            int segs = pts.Length - 1;
            var packed = new double[segs * 4][];
            for (int i = 0; i < segs; i++)
            {
                packed[4 * i + 0] = pts[i];
                packed[4 * i + 1] = pts[i + 1];
                packed[4 * i + 2] = tans[i];
                packed[4 * i + 3] = tans[i + 1];
            }
            return packed;
        }

        private static double[,] ConstructBasis()
        {
            // Right-multiply basis for p = G · B · [1,u,u^2,u^3]^T with G=[P_i, P_{i+1}, M_i, M_{i+1}]
            return new double[,]
            {
                { 1,  0, -3,  2 },  // h00
                { 0,  0,  3, -2 },  // h01
                { 0,  1, -2,  1 },  // h10
                { 0,  0, -1,  1 }   // h11
            };
        }
    }
}
