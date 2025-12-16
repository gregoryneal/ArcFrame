using System;

namespace ArcFrame.Core.Math.Geometry.Splines
{
    /// <summary>
    /// Create a Cubic Bezier curve out of control points.
    /// </summary>
    public class BezierSpline : Spline
    {
        /// <inheritdoc/>
        public override int ComputeSegmentCount()
        {
            int cps = ControlPoints.Length;
            if ((cps - 1) % 3 != 0)
                throw new ArgumentException("BezierSpline requires 3n+1 control points.");

            return (cps - 1) / 3;
        }

        /// <inheritdoc/>
        public override int GetControlPointStartForSegment(int segIndex)
        {
            // in precomputed _gbCache, use blocks 0, 3, 6, ..
            return 3 * segIndex;
        }

        /// <summary>
        /// Create the cubic bezier curve with a set of control points.
        /// </summary>
        /// <param name="controlPoints"></param>
        /// <param name="frame"></param>
        /// <param name="fastMode"></param>
        /// <param name="cacheSamplesOverride"></param>
        /// <exception cref="ArgumentException"></exception>
        public BezierSpline(double[][] controlPoints, 
                            FrameModel frame = FrameModel.Frenet,
                            bool fastMode = false,
                            int cacheSamplesOverride = 0) : base(controlPoints, ConstructBasis(), frame, fastMode, cacheSamplesOverride)
        {
            if ((controlPoints.Length - 1) % 3 != 0) throw new ArgumentException("BezierSpline requires 3n+1 control points.");
        }

        private static double[,] ConstructBasis()
        {
            return new double[,]
            {
                {1, -3, 3, -1},
                {0, 3, -6, 3 },
                {0, 0, 3, -3 },
                {0, 0, 0,  1 }
            };
        }

        /// <summary>
        /// Given a parameter t in [0, 1] return the start segment index for a 
        /// stride 3 control point window [pi, pi+1, pi+2, pi+3].
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        /*
        protected override (int k, double u) Locate(double t)
        {
            t = System.Math.Clamp(t, 0, 1);
            int segCount = ComputeSegmentCount();
            double seg = t * segCount;
            int i = System.Math.Min(segCount - 1, (int)System.Math.Floor(seg));
            return (3 * i, seg - i);
        }*/

        /// <summary>
        /// Create a stride 3 window of control points
        /// [pi, pi+1, pi+2, pi+3]
        /// </summary>
        /// <param name="cpi"></param>
        /// <returns></returns>
        protected override double[,] CreateG(int cpi)
        {
            // cpi is already 3*i
            double[,] G = new double[Dimension, 4];
            for (int col = 0; col < 4; col++)
            {
                var P = ControlPoints[cpi + col];
                for (int row = 0; row < Dimension; row++) G[row, col] = P[row];
            }
            return G;
        }
    }
}
