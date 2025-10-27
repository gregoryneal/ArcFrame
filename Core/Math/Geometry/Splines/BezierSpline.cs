namespace ArcFrame.Core.Math.Geometry.Splines
{
    /// <summary>
    /// Create a Cubic Bezier curve out of control points.
    /// </summary>
    public class BezierSpline : Spline
    {
        public BezierSpline(double[][] controlPoints, FrameModel frame = FrameModel.Frenet) : base(controlPoints, ConstructBasis(), frame)
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

        // stride 3 segmentation vs stride 1 in original implementation
        protected override (int k, double u) Locate(double t)
        {
            t = System.Math.Clamp(t, 0, 1);
            int segCount = (ControlPoints.Length - 1) / 3;
            double seg = t * segCount;
            int i = System.Math.Min(segCount - 1, (int)System.Math.Floor(seg));
            return (3 * i, seg - i);
        }

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
