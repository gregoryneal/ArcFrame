namespace ArcFrame.Core.Math.Geometry.Splines
{
    /// <summary>
    /// A cubic B-spline
    /// </summary>
    public class BSpline : Spline
    {
        public BSpline(double[][] controlPoints, FrameModel frame = FrameModel.Frenet) : base(controlPoints, ConstructBasis(), frame) { }


        private static double[,] ConstructBasis()
        {
            return Helpers.Multiply(new double[,]
            {
                {1, -3, 3, -1},
                {4, 0, -6, 3 },
                {1, 3, 3, -3 },
                {0, 0, 0,  1 }
            }, 1.0/6.0);
        }
    }
}
