namespace ArcFrame.Core.Math.Geometry.Splines
{
    /// <summary>
    /// A cubic B-spline
    /// </summary>
    public class BSpline : Spline
    {
        /// <summary>
        /// Create a cubic B-Spline
        /// </summary>
        /// <param name="controlPoints"></param>
        /// <param name="frame"></param>
        /// <param name="fastMode"></param>
        /// <param name="cacheSamplesOverride"></param>
        public BSpline(double[][] controlPoints, 
            FrameModel frame = FrameModel.Frenet,
            bool fastMode = false,
            int cacheSamplesOverride = 0) : base(controlPoints, ConstructBasis(), frame, fastMode, cacheSamplesOverride) { }


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
