using ArcFrame.Core.Math;

namespace ArcFrame.Core.Params
{
    /// <summary>
    /// Data class to hold intrisic curve data.
    /// </summary>
    public class CurveSpec
    {
        /// <summary>
        /// Dimension of the curve.
        /// </summary>
        public int N { get; }
        /// <summary>
        /// Total length of the curve.
        /// </summary>
        public double Length { get; }
        /// <summary>
        /// Initial position of the curve.
        /// </summary>
        public double[] P0 { get; }
        /// <summary>
        /// Initial orthonormal basis frame.
        /// </summary>
        public double[,] R0 { get; }
        /// <summary>
        /// The curvature model (how curvature will evolve with arc length)
        /// </summary>
        public ICurvatureLaw Kappa { get; }
        /// <summary>
        /// The type of frame (Frenet: curvature guided ONFrame | Bishop: tangent guided ONFrame)
        /// </summary>
        public FrameModel Frame { get; }

        public CurveSpec(int n, double length, double[] p0, double[,] r0, ICurvatureLaw kappa, FrameModel frame)
        {
            N = n;
            Length = length;
            P0 = (double[])p0.Clone();
            R0 = (double[,])r0.Clone();
            Kappa = kappa;
            Frame = frame;
        }

        /// <summary>
        /// Get a column from the R0 matrix.
        /// </summary>
        /// <param name="columnNumber"></param>
        /// <returns></returns>
        public double[] GetONAxis(int columnNumber)
        {
            int n = R0.GetLength(0);
            double[] v = new double[n];
            for(int i = 0; i < n; i++)
            {
                v[i] = R0[i, columnNumber];
            }
            return v;
        }
    }
}
