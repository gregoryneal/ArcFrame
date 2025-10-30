using ArcFrame.Core.Math;

namespace ArcFrame.Core.Results
{
    /// <summary>
    /// A sample from a generalized ND arc length parameterized curve
    /// </summary>
    public struct Sample
    {
        /// <summary>
        /// Position
        /// </summary>
        public double[] P;
        /// <summary>
        /// SO(N)
        /// </summary>
        public double[,] R;
        /// <summary>
        /// Curvatures
        /// </summary>
        public double[] k;
        /// <summary>
        /// Arc length
        /// </summary>
        public double s;

        /// <param name="P">N dimensional position at s</param>
        /// <param name="R">NxN matrix of orthonormal basis vectors at s. Columns are [T, N, B, ...]</param>
        /// <param name="s">Arc length</param>
        /// <param name="k">N-1 dimensional generalized curvatures at s.
        /// Example: let N = 2, then k = [k1] with k1 the normal curvature. 
        /// Let N = 3, then k = [k1, k2] with k2 the normal torsion.</param>
        public Sample(double[] P, double[,] R, double s, double[] k)
        {
            this.P = P;
            this.R = R;
            this.k = k;
            this.s = s;
        }

        /// <summary>
        /// Tangent vector
        /// </summary>
        public double[] T
        {
            get
            {
                return ONFrame.GetCol(R, 0);
            }
        }
    }
}
