namespace ArcFrame.Core.Results
{
    /// <summary>
    /// A sample from a generalized ND arc length parameterized curve
    /// </summary>
    /// <param name="P">N dimensional position at s</param>
    /// <param name="R">NxN matrix of orthonormal basis vectors at s. Columns are [T, N, B, ...]</param>
    /// <param name="s">Arc length</param>
    /// <param name="k">N-1 dimensional generalized curvatures at s.
    /// Example: let N = 2, then k = [k1] with k1 the normal curvature. 
    /// Let N = 3, then k = [k1, k2] with k2 the normal torsion.</param>
    public struct Sample(double[] P, double[,] R, double s, double[] k)
    {
        /// <summary>
        /// Position
        /// </summary>
        public double[] P = P;
        /// <summary>
        /// SO(N)
        /// </summary>
        public double[,] R = R;
        /// <summary>
        /// Curvatures
        /// </summary>
        public double[] k = k;
        /// <summary>
        /// Arc length
        /// </summary>
        public double s = s;

        private double[]? _T;
        /// <summary>
        /// Tangent vector
        /// </summary>
        public double[] T
        {
            get
            {
                if (_T != null) return _T;
                int n = R.GetLength(0);
                double[] T = new double[n];
                for (int i = 0; i < n; i++)
                {
                    T[i] = R[i, 0];
                }
                _T = T;
                return _T;
            }
        }
    }
}
