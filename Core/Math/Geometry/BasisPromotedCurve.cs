using ArcFrame.Core.Math;
using ArcFrame.Core.Results;

namespace ArcFrame.Core.Geometry
{
    /// <summary>
    /// Embed curves into another larger dimensional subspace with arbitrary basis vectors.
    /// </summary>
    public sealed class BasisPromotedCurve : IArcLengthCurve
    {
        private readonly IArcLengthCurve _inner; // dimension d
        private readonly int _N;                  // target ambient dimension
        private readonly int _d;                  // inner curve dimension
        /// <summary>
        /// Nxd, ON columns spanning embedding subspace.
        /// Example: if our inner curve is embedded in 2D and we do canonical embed in 3D, it might look like:
        /// _E = [
        ///         [0, 1],
        ///         [0, 0],
        ///         [1, 0]
        ///      ]
        /// This means our new _e1 basis vector is along the Z axis, and
        /// </summary>
        private readonly double[,] _E;            // N x d, ON columns spanning the embedding subspace
        private readonly double[] _p0;            // optional offset (N)
        public double[] Position(double s) => Evaluate(s).P;
        public double[] Tangent(double s) => Evaluate(s).T;

        /// <summary>
        /// Embed a d dimensional curve in an arbitrary N dimensional subspace. where d < N.
        /// </summary>
        /// <param name="inner">Source curve in R^d.</param>
        /// <param name="E">N x d matrix, columns are (ideally) orthonormal basis vectors of the target subspace.</param>
        /// <param name="p0">Optional translation in R^N (defaults to 0).</param>
        /// <param name="orthonormalizeColumns">If true, run a light Gram–Schmidt to clean E's columns.</param>
        /// <param name="tol">for orthornomalization of E column vectors</param>
        /// <exception cref="ArgumentException"></exception>
        public BasisPromotedCurve(IArcLengthCurve inner, double[,] E, double[]? p0 = null, bool orthonormalizeColumns = true, double tol = 1e-10)
        {
            if (E == null) throw new ArgumentException("Basis vectors matrix is null");
            int n0 = E.GetLength(0);
            int n1 = E.GetLength(1);
            if (n1 != inner.Dimension) throw new ArgumentException("E must have d columns (d = inner.Dimension).");

            _inner = inner;
            _d = inner.Dimension;
            _N = n0;
            _E = (double[,])E.Clone();

            if (orthonormalizeColumns) OrthonormalizeColumnsInPlace(_E, tol);
            ONFrame.ValidateONColumns(_E, tol);

            _p0 = p0 != null ? CheckLen(p0, _N) : new double[_N];
        }

        public int Dimension => _N;

        public double Length => _inner.Length;

        public Sample Evaluate(double s)
        {
            var q = _inner.Evaluate(s); // q.P, q.T in R^d ; q.K length d-1 (or null)

            // P_out = E * P_in + p0, R_out = E * R_in, T_out = E * T_in
            double[] mult = Helpers.Multiply(_E, q.P);
            double[] tempP = Helpers.Add(mult, _p0);
            // double[] tempT = Helpers.Multiply(_E, q.T);
            // _E is Nxd -> q.R is dxd 
            // we do _EqR = _E * q.R to rotate inner inside _E
            // we find the complement of _E, _Ec, where _E concat _Ec is SO(N) 
            // we do _EqR concat _Ec to build R final (combine last two steps)
            double[,] _EqR = Helpers.Multiply(_E, q.R);
            double[,] _Ecomplete = ONFrame.R0_FromR_Complete(_E);
            double[,] _Ef = new double[_N, _N];
            for (int i = 0; i < _d; i++)
            {
                for (int j = 0; j < _N; j++)
                {
                    _Ef[j, i] = _EqR[j, i];
                }
            }
            for (int i = _d; i < _N; i++)
            {
                for (int j = 0; j < _N; j++)
                {
                    _Ef[j, i] = _Ecomplete[j, i];
                }
            }

            // Curvatures: copy inner generalized curvatures and pad zeros.
            // (Embedding is an isometry on the subspace → intrinsic κ's preserved.)
            int outK = System.Math.Max(0, _N - 1);
            var K = new double[outK];
            if (q.k != null)
            {
                int m = System.Math.Min(q.k.Length, outK);
                for (int i = 0; i < m; i++) K[i] = q.k[i];
            }

            return new Sample(tempP, _Ef, q.s, K);
        }

        /// <summary>
        /// Get a number of evenly spaced samples along the arc length
        /// </summary>
        /// <param name="count"></param>
        /// <returns></returns>
        public Sample[] GetSamples(int count)
        {
            Sample[] samples = new Sample[count];
            if (count <= 0) count = 1;
            double s;
            for (int i = 0; i < count; i++)
            {
                s = (i == count - 1) ? Length : i * (Length / System.Math.Max(1, count - 1));
                samples[i] = Evaluate(s);
            }
            return samples;
        }

        private static double[] CheckLen(double[] v, int n)
        {
            if (v.Length != n) throw new ArgumentException("Vector length mismatch.");
            return (double[])v.Clone();
        }

        private static void OrthonormalizeColumnsInPlace(double[,] E, double tol)
        {
            int n = E.GetLength(0), d = E.GetLength(1);
            for (int j = 0; j < d; j++)
            {
                // subtract projections onto previous columns
                for (int k = 0; k < j; k++)
                {
                    double proj = 0;
                    for (int i = 0; i < n; i++) proj += E[i, j] * E[i, k];
                    for (int i = 0; i < n; i++) E[i, j] -= proj * E[i, k];
                }
                // normalize
                double norm2 = 0; for (int i = 0; i < n; i++) norm2 += E[i, j] * E[i, j];
                double norm = System.Math.Sqrt(norm2);
                if (norm < tol) throw new ArgumentException("E has dependent/degenerate columns.");
                for (int i = 0; i < n; i++) E[i, j] /= norm;
            }
        }

        /// <summary>
        /// Build E from column vectors (each length N). Columns will be orthonormalized unless disabled.
        /// </summary>
        public static BasisPromotedCurve FromColumns(IArcLengthCurve inner, double[][] columns, double[]? p0 = null, bool orthonormalizeColumns = true, double tol = 1e-10)
        {
            if (columns == null || columns.Length != inner.Dimension) throw new ArgumentException("Need exactly d column vectors (d = inner.Dimension).");
            int N = columns[0].Length;
            var E = new double[N, columns.Length];
            for (int j = 0; j < columns.Length; j++)
            {
                if (columns[j].Length != N) throw new ArgumentException("Column length mismatch.");
                for (int i = 0; i < N; i++) E[i, j] = columns[j][i];
            }
            return new BasisPromotedCurve(inner, E, p0, orthonormalizeColumns, tol);
        }
    }
}
