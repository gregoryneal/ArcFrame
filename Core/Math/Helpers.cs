using System;

namespace ArcFrame.Core.Math
{
    /// <summary>
    /// Other math helpers.
    /// </summary>
    public static class Helpers
    {
        /// <summary>
        /// Normalize an N dimensional vector
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static double[] Normalize(double[] v)
        {
            double s = 0;
            foreach (var x in v) s += x * x;
            s = System.Math.Sqrt(s);
            var r = (double[])v.Clone();
            if (s > 0) for (int i = 0; i < r.Length; i++) r[i] /= s;
            return r;
        }

        /// <summary>
        /// 3D determinant on square matrix
        /// </summary>
        /// <param name="M"></param>
        /// <returns></returns>
        public static double Det3(double[,] M)
        {
            return
                M[0, 0] * (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]) -
                M[0, 1] * (M[1, 0] * M[2, 2] - M[1, 2] * M[2, 0]) +
                M[0, 2] * (M[1, 0] * M[2, 1] - M[1, 1] * M[2, 0]);
        }

        /// <summary>
        /// 3D cross product
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static double[] Cross3(double[] a, double[] b) => new[] { (a[1] * b[2]) - (a[2] * b[1]), (a[2] * b[0]) - (a[0] * b[2]), (a[0] * b[1]) - (a[1] * b[0]) };

        /// <summary>
        /// N-D Dot product
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        /// <exception cref="RankException"></exception>
        public static double Dot(double[] a, double[] b)
        {
            if (a.Length != b.Length) throw new RankException("Array size mismatch");
            double s = 0;
            for (int i = 0; i < a.Length; i++)
            {
                s += a[i] * b[i];
            }
            return s;
        }

        /// <summary>
        /// Lerp one vector to another
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="u"></param>
        /// <returns></returns>
        public static double[] Lerp(double[] a, double[] b, double u)
        {
            var r = new double[a.Length];
            for (int i = 0; i < r.Length; i++)
            {
                r[i] = a[i] + u * (b[i] - a[i]);
            }
            return r;
        }

        /// <summary>
        /// Squared length of a vector.
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static double Len2(double[] a)
        {
            return Dot(a, a);
        }

        /// <summary>
        /// Equal to the sum of the largest column vector of M
        /// </summary>
        /// <param name="M"></param>
        /// <returns></returns>
        public static double OneNorm(double[,] M)
        {
            int n = M.GetLength(0);
            int m = M.GetLength(1);
            double mx = 0;
            for (int j = 0; j < m; j++)
            {
                double col = 0;
                for (int i = 0; i < n; i++)
                    col += System.Math.Abs(M[i, j]);
                mx = System.Math.Max(mx, col);
            }
            return mx;
        }
        /// <summary>
        /// Add two N dimensional vectors
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        /// <exception cref="RankException"></exception>
        public static double[] Add(double[] a, double[] b)
        {
            if (a.Length != b.Length) throw new RankException($"Array size mismatch ({a.Length} != {b.Length})");
            double[] s = new double[a.Length];
            for (int i = 0; i < a.Length; i++)
            {
                s[i] = a[i] + b[i];
            }
            return s;
        }

        /// <summary>
        /// Add two NxM dimensional matrices
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        /// <exception cref="RankException"></exception>
        public static double[,] Add(double[,] a, double[,] b)
        {
            int n = a.GetLength(0);
            int m = a.GetLength(1);
            int o = b.GetLength(0);
            int p = b.GetLength(1);
            if (n != o || m != p) throw new ArgumentException("Dimension mismatch!");
            double[,] s = new double[n, m];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    s[i, j] = a[i, j] + b[i, j];
                }
            }
            return s;
        }

        /// <summary>
        /// Subtract two N dimensional vectors, a - b.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        /// <exception cref="RankException"></exception>
        public static double[] Subtract(double[] a, double[] b)
        {
            if (a.Length != b.Length) throw new RankException("Array size mismatch");
            double[] s = new double[a.Length];
            for (int i = 0; i < a.Length; i++)
            {
                s[i] = a[i] - b[i];
            }
            return s;
        }

        /// <summary>
        /// Add two NxM dimensional matrices
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        /// <exception cref="RankException"></exception>
        public static double[,] Subtract(double[,] a, double[,] b)
        {
            int n = a.GetLength(0);
            int m = a.GetLength(1);
            int o = b.GetLength(0);
            int p = b.GetLength(1);
            if (n != o || m != p) throw new ArgumentException("Dimension mismatch!");
            double[,] s = new double[n, m];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    s[i, j] = a[i, j] - b[i, j];
                }
            }
            return s;
        }

        /// <summary>
        /// Multiply an N dimensional vector by a scalar
        /// </summary>
        /// <param name="a">Scalar</param>
        /// <param name="b">Vector</param>
        /// <returns></returns>
        /// <exception cref="RankException"></exception>
        public static double[] Multiply(double a, double[] b)
        {
            double[] s = new double[b.Length];
            for (int i = 0; i < b.Length; i++)
            {
                s[i] = a * b[i];
            }
            return s;
        }

        /// <summary>
        /// Multiply an N dimensional vector by a scalar
        /// </summary>
        /// <param name="a">Scalar</param>
        /// <param name="b">Vector</param>
        /// <returns></returns>
        /// <exception cref="RankException"></exception>
        public static double[] Multiply(double[] a, double b)
        {
            return Multiply(b, a);
        }

        /// <summary>
        /// Multiply and MxN matrix a by an N row vector b.
        /// </summary>
        /// <param name="a">MxN matrix</param>
        /// <param name="b">N row vector</param>
        /// <returns>An N row vector</returns>
        public static double[] Multiply(double[,] a, double[] b)
        {
            int n = b.Length;
            int n0 = a.GetLength(0);
            int n1 = a.GetLength(1);
            if (n1 != n) throw new ArgumentException($"Dimension mismatch. a is [{n0},{n1}] and b is length {n}. You can only multiply an MxN matrix with a length N vector");

            double[] c = new double[n0];
            double dot;
            //rows of a, rows of c
            for (int i = 0; i < n0; i++)
            {
                dot = 0;
                //cols of a, index (row) of b
                for (int j = 0; j < n; j++)
                {
                    dot += a[i, j] * b[j];
                }
                c[i] = dot;
            }
            return c;
        }

        /// <summary>
        /// Multiply two matrices.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static double[,] Multiply(double[,] a, double[,] b)
        {
            //a is n x k, b is l x m. In order to multiply we need k = l
            //we will end up with and n x m matrix.
            int n = a.GetLength(0);
            int l = b.GetLength(0);
            int m = b.GetLength(1);
            int k = a.GetLength(1);
            if (k != l) throw new ArgumentException("Dimension mismatch");
            var Z = new double[n, m];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    double s = 0;
                    for (int t = 0; t < k; t++)
                    {
                        s += a[i, t] * b[t, j];
                    }
                    Z[i, j] = s;
                }
            }
            return Z;
        }
        /// <summary>
        /// Multiply a matrix and a scalar.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static double[,] Multiply(double a, double[,] b)
        {
            int n = b.GetLength(0);
            int m = b.GetLength(1);
            var Z = new double[n, m];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    Z[i, j] = b[i, j] * a;
                }
            }
            return Z;
        }

        /// <summary>
        /// Multiply a matrix and a scalar.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static double[,] Multiply(double[,] a, double b)
        {
            return Multiply(b, a);
        }
        /// <summary>
        /// Length of a vector. Same as the norm.
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static double Len(double[] a)
        {
            return System.Math.Sqrt(Len2(a));
        }

        //Gram-Schmidt style orthonormalization for building an N dimensional frame from
        //lower dimensional starting vectors. 

        /// <summary>
        /// 
        ///       v   ↗^
        ///          / |
        ///         /  |
        ///        /   |
        ///       /    |
        ///      /     | Reject_u(v)
        ///     /      |
        ///    /       |
        ///   /        |
        ///  /         |
        ///  --------->----------------> u
        ///  Project_u(v)
        /// 
        /// Reject_u(v) = v - Dot(u, v)u
        /// </summary>
        /// <param name="v">Normalized ND vector</param>
        /// <param name="u">Normalized ND vector</param>
        /// <returns></returns>
        public static double[] Reject(double[] v, double[] u)
        {
            return Subtract(v, Multiply(Dot(v, u), u));
        }

        /// <summary>
        /// Create the kth basis vector in N dimensional space
        /// </summary>
        /// <param name="N">Dimensions</param>
        /// <param name="k">Basis index</param>
        /// <returns></returns>
        public static double[] StdBasis(int N, int k)
        {
            double[] v = new double[N];
            v[System.Math.Min(k, N - 1)] = 1;
            return v;
        }

        /// <summary>
        /// Create a basis vector least aligned with t
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        public static double[] StdSeed(double[] t)
        {
            double min = double.MaxValue;
            int best = -1;

            double a;
            for (int i = 0; i < t.Length; i++)
            {
                a = System.Math.Abs(t[i]);
                if (a < min)
                {
                    min = a;
                    best = i;
                }
            }

            double[] v = new double[t.Length];
            v[best] = 1;
            return v;
        }

        //===== End Gram-Schmidt style orthonormalization functions

        /// <summary>
        /// Performs the matrix calculation:
        /// M += alpha * (u * v^T)
        /// Note: M is an NxN square matrix.
        /// Note: u * v^T is the outer product of two N dimensional vectors.
        /// </summary>
        /// <param name="M">Matrix to operate on</param>
        /// <param name="u">vector 1</param>
        /// <param name="v">vector 2</param>
        /// <param name="alpha">scalar</param>
        public static void AddOuterScaled(double[,] M, double[] u, double[] v, double alpha)
        {
            int n = u.Length;
            for (int i = 0; i < n; i++)
            {
                double ui = u[i] * alpha;
                for (int j = 0; j < n; j++) M[i, j] += ui * v[j];
            }
        }

        /// <summary>
        /// Invert an arbitrary square matrix A using Gauss–Jordan elimination with
        /// scaled partial pivoting on the augmented system [A|I] -> [I|A^{-1}].
        /// Throws InvalidOperationException if A is (numerically) singular.
        /// </summary>
        public static double[,] Invert(double[,] A, double tol = 1e-14)
        {
            if (A is null) throw new ArgumentNullException(nameof(A));
            int n = A.GetLength(0);
            if (n != A.GetLength(1)) throw new ArgumentException("Matrix must be square.", nameof(A));

            // Build augmented matrix [A | I]
            int m = 2 * n;
            var aug = new double[n, m];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++) aug[i, j] = A[i, j];
                aug[i, n + i] = 1.0; // identity on the right
            }

            // Scaled partial pivoting: s[i] = max_j |A[i,j]| to scale pivot tests
            var s = new double[n];
            for (int i = 0; i < n; i++)
            {
                double max = 0.0;
                for (int j = 0; j < n; j++) max = System.Math.Max(max, System.Math.Abs(aug[i, j]));
                if (max <= tol) throw new InvalidOperationException("Matrix is (numerically) singular: zero scaling row.");
                s[i] = max;
            }

            // Forward pass: make left half upper-triangular with 1.0 pivots
            for (int col = 0; col < n; col++)
            {
                // Pick pivot row r >= col maximizing |A[r,col]|/s[r]
                int r = col;
                double best = System.Math.Abs(aug[col, col]) / s[col];
                for (int i = col + 1; i < n; i++)
                {
                    double val = System.Math.Abs(aug[i, col]) / s[i];
                    if (val > best) { best = val; r = i; }
                }

                // Singular check
                if (System.Math.Abs(aug[r, col]) <= tol) throw new InvalidOperationException("Matrix is (numerically) singular: zero pivot.");

                // Swap rows if needed (both halves)
                if (r != col) SwapRows(aug, r, col);

                // Scale pivot row to make pivot = 1
                double piv = aug[col, col];
                double invPiv = 1.0 / piv;
                for (int j = col; j < m; j++) aug[col, j] *= invPiv; // left entries before 'col' are zeros by construction

                // Eliminate this column from all OTHER rows
                for (int i = 0; i < n; i++)
                {
                    if (i == col) continue;
                    double factor = aug[i, col];
                    if (System.Math.Abs(factor) <= 0.0) continue;
                    for (int j = col; j < m; j++) aug[i, j] -= factor * aug[col, j];
                    // Explicitly zero the eliminated entry to control roundoff
                    aug[i, col] = 0.0;
                }
            }

            // At this point, left half should be identity; extract inverse from right half
            var inv = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    inv[i, j] = aug[i, n + j];

            return inv;
        }

        /// <summary>
        /// Try-invert variant that returns false instead of throwing on singularity.
        /// </summary>
        public static bool TryInvert(double[,] A, out double[,]? Ainv, double tol = 1e-14)
        {
            try
            {
                Ainv = Invert(A, tol);
                return true;
            }
            catch
            {
                Ainv = null;
                return false;
            }
        }
        /// <summary>
        /// Take the transpose of a matrix.
        /// </summary>
        /// <param name="matrix"></param>
        /// <returns></returns>
        public static double[,] Transpose(double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            double[,] result = new double[cols, rows];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    result[j, i] = matrix[i, j];
                }
            }
            return result;
        }

        private static void SwapRows(double[,] M, int r1, int r2)
        {
            if (r1 == r2) return;
            int cols = M.GetLength(1);
            for (int j = 0; j < cols; j++)
            {
                double tmp = M[r1, j];
                M[r1, j] = M[r2, j];
                M[r2, j] = tmp;
            }
        }// A function to calculate binomial coefficients "n choose k"
        private static long BinomialCoefficient(int n, int k)
        {
            if (k < 0 || k > n)
            {
                return 0;
            }
            if (k == 0 || k == n)
            {
                return 1;
            }
            if (k > n / 2)
            {
                k = n - k;
            }

            long result = 1;
            for (int i = 1; i <= k; i++)
            {
                result = result * (n - i + 1) / i;
            }
            return result;
        }

        /// <summary>
        /// Computes the nth derivative of a function using the central finite difference method.
        /// </summary>
        /// <param name="f">The function to differentiate.</param>
        /// <param name="x">The point at which to evaluate the derivative.</param>
        /// <param name="n">The order of the derivative.</param>
        /// <param name="len">The expected length of the output array</param>
        /// <param name="h">The step size.</param>
        /// <returns>The approximate value of the nth derivative.</returns>
        public static double[] Differentiate(Func<double, double[]> f, double x, int n, int len, double h = 1E-6)
        {
            if (n < 0)
            {
                throw new ArgumentException("The derivative order 'n' must be non-negative.");
            }
            if (h <= 0)
            {
                throw new ArgumentException("The step size 'h' must be positive.");
            }

            double[] sum = new double[len];
            for (int k = 0; k <= n; k++)
            {
                double coefficient = BinomialCoefficient(n, k);
                double sign = (n - k) % 2 == 0 ? 1.0 : -1.0;
                double[] functionValue = f(x + (k - (double)n / 2.0) * h);
                sum = Helpers.Add(sum, Helpers.Multiply(sign * coefficient, functionValue));
            }

            return Helpers.Multiply(sum, 1 / System.Math.Pow(h, n));
        }

        /// <summary>
        /// Helper to print a matrix to the console.
        /// </summary>
        /// <param name="R"></param>
        public static void PrintMat(double[,] R)
        {
            int N = R.GetLength(0);
            int M = R.GetLength(1);
            Console.WriteLine("[");
            for (int i = 0; i < N; i++)
            {
                Console.Write("[");
                for (int j = 0; j < M; j++)
                {
                    Console.Write(R[i, j]);
                    if (j < M - 1) Console.Write(", ");
                }
                if (i < N - 1) Console.WriteLine("],");
                else Console.WriteLine("]");
            }
            Console.WriteLine("]");
        }

        /// <summary>
        /// Helper to print a vector to the console.
        /// </summary>
        /// <param name="T"></param>
        public static void PrintVector(double[] T)
        {
            int N = T.Length;
            Console.Write("[");
            for (int i = 0; i < N; i++)
            {
                Console.Write(T[i]);
                if (i < N - 1) Console.Write(", ");
            }
            Console.WriteLine("]");
        }

    }
}
