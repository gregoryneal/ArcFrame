namespace ArcFrame.Core.Math
{
    /// <summary>
    /// Linear algebra helper classes
    /// </summary>
    public static class LA
    {
        /// <summary>
        /// LU factorization with partial pivoting: P*A = L*U.
        /// Returns (LU, piv, pivSign). L has unit diagonal and is stored in the strict lower triangle of LU.
        /// U is stored in the upper triangle (incl. diagonal). piv is the row permutation.
        /// Throws InvalidOperationException if matrix is (numerically) singular.
        /// </summary>
        public static (double[,] LU, int[] piv, int pivSign) LUFactor(double[,] A, double tol = 1e-14)
        {
            if (A is null) throw new ArgumentNullException(nameof(A));
            int n = A.GetLength(0);
            if (n != A.GetLength(1)) throw new ArgumentException("Matrix must be square.");

            // Make a working copy
            var LU = (double[,])A.Clone();
            var piv = new int[n];
            for (int i = 0; i < n; i++) piv[i] = i;
            int pivSign = 1;

            for (int k = 0; k < n; k++)
            {
                // Find pivot row
                int p = k;
                double max = System.Math.Abs(LU[k, k]);
                for (int i = k + 1; i < n; i++)
                {
                    double v = System.Math.Abs(LU[i, k]);
                    if (v > max) { max = v; p = i; }
                }
                if (max <= tol) throw new InvalidOperationException("Matrix is (numerically) singular.");

                // Swap rows if needed
                if (p != k)
                {
                    SwapRows(LU, p, k);
                    (piv[p], piv[k]) = (piv[k], piv[p]);
                    pivSign = -pivSign;
                }

                // Factorize column k
                double pivot = LU[k, k];
                // Compute multipliers below the pivot
                for (int i = k + 1; i < n; i++)
                    LU[i, k] /= pivot;

                // Schur complement update
                for (int i = k + 1; i < n; i++)
                {
                    double lik = LU[i, k];
                    if (lik == 0.0) continue;
                    for (int j = k + 1; j < n; j++)
                        LU[i, j] -= lik * LU[k, j];
                }
            }

            return (LU, piv, pivSign);
        }

        /// <summary>
        /// Solve A x = b given LU and piv from LUFactor. Returns the solution x.
        /// </summary>
        public static double[] LUSolve(double[,] LU, int[] piv, double[] b)
        {
            int n = LU.GetLength(0);
            if (LU.GetLength(1) != n) throw new ArgumentException("LU not square.");
            if (b.Length != n) throw new ArgumentException("b length mismatch.");
            if (piv.Length != n) throw new ArgumentException("piv length mismatch.");

            var x = new double[n];

            // Apply permutation: x = P*b
            for (int i = 0; i < n; i++) x[i] = b[piv[i]];

            // Forward solve L y = P*b  (L has unit diagonal)
            for (int i = 0; i < n; i++)
            {
                double sum = x[i];
                for (int j = 0; j < i; j++) sum -= LU[i, j] * x[j];
                x[i] = sum; // since L_ii = 1
            }

            // Back solve U x = y
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = x[i];
                for (int j = i + 1; j < n; j++) sum -= LU[i, j] * x[j];
                x[i] = sum / LU[i, i];
            }

            return x;
        }

        /// <summary>
        /// Solve A X = B (many RHS). B is (n x m) and overwritten with X.
        /// </summary>
        public static void LUSolveInPlace(double[,] LU, int[] piv, double[,] B)
        {
            int n = LU.GetLength(0);
            int m = B.GetLength(1);
            if (LU.GetLength(1) != n) throw new ArgumentException("LU not square.");
            if (B.GetLength(0) != n) throw new ArgumentException("B row mismatch.");
            if (piv.Length != n) throw new ArgumentException("piv length mismatch.");

            // Apply permutation to B: B <- P*B
            var tmp = new double[m];
            for (int i = 0; i < n; i++)
            {
                int pi = piv[i];
                if (pi == i) continue;
                // swap row i with row pi in B
                for (int j = 0; j < m; j++) tmp[j] = B[i, j];
                for (int j = 0; j < m; j++) B[i, j] = B[pi, j];
                for (int j = 0; j < m; j++) B[pi, j] = tmp[j];
            }

            // Forward solve: L Y = P*B
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    double sum = B[i, j];
                    for (int k = 0; k < i; k++) sum -= LU[i, k] * B[k, j];
                    B[i, j] = sum; // L_ii = 1
                }
            }

            // Back solve: U X = Y
            for (int i = n - 1; i >= 0; i--)
            {
                double uii = LU[i, i];
                for (int j = 0; j < m; j++)
                {
                    double sum = B[i, j];
                    for (int k = i + 1; k < n; k++) sum -= LU[i, k] * B[k, j];
                    B[i, j] = sum / uii;
                }
            }
        }

        /// <summary> Compute det(A) from LU (det = pivSign * product(diag(U))). </summary>
        public static double DeterminantFromLU(double[,] LU, int pivSign)
        {
            int n = LU.GetLength(0);
            double det = pivSign;
            for (int i = 0; i < n; i++) det *= LU[i, i];
            return det;
        }

        /// <summary> Build A^{-1} from LU by solving A X = I. </summary>
        public static double[,] InvertFromLU(double[,] LU, int[] piv)
        {
            int n = LU.GetLength(0);
            var I = new double[n, n];
            for (int i = 0; i < n; i++) I[i, i] = 1.0;

            // Solve for all columns (in place): result overwrites I
            LUSolveInPlace(LU, piv, I);
            return I;
        }

        // --- helpers ---
        private static void SwapRows(double[,] M, int r1, int r2)
        {
            if (r1 == r2) return;
            int cols = M.GetLength(1);
            for (int j = 0; j < cols; j++)
            {
                double tmp = M[r1, j]; M[r1, j] = M[r2, j]; M[r2, j] = tmp;
            }
        }
    }

}
