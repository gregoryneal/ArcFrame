using System;

namespace ArcFrame.Core.Math
{
    internal class SVDJacobiProgram
    {
        public static bool SVD_Jacobi(double[,] M, out double[,] U, out double[,] Vh, out double[] s)
        {
            int nRows = M.GetLength(0), nCols = M.GetLength(1);
            double[][] M1 = new double[nRows][];
            double[] valS = new double[nRows];

            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    M1[i][j] = M[i, j];
                }
            }

            bool val = SVD_Jacobi(M1, out double[][] U1, out double[][] Vh1, out s);

            U = new double[nRows, nCols];
            Vh = new double[nRows, nCols];
            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    U[i, j] = U1[i][j];
                    Vh[i, j] = Vh1[i][j];
                }
            }

            return val;
        }

        /// <summary>
        /// Calculate single value decomposition of M
        /// </summary>
        /// <param name="M"></param>
        /// <param name="U"></param>
        /// <param name="Vh"></param>
        /// <param name="s"></param>
        /// <returns></returns>
        public static bool SVD_Jacobi(double[][] M, out double[][] U, out double[][] Vh, out double[] s)
        {
            double DBL_EPSILON = 1.0e-15;

            double[][] A = MatCopy(M); // working U
            int m = A.Length;
            int n = A[0].Length;
            double[][] Q = MatIdentity(n); // working V
            double[] t = new double[n];    // working s

            // initialize counters
            int count = 1;
            int sweep = 0;
            //int sweepmax = 5 * n;

            double tolerance = 10 * m * DBL_EPSILON; // heuristic

            // Always do at least 12 sweeps
            int sweepmax = System.Math.Max(5 * n, 12); // heuristic

            // store the column error estimates in St for use
            // during orthogonalization

            for (int j = 0; j < n; ++j)
            {
                double[] cj = MatGetColumn(A, j);
                double sj = VecNorm(cj);
                t[j] = DBL_EPSILON * sj;
            }

            // orthogonalize A by plane rotations
            while (count > 0 && sweep <= sweepmax)
            {
                // initialize rotation counter
                count = n * (n - 1) / 2;

                for (int j = 0; j < n - 1; ++j)
                {
                    for (int k = j + 1; k < n; ++k)
                    {
                        double cosine, sine;

                        double[] cj = MatGetColumn(A, j);
                        double[] ck = MatGetColumn(A, k);

                        double p = 2.0 * VecDot(cj, ck);
                        double a = VecNorm(cj);
                        double b = VecNorm(ck);

                        double q = a * a - b * b;
                        double v = Hypot(p, q);

                        // test for columns j,k orthogonal,
                        // or dominant errors 
                        double abserr_a = t[j];
                        double abserr_b = t[k];

                        bool sorted = (a >= b);
                        bool orthog = (System.Math.Abs(p) <=
                            tolerance * (a * b));
                        bool noisya = (a < abserr_a);
                        bool noisyb = (b < abserr_b);

                        if (sorted && (orthog ||
                            noisya || noisyb))
                        {
                            --count;
                            continue;
                        }

                        // calculate rotation angles
                        if (v == 0 || !sorted)
                        {
                            cosine = 0.0;
                            sine = 1.0;
                        }
                        else
                        {
                            cosine = System.Math.Sqrt((v + q) / (2.0 * v));
                            sine = p / (2.0 * v * cosine);
                        }

                        // apply rotation to A (U)
                        for (int i = 0; i < m; ++i)
                        {
                            double Aik = A[i][k];
                            double Aij = A[i][j];
                            A[i][j] = Aij * cosine + Aik * sine;
                            A[i][k] = -Aij * sine + Aik * cosine;
                        }

                        // update singular values
                        t[j] = System.Math.Abs(cosine) * abserr_a +
                            System.Math.Abs(sine) * abserr_b;
                        t[k] = System.Math.Abs(sine) * abserr_a +
                            System.Math.Abs(cosine) * abserr_b;

                        // apply rotation to Q (V)
                        for (int i = 0; i < n; ++i)
                        {
                            double Qij = Q[i][j];
                            double Qik = Q[i][k];
                            Q[i][j] = Qij * cosine + Qik * sine;
                            Q[i][k] = -Qij * sine + Qik * cosine;
                        } // i
                    } // k
                } // j

                ++sweep;
            } // while

            //    compute singular values
            double prevNorm = -1.0;

            for (int j = 0; j < n; ++j)
            {
                double[] column = MatGetColumn(A, j);
                double norm = VecNorm(column);

                // determine if singular value is zero
                if (norm == 0.0 || prevNorm == 0.0
                    || (j > 0 &&
                        norm <= tolerance * prevNorm))
                {
                    t[j] = 0.0;
                    for (int i = 0; i < m; ++i)
                        A[i][j] = 0.0;
                    prevNorm = 0.0;
                }
                else
                {
                    t[j] = norm;
                    for (int i = 0; i < m; ++i)
                        A[i][j] = A[i][j] * 1.0 / norm;
                    prevNorm = norm;
                }
            }

            U = A;
            Vh = MatTranspose(Q);
            s = t;

            // to sync with default np.linalg.svd() shapes:
            // if m "lt" n, extract 1st m columns of U
            //     extract 1st m values of s
            //     extract 1st m rows of Vh

            if (m < n)
            {
                U = MatExtractFirstColumns(U, m);
                s = VecExtractFirst(s, m);
                Vh = MatExtractFirstRows(Vh, m);
            }

            if (count > 0)
            {
                //Console.WriteLine("Jacobi iterations did not converge");
                return false;
            }

            return true;

        } // SVD_Jacobi()

        // === helper functions =================================
        //
        // MatMake, MatCopy, MatIdentity, MatGetColumn,
        // MatExtractFirstColumns, MatExtractFirstRows,
        // MatTranspose, MatDiag, MatProduct, VecNorm, VecDot,
        // Hypot, VecExtractFirst, MatShow, VecShow
        //
        // ======================================================

        public static double[][] MatMake(int r, int c)
        {
            double[][] result = new double[r][];
            for (int i = 0; i < r; ++i)
                result[i] = new double[c];
            return result;
        }

        static double[][] MatCopy(double[][] m)
        {
            int r = m.Length; int c = m[0].Length;
            double[][] result = MatMake(r, c);
            for (int i = 0; i < r; ++i)
                for (int j = 0; j < c; ++j)
                    result[i][j] = m[i][j];
            return result;
        }

        static double[][] MatIdentity(int n)
        {
            double[][] result = MatMake(n, n);
            for (int i = 0; i < n; ++i)
                result[i][i] = 1.0;
            return result;
        }

        static double[] MatGetColumn(double[][] m, int j)
        {
            int rows = m.Length;
            double[] result = new double[rows];
            for (int i = 0; i < rows; ++i)
                result[i] = m[i][j];
            return result;
        }

        static double[][] MatExtractFirstColumns(double[][] src,
            int n)
        {
            int nRows = src.Length;
            // int nCols = src[0].Length;
            double[][] result = MatMake(nRows, n);
            for (int i = 0; i < nRows; ++i)
                for (int j = 0; j < n; ++j)
                    result[i][j] = src[i][j];
            return result;
        }

        static double[][] MatExtractFirstRows(double[][] src,
            int n)
        {
            // int nRows = src.Length;
            int nCols = src[0].Length;
            double[][] result = MatMake(n, nCols);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < nCols; ++j)
                    result[i][j] = src[i][j];
            return result;
        }

        public static double[][] MatTranspose(double[][] m)
        {
            int r = m.Length;
            int c = m[0].Length;
            double[][] result = MatMake(c, r);
            for (int i = 0; i < r; ++i)
                for (int j = 0; j < c; ++j)
                    result[j][i] = m[i][j];
            return result;
        }

        static double[][] MatDiag(double[] vec)
        {
            int n = vec.Length;
            double[][] result = MatMake(n, n);
            for (int i = 0; i < n; ++i)
                result[i][i] = vec[i];
            return result;
        }

        public static double[][] MatProduct(double[][] matA, double[][] matB)
        {
            int aRows = matA.Length;
            int aCols = matA[0].Length;
            int bRows = matB.Length;
            int bCols = matB[0].Length;
            if (aCols != bRows)
            {
                //MatShow(matA);
                //MatShow(matB);
                throw new Exception("Non-conformable matrices");
            }

            double[][] result = MatMake(aRows, bCols);

            for (int i = 0; i < aRows; ++i)
                for (int j = 0; j < bCols; ++j)
                    for (int k = 0; k < aCols; ++k)
                        result[i][j] += matA[i][k] * matB[k][j];

            return result;
        }

        //Can't believe this wasn't already here.
        public static double[][] MatAdd(double[][] matA, double[][] matB)
        {
            int aRows = matA.Length;
            int aCols = matA[0].Length;
            int bRows = matB.Length;
            int bCols = matB[0].Length;
            if (aRows != bRows || aCols != bCols)
            {
                throw new Exception("Non-conformable matrices");
            }
            double[][] result = MatMake(aRows, aCols);
            for (int i = 0; i < aRows; ++i)
                for (int j = 0; j < aCols; ++j)
                    result[i][j] = matA[i][j] + matB[i][j];
            return result;
        }

        static double VecNorm(double[] vec)
        {
            double sum = 0.0;
            int n = vec.Length;
            for (int i = 0; i < n; ++i)
                sum += vec[i] * vec[i];
            return System.Math.Sqrt(sum);
        }

        static double VecDot(double[] v1, double[] v2)
        {
            int n = v1.Length;
            double sum = 0.0;
            for (int i = 0; i < n; ++i)
                sum += v1[i] * v2[i];
            return sum;
        }

        static double Hypot(double x, double y)
        {
            // fancy sqrt(x^2 + y^2)
            double xabs = System.Math.Abs(x);
            double yabs = System.Math.Abs(y);
            double min, max;

            if (xabs < yabs)
            {
                min = xabs; max = yabs;
            }
            else
            {
                min = yabs; max = xabs;
            }

            if (min == 0)
                return max;

            double u = min / max;
            return max * System.Math.Sqrt(1 + u * u);
        }

        static double[] VecExtractFirst(double[] vec, int n)
        {
            double[] result = new double[n];
            for (int i = 0; i < n; ++i)
                result[i] = vec[i];
            return result;
        }
    }
}
