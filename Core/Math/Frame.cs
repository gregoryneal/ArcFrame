using System;
using System.Collections.Generic;

namespace ArcFrame.Core.Math
{
    /// <summary>
    /// Type of frame model, in the Frenet case, the frame evolves as a function of curvature. In the degenerate case of 0 curvature it gets all wonky.
    /// The Bishop frame solves these issues by evolving the frame with minimal rotation. It will basically follow the tangent vector in this case.
    /// </summary>
    public enum FrameModel
    {
        Frenet,
        Bishop
    }

    /// <summary>
    /// Orthonormal frame in D dimensions
    /// </summary>
    public struct Frame
    {
        /// <summary>
        /// Initial position of the frame
        /// </summary>
        public double[] P;
        /// <summary>
        /// Orthonormal column vectors in N dimension
        /// [T, N, B, ... ]
        /// </summary>
        public double[,] R0;
    }

    /// <summary>
    /// Helper methods for the creation of orthonormal frames in N dimensions
    /// </summary>
    public static class ONFrame
    {
        private static double eps = 1E-12;

        public static double[,] R0_FromTNB_Complete(double[] T, double[] N, double[] B)
        {
            return R0_FromR_Complete(R0_FromTNB(T, N, B));
        }

        public static double[,] R0_FromTNB(double[] T, double[] N, double[] B)
        {
            double[,] R = new double[T.Length, 3];
            SetCol(R, 0, T);
            SetCol(R, 1, N);
            SetCol(R, 2, B);
            return R;
        }

        public static double[,] R0_FromTN_Complete(double[] T, double[] N)
        {
            return R0_FromR_Complete(R0_FromTN(T, N));
        }

        public static double[,] R0_FromTN(double[] T, double[] N)
        {
            double[,] R = new double[T.Length, 2];
            SetCol(R, 0, T);
            SetCol(R, 1, N);
            return R;
        }

        /// <summary>
        /// Create an orthonormal vector matrix from a starting direction, your curve tangent, and an optional up vector to bias the second orthonormal basis vector N1.
        /// </summary>
        /// <param name="T">D Dimensional tangent vector</param>
        /// <param name="up">D Dimensional "up" vector, optional.</param>
        /// <returns>a DxD matrix with column vectors C = [T=N0, N1, N2, N3, ... ND-1] being orthonormal to all other columns</returns>
        public static double[,] R0_FromT_Complete(double[] T, double[]? up = null)
        {
            int n = T.Length;
            double[][] cols = new double[T.Length][];
            cols[0] = Helpers.Normalize(T);
            double[] n1;
            if (up == null)
            {
                //Find the basis least aligned with norm(T), find the rejection vector
                n1 = Helpers.Reject(Helpers.StdSeed(cols[0]), cols[0]);
            }
            else
            {
                n1 = Helpers.Reject(up, cols[0]);
                if (Helpers.Len2(n1) < eps) n1 = Helpers.Reject(Helpers.StdSeed(cols[0]), cols[0]);
            }

            cols[1] = Helpers.Normalize(n1);

            double[] seed;
            double[] v;
            bool ok;
            //Graham schmidt style orthonormalization
            int c = 2;
            for (int k = 0; k < n && c < n; k++)
            {
                seed = Helpers.StdBasis(n, k);
                ok = true;

                for (int j = 0; j < c; j++)
                {
                    //check if any of our basis vectors are close to the seed
                    if (System.Math.Abs(Helpers.Dot(seed, cols[j])) > 0.9)
                    {
                        ok = false;
                        break;
                    }
                }

                if (!ok) continue;

                v = (double[])seed.Clone();
                for (int j = 0; j < c; j++) v = Helpers.Reject(v, cols[j]);
                if (Helpers.Len2(v) < eps) continue;
                cols[c++] = Helpers.Normalize(v);
            }

            while (c < n) cols[c++] = Helpers.StdBasis(n, c - 1);

            double[,] R = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    R[i, j] = cols[j][i];
                }
            }
            ValidateONColumns(R, 1E-7);
            return R;
        }

        /// <summary>
        /// Given an Nxd matrix where the columns are orthonormal basis vectors in R^N,
        /// and N < d, Calculate the other N-d basis vectors and add them to the matrix.
        /// 
        /// Example:
        ///     in -> [
        ///             [0, 1]
        ///             [0, 0]
        ///             [1, 0]
        ///           ]
        ///           
        /// N = 3, d = 2. We have 2 (d) orthonormal basis vectors in R^3 (R^N). We want
        /// the other basis vector. This outputs:
        ///     out -> [
        ///              [0, 1, 0]
        ///              [0, 0, 1]
        ///              [1, 0, 0]
        ///            ]
        /// </summary>
        /// <param name="A">The incomplete basis matrix</param>
        /// <returns></returns>
        public static double[,] R0_FromR_Complete(double[,] A, double tol = 1E-10)
        {
            int N = A.GetLength(0);
            int d = A.GetLength(1);
            if (d > N) throw new ArgumentException("d must be <= N.");
            if (d == N) return (double[,])A.Clone();

            // Copy A into a working list of columns (as double[] for convenience)
            var cols = new List<double[]>(d);
            for (int j = 0; j < d; j++)
                cols.Add(GetCol(A, j));

            // Verify approximate orthonormality
            for (int i = 0; i < d; i++)
                for (int j = i; j < d; j++)
                {
                    double dot = Helpers.Dot(cols[i], cols[j]);
                    if ((i == j && System.Math.Abs(dot - 1.0) > 1e-6) ||
                        (i != j && System.Math.Abs(dot) > 1e-6))
                    {
                        // Not strictly necessary to throw; you could re-orthonormalize A first.
                        throw new ArgumentException("Input columns are not orthonormal within tolerance.");
                    }
                }

            // Candidates: start with standard basis e1..eN; you can also add randoms if needed.
            for (int k = 0; k < N && cols.Count < N; k++)
            {
                var v = new double[N];
                v[k] = 1.0; // e_k

                // Project out components along existing columns
                OrthoProjectOut(v, cols);

                double nrm = Helpers.Norm(v);
                if (nrm > tol)
                {
                    Helpers.Multiply(v, 1.0 / nrm);
                    cols.Add(v);
                }
            }

            // If still short (pathological alignment), add a few randoms
            var rng = new Random(1234);
            int guard = 0;
            while (cols.Count < N && guard < 10 * N)
            {
                var v = new double[N];
                for (int i = 0; i < N; i++) v[i] = rng.NextDouble() - 0.5;
                OrthoProjectOut(v, cols);
                double nrm = Helpers.Norm(v);
                if (nrm > tol)
                {
                    Helpers.Multiply(v, 1.0 / nrm);
                    cols.Add(v);
                }
                guard++;
            }

            if (cols.Count != N)
                throw new InvalidOperationException("Failed to complete orthonormal basis; try loosening tol.");

            // Pack into N x N matrix [A | B]
            var M = new double[N, N];
            for (int j = 0; j < N; j++)
                SetCol(M, j, cols[j]);

            return M;
        }

        /// <summary>
        /// Frenet: tridiagonal skew with κ on super/sub-diagonal
        /// </summary>
        /// <param name="k">Principal curvatures</param>
        /// <returns></returns>
        public static double[,] BuildFrenetSkew(double[] k)
        {
            int m = k.Length + 1;                 // m = N
            var A = new double[m, m];              // zeros
            for (int i = 0; i < k.Length; i++)
            {
                A[i + 1, i] = k[i];
                A[i, i + 1] = -k[i];
            }
            return A;
        }

        /// <summary>
        /// Bishop (parallel transport): rotate T minimally; k is the curvature vector in current normal basis.
        /// Here we interpret k[0]..k[N-2] as components along the current normal columns N1..N_{N-1}.
        /// </summary>
        /// <param name="k"></param>
        /// <returns></returns>
        public static double[,] BuildBishopSkew(double[] k)
        {
            int N = k.Length + 1;
            var A = new double[N, N]; // zeros
                                      // Couple T (col 0) with normal columns via k components; keep normal plane rotation 0.
                                      // A_{0,j} = -k_{j-1}, A_{j,0} = k_{j-1} for j=1..N-1
            for (int j = 1; j < N; j++) 
            { 
                A[0, j] = -k[j - 1]; 
                A[j, 0] = k[j - 1]; 
            }
            return A;
        }
        /*
        public static double[,] BuildBishopSkew(double[] k)
        {
            int N = k.Length + 1;
            var A = new double[N, N]; // zeros
                                      // Couple T (col 0) with normal columns via k components; keep normal plane rotation 0.
                                      // A_{0,j} = -k_{j-1}, A_{j,0} = k_{j-1} for j=1..N-1
            for (int j = 1; j < N; j++)
            {
                A[0, j] = -k[j - 1];
                A[j, 0] = k[j - 1];
            }
            return A;
        }*/

        //helpers
        public static double[] GetCol(double[,] M, int j)
        {
            int N = M.GetLength(0);
            var c = new double[N];
            for (int i = 0; i < N; i++) c[i] = M[i, j];
            return c;
        }
        public static void SetCol(double[,] M, int j, double[] c)
        {
            for (int i = 0; i < M.GetLength(0); i++) M[i, j] = c[i];
        }
        private static void Axpy(double[] y, double alpha, double[] x)
        {
            for (int i = 0; i < y.Length; i++) y[i] += alpha * x[i];
        }
        private static void OrthoProjectOut(double[] v, List<double[]> basis)
        {
            // Modified Gram–Schmidt step: v := v - sum_j (q_j^T v) q_j
            for (int j = 0; j < basis.Count; j++)
            {
                double coeff = Helpers.Dot(basis[j], v);
                if (System.Math.Abs(coeff) > 0) Axpy(v, -coeff, basis[j]);
            }
            // (Optional) one re-orthogonalization pass for stability
            for (int j = 0; j < basis.Count; j++)
            {
                double coeff = Helpers.Dot(basis[j], v);
                if (System.Math.Abs(coeff) > 0) Axpy(v, -coeff, basis[j]);
            }
        }
        /// <summary>
        /// Ensure columns are orthonormal
        /// </summary>
        /// <param name="E"></param>
        /// <param name="tol"></param>
        /// <exception cref="ArgumentException"></exception>
        internal static void ValidateONColumns(double[,] E, double tol)
        {
            int n = E.GetLength(0), d = E.GetLength(1);
            // Check E^T E ≈ I_d
            for (int a = 0; a < d; a++)
                for (int b = 0; b < d; b++)
                {
                    double s = 0;
                    for (int i = 0; i < n; i++) s += E[i, a] * E[i, b];
                    if (a == b)
                    {
                        if (System.Math.Abs(s - 1.0) > 10 * tol)
                            throw new ArgumentException("E columns are not unit length.");
                    }
                    else
                    {
                        if (System.Math.Abs(s) > 10 * tol)
                            throw new ArgumentException("E columns are not orthogonal.");
                    }
                }
        }

        internal static double[,] ExtractNormalBlock(double[,] aF)
        {
            int n = aF.GetLength(0) - 1;
            int m = aF.GetLength(1) - 1;
            double[,] a = new double[n, m];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    a[i, j] = aF[i + 1, j + 1];
                }
            }
            return a;
        }
    }
}
