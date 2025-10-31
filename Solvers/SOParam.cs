using ArcFrame.Core.Math;

namespace ArcFrame.Solvers
{
    /// Minimal parameterization of R in SO(N) via planar Givens rotations.
    /// We map a vector w of length N*(N-1)/2 to a left-multiplied rotation:
    ///    R(w) = Π_{i<j} G(i,j, θ_ij)
    /// Zero vector => Identity. Stable for LM with FD Jacobians.
    internal static class SOParam
    {
        /// <summary>
        /// The parameter count of the NxN matrix.
        /// </summary>
        /// <param name="N"></param>
        /// <returns></returns>
        public static int ParamCount(int N) => (N * (N - 1)) / 2;

        /// Build rotation from w (applied on the LEFT of the seed).
        public static double[,] BuildRotation(int N, double[] omega)
        {
            double[,] R = RigidTransform.Identity(N).R;
            int k = 0;
            for (int i = 0; i < N - 1; i++)
                for (int j = i + 1; j < N; j++)
                {
                    double theta = omega[k++];
                    LeftMultiplyGivens(R, N, i, j, theta);
                }
            return R;
        }

        /// R <- G(i,j,θ) * R   where G is identity except 2x2 block [[c,-s],[s,c]] in (i,j)-plane.
        private static void LeftMultiplyGivens(double[,] R, int N, int i, int j, double theta)
        {
            double c = System.Math.Cos(theta);
            double s = System.Math.Sin(theta);
            for (int col = 0; col < N; col++)
            {
                double Ri = R[i, col];
                double Rj = R[j, col];
                R[i, col] = c * Ri - s * Rj;
                R[j, col] = s * Ri + c * Rj;
            }
        }
    }
}
