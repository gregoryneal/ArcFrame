using ArcFrame.Core.Math;
using System;

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

    /// <summary>
    /// Smart parameterization, converts the N-dimensional rotation matrix
    /// into a parameter vector in N(N-1)/2 dimensions. 
    /// 
    /// N = 2: planar rotation parameterized as [θ]
    /// N = 3: rodriguez formula
    /// N = 4: cayley factorization
    /// </summary>
    internal static class SOParam2
    {
        /// <summary>
        /// The parameter count of the NxN matrix.
        /// </summary>
        /// <param name="N"></param>
        /// <returns></returns>
        public static int ParamCount(int N) => (N * (N - 1)) / 2;

        /// <summary>
        /// Build a rotation matrix dR with theta, then
        /// multiply dR by R to update the rotation matrix.
        /// 
        /// R_f = dR * R_i
        /// 
        /// where dR = exp(A)
        /// where A is skew symmetric
        /// </summary>
        /// <param name="R"></param>
        /// <param name="theta"></param>
        public static void UpdateRotation(double[,] R, double[] theta)
        {
            int N = R.GetLength(0);
            double[,] dR;
            // Optimize for 2, 3 and 4d
            // left multiply
            switch (N)
            {
                //      | cosθ  -sinθ | | R00  R01 |   | R00cosθ - R10sinθ  R01cosθ - R11sinθ |
                // dR = |             | |          | = |                                      |
                //      | sinθ   cosθ | | R10  R11 |   | R00sinθ + R10cosθ  R01sinθ + R11cosθ |
                case 2:
                    double t = theta[0];
                    double c = Math.Cos(t);
                    double s = Math.Sin(t);
                    dR = new double[,] { { c, -s },
                                                   { s,  c } };
                    R = Helpers.Multiply(dR, R);
                    return;
                case 3:
                    // create skew matrix and use MatrixExp
                    double wx = theta[0], wy = theta[1], wz = theta[2];
                    double[,] A = new double[,] { { 0, -wz, wy },
                                                  { wz, 0, -wx },
                                                  { -wy, wx, 0 } };
                    double[,] exp = new double[3, 3];
                    MatrixExp.Rodrigues3(A, 1, ref exp);
                    R = Helpers.Multiply(exp, R);
                    return;
                case 4:
                    double a01 = theta[0], a02 = theta[1], a03 = theta[2],
                    a12 = theta[3], a13 = theta[4], a23 = theta[5];

                    // Build A and B = (I - 0.5 A)^{-1}(I + 0.5 A)

                    // I + 0.5A
                    double h = 0.5;

                    double[,] Up = new double[4, 4] {
                        {1,    h*a01, h*a02, h*a03},
                        {-h*a01,1,    h*a12, h*a13},
                        {-h*a02,-h*a12,1,    h*a23},
                        {-h*a03,-h*a13,-h*a23,1}
                    };

                    // I - 0.5A
                    double[,] Um = new double[4, 4] {
                        {1,    -h*a01, -h*a02, -h*a03},
                        {h*a01, 1,     -h*a12, -h*a13},
                        {h*a02, h*a12, 1,      -h*a23},
                        {h*a03, h*a13, h*a23,  1}
                    };

                    // Compute exp(A) ≈ Um^{-1} * Up
                    double[,] Umi = new double[4, 4];
                    if (Helpers.TryInvert(Um, out Umi!))
                    {
                        dR = Helpers.Multiply(Umi!, Up);
                        R = Helpers.Multiply(dR, R);
                        return;
                    }
                    goto default;
                default:
                    dR = SOParam.BuildRotation(N, theta);
                    R = Helpers.Multiply(dR, R);
                    return;
            }
        }
    }
}
