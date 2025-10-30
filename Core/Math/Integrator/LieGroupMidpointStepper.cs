
using System;

namespace ArcFrame.Core.Math
{
    /// <summary>
    /// Minimal Lie group stepper using exp of small skew matrices
    /// </summary>
    public sealed class LieGroupMidpointStepper : IFrameStepper
    {
        /// <summary>
        /// Integrate R′ = RA(s), P′=Re1​
        /// where A(s) is skew-symmetric (so R∈SO(N) stays orthonormal)
        /// </summary>
        /// <param name="P">Current position at arc length s</param>
        /// <param name="R">Current orthonormal basis vectors SO(n) at arc length s</param>
        /// <param name="kappa">Function that calculates curvature at arc length s</param>
        /// <param name="s">The current arc length</param>
        /// <param name="h">The step size, new arc length = s + h</param>
        /// <param name="frame">Type of frame system we are using</param>
        /// <exception cref="NotImplementedException"></exception>
        public void Step(ref double[] P, ref double[,] R, Func<double, double[]> kappa, double s, double h, FrameModel frame)
        {
            /*
            int N = R.GetLength(0);
            var k = kappa(s + 0.5 * h);

            var AF = ONFrame.BuildFrenetSkew(k);
            var Ehalf = MatrixExp.ExpSkew(AF, 0.5 * h);
            var Efull = MatrixExp.ExpSkew(AF, h);

            // position update with midpoint tangent
            var Tmid = MultiplyRe1(R, Ehalf);
            for (int i = 0; i < N; i++) P[i] += h * Tmid[i];

            if (frame == FrameModel.Frenet)
            {
                // advance frame
                R = Helpers.Multiply(R, Efull);
                return;
            }

            // Bishop (RMF), shape-preserving:
            else {
                var RF = Helpers.Multiply(R, Efull);

                // 2) Extract the normal-normal block O_F (indices 1..N-1)
                double[,] Omega = ONFrame.ExtractNormalBlock(AF);        // (N-1)x(N-1), skew

                // 3) Untwist normal space: Q = exp(-O_F h), then R = R_F * blockdiag(1, Q)
                double[,] Q = MatrixExp.ExpSkew(Omega, -h);
                double[,] B = BlockDiag1Q(Q);
                R = Helpers.Multiply(RF, B);
                return;
            }

            // local helpers
            double[] MultiplyRe1(double[,] Rcur, double[,] EhalfCur)
            {
                var v = new double[N];
                for (int i = 0; i < N; i++)
                {
                    double srow = 0;
                    for (int j = 0; j < N; j++) srow += Rcur[i, j] * EhalfCur[j, 0];
                    v[i] = srow;
                }
                return v;
            }

            double[,] BlockDiag1Q(double[,] Q)
            {
                int m = Q.GetLength(0);
                var B = new double[m + 1, m + 1];
                B[0, 0] = 1.0;
                for (int i = 0; i < m; i++)
                    for (int j = 0; j < m; j++)
                        B[i + 1, j + 1] = Q[i, j];
                return B;
            }*/
            
            int N = R.GetLength(0);
            var k = kappa(s + 0.5 * h);                  // midpoint sample (size N-1)

            // Build skew A at the midpoint
            var A = (frame == FrameModel.Frenet)
                  ? ONFrame.BuildFrenetSkew(k)                   // uses FS invariants
                  : ONFrame.BuildBishopSkew(k);                  // uses Bishop components

            // Exponentials on SO(N)
            double[,] Ehalf = MatrixExp.ExpSkew(A, 0.5 * h);
            double[,] Efull = MatrixExp.ExpSkew(A, h);

            // T_mid = R * Ehalf * e1 
            var Tmid = new double[N];
            for (int i = 0; i < N; i++)
            {
                double srow = 0.0;
                for (int j = 0; j < N; j++) srow += R[i, j] * Ehalf[j, 0];
                Tmid[i] = srow;
            }

            // P <- P + h * T_mid
            for (int i = 0; i < N; i++) P[i] += h * Tmid[i];

            // R <- R * Efull
            R = Helpers.Multiply(R, Efull);
        }
    }
}
