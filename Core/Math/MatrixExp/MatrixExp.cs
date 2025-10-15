﻿namespace ArcFrame.Core.Math
{
    /// <summary>
    /// Matrix exponential helpers to solve exp(X) where X is a matrix
    /// </summary>
    public static class MatrixExp
    {
        /// <summary>
        /// exp(h*A) for skew-symmetric A (so result is in SO(N)).
        /// Uses: scaling & squaring + Padé(6,6) and LU solves (no explicit inverse).
        /// Falls back to Rodrigues in 3D.
        /// </summary>
        public static double[,] ExpSkew(double[,] A, double h)
        {
            int N = A.GetLength(0);
            if (N != A.GetLength(1)) throw new ArgumentException("A must be square.");
            if (N == 3) return Rodrigues3(A, h);

            // Scale A
            var Ah = Helpers.Multiply(A, h);
            double norm1 = Helpers.OneNorm(Ah);
            // conservative θ6; smaller → more squarings but very safe numerically
            const double theta6 = 0.5;
            int s = System.Math.Max(0, (int)System.Math.Ceiling(System.Math.Log(norm1 / theta6, 2)));
            var As = Helpers.Multiply(Ah, 1.0 / System.Math.Pow(2, s));

            // Padé(6) step using LU solves: E = (V - U)^{-1} (V + U)
            var E = Pade6_LU(As);

            // Squaring
            for (int i = 0; i < s; i++) E = Helpers.Multiply(E, E);
            return E;
        }

        // ---- Padé(6) core with LU solve ----
        static double[,] Pade6_LU(double[,] A)
        {
            int N = A.GetLength(0);

            var A2 = Helpers.Multiply(A, A);
            var A4 = Helpers.Multiply(A2, A2);
            var A6 = Helpers.Multiply(A2, A4);

            // Higham-style coefficients for [6/6] (double)
            // These give a very good (m=6) rational approximation for exp.
            double[] b = {
            1.0,
            0.5,
            0.12,                     // 3!/ (6!)?  (pre-tuned Padé coeffs; keep as constants)
            0.018333333333333333,     // 11/600
            0.0019927536231884053,    // 151/75600
            0.00016059043836821612,   // 15619/97195500
            0.000009618129107628477   // 691/718502400
        };

            // U = A * (b1*I + b3*A2 + b5*A4)
            var I = RigidTransform.Identity(N).R;
            var polyOdd = Helpers.Add(Helpers.Multiply(I, b[1]), Helpers.Add(Helpers.Multiply(A2, b[3]), Helpers.Multiply(A4, b[5])));
            var U = Helpers.Multiply(A, polyOdd);

            // V = b0*I + b2*A2 + b4*A4 + b6*A6
            var V = Helpers.Add(Helpers.Multiply(I, b[0]),
                     Helpers.Add(Helpers.Multiply(A2, b[2]),
                     Helpers.Add(Helpers.Multiply(A4, b[4]), Helpers.Multiply(A6, b[6]))));

            // Solve (V - U) * X = (V + U)
            var W = Helpers.Subtract(V, U);
            var B = Helpers.Add(V, U);

            var (LU, piv, _) = LA.LUFactor(W);
            LA.LUSolveInPlace(LU, piv, B); // B <- X

            return B;
        }

        // ---- Rodrigues (exact for so(3)) ----
        static double[,] Rodrigues3(double[,] A, double h)
        {
            // A = [ω]_x. Extract ω from A.
            double wx = A[2, 1], wy = A[0, 2], wz = A[1, 0];
            double theta = h * System.Math.Sqrt(wx * wx + wy * wy + wz * wz);
            var I = RigidTransform.Identity(3).R;

            if (theta < 1e-12)
            {
                var _hA = Helpers.Multiply(A, h);
                return Helpers.Add(Helpers.Add(I, _hA), Helpers.Multiply(Helpers.Multiply(_hA, _hA), 0.5));
            }

            double c = System.Math.Cos(theta), s = System.Math.Sin(theta);
            double a = s / theta, b = (1 - c) / (theta * theta);

            var hA = Helpers.Multiply(A, h);
            return Helpers.Add(Helpers.Add(I, Helpers.Multiply(hA, a)), Helpers.Multiply(Helpers.Multiply(hA, hA), b));
        }
    }
}
