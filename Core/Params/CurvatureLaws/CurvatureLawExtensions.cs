using ArcFrame.Core.Math;
using System;

namespace ArcFrame.Core.Params.CurvatureLaws
{
    /// <summary>
    /// Helper methods for Curvature laws
    /// </summary>
    public static class CurvatureLawExtensions
    {
        /// Returns a jet up to `order` (>=0). outJet[o][j] = d^o κ_j(s)/ds^o.
        /// <summary>
        /// Returns a jet up to `order` (>=0). outJet[o][j] = d^o κ_j(s)/ds^o.
        /// </summary>
        /// <param name="law"></param>
        /// <param name="s"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double[][] EvalJet(this ICurvatureLaw law, double s, int order = 1)
        {
            if (order < 0) order = 0;
            var jet = AllocateJet(law.Order, order);

            // 1) Exact path if available
            if (law is ICurvatureJetProvider jp && order <= jp.MaxSupportedJetOrder)
            {
                jet = jp.EvalJet(s, order);
                return jet;
            }

            if (law is LinearCurvatureLaw lin) return LinearJet(lin, s, order);
            //if (law is SplineCurvatureLaw spline) return spline.EvalJet(s, order);

            return NumericJet(law, s, order);
        }

        static double[][] AllocateJet(int comps, int order)
        {
            var jet = new double[order + 1][];
            for (int o = 0; o <= order; o++) jet[o] = new double[comps];
            return jet;
        }

        static double[][] LinearJet(LinearCurvatureLaw lin, double s, int order)
        {
            var jet = AllocateJet(lin.Order, order);
            // κ(s) = k0 + a*s
            jet[0] = lin.Eval(s);
            var a = lin.Slope; // double[] of length Order
            if (order >= 1) for (int j = 0; j < lin.Order; j++) jet[1][j] = a[j];
            if (order >= 2) { /* zeros already */ }
            return jet;
        }

        static double[][] NumericJet(ICurvatureLaw law, double s, int order)
        {
            int m = law.Order;
            var jet = AllocateJet(m, order);

            // base value
            jet[0] = law.Eval(s);

            // step heuristic
            const double eps = 2.22e-16;                // machine eps
            double scale = System.Math.Max(1.0, System.Math.Abs(s));
            double h1 = System.Math.Sqrt(eps) * scale;         // good for first derivative
            double h = System.Math.Max(1e-9, h1);

            // 5-point central for 1st & 2nd; Richardson refine once
            if (order >= 1)
            {
                var fm2 = law.Eval(s - 2 * h);
                var fm1 = law.Eval(s - h);
                var fp1 = law.Eval(s + h);
                var fp2 = law.Eval(s + 2 * h);

                // f'(s) ≈ ( -f(s+2h) + 8 f(s+h) - 8 f(s-h) + f(s-2h) ) / (12h)
                for (int j = 0; j < m; j++)
                    jet[1][j] = (-fp2[j] + 8 * fp1[j] - 8 * fm1[j] + fm2[j]) / (12 * h);

                if (order >= 2)
                {
                    // f''(s) ≈ ( -f(s+2h) + 16 f(s+h) - 30 f(s) + 16 f(s-h) - f(s-2h) ) / (12 h^2)
                    for (int j = 0; j < m; j++)
                        jet[2][j] = (-fp2[j] + 16 * fp1[j] - 30 * jet[0][j] + 16 * fm1[j] - fm2[j]) / (12 * h * h);
                }
            }

            if (order >= 3)
            {
                // simple higher order via smaller step + Richardson (kept compact)
                double h2 = 0.5 * h;
                var j1 = NumericJet(law, s, 2);        // base (h)
                var j2 = NumericJetWithStep(law, s, 2, h2);

                // Richardson for first derivative (order-4 vs order-2): D ≈ (16 D(h/2) - D(h)) / 15
                for (int j = 0; j < m; j++)
                {
                    jet[1][j] = (16 * j2[1][j] - j1[1][j]) / 15.0;
                    jet[2][j] = (16 * j2[2][j] - j1[2][j]) / 15.0;
                }

                // 3rd derivative by finite difference of 2nd (reuse values)
                double hh = System.Math.Max(1e-9, System.Math.Cbrt(eps) * scale);
                var k_m = new double[m]; var k_0 = new double[m]; var k_p = new double[m];
                // use NumericJetWithStep to get f'' at s±hh
                var jm = NumericJetWithStep(law, s - hh, 2, h);
                var j0 = NumericJetWithStep(law, s, 2, h);
                var jp = NumericJetWithStep(law, s + hh, 2, h);
                for (int j = 0; j < m; j++)
                    jet[3][j] = (jp[2][j] - jm[2][j]) / (2 * hh);
            }

            return jet;
        }

        static double[][] NumericJetWithStep(ICurvatureLaw law, double s, int order, double h)
        {
            int m = law.Order;
            var jet = AllocateJet(m, order);
            jet[0] = law.Eval(s);

            var fm2 = law.Eval(s - 2 * h);
            var fm1 = law.Eval(s - h);
            var fp1 = law.Eval(s + h);
            var fp2 = law.Eval(s + 2 * h);


            if (order >= 1)
                for (int j = 0; j < m; j++)
                    jet[1][j] = (-fp2[j] + 8 * fp1[j] - 8 * fm1[j] + fm2[j]) / (12 * h);
            if (order >= 2)
                for (int j = 0; j < m; j++)
                    jet[2][j] = (-fp2[j] + 16 * fp1[j] - 30 * jet[0][j] + 16 * fm1[j] - fm2[j]) / (12 * h * h);
            return jet;
        }

        /// <summary>
        /// Convert a Frenet ICurvatureLaw (k1,k2,...,k_{N-1}) into a Bishop/RMF FunctionCurvatureLaw u(s).
        /// Assumes s is queried monotonically nondecreasing.
        /// </summary>
        public static FunctionCurvatureLaw FromFrenetToBishopFunction(
            this ICurvatureLaw frenetLaw,  // returns [k1, k2, ..., k_{N-1}]
            double s0 = 0.0,
            double[,]? Q0 = null,     // optional initial normal rotation (SO(N-1)); null => Identity
            double maxStep = 0.05)    // should be ≤ your IntegratorOptions.MaxStep
        {
            int N = frenetLaw.Order + 1;
            int m = N - 1; // normal-space size
            if (m <= 0) throw new ArgumentException("N must be ≥ 2.");

            // normal rotation state Q(s) in SO(m)
            double[,] Q = Q0 ?? RigidTransform.Identity(m).R;
            double sLast = s0;

            // pre-allocate scratch
            double[] k = new double[m];     // holds kappa[0..m-1]
            double[,] Omega = new double[m, m];

            return new FunctionCurvatureLaw(m, s =>
            {
                // advance Q from sLast to s with Lie midpoint on SO(m)
                if (s < sLast - 1e-12)
                {
                    // if someone goes backwards, restart from s0
                    Q = Q0 ?? RigidTransform.Identity(m).R;
                    sLast = s0;
                }

                while (sLast < s - 1e-16)
                {
                    double h = System.Math.Min(maxStep, s - sLast);
                    double smid = sLast + 0.5 * h;

                    var kMid = frenetLaw.Eval(smid); // [k1,k2,...]
                    Omega = ONFrame.ExtractNormalBlock(ONFrame.BuildFrenetSkew(kMid)); // Ω_F from torsions (k2..)

                    // Q <- exp(Ω_F * -h) * Q 
                    var E = MatrixExp.ExpSkew(Omega, -h);
                    Q = Helpers.Multiply(E, Q);

                    sLast += h;
                }

                // Now emit u(s) = k1(s) * Q^T e1  (size m)
                var kNow = frenetLaw.Eval(s);
                double k1 = kNow[0];
                var u = new double[m];

                // v = Q^T e1  ⇒ v_i = Q_{0,i} (first row of Q)
                for (int i = 0; i < m; i++) u[i] = k1 * Q[0, i];
                return u;
            });
        }
    }
}
