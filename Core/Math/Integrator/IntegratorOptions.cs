using System;

namespace ArcFrame.Core.Math
{
    /// <summary>
    /// Settings for integrating the moving frame along the curve.
    /// </summary>
    public sealed class IntegratorOptions
    {
        /// <summary>
        /// Maximum step size (absolute upper bound).
        /// </summary>
        public double MaxStep { get; set; } = 0.05;

        /// <summary>
        /// Minimum step size (guards against shrinking to zero).
        /// </summary>
        public double MinStep { get; set; } = 1e-5;

        /// <summary>
        /// Limit on frame rotation per step: ||kappa(s)|| * h <= ThetaPerStep.        
        /// Typical: 0.05..0.2 rad (≈3°..12°)
        /// </summary>
        public double ThetaPerStep { get; set; } = 0.1;

        /// <summary>
        /// Optional: bound how much curvature is allowed to *change* in one step
        /// (measured as L2 norm of delta-kappa). Set <= 0 to disable.
        /// </summary>
        public double CurvDeltaPerStep { get; set; } = 0.0;

        /// <summary>
        /// Probe fraction for curvature-change estimate (relative to proposed h).
        /// </summary>
        public double CurvProbeFrac { get; set; } = 0.5;

        /// <summary>
        /// Suggest a step h at arclength s0 toward sTarget using kappa accessor.
        /// </summary>
        /// <param name="kappa">Principal curvatures at arc length s. Length N-1 where N is the intrinsic curve dimension</param>
        /// <param name="s0">start arc length</param>
        /// <param name="sTarget">target arc length</param>
        /// <returns></returns>
        public double SuggestStep(Func<double, double[]> kappa, double s0, double sTarget)
        {
            double remaining = System.Math.Max(0, sTarget - s0);
            if (remaining <= MinStep) return remaining;

            // 1) rotation-based bound: h <= ThetaPerStep / ||k(s0)||
            var k0 = kappa(s0);
            double k0norm = Helpers.Len(k0);
            double hRot = (k0norm > 1e-12) ? (ThetaPerStep / k0norm) : MaxStep;

            double h = System.Math.Min(MaxStep, System.Math.Min(hRot, remaining));

            // 2) curvature-change bound (optional)
            if (CurvDeltaPerStep > 0.0)
            {
                double probe = System.Math.Max(MinStep, CurvProbeFrac * h);
                probe = System.Math.Min(probe, remaining);
                var k1 = kappa(s0 + probe);
                double dk = Norm2Diff(k1, k0);
                if (dk > 1e-12)
                {
                    // scale h so that (dk / probe) * h <= CurvDeltaPerStep
                    double dk_ds = dk / probe;
                    double hCurv = CurvDeltaPerStep / dk_ds;
                    h = System.Math.Min(h, System.Math.Max(MinStep, hCurv));
                }
            }

            if (h < MinStep && remaining > 0) h = System.Math.Min(remaining, MinStep);
            return h;

            static double Norm2Diff(double[] a, double[] b) { double s = 0; for (int i = 0; i < a.Length; i++) { double d = a[i] - b[i]; s += d * d; } return System.Math.Sqrt(s); }
        }

        public static IntegratorOptions Default => new IntegratorOptions();
    }

}
