namespace ArcFrame.Core.Geometry
{
    using System;
    using ArcFrame.Core.Math;
    using ArcFrame.Core.Results;
    using ArcFrame.Solvers.G1;

    /// <summary>
    /// 3D arc-length curve implementing the Makino / Lan / Miura triple-clothoid
    /// construction via piecewise-quadratic pitch α(S) and yaw β(S) over
    /// normalized arc-length S in [0,1], with fixed internal breakpoints.
    /// 
    /// The curve is defined in a local frame where:
    ///   - P0 is at the origin,
    ///   - the chord (P1 - P0) is aligned with +X,
    ///   - and the initial tangent T0 lies in the X-Y plane as much as possible.
    /// 
    /// For now, the curvature vector k(s) is [κ(s), 0], where κ(s) is computed
    /// from α,β via ||d u / d s||, and torsion is left as 0. This is sufficient
    /// for many track-building use-cases where tangent and position dominate.
    /// </summary>
    public sealed class Makino3DTripleClothoidCurve : IArcLengthCurve
    {
        private readonly double[] _P0;
        private readonly double[] _ex;
        private readonly double[] _ey;
        private readonly double[] _ez;

        private readonly Makino3DTripleClothoidProblem.AngleCoeffs _alpha;
        private readonly Makino3DTripleClothoidProblem.AngleCoeffs _beta;
        private readonly double _h;

        private readonly double _alpha0;
        private readonly double _alpha1;
        private readonly double _beta0;
        private readonly double _beta1;

        private const int IntegrationSteps = 64;

        private Makino3DTripleClothoidCurve(
            double[] P0,
            double[] ex,
            double[] ey,
            double[] ez,
            Makino3DTripleClothoidProblem.AngleCoeffs alphaCoeffs,
            Makino3DTripleClothoidProblem.AngleCoeffs betaCoeffs,
            double h,
            double alpha0,
            double alpha1,
            double beta0,
            double beta1)
        {
            _P0 = (double[])P0.Clone();
            _ex = (double[])ex.Clone();
            _ey = (double[])ey.Clone();
            _ez = (double[])ez.Clone();
            _alpha = alphaCoeffs;
            _beta = betaCoeffs;
            _h = h;
            _alpha0 = alpha0;
            _alpha1 = alpha1;
            _beta0 = beta0;
            _beta1 = beta1;
        }

        /// <summary>
        /// Factory: rebuild the curve from endpoints and fitted parameters p = [qAlpha, qBeta, h].
        /// Recomputes the local frame and endpoint angles using the same logic as the problem.
        /// </summary>
        public static Makino3DTripleClothoidCurve CreateFromFit(
            double[] P0, double[] T0,
            double[] P1, double[] T1,
            double[] parameters)
        {
            if (P0 == null || T0 == null || P1 == null || T1 == null)
                throw new ArgumentNullException("CreateFromFit: null endpoints.");
            if (P0.Length != 3 || T0.Length != 3 || P1.Length != 3 || T1.Length != 3)
                throw new ArgumentException("CreateFromFit: only 3D data is supported.");
            if (parameters == null || parameters.Length != 3)
                throw new ArgumentException("CreateFromFit: parameters must have length 3.");

            var p0 = (double[])P0.Clone();
            var t0 = Helpers.Normalize((double[])T0.Clone());
            var p1 = (double[])P1.Clone();
            var t1 = Helpers.Normalize((double[])T1.Clone());

            // Local chord and frame (same as in the problem).
            var chord = Helpers.Subtract(p1, p0);
            double r2 = Helpers.Len2(chord);
            if (r2 < 1e-12)
            {
                // Degenerate: use tangent as chord.
                chord[0] = t0[0];
                chord[1] = t0[1];
                chord[2] = t0[2];
                r2 = Helpers.Len2(chord);
            }
            double r = Math.Sqrt(r2);

            double[] ex = Helpers.Multiply(1.0 / r, chord);

            double[] yCand = Helpers.Reject(t0, ex);
            if (Helpers.Len2(yCand) < 1e-8)
            {
                double[] up = new double[] { 0.0, 0.0, 1.0 };
                yCand = Helpers.Reject(up, ex);
                if (Helpers.Len2(yCand) < 1e-8)
                {
                    up = new double[] { 0.0, 1.0, 0.0 };
                    yCand = Helpers.Reject(up, ex);
                }
            }
            double[] ey = Helpers.Normalize(yCand);
            double[] ez = new double[3]
            {
                ex[1] * ey[2] - ex[2] * ey[1],
                ex[2] * ey[0] - ex[0] * ey[2],
                ex[0] * ey[1] - ex[1] * ey[0],
            };

            // Local tangents
            double[] T0loc = new double[3]
            {
                Helpers.Dot(t0, ex),
                Helpers.Dot(t0, ey),
                Helpers.Dot(t0, ez)
            };
            double[] T1loc = new double[3]
            {
                Helpers.Dot(t1, ex),
                Helpers.Dot(t1, ey),
                Helpers.Dot(t1, ez)
            };

            // Endpoint angles
            (double alpha0, double beta0) = ToAngles(T0loc);
            (double alpha1, double beta1) = ToAngles(T1loc);

            double qAlpha = parameters[0];
            double qBeta = parameters[1];
            double h = parameters[2];

            // Clamp h similarly to the fitter
            double minH = 0.25 * r;
            double maxH = 4.0 * Math.Max(1e-3, r);
            if (double.IsNaN(h) || double.IsInfinity(h) || h <= 0.0)
                h = r;
            if (h < minH) h = minH;
            if (h > maxH) h = maxH;

            // Precompute angle coefficients
            var alphaCoeffs = ComputeAngleCoeffs(alpha0, alpha1, qAlpha);
            var betaCoeffs = ComputeAngleCoeffs(beta0, beta1, qBeta);

            return new Makino3DTripleClothoidCurve(
                p0,
                ex,
                ey,
                ez,
                alphaCoeffs,
                betaCoeffs,
                h,
                alpha0,
                alpha1,
                beta0,
                beta1);
        }

        private static (double alpha, double beta) ToAngles(double[] u)
        {
            double ux = u[0];
            double uy = u[1];
            double uz = u[2];

            double beta = Math.Atan2(uy, ux);
            double cosAlpha = Math.Sqrt(ux * ux + uy * uy);
            if (cosAlpha < 1e-12)
            {
                cosAlpha = 0.0;
            }

            double alpha = Math.Atan2(-uz, cosAlpha);
            return (alpha, beta);
        }

        private static Makino3DTripleClothoidProblem.AngleCoeffs ComputeAngleCoeffs(
            double theta0, double theta1, double q)
        {
            double a10 = theta0;
            double a11 = 0.0;
            double a12 = q;

            double a20 = theta0 + q * (1.0 / 16.0);
            double a21 = q * 0.5;
            double a22 = -q - (8.0 / 3.0) * theta0 + (8.0 / 3.0) * theta1;

            double a30 = q * (1.0 / 16.0) + theta0 / 3.0 + 2.0 * theta1 / 3.0;
            double a31 = -q * 0.5 - (8.0 / 3.0) * theta0 + (8.0 / 3.0) * theta1;
            double a32 = q + (16.0 / 3.0) * theta0 - (16.0 / 3.0) * theta1;

            return new Makino3DTripleClothoidProblem.AngleCoeffs(
                a10, a11, a12,
                a20, a21, a22,
                a30, a31, a32);
        }

        public int Dimension => 3;

        public double Length => _h;

        public Sample Evaluate(double s)
        {
            s = Math.Clamp(s, 0.0, _h);

            if (s <= 0.0)
            {
                double[,] Rstart = ONFrame.R0_FromT_Complete(Tangent(0.0));
                double[] kstart = new double[2] { 0.0, 0.0 };
                return new Sample((double[])_P0.Clone(), Rstart, 0.0, kstart);
            }

            double S_end = s / _h;
            if (S_end < 0.0) S_end = 0.0;
            if (S_end > 1.0) S_end = 1.0;

            double xLoc, yLoc, zLoc;
            IntegrateLocal(_alpha, _beta, _h, S_end, IntegrationSteps,
                           out xLoc, out yLoc, out zLoc);

            double alpha = _alpha.Eval(S_end);
            double beta = _beta.Eval(S_end);

            double cosA = Math.Cos(alpha);
            double sinA = Math.Sin(alpha);
            double cosB = Math.Cos(beta);
            double sinB = Math.Sin(beta);

            double ux = cosA * cosB;
            double uy = cosA * sinB;
            double uz = -sinA;

            double[] Pglob = new double[3]
            {
                _P0[0] + xLoc * _ex[0] + yLoc * _ey[0] + zLoc * _ez[0],
                _P0[1] + xLoc * _ex[1] + yLoc * _ey[1] + zLoc * _ez[1],
                _P0[2] + xLoc * _ex[2] + yLoc * _ey[2] + zLoc * _ez[2],
            };

            double[] Tglob = new double[3]
            {
                ux * _ex[0] + uy * _ey[0] + uz * _ez[0],
                ux * _ex[1] + uy * _ey[1] + uz * _ez[1],
                ux * _ex[2] + uy * _ey[2] + uz * _ez[2],
            };

            // curvature magnitude from α'(s), β'(s); torsion left as 0
            double dAlpha_dS = _alpha.EvalPrimeWrtS(S_end);
            double dBeta_dS = _beta.EvalPrimeWrtS(S_end);
            double invH = 1.0 / _h;
            double dAlpha_ds = dAlpha_dS * invH;
            double dBeta_ds = dBeta_dS * invH;

            CurvatureTorsion(s, out double k, out double t);

            double[] kappa = new double[2] { k, t };

            double[,] R = ONFrame.R0_FromT_Complete(Tglob);

            return new Sample(Pglob, R, s, kappa);
        }

        public double[] Tangent(double s)
        {
            s = Math.Clamp(s, 0.0, _h);
            double S = (s <= 0.0) ? 0.0 : s / _h;

            double alpha = _alpha.Eval(S);
            double beta = _beta.Eval(S);

            double cosA = Math.Cos(alpha);
            double sinA = Math.Sin(alpha);
            double cosB = Math.Cos(beta);
            double sinB = Math.Sin(beta);

            double ux = cosA * cosB;
            double uy = cosA * sinB;
            double uz = -sinA;

            return new double[3]
            {
                ux * _ex[0] + uy * _ey[0] + uz * _ez[0],
                ux * _ex[1] + uy * _ey[1] + uz * _ez[1],
                ux * _ex[2] + uy * _ey[2] + uz * _ez[2],
            };
        }

        private static void IntegrateLocal(
            Makino3DTripleClothoidProblem.AngleCoeffs alpha,
            Makino3DTripleClothoidProblem.AngleCoeffs beta,
            double h,
            double S_end,
            int steps,
            out double xEnd,
            out double yEnd,
            out double zEnd)
        {
            if (steps < 1) steps = 1;

            double Smax = Math.Clamp(S_end, 0.0, 1.0);
            double dS = Smax / steps;
            double ds = h * dS;

            double px = 0.0;
            double py = 0.0;
            double pz = 0.0;

            for (int i = 0; i < steps; ++i)
            {
                double S = (i + 0.5) * dS;
                double a = alpha.Eval(S);
                double b = beta.Eval(S);

                double cosA = Math.Cos(a);
                double sinA = Math.Sin(a);
                double cosB = Math.Cos(b);
                double sinB = Math.Sin(b);

                double ux = cosA * cosB;
                double uy = cosA * sinB;
                double uz = -sinA;

                px += ux * ds;
                py += uy * ds;
                pz += uz * ds;
            }

            xEnd = px;
            yEnd = py;
            zEnd = pz;
        }

        public void CurvatureTorsion(double s, out double kappa, out double tau)
        {
            if (s <= 0.0 || s >= _h)
            {
                // endpoints: we designed it so curvature → 0
                kappa = 0.0;
                tau = 0.0;
                return;
            }

            double S = Math.Clamp(s / _h, 0.0, 1.0);

            double alpha = _alpha.Eval(S);
            double beta = _beta.Eval(S);

            double dAlpha_dS = _alpha.EvalPrimeWrtS(S);
            double dBeta_dS = _beta.EvalPrimeWrtS(S);
            double d2Alpha_dS2 = _alpha.EvalSecondPrimeWrtS(S);
            double d2Beta_dS2 = _beta.EvalSecondPrimeWrtS(S);

            double invH = 1.0 / _h;
            double dAlpha_ds = dAlpha_dS * invH;
            double dBeta_ds = dBeta_dS * invH;
            double d2Alpha_ds2 = d2Alpha_dS2 * invH * invH;
            double d2Beta_ds2 = d2Beta_dS2 * invH * invH;

            double[] u, up, upp;
            TangentDerivativesLocal(
                alpha, beta,
                dAlpha_ds, dBeta_ds,
                d2Alpha_ds2, d2Beta_ds2,
                out u, out up, out upp);

            // curvature: κ = ||u'||
            double upNorm2 = up[0] * up[0] + up[1] * up[1] + up[2] * up[2];
            double eps = 1e-12;
            if (upNorm2 < eps)
            {
                kappa = 0.0;
                tau = 0.0;
                return;
            }

            kappa = Math.Sqrt(upNorm2);

            // torsion: τ = det(u, u', u'') / ||u × u'||^2
            double cx = u[1] * up[2] - u[2] * up[1];
            double cy = u[2] * up[0] - u[0] * up[2];
            double cz = u[0] * up[1] - u[1] * up[0];
            double crossNorm2 = cx * cx + cy * cy + cz * cz;
            if (crossNorm2 < eps)
            {
                tau = 0.0;
                return;
            }

            double det = u[0] * (up[1] * upp[2] - up[2] * upp[1])
                       - u[1] * (up[0] * upp[2] - up[2] * upp[0])
                       + u[2] * (up[0] * upp[1] - up[1] * upp[0]);

            tau = det / crossNorm2;
        }

        private static void TangentDerivativesLocal(double alpha, double beta,
                                                    double dAlpha_ds, double dBeta_ds,
                                                    double d2Alpha_ds2, double d2Beta_ds2,
                                                    out double[] u, out double[] up, out double[] upp)
        {
            double cA = Math.Cos(alpha);
            double sA = Math.Sin(alpha);
            double cB = Math.Cos(beta);
            double sB = Math.Sin(beta);

            // u
            double ux = cA * cB;
            double uy = cA * sB;
            double uz = -sA;

            // ∂u/∂α and ∂u/∂β
            double ux_a = -sA * cB;
            double uy_a = -sA * sB;
            double uz_a = -cA;

            double ux_b = -cA * sB;
            double uy_b = cA * cB;
            double uz_b = 0.0;

            // second partials
            double ux_aa = -cA * cB;
            double uy_aa = -cA * sB;
            double uz_aa = sA;

            double ux_bb = -cA * cB;
            double uy_bb = -cA * sB;
            double uz_bb = 0.0;

            double ux_ab = sA * sB;
            double uy_ab = -sA * cB;
            double uz_ab = 0.0;

            double a1 = dAlpha_ds;
            double b1 = dBeta_ds;
            double a2 = d2Alpha_ds2;
            double b2 = d2Beta_ds2;

            // u' = a1*u_a + b1*u_b
            double upx = a1 * ux_a + b1 * ux_b;
            double upy = a1 * uy_a + b1 * uy_b;
            double upz = a1 * uz_a + b1 * uz_b;

            // u'' = a2*u_a + b2*u_b + a1²*u_aa + 2 a1 b1*u_ab + b1²*u_bb
            double a1a1 = a1 * a1;
            double b1b1 = b1 * b1;
            double two_a1b1 = 2.0 * a1 * b1;

            double uppx = a2 * ux_a + b2 * ux_b +
                          a1a1 * ux_aa + two_a1b1 * ux_ab + b1b1 * ux_bb;
            double uppy = a2 * uy_a + b2 * uy_b +
                          a1a1 * uy_aa + two_a1b1 * uy_ab + b1b1 * uy_bb;
            double uppz = a2 * uz_a + b2 * uz_b +
                          a1a1 * uz_aa + two_a1b1 * uz_ab + b1b1 * uz_bb;

            u = new double[3] { ux, uy, uz };
            up = new double[3] { upx, upy, upz };
            upp = new double[3] { uppx, uppy, uppz };
        }
    }
}