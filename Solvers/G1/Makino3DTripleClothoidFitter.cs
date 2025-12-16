using ArcFrame.Core.Math;
using ArcFrame.Core.Geometry;
using ArcFrame.Solvers.Core;
using System;

namespace ArcFrame.Solvers.G1
{
    /// <summary>
    /// Internal nonlinear problem for fitting a 3D triple "Makino" clothoid between
    /// two poses [P0, T0] and [P1, T1], following the construction in the
    /// Makino / Lan / Miura 3D triple clothoid paper.
    /// 
    /// Parameter vector p has length 3:
    ///  p[0] = qAlpha = a12  (second derivative coefficient for pitch angle alpha(S) on [0, S1])
    ///  p[1] = qBeta  = b12  (second derivative coefficient for yaw   angle beta(S)  on [0, S1])
    ///  p[2] = h      = total arc length of the triple segment
    /// 
    /// Residual vector r has length 3 (local chord constraints):
    ///  r[0] = (x(h) - r) / r
    ///  r[1] =  y(h)      / r
    ///  r[2] =  z(h)      / r
    /// where (x(h), y(h), z(h)) is the endpoint of the integrated curve in the
    /// local frame where the chord is aligned with +X and P0 is at the origin.
    /// </summary>
    internal sealed class Makino3DTripleClothoidProblem : ISmallNonlinearLeastSquaresProblem
    {
        public int ParameterCount => 3;
        public int ResidualCount => 3;

        // Endpoints in global coordinates
        private readonly double[] _P0;
        private readonly double[] _T0;
        private readonly double[] _P1;
        private readonly double[] _T1;

        // Local chord length
        private readonly double _r;
        private readonly double _invR;

        // Local basis vectors (columns of R0)
        // ex: along chord, ey: approximately along T0 projected onto chord-orthogonal plane, ez: ex x ey.
        private readonly double[] _ex;
        private readonly double[] _ey;
        private readonly double[] _ez;

        // Endpoint angles in local frame
        internal readonly double _alpha0;
        internal readonly double _alpha1;
        internal readonly double _beta0;
        internal readonly double _beta1;

        // Fixed normalized breakpoints, following the triple-clothoid papers.
        private const double S1 = 0.25;
        private const double S2 = 0.75;

        // Integration resolution for the chord constraint integrals.
        private const int IntegrationSteps = 64;

        public Makino3DTripleClothoidProblem(double[] P0, double[] T0, double[] P1, double[] T1)
        {
            if (P0 == null || T0 == null || P1 == null || T1 == null)
                throw new ArgumentNullException("Makino3DTripleClothoidProblem: null endpoint data.");
            if (P0.Length != 3 || T0.Length != 3 || P1.Length != 3 || T1.Length != 3)
                throw new ArgumentException("Makino3DTripleClothoidProblem: only 3D endpoints are supported.");

            _P0 = (double[])P0.Clone();
            _T0 = Helpers.Normalize((double[])T0.Clone());
            _P1 = (double[])P1.Clone();
            _T1 = Helpers.Normalize((double[])T1.Clone());

            // chord and local frame
            var chord = Helpers.Subtract(_P1, _P0);
            double r2 = Helpers.Len2(chord);
            if (r2 < 1e-12)
            {
                // Degenerate: treat as tiny chord
                chord[0] = _T0[0];
                chord[1] = _T0[1];
                chord[2] = _T0[2];
                r2 = Helpers.Len2(chord);
            }

            _r = Math.Sqrt(r2);
            _invR = 1.0 / Math.Max(1e-6, _r);

            _ex = Helpers.Multiply(1.0 / _r, chord);

            // Build ey from T0 projected onto chord-orthogonal plane,
            // with fallbacks if T0 ~ parallel to ex.
            double[] yCand = Helpers.Reject(_T0, _ex);
            if (Helpers.Len2(yCand) < 1e-8)
            {
                // Try global up
                double[] up = new double[] { 0.0, 0.0, 1.0 };
                yCand = Helpers.Reject(up, _ex);
                if (Helpers.Len2(yCand) < 1e-8)
                {
                    // Try another axis
                    up = new double[] { 0.0, 1.0, 0.0 };
                    yCand = Helpers.Reject(up, _ex);
                }
            }
            _ey = Helpers.Normalize(yCand);
            _ez = Cross(_ex, _ey);

            // Local tangents
            var T0loc = ToLocal(_T0);
            var T1loc = ToLocal(_T1);

            // Endpoint pitch/yaw angles in local frame, matching the
            // paper's convention:
            //  u = [cosα cosβ, cosα sinβ, -sinα]^T
            (_alpha0, _beta0) = ToAngles(T0loc);
            (_alpha1, _beta1) = ToAngles(T1loc);
        }

        /// <summary>
        /// Evaluate the residual for a given parameter vector p = [qAlpha, qBeta, h].
        /// </summary>
        public void Evaluate(double[] parameters, double[] residuals)
        {
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));
            if (residuals == null) throw new ArgumentNullException(nameof(residuals));
            if (parameters.Length != ParameterCount) throw new ArgumentException("parameter length mismatch");
            if (residuals.Length != ResidualCount) throw new ArgumentException("residual length mismatch");

            double qAlpha = parameters[0];
            double qBeta = parameters[1];
            double h = parameters[2];

            // Clamp the parameter range modestly to avoid wild oscillations.
            const double QMax = 8.0;
            if (qAlpha > QMax) qAlpha = QMax;
            if (qAlpha < -QMax) qAlpha = -QMax;
            if (qBeta > QMax) qBeta = QMax;
            if (qBeta < -QMax) qBeta = -QMax;

            // Total length bounded relative to chord length.
            double minH = 0.25 * _r;
            double maxH = 4.0 * Math.Max(1e-3, _r);
            if (double.IsNaN(h) || double.IsInfinity(h) || h <= 0.0)
                h = _r;
            if (h < minH) h = minH;
            if (h > maxH) h = maxH;

            // Compute angle coefficients for alpha(S) and beta(S).
            AngleCoeffs alphaCoeffs = ComputeAngleCoeffs(_alpha0, _alpha1, qAlpha);
            AngleCoeffs betaCoeffs = ComputeAngleCoeffs(_beta0, _beta1, qBeta);

            // Integrate local position up to s = h (S_end = 1).
            IntegrateLocal(alphaCoeffs, betaCoeffs, h, 1.0, IntegrationSteps,
                           out double xEnd, out double yEnd, out double zEnd);

            // Residual: chord constraints in local frame, scaled by chord length.
            residuals[0] = (xEnd - _r) * _invR;
            residuals[1] = yEnd * _invR;
            residuals[2] = zEnd * _invR;
        }

        /// <summary>
        /// Build a reasonable initial guess for p = [qAlpha, qBeta, h].
        /// </summary>
        public void BuildInitialGuess(double[] p)
        {
            if (p == null) throw new ArgumentNullException(nameof(p));
            if (p.Length != ParameterCount) throw new ArgumentException("Makino3DTripleClothoidProblem.BuildInitialGuess: length mismatch.");

            // Angle between endpoint tangents (in global space)
            double dotT = Helpers.Dot(_T0, _T1);
            dotT = Math.Clamp(dotT, -1.0, 1.0);
            double theta = Math.Acos(dotT);

            // Total length: chord plus some slack depending on bend.
            double h0 = _r * (1.0 + 0.5 * (theta / Math.PI));
            double minH = 0.5 * _r;
            double maxH = 4.0 * Math.Max(1e-3, _r);
            if (h0 < minH) h0 = minH;
            if (h0 > maxH) h0 = maxH;

            p[0] = 0.0; // qAlpha
            p[1] = 0.0; // qBeta
            p[2] = h0;  // h
        }

        #region Helper types / methods

        /// <summary>
        /// Simple struct holding the piecewise quadratic coefficients for α(S) or β(S).
        /// </summary>
        internal readonly struct AngleCoeffs
        {
            public readonly double a10, a11, a12;
            public readonly double a20, a21, a22;
            public readonly double a30, a31, a32;

            public AngleCoeffs(
                double a10, double a11, double a12,
                double a20, double a21, double a22,
                double a30, double a31, double a32)
            {
                this.a10 = a10; this.a11 = a11; this.a12 = a12;
                this.a20 = a20; this.a21 = a21; this.a22 = a22;
                this.a30 = a30; this.a31 = a31; this.a32 = a32;
            }

            public double Eval(double S)
            {
                if (S <= S1)
                {
                    return a10 + a11 * S + a12 * S * S;
                }
                else if (S <= S2)
                {
                    double u = S - S1;
                    return a20 + a21 * u + a22 * u * u;
                }
                else
                {
                    double u = S - S2;
                    return a30 + a31 * u + a32 * u * u;
                }
            }

            public double EvalPrimeWrtS(double S)
            {
                if (S <= S1)
                {
                    return a11 + 2.0 * a12 * S;
                }
                else if (S <= S2)
                {
                    double u = S - S1;
                    return a21 + 2.0 * a22 * u;
                }
                else
                {
                    double u = S - S2;
                    return a31 + 2.0 * a32 * u;
                }
            }

            public double EvalSecondPrimeWrtS(double S)
            {
                if (S <= S1)
                {
                    // a10 + a11 S + a12 S^2 -> derivative a11 + 2a12 S -> second derivative 2a12
                    return 2.0 * a12;
                }
                else if (S <= S2)
                {
                    // a20 + a21 u + a22 u^2, u = S - S1
                    return 2.0 * a22;
                }
                else
                {
                    return 2.0 * a32;
                }
            }
        }

        /// <summary>
        /// Compute the piecewise quadratic coefficients for an angle function
        /// θ(S) on [0,1] with:
        ///  - θ(0) = θ0
        ///  - θ(1) = θ1
        ///  - θ'(0) = 0, θ'(1) = 0 (zero curvature at endpoints)
        ///  - C¹ continuity at S1 and S2.
        /// 
        /// The only free "shape" parameter is q = θ''(0) = 2*a12.
        /// The rest of the coefficients are determined analytically, matching
        /// the derivation in the Makino triple-clothoid paper for S1=1/4, S2=3/4.
        /// </summary>
        private static AngleCoeffs ComputeAngleCoeffs(double theta0, double theta1, double q)
        {
            // Closed-form solution for S1 = 1/4 and S2 = 3/4:
            //  a10 = θ0
            //  a11 = 0
            //  a12 = q
            //
            //  a20 = θ0 + q / 16
            //  a21 = q / 2
            //  a22 = -q - (8/3) θ0 + (8/3) θ1
            //
            //  a30 = q / 16 + θ0 / 3 + 2 θ1 / 3
            //  a31 = -q / 2 - (8/3) θ0 + (8/3) θ1
            //  a32 =  q + (16/3) θ0 - (16/3) θ1

            double a10 = theta0;
            double a11 = 0.0;
            double a12 = q;

            double a20 = theta0 + q * (1.0 / 16.0);
            double a21 = q * 0.5;
            double a22 = -q - (8.0 / 3.0) * theta0 + (8.0 / 3.0) * theta1;

            double a30 = q * (1.0 / 16.0) + theta0 / 3.0 + 2.0 * theta1 / 3.0;
            double a31 = -q * 0.5 - (8.0 / 3.0) * theta0 + (8.0 / 3.0) * theta1;
            double a32 = q + (16.0 / 3.0) * theta0 - (16.0 / 3.0) * theta1;

            return new AngleCoeffs(a10, a11, a12, a20, a21, a22, a30, a31, a32);
        }

        /// <summary>
        /// Integrate the local position from s=0 to s=h*S_end using a midpoint rule
        /// over the unit interval in S = s/h. The tangent u(S) is given by the
        /// piecewise angle functions alpha(S) and beta(S).
        /// </summary>
        private static void IntegrateLocal(
            AngleCoeffs alpha,
            AngleCoeffs beta,
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

        private double[] ToLocal(double[] v)
        {
            return new double[3]
            {
                Helpers.Dot(v, _ex),
                Helpers.Dot(v, _ey),
                Helpers.Dot(v, _ez)
            };
        }

        /// <summary>
        /// Convert a local 3D vector u = (ux,uy,uz) to pitch/yaw angles (α,β)
        /// matching the paper's convention:
        ///   u = [cosα cosβ, cosα sinβ, -sinα]^T.
        /// </summary>
        private static (double alpha, double beta) ToAngles(double[] u)
        {
            double ux = u[0];
            double uy = u[1];
            double uz = u[2];

            double beta = Math.Atan2(uy, ux);
            double cosAlpha = Math.Sqrt(ux * ux + uy * uy);
            if (cosAlpha < 1e-12)
            {
                // Nearly vertical: treat as cosα≈0 and recover α from uz
                cosAlpha = 0.0;
            }

            double alpha = Math.Atan2(-uz, cosAlpha);
            return (alpha, beta);
        }

        private static double[] Cross(double[] a, double[] b)
        {
            return new double[3]
            {
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0]
            };
        }
        #endregion
    }

    /// <summary>
    /// Public fitter that wraps Makino3DTripleClothoidProblem and uses the
    /// small Levenberg-Marquardt solver with a Nelder–Mead fallback to find
    /// the best parameters p = [qAlpha, qBeta, h].
    /// 
    /// On success it constructs a Makino3DTripleClothoidCurve that implements
    /// IArcLengthCurve and can be sampled like any other ArcFrame curve.
    /// </summary>
    public static class Makino3DTripleClothoidFitter
    {
        public sealed class Result
        {
            public bool Succeeded;
            public double[] Parameters;
            public double FinalCost;
            public int LmIterations;
            public bool LmConverged;
            public bool NelderMeadUsed;
            public int NelderMeadIterations;
            public Makino3DTripleClothoidCurve? Curve;
        }

        /// <summary>
        /// Fit a 3D triple Makino clothoid between [P0, T0] and [P1, T1].
        /// </summary>
        public static Result Solve(double[] P0, double[] T0, double[] P1, double[] T1)
        {
            var problem = new Makino3DTripleClothoidProblem(P0, T0, P1, T1);

            double[] p0 = new double[problem.ParameterCount];
            problem.BuildInitialGuess(p0);

            // 1) LM phase
            var lmOpts = SmallLevenbergMarquardt.Options.Default;
            var lmRes = SmallLevenbergMarquardt.Solve(problem, p0, lmOpts);

            double[] bestP = lmRes.Parameters;
            double bestCost = lmRes.FinalCost;

            bool usedNM = false;
            int nmIter = 0;

            const double targetCost = 1e-6;

            // 2) Nelder–Mead fallback if needed
            if (!lmRes.Converged || bestCost > targetCost)
            {
                var nmOpts = SmallNelderMead.Options.Default;
                var nmRes = SmallNelderMead.Solve(problem, bestP, nmOpts);

                usedNM = true;
                nmIter = nmRes.Iterations;

                if (nmRes.Cost < bestCost)
                {
                    bestP = nmRes.Parameters;
                    bestCost = nmRes.Cost;
                }
            }

            bool succeeded = bestCost <= targetCost;

            Makino3DTripleClothoidCurve? curve = null; 
            curve = Makino3DTripleClothoidCurve.CreateFromFit(P0, T0, P1, T1, bestP);

            return new Result
            {
                Succeeded = succeeded,
                Parameters = bestP,
                FinalCost = bestCost,
                LmIterations = lmRes.Iterations,
                LmConverged = lmRes.Converged,
                NelderMeadUsed = usedNM,
                NelderMeadIterations = nmIter,
                Curve = curve
            };
        }
    }
}

