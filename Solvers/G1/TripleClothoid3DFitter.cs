using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using ArcFrame.Core.Results;
using ArcFrame.Solvers.Core;
using System;

namespace ArcFrame.Solvers.G1
{
    /*
    internal sealed class TripleClothoidProblem : ISmallNonlinearLeastSquaresProblem
    {
        public int ParameterCount => 4; // e.g. Ltot, a, kMid, tMid
        public int ResidualCount => 6; // 3 pos + 3 tangent

        private readonly double[] _p0;
        private readonly double[] _t0;
        private readonly double[] _p1;
        private readonly double[] _t1;

        public TripleClothoidProblem(double[] p0, double[] t0, double[] p1, double[] t1)
        {
            _p0 = p0;
            _t0 = t0;
            _p1 = p1;
            _t1 = t1;
        }

        public void Evaluate(double[] parameters, double[] residuals)
        {
            // BuildTripleSpecs(...) -> evaluate end pose -> fill residuals
            // residuals[0..2] = Pend - P1
            // residuals[3..5] = Tend - T1

            var c = BuildTripleSpecs(_p0, _t0, parameters);
            var lastSample = c.Evaluate(c.Length);
            var Pf = lastSample.P;
            var Tf = lastSample.T;

            // [dx, dy, dz, dtx, dty, dtz]
            for (int i = 0; i < 3; i++)
            {
                residuals[i] = Pf[i] - _p1[i];
                residuals[i + 3] = Tf[i] - _t1[i];
            }
        }

        /// <summary>
        /// Build the CompositeCurve of the triple clothoid from the parameter list
        /// We constrain the clothoid shape to be symmetrical to clamp down on DOF.
        /// Use the compositecurve for residuals
        /// </summary>
        /// <param name="P0"></param>
        /// <param name="T0"></param>
        /// <param name="p">Parameter vector: [Ltot, a, k_mid, tau_mid]</param>
        /// <returns></returns>
        public static CompositeCurve BuildTripleSpecs(double[] P0, double[] T0, double[] p)
        {
            // decode p
            var Ltot = p[0];
            // ensure a is less than 0.5
            // because if it is longer it causes the middle segment to be negative length
            // Ltot = L0 + L1 + L2
            // L0 = L2 = a * Ltot
            // L1 = (1 - 2a) * Ltot
            // so we need to ensure 2a < 1
            var a = Math.Clamp(p[1], 0, 0.49999);
            var kMid = p[2];
            var tMid = p[3];

            var L0 = a * Ltot;
            var L1 = (1 - (2 * a)) * Ltot;
            var L2 = a * Ltot;

            var k0 = 0.0;
            var k1 = 0.0;
            var t0 = 0.0;
            var t1 = 0.0;

            var dk0 = (kMid - k0) / L0;
            var dk1 = (k1 - kMid) / L1;
            var dk2 = 0.0; // (k1 - k1) / L2;

            var dt0 = (tMid - t0) / L0;
            var dt1 = (t1 - tMid) / L1;
            var dt2 = 0.0;

            var R0 = ONFrame.R0_FromT_Complete(T0);

            var seg0k = new double[] { k0, t0 };
            var seg0dk = new double[] { dk0, dt0 };
            var seg1k = new double[] { kMid, tMid };
            var seg1dk = new double[] { dk1, dt1 };
            var seg2k = new double[] { k1, t1 };
            var seg2dk = new double[] { dk2, dt2 };
            var law0 = new LinearCurvatureLaw(seg0k, seg0dk);
            var law1 = new LinearCurvatureLaw(seg1k, seg1dk);
            var law2 = new LinearCurvatureLaw(seg2k, seg2dk);
            // enforce G1 via construction
            CompositeCurve c = new CompositeCurve();
            CurveSpec spec0 = new CurveSpec(3, L0, P0, R0, law0, FrameModel.Frenet);
            // these two in local space and add all to composite curve
            var R00 = RigidTransform.Identity(3).R;
            CurveSpec spec1 = new CurveSpec(3, L1, new double[3], R00, law0, FrameModel.Frenet);
            CurveSpec spec2 = new CurveSpec(3, L2, new double[3], (double[,])R00.Clone(), law0, FrameModel.Frenet);

            spec0.ShowInfo();
            spec1.ShowInfo();
            spec2.ShowInfo();


            c.Add(new CachedIntrinsicCurve(spec0))
            .AddG1FullFrame(new CachedIntrinsicCurve(spec1), out _)
            .AddG1FullFrame(new CachedIntrinsicCurve(spec2), out _);

            return c;
        }
    }


    /// <summary>
    /// 3D G1 triple clothoid segment fitter between [P0, T0] and [P1, T1].
    /// Uses a symmetric 3d clothoid to constrain DoF for solver. 
    /// </summary>
    public class TripleClothoid3DFitter
    {
        /// <summary>
        /// Solve a TripleClothoidProblem and build a 3d composite clothoid with desired endpoint position and tangents.
        /// </summary>
        /// <param name="P0"></param>
        /// <param name="T0"></param>
        /// <param name="P1"></param>
        /// <param name="T1"></param>
        /// <param name="normalHint"></param>
        /// <returns></returns>
        public static Result Solve(double[] P0, double[] T0, double[] P1, double[] T1, double[]? normalHint = null)
        {
            // Build seed curve
            var eps = 1E-4;
            var chord = Helpers.Subtract(P1, P0);
            var D = Helpers.Len(chord);
            if (D < eps)
            {
                // degen case
                return new Result(true, 0, 0, 0, new double[] { 0, 0, 0, 0 });
            }
            var Ltot = 1.5 * D;
            // scale the outer segments by total len
            var a = 0.25;

            // curvature and torsion midpoints
            var u = Helpers.Cross3(T0, T1);
            double kMid;
            double tMid = 0;
            if (Helpers.Len2(u) < eps) 
            {
                // nearly parallel 
                kMid = 0;
            }
            else
            {
                double[,] frame = ONFrame.R0_FromT_Complete(u, normalHint);
                var B0 = ONFrame.GetCol(frame, 2);
                var sign = Math.Sign(Helpers.Dot(u, B0));
                kMid = sign * (0.1 / D);
            }

            // Pack [Ltot, a, kMid, tMid]
            double[] p0 = new double[] { Ltot, a, kMid, tMid };

            var problem = new TripleClothoidProblem(P0, T0, P1, T1);
            var opts = SmallLevenbergMarquardt.Options.Default;
            var result = SmallLevenbergMarquardt.Solve(problem, p0, opts);

            // build the final curve and attach to problem
            result.Curve = TripleClothoidProblem.BuildTripleSpecs(P0, T0, result.Parameters);            

            return result;
        }
    }*/

    /// <summary>
    /// G1 3D triple-clothoid fitter between two poses [P0, T0] and [P1, T1].
    /// 
    /// Model:
    ///  - 3 segments, Frenet frame.
    ///  - Each segment has linear curvature and torsion (kappa, tau).
    ///  - kappa(s) and tau(s) are piecewise linear over the full length with nodes at
    ///    normalized arc-length S = 0, 0.25, 0.75, 1.0.
    ///  - At each node we specify kappa and tau (except start where we take kappa0=tau0=0).
    /// 
    /// Parameters p (length 7):
    ///  p[0] = Ltot  (total length of triple segment)
    ///  p[1] = kappa(S = 0.25)
    ///  p[2] = kappa(S = 0.75)
    ///  p[3] = kappa(S = 1.0)
    ///  p[4] = tau(S = 0.25)
    ///  p[5] = tau(S = 0.75)
    ///  p[6] = tau(S = 1.0)
    /// 
    /// So there are 6 shape DOFs (as in the triple-clothoid papers) plus the
    /// total length Ltot. This is the 3D extension with linear curvature and
    /// torsion segments.
    /// 
    /// The solver:
    ///  - Runs a small LM on these 7 parameters.
    ///  - If LM doesn't converge or cost is still too large, it runs a small
    ///    Nelder–Mead starting from the LM result.
    /// </summary>
    public static class TripleClothoid3DFitter
    {
        public sealed class Result
        {
            public bool Solved;
            public string Message;
            public CompositeCurve Curve;
            public double[] Parameters;
            public double FinalCost;
            public int LmIterations;
            public bool LmConverged;
            public bool NelderMeadUsed;
            public int NelderMeadIterations;
        }

        /// <summary>
        /// Solve for a 3D triple-clothoid between [P0,T0] and [P1,T1].
        /// All vectors must be length-3 and T0, T1 must be unit or close.
        /// </summary>
        public static Result Solve(double[] P0, double[] T0, double[] P1, double[] T1)
        {
            if (P0 == null || T0 == null || P1 == null || T1 == null)
                throw new ArgumentNullException("P0, T0, P1, T1 must be non-null.");
            if (P0.Length != 3 || T0.Length != 3 || P1.Length != 3 || T1.Length != 3)
                throw new ArgumentException("TripleClothoid3DFitter.Solve expects 3D vectors.");

            var problem = new TripleClothoid3DProblem(P0, T0, P1, T1);

            int n = problem.ParameterCount;
            double[] p0 = new double[n];
            problem.BuildInitialGuess(p0);

            // --- 1) LM phase ---
            var lmOpts = SmallLevenbergMarquardt.Options.Default;
            var lmRes = SmallLevenbergMarquardt.Solve(problem, p0, lmOpts);

            double[] bestP = lmRes.Parameters;
            double bestCost = lmRes.FinalCost;
            bool usedNelder = false;
            int nmIter = 0;
            bool solved = lmRes.Converged;
            const double targetCost = 1e-6;

            // --- 2) Nelder–Mead fallback if needed ---
            if (!lmRes.Converged || bestCost > targetCost)
            {
                var nmOpts = SmallNelderMead.Options.Default;
                var nmRes = SmallNelderMead.Solve(problem, bestP, nmOpts);

                usedNelder = true;
                nmIter = nmRes.Iterations;

                if (nmRes.Cost < bestCost)
                {
                    bestP = nmRes.Parameters;
                    bestCost = nmRes.Cost;
                    solved = nmRes.Converged || bestCost <= targetCost;
                }
            }

            // --- Build final segments from best parameters ---
            var curve = TripleClothoid3DProblem.BuildSegments(P0, T0, bestP);

            string msg = solved
                ? $"TripleClothoid3DFitter: OK, cost={bestCost:F6}, LM iters={lmRes.Iterations}, NM used={usedNelder}"
                : $"TripleClothoid3DFitter: NOT CONVERGED, cost={bestCost:F6}, LM iters={lmRes.Iterations}, NM used={usedNelder}";

            return new Result
            {
                Solved = solved,
                Message = msg,
                Curve = curve,
                Parameters = bestP,
                FinalCost = bestCost,
                LmIterations = lmRes.Iterations,
                LmConverged = lmRes.Converged,
                NelderMeadUsed = usedNelder,
                NelderMeadIterations = nmIter
            };
        }

        /// <summary>
        /// Internal problem wrapper: maps parameter vector p to residuals
        /// (position + tangent error).
        /// </summary>
        internal sealed class TripleClothoid3DProblem : ISmallNonlinearLeastSquaresProblem
        {
            public int ParameterCount => 7; // Ltot + 6 shape params
            public int ResidualCount => 6; // 3 pos + 3 tangent

            private readonly double[] _P0;
            private readonly double[] _T0;
            private readonly double[] _P1;
            private readonly double[] _T1;
            private readonly double[,] _R0;

            private readonly double _chordLen;
            private readonly double _invChordLen;

            private static readonly IFrameStepper Stepper = new LieGroupMidpointStepper();
            private static readonly IntegratorOptions Opt = IntegratorOptions.Default;

            private const double S1 = 0.25;
            private const double S2 = 0.75;

            public TripleClothoid3DProblem(double[] P0, double[] T0, double[] P1, double[] T1)
            {
                _P0 = (double[])P0.Clone();
                _T0 = Helpers.Normalize((double[])T0.Clone());
                _P1 = (double[])P1.Clone();
                _T1 = Helpers.Normalize((double[])T1.Clone());

                _R0 = ONFrame.R0_FromT_Complete(_T0);

                var chord = Helpers.Subtract(_P1, _P0);
                _chordLen = Helpers.Len(chord);
                _invChordLen = 1.0 / Math.Max(1e-6, _chordLen);
            }

            public void Evaluate(double[] parameters, double[] residuals)
            {
                // Compute end pose from current parameters
                double[] Pend = new double[3];
                double[] Tend = new double[3];
                ComputeEndPose(parameters, Pend, Tend);

                // Residual: normalized position error + tangent difference
                for (int i = 0; i < 3; ++i)
                {
                    residuals[i] = (Pend[i] - _P1[i]) * _invChordLen;
                    residuals[3 + i] = Tend[i] - _T1[i];
                }
            }

            /// <summary>
            /// Build a reasonable initial guess for p
            /// based on chord length and the angle between T0 and T1.
            /// </summary>
            public void BuildInitialGuess(double[] p)
            {
                if (p == null) throw new ArgumentNullException(nameof(p));
                if (p.Length != ParameterCount) throw new ArgumentException("Initial guess array has wrong length.");

                // chord direction and angle between tangents
                var chord = Helpers.Subtract(_P1, _P0);
                double Lc = _chordLen;
                double[] tChord = Lc > 1e-8 ? Helpers.Multiply(1.0 / Lc, chord) : new double[] { _T0[0], _T0[1], _T0[2] };

                double dotT = Helpers.Dot(_T0, _T1);
                if (dotT < -1.0) dotT = -1.0;
                if (dotT > 1.0) dotT = 1.0;
                double theta = Math.Acos(dotT); // angle between tangents

                // Total length: chord length with a little slack depending on bend
                double Ltot = Math.Max(1e-3, Lc * (1.0 + 0.5 * (theta / Math.PI)));
                double minL = 0.5 * Lc;
                double maxL = 5.0 * Math.Max(1e-3, Lc);
                if (Ltot < minL) Ltot = minL;
                if (Ltot > maxL) Ltot = maxL;

                // Curvature magnitude to roughly produce the turning angle: theta ≈ k_avg * Ltot
                double kBase = (Ltot > 1e-6) ? theta / Ltot : 0.0;

                // Choose the sign of curvature based on cross(T0, T1)·B0
                double[] cross = Cross(_T0, _T1);
                double[] B0 = new double[3];
                for (int i = 0; i < 3; ++i) B0[i] = _R0[i, 2]; // binormal = column 2
                double signK = 1.0;
                double crossLen2 = Helpers.Len2(cross);
                if (crossLen2 > 1e-12)
                {
                    double dot = Helpers.Dot(cross, B0);
                    signK = Math.Sign(dot);
                    if (signK == 0.0) signK = 1.0;
                }

                double kAmp = signK * kBase;

                // Simple shape: kappa ramps up then down.
                // We keep endpoints (S=0) zero curvature and allow kappa(1.0) small.
                double k1 = 0.5 * kAmp;
                double k2 = 1.0 * kAmp;
                double k3 = 0.25 * kAmp;

                // Torsion initial guess: 0 (planar); solver can introduce torsion if needed.
                double t1 = 0.0;
                double t2 = 0.0;
                double t3 = 0.0;

                p[0] = Ltot;
                p[1] = k1;
                p[2] = k2;
                p[3] = k3;
                p[4] = t1;
                p[5] = t2;
                p[6] = t3;
            }

            /// <summary>
            /// Build the three CurveSpecs for a given parameter vector,
            /// starting from [P0, T0].
            /// </summary>
            public static CompositeCurve BuildSegments(double[] P0, double[] T0, double[] p)
            {
                // We construct a temporary problem just to reuse the pose-building logic.
                var tmp = new TripleClothoid3DProblem(P0, T0, P0, T0); // P1,T1 dummy
                return tmp.BuildSegmentsInternal(p);
            }

            private CompositeCurve BuildSegmentsInternal(double[] parameters)
            {
                // clamp total length to reasonable range
                double Lc = _chordLen;
                double minL = 0.5 * Lc;
                double maxL = 5.0 * Math.Max(1e-3, Lc);

                double Ltot = parameters[0];
                if (double.IsNaN(Ltot) || double.IsInfinity(Ltot)) Ltot = Lc;
                if (Ltot < minL) Ltot = minL;
                if (Ltot > maxL) Ltot = maxL;

                // segment lengths
                double L0 = S1 * Ltot;
                double L1 = (S2 - S1) * Ltot;
                double L2 = (1.0 - S2) * Ltot;

                if (L0 <= 0) L0 = 1e-6;
                if (L1 <= 0) L1 = 1e-6;
                if (L2 <= 0) L2 = 1e-6;

                double k0 = 0.0;
                double k1 = parameters[1];
                double k2 = parameters[2];
                double k3 = parameters[3];

                double t0 = 0.0;
                double t1 = parameters[4];
                double t2 = parameters[5];
                double t3 = parameters[6];

                // Segment 0: s in [0,L0], kappa, tau linear from (k0,t0) to (k1,t1)
                double[] k0Vec0 = new double[2] { k0, t0 };
                double[] dkVec0 = new double[2]
                {
                    (k1 - k0) / L0,
                    (t1 - t0) / L0
                };

                // Segment 1: s in [0,L1], local start at (k1,t1) to (k2,t2)
                double[] k0Vec1 = new double[2] { k1, t1 };
                double[] dkVec1 = new double[2]
                {
                    (k2 - k1) / L1,
                    (t2 - t1) / L1
                };

                // Segment 2: s in [0,L2], local start at (k2,t2) to (k3,t3)
                double[] k0Vec2 = new double[2] { k2, t2 };
                double[] dkVec2 = new double[2]
                {
                    (k3 - k2) / L2,
                    (t3 - t2) / L2
                };

                var law0 = new LinearCurvatureLaw(k0Vec0, dkVec0);
                var law1 = new LinearCurvatureLaw(k0Vec1, dkVec1);
                var law2 = new LinearCurvatureLaw(k0Vec2, dkVec2);

                // Build three specs in sequence
                //var specs = new CurveSpec[3];

                double[] Pstart = (double[])_P0.Clone();
                double[,] Rstart = (double[,])_R0.Clone();

                // Segment 0
                var spec0 = new CurveSpec(3, L0, Pstart, Rstart, law0, FrameModel.Frenet);
                //var c0 = new IntrinsicCurve(seg0, Stepper, Opt);

                //var s0 = c0.Evaluate(L0);
                //specs[0] = seg0;


                // these two in local space and add all to composite curve
                var R00 = RigidTransform.Identity(3).R;
                CurveSpec spec1 = new CurveSpec(3, L1, new double[3], R00, law0, FrameModel.Frenet);
                CurveSpec spec2 = new CurveSpec(3, L2, new double[3], (double[,])R00.Clone(), law0, FrameModel.Frenet);
                //-------

                // Segment 1
                /*
                var seg1 = new CurveSpec(3, L1, s0.P, s0.R, law1, FrameModel.Frenet);
                var c1 = new IntrinsicCurve(seg1, Stepper, Opt);
                var s1 = c1.Evaluate(L1);
                specs[1] = seg1;

                // Segment 2
                var seg2 = new CurveSpec(3, L2, s1.P, s1.R, law2, FrameModel.Frenet);
                specs[2] = seg2;
                */

                CompositeCurve c = new CompositeCurve();
                c.Add(new CachedIntrinsicCurve(spec0))
                .AddG1FullFrame(new CachedIntrinsicCurve(spec1), out _)
                .AddG1FullFrame(new CachedIntrinsicCurve(spec2), out _);


                return c;
            }

            private void ComputeEndPose(double[] parameters, double[] Pout, double[] Tout)
            {
                //var specs = BuildSegmentsInternal(parameters);

                // Integrate the last segment to its end to get final pose
                //var seg2 = specs[2];
                //var c2 = new IntrinsicCurve(seg2, Stepper, Opt);
                //var s2 = c2.Evaluate(seg2.Length);

                var curve = BuildSegmentsInternal(parameters);
                var s2 = curve.Evaluate(curve.Length);

                // Copy position
                var P = s2.P;
                Pout[0] = P[0];
                Pout[1] = P[1];
                Pout[2] = P[2];

                var T = s2.T;
                Tout[0] = T[0];
                Tout[1] = T[1];
                Tout[2] = T[2];
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
        }
    }
}
