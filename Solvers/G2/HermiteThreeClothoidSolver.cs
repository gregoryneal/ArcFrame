using ArcFrame.Core.Constraints;
using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;
using ArcFrame.Core.Params;

namespace ArcFrame.Solvers.G2
{
    /// <summary>
    /// G²-Hermite interpolation in the plane with 3 clothoid segments.
    /// Inputs: (P0, θ0, k0)  and  (P1, θ1, k1).
    /// Output: three CurveSpecs that chain with C0/C1/C2 at the two joints.
    ///
    /// Method: numeric LM on 3 segments' (P0,R0,L, k0,dk) parameters with constraints:
    ///  - StartPose( seg0, P0, T(θ0) )   [hard]
    ///  - EndPose(   seg2, P1, T(θ1) )   [hard]
    ///  - CurvBoundary( seg0, k0, Start) [hard]
    ///  - CurvBoundary( seg2, k1, End  ) [hard]
    ///  - Joint C0/C1/C2 between (seg0,seg1) and (seg1,seg2)   [hard]
    /// Optionally add mild regularizers on slopes/lengths to pick a unique solution.
    /// </summary>
    public sealed class HermiteThreeClothoidSolver
    {
        public sealed record Result(CurveSpec[] Segments, CompositeCurve Composite, bool Converged, int Iterations);

        public Result Solve(
            double[] P0, double theta0, double k0,
            double[] P1, double theta1, double k1,
            double fairnessWeight = 0.0,  // set >0 to discourage large |dk| or overlong segments
            int maxIter = 80)
        {
            if (P0.Length != 2 || P1.Length != 2)
                throw new ArgumentException("Planar implementation: P0 and P1 must be length-2.");

            // --- 0) Seeds ----------------------------------------------------
            var (seg0, seg1, seg2) = BuildInitialSeeds(P0, theta0, k0, P1, theta1, k1);

            // --- 1) Problem & constraints -----------------------------------
            var prob = new CompositeCurveProblem(new[] { seg0, seg1, seg2 });

            // Start/end pose
            var T0 = new[] { Math.Cos(theta0), Math.Sin(theta0) };
            var T1 = new[] { Math.Cos(theta1), Math.Sin(theta1) };
            prob.Add(new StartPoseConstraint(0, P0, T0, ConstraintType.Hard));
            prob.Add(new EndPoseConstraint(2, P1, T1, ConstraintType.Hard));

            // Endpoint curvature values
            prob.Add(new CurvatureBoundaryConstraint(0, k0, CurvatureBoundaryConstraint.Where.Start, ConstraintType.Hard));
            prob.Add(new CurvatureBoundaryConstraint(2, k1, CurvatureBoundaryConstraint.Where.End, ConstraintType.Hard));

            // Joint G2 at both internal joints
            prob.Add(new CompositeCurveJointConstraint(0, 1, c0: true, c1: true, c2: true, type: ConstraintType.Hard));
            prob.Add(new CompositeCurveJointConstraint(1, 2, c0: true, c1: true, c2: true, type: ConstraintType.Hard));

            // Optional mild fairness: discourage big slopes |dk| and excessive lengths
            if (fairnessWeight > 0)
            {
                prob.Add(new SlopeRegularizer(0, fairnessWeight));
                prob.Add(new SlopeRegularizer(1, fairnessWeight));
                prob.Add(new SlopeRegularizer(2, fairnessWeight));
                prob.Add(new LengthRegularizer(0, fairnessWeight * 0.1));
                prob.Add(new LengthRegularizer(1, fairnessWeight * 0.1));
                prob.Add(new LengthRegularizer(2, fairnessWeight * 0.1));
            }

            // --- 2) Solve ----------------------------------------------------
            var solver = new CompositeCurveSolver
            {
                HardWeight = 1e7,
                EpsFD = 1e-6,
                MaxIter = maxIter,
                RelTol = 1e-8
            };
            var specs = solver.Solve(prob, optimizeP0: true, optimizeLength: true, optimizeR0: true);

            // --- 3) Rebuild a composite curve for convenience ----------------
            // TODO: Change this to add Clothoids with BF evaluator for speed.
            var comp = new CompositeCurve();
            comp.AddG1(new IntrinsicCurve(specs[0]), out _);
            comp.AddG1(new IntrinsicCurve(specs[1]), out _);
            comp.AddG1(new IntrinsicCurve(specs[2]), out _);

            return new Result(specs, comp, Converged: true, Iterations: maxIter);
        }

        // ------------------ seed builder ------------------
        private static (CurveSpec, CurveSpec, CurveSpec) BuildInitialSeeds(
            double[] P0, double theta0, double k0,
            double[] P1, double theta1, double k1)
        {
            int N = 2;
            var R0 = FrameFromTheta(theta0);
            var R2 = FrameFromTheta(theta1);
            var D = Helpers.Subtract(P1, P0);
            double Ltot = Helpers.Len(D);
            if (Ltot < 1e-9) Ltot = 1.0;

            // Split total length roughly 1/3 each (tweak with cosine of turn)
            double L0 = 0.33 * Ltot;
            double L1 = 0.34 * Ltot;
            double L2 = 0.33 * Ltot;

            // Heuristic slopes: linearly transition k0 -> k1 across whole chain
            double dkTotal = (k1 - k0) / Math.Max(1e-6, (L0 + L1 + L2));
            double dk0 = dkTotal * 1.0;
            double kJoin1 = k0 + dk0 * L0;
            double dk1 = dkTotal * 1.0;
            double kJoin2 = kJoin1 + dk1 * L1;
            double dk2 = (k1 - kJoin2) / Math.Max(1e-6, L2);

            // Build parameterizable curvature laws
            var law0 = new LinearCurvatureLawParamAdapter(new LinearCurvatureLaw(new[] { k0 }, new[] { dk0 }));
            var law1 = new LinearCurvatureLawParamAdapter(new LinearCurvatureLaw(new[] { kJoin1 }, new[] { dk1 }));
            var law2 = new LinearCurvatureLawParamAdapter(new LinearCurvatureLaw(new[] { kJoin2 }, new[] { dk2 }));

            // Place seed segment 0 exactly at P0, frame θ0
            var spec0 = new CurveSpec(N, L0, (double[])P0.Clone(), R0, law0, FrameModel.Bishop);

            // Predict intermediate frames/anchors by forward marching the previous ones:
            var E0 = new CachedIntrinsicCurve(spec0).Evaluate(L0);
            var spec1 = new CurveSpec(N, L1, E0.P, FrameFromT(E0.T), law1, FrameModel.Bishop);

            var E1 = new CachedIntrinsicCurve(spec1).Evaluate(L1);
            var spec2 = new CurveSpec(N, L2, E1.P, FrameFromTheta(theta1), law2, FrameModel.Bishop);

            return (spec0, spec1, spec2);
        }

        private static double[,] FrameFromTheta(double theta)
        {
            // 2D ONB with X along the tangent
            double c = Math.Cos(theta);
            double s = Math.Sin(theta);
            var R = new double[2, 2];
            R[0, 0] = c; R[0, 1] = -s;
            R[1, 0] = s; R[1, 1] = c;
            return R;
        }
        private static double[,] FrameFromT(double[] T)
        {
            double theta = Math.Atan2(T[1], T[0]);
            return FrameFromTheta(theta);
        }

        // --------- simple regularizers to pick a "nice" solution ----------
        private sealed class SlopeRegularizer : ICurveConstraint
        {
            public ConstraintType Type { get; } = ConstraintType.Soft;
            public double Weight { get; }
            public int SegmentIndex { get; }

            public SlopeRegularizer(int segIdx, double weight) { SegmentIndex = segIdx; Weight = weight; }

            public double[] Residual(CurveSpec spec)
            {
                // Penalize dk magnitude: weight * dk
                var p = (spec.Kappa as IParamCurvatureLaw)?.GetParams();
                // For LinearCurvatureLawParamAdapter: p = [k0, dk]
                double dk = (p != null && p.Length >= 2) ? p[1] : 0.0;
                return new[] { Weight * dk };
            }
        }

        private sealed class LengthRegularizer : ICompositeConstraint
        {
            public ConstraintType Type { get; } = ConstraintType.Soft;
            public double Weight { get; }
            private readonly int _segIdx;

            public LengthRegularizer(int segIdx, double weight) { _segIdx = segIdx; Weight = weight; }

            public double[] Residual(IReadOnlyList<CurveSpec> specs)
            {
                return new[] { Weight * specs[_segIdx].Length };
            }
        }
    }
}
