using ArcFrame.Core.Constraints;
using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using System;

namespace ArcFrame.Solvers.G1
{
    /// <summary>
    /// Solve for a 3D composite clothoid with three arcs, with start pose [P0, T0] and end pose [P1, T1].
    /// </summary>
    public static class Hermite3DThreeClothoidSolver
    {
        /// <summary>
        /// Idea: Fallback: 3 segments in 3D (G1 joints) – more DOF, very robust
        /// </summary>
        public static CompositeCurveSolverResult Solve(double[] P0, double[] T0, double[] P1, double[] T1)
        {
            // Initial seed curve is a straight line from P0 to P1
            var chord = Helpers.Subtract(P1, P0);
            var t0 = Helpers.Normalize(chord);
            var R0 = ONFrame.R0_FromT_Complete(t0);
            double Ltot = Math.Max(1.0, Helpers.Len(chord));
            double L0 = 0.33 * Ltot, L1 = 0.34 * Ltot, L2 = 0.33 * Ltot;

            // Seed three vector-curvature laws (order=2) with small slopes
            var law0 = new LinearCurvatureLaw(new[] { 0.0, 0.0 }, new[] { 0.01, 0.01 });
            var spec0 = new CurveSpec(3, L0, (double[])P0.Clone(), R0, law0, FrameModel.Frenet);

            var E0 = new CachedIntrinsicCurve(spec0).Evaluate(L0);
            var law1 = new LinearCurvatureLaw(new[] { 0.01, 0.0 }, new[] { 0.01, 0.01 });
            var spec1 = new CurveSpec(3, L1, E0.P, ONFrame.R0_FromT_Complete(E0.T), law1, FrameModel.Frenet);

            var E1 = new CachedIntrinsicCurve(spec1).Evaluate(L1);
            var law2 = new LinearCurvatureLaw(new[] { 0.0, 0.01 }, new[] { 0.01, 0.01 });
            var spec2 = new CurveSpec(3, L2, E1.P, ONFrame.R0_FromT_Complete(T1), law2, FrameModel.Frenet);

            var prob = new CompositeCurveProblem(new[] { spec0, spec1, spec2 });
            // Hard start/end pose
            prob.Add(new StartPoseConstraint(0, P0, T0, ConstraintType.Hard));
            prob.Add(new EndPoseConstraint(2, P1, T1, ConstraintType.Hard));
            // G1 joints: C0 & tangent match at internal junctions
            prob.Add(new CompositeCurveJointConstraint(0, 1, c0: true, c1: true, c2: false, type: ConstraintType.Hard));
            prob.Add(new CompositeCurveJointConstraint(1, 2, c0: true, c1: true, c2: false, type: ConstraintType.Hard));

            // Mild fairness to avoid wild curvature
            prob.Add(new CurvatureMagnitudeRegularizer(0, 5e-4));
            prob.Add(new CurvatureMagnitudeRegularizer(1, 5e-4));
            prob.Add(new CurvatureMagnitudeRegularizer(2, 5e-4));

            var solver = new CompositeCurveSolver { HardWeight = 1e7, EpsFD = 1e-6, MaxIter = 100, RelTol = 1e-8 };
            var result = solver.Solve(prob,
                                     optimizeP0: false,   // let internal anchors settle
                                     optimizeLength: true, // try to get shorter curve length
                                     optimizeR0: true);  // OK here; joints enforce G1

            var specs = result.FinalSolution;

            var comp = new CompositeCurve();
            comp.AddG1(specs[0].GetOptimizedCurve(), out _);
            comp.AddG1(specs[1].GetOptimizedCurve(), out _);
            comp.AddG1(specs[2].GetOptimizedCurve(), out _);
            return result;
        }
    }
}
