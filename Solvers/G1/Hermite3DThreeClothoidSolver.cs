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
            var R0 = ONFrame.R0_FromT_Complete(T0);
            double Ltot = Math.Max(1.0, Helpers.Len(chord));
            double L0 = 0.33 * Ltot, L1 = 0.34 * Ltot, L2 = 0.33 * Ltot;

            // actual tolerance scaled by the chord length
            // ensures when we draw our resulting line we
            // can control how close we are to P1 since it
            // is constructed. 
            double alpha = 1E-5; // max error at Ltot = 1
            double eps = 1 - (alpha / System.Math.Max(1, Ltot));
            if (Helpers.Dot(t0, T0) > eps && Helpers.Dot(t0, T1) > eps)
            {
                var law = new ConstantCurvatureLaw(new double[] { 0, 0 });
                var sp1 = new double[3];
                var sp2 = new double[3];
                var s1 = new CurveSpec(3, L0, P0, R0, law, FrameModel.Frenet);
                var s2 = new CurveSpec(3, L1, sp1, RigidTransform.Identity(3).R, law.CloneWithParams(law.GetParams()), FrameModel.Frenet);
                var s3 = new CurveSpec(3, L2, sp2, RigidTransform.Identity(3).R, law.CloneWithParams(law.GetParams()), FrameModel.Frenet);
                CurveSpec[][] trials = new CurveSpec[][] { new CurveSpec[] { s1, s2, s3 }, };
                return new CompositeCurveSolverResult(trials, trials[0], true, "Degenerate Line Result", 0, 0);
            }
            /*
            var law1 = new LinearCurvatureLaw(new double[] { 0, 0 }, new double[] { 0, 0 });
            var sp11 = new double[3];
            var sp21 = new double[3];
            var s11 = new CurveSpec(3, L0, P0, R0, law1, FrameModel.Frenet);
            var s21 = new CurveSpec(3, L1, sp11, RigidTransform.Identity(3).R, law1.CloneWithParams(law1.GetParams()), FrameModel.Frenet);
            var s31 = new CurveSpec(3, L2, sp21, RigidTransform.Identity(3).R, law1.CloneWithParams(law1.GetParams()), FrameModel.Frenet);*/

            
            // Seed three vector-curvature laws (order=2) with small slopes
            // first one is posed, second two are in local space
            var law0 = new LinearCurvatureLaw(new[] { 0.0, 0.0 }, new[] { 0.01, 0.01 });
            var spec0 = new CurveSpec(3, L0, (double[])P0.Clone(), R0, law0, FrameModel.Frenet);

            var R1 = RigidTransform.Identity(3).R;

            var p1 = new double[P0.Length];
            var law1 = new LinearCurvatureLaw(new[] { 0.01, 0.0 }, new[] { 0.01, 0.01 });
            var spec1 = new CurveSpec(3, L1, p1, (double[,])R1.Clone(), law1, FrameModel.Frenet);

            var p2 = new double[P0.Length];
            var law2 = new LinearCurvatureLaw(new[] { 0.0, 0.01 }, new[] { 0.01, 0.01 });
            var spec2 = new CurveSpec(3, L2, p2, (double[,])R1.Clone(), law2, FrameModel.Frenet);

            var prob = new CompositeCurveProblem(new[] { spec0, spec1, spec2 });
            // Hard start/end pose
            // leave out start pose when optimize P0 and R0 is off, they will stay fixed.
            //prob.Add(new StartPoseConstraint(0, P0, T0, ConstraintType.Hard));
            prob.Add(new EndPoseConstraint(2, P1, T1, ConstraintType.Hard));
            // G1 joints: C0 & tangent match at internal junctions
            // don't use these for chained solver mode
            //prob.Add(new CompositeCurveJointConstraint(0, 1, c0: true, c1: true, c2: false, type: ConstraintType.Hard));
            //prob.Add(new CompositeCurveJointConstraint(1, 2, c0: true, c1: true, c2: false, type: ConstraintType.Hard));

            // Mild fairness to avoid wild curvature
            prob.Add(new CurvatureMagnitudeRegularizer(0, 1e-2));
            prob.Add(new CurvatureMagnitudeRegularizer(1, 1e-2));
            prob.Add(new CurvatureMagnitudeRegularizer(2, 1e-2));

            var solver = new CompositeCurveSolver { HardWeight = 1e4, EpsFD = 1e-5, MaxIter = 500, RelTol = 1e-8, NumSamples = 20, Mode = CompositeCurveSolver.CompositeSolverMode.Chained };
            var result = solver.Solve(prob,
                                     optimizeP0: false,   // don't optimize
                                     optimizeLength: true, // try to get shorter curve length
                                     optimizeR0: false);  // OK here; joints enforce G1

            //var specs = result.FinalSolution;

            //var comp = new CompositeCurve();
            //comp.AddG1FullFrame(specs[0].GetOptimizedCurve(), out _);
            //comp.AddG1FullFrame(specs[1].GetOptimizedCurve(), out _);
            //comp.AddG1FullFrame(specs[2].GetOptimizedCurve(), out _);
            return result;
        }
    }
}
