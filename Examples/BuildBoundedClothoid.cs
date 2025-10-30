using ArcFrame.Core.Constraints;
using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;
using ArcFrame.Solvers;
using ArcFrame.Solvers.G2;
using System;

namespace ArcFrame.Examples
{
    /// <summary>
    /// Build a clothoid curve bounded inside a track defined by two IArcLengthCurves.
    /// 
    /// 
    ///                                                 |            | 
    ///                                                 |            | 
    ///                                                 |     |      | 
    ///                                                 |     |      | 
    ///                                         line    |            | line
    ///                                                 |            | 
    ///                                                 |     |      | 
    ///                                                 |     |      | 
    ///                                                 |            | 
    ///                                                 |            | 
    ///                                                 |     |      | 
    ///                                                 |     |      |
    ///                                           arc  /             |
    /// curve1 =    line                              /              |
    /// ----------------------------------------------      /       /
    ///                                                    /       /
    ///                                                           /
    ///  ---      ---      ---      ---      ---      ---        /  arc
    ///                                                         /
    ///                                                        /
    /// -------------------------------------------------------
    /// curve2 =    line
    /// 
    /// 
    /// 
    /// 
    /// 
    /// </summary>
    public class BuildBoundedClothoid
    {
        public static (CompositeCurve[] trialSolutions, CompositeCurve finalCurve, CompositeCurve leftBound, CompositeCurve rightBound) Main()
        {
            double laneWidth = 2;
            double laneStartX = -10;
            CompositeCurve upperLaneCurve = new CompositeCurve();
            CompositeCurve lowerLaneCurve = new CompositeCurve();

            // Start, tangents, specs of curve 1, the upper lane bound:
            double[] P01 = new double[] { laneStartX, laneWidth / 2 };
            double[] T01 = new double[] { 1, 0 };
            double length1 = -laneStartX - (laneWidth / 2) - 0.5;
            Line line10 = new Line(P01, Helpers.Add(P01, Helpers.Multiply(length1, T01)));
            Line line11 = new Line(P01, Helpers.Add(P01, Helpers.Multiply(length1, T01)));
            Arc arc10 = Arc.From2D(0, 0, 0.5, 0, Math.PI / 2);
            upperLaneCurve.Add(line10).AddG1(arc10, out _).AddG1(line11, out _);

            // Start, tangents, specs of curve 2, the lower lane bound:
            double[] P02 = new double[] { laneStartX, -laneWidth / 2 };
            double[] T02 = new double[] { 1, 0 };
            double length2 = -laneStartX;
            Line line20 = new Line(P02, Helpers.Add(P02, Helpers.Multiply(length2, T02)));
            Line line21 = new Line(P02, Helpers.Add(P02, Helpers.Multiply(length2, T02)));
            Arc arc20 = Arc.From2D(0, 0, laneWidth / 2, 0, Math.PI / 2);
            lowerLaneCurve.Add(line20).AddG1(arc20, out _).AddG1(line21, out _);

            // Initial frame of the clothoid
            double[] pC0 = new double[] { 5 * laneStartX / 6, 0 }; // just inside the start of the road
            double tC0 = 0; //radians
            double kC0 = 0;
            // Final frame of the clothoid
            double[] pC1 = new double[] { 0, -5 * laneStartX / 6 };
            double tC1 = Math.PI / 2;
            double kC1 = 0;

            // Create the seed curve and problem.
            var seedCurve = new HermiteThreeClothoidSolver();
            HermiteThreeClothoidSolver.Result curveResult = seedCurve.Solve(pC0, tC0, kC0, pC1, tC1, kC1, 0, 800);
            var problem = new CompositeCurveProblem(curveResult.Segments);

            // Now we need to add constraints to the problem:
            // Start / End pose
            problem.Add(new StartPoseConstraint(0, pC0, T01));
            problem.Add(new EndPoseConstraint(problem.Seeds.Length - 1, pC1, T02));
            // Start / End curvature at boundary
            problem.Add(new CurvatureBoundaryConstraint(0, kC0, CurvatureBoundaryConstraint.Where.Start));
            problem.Add(new CurvatureBoundaryConstraint(problem.Seeds.Length - 1, kC1, CurvatureBoundaryConstraint.Where.End));
            // First reimpose the G2 boundaries across all segments
            for (int i = 0; i + 1 < problem.Seeds.Length; i++)
            {
                problem.Add(new CompositeCurveJointConstraint(i, i + 1, true, true, true));
            }

            // Add the bounded lane constraint
            problem.Add(new BoundedCurveConstraint(upperLaneCurve, lowerLaneCurve, 0.15, 4, 0.5, 512, 20));

            // Penalize hard steering
            problem.Add(new SlopeJumpRegularizer(5E-3));
            //problem.Add(new LocalDkPenalty(0, 1E-2));
            //problem.Add(new LocalDkPenalty(1, 1E-2));
            //problem.Add(new LocalDkPenalty(2, 1E-2));

            // Setup the solver with non default parameters
            var solver = new CompositeCurveSolver() { HardWeight = 1E-7, EpsFD = 1E-6, MaxIter = 80, RelTol = 1E-8 };
            var solverResult = solver.Solve(problem, true, true, true);
            var solvedSpecs = solverResult.FinalSolution;

            // Draw the final curve
            CompositeCurve constrainedCurve = new CompositeCurve();
            for (int i = 0; i < solvedSpecs.Length; i++)
            {
                Console.WriteLine($"Adding curve to final solution: {i} : {solvedSpecs.Length - 1}");
                constrainedCurve.AddG1(new CachedIntrinsicCurve(solvedSpecs[i]), out _);
            }

            CompositeCurve[] trialCurves = new CompositeCurve[solverResult.TrialSolutions.Length];
            for (int i = 0; i < solverResult.TrialSolutions.Length; i++)
            {
                CompositeCurve trialCurve = new CompositeCurve();
                for (int j = 0; j < solverResult.TrialSolutions[i].Length; j++)
                {
                    trialCurve.AddG1(new CachedIntrinsicCurve(solverResult.TrialSolutions[i][j]), out _);
                }
                trialCurves[i] = trialCurve;
            }

            return (trialCurves, constrainedCurve, upperLaneCurve, lowerLaneCurve);
        }
    }
}
