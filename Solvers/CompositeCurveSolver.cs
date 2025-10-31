using ArcFrame.Core.Constraints;
using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using System.Collections.Generic;
using System;
using System.IO.Compression;
using ArcFrame.Solvers.Core;

namespace ArcFrame.Solvers
{
    /// <summary>
    /// Describe a constrained composite curve problem to be solved.
    /// </summary>
    public sealed class CompositeCurveProblem
    {
        /// <summary>
        /// The seed specs, these will be altered to fit the constraints
        /// </summary>
        public CurveSpec[] Seeds { get; }
        /// <summary>
        /// Constraints on the problem to be solved.
        /// </summary>
        public List<ICompositeConstraint> Constraints { get; }
        /// <summary>
        /// Create a new composite curve problem with input seeds.
        /// </summary>
        /// <param name="seeds"></param>
        public CompositeCurveProblem(CurveSpec[] seeds) 
        {
            Constraints = new List<ICompositeConstraint>();
            Seeds = seeds; 
        }

        /// <summary>
        /// Wrap a single curve constraint in a SegmentConstraintWrapper for ICompositeConstraint
        /// </summary>
        /// <param name="single"></param>
        public void Add(ICurveConstraint single) => Constraints.Add(new SegmentConstraintWrapper(single));
        /// <summary>
        /// Add a single composite curve constraint.
        /// </summary>
        /// <param name="joint"></param>
        public void Add(ICompositeConstraint joint) => Constraints.Add(joint);
    }

    /// <summary>
    /// Class to hold information about the solver results.
    /// </summary>
    public sealed class CompositeCurveSolverResult
    {
        /// <summary>
        /// All of the trial solutions on the way to finding the real solution.
        /// </summary>
        public CurveSpec[][] TrialSolutions { get; set; }
        /// <summary>
        /// The final curve spec solution list.
        /// </summary>
        public CurveSpec[] FinalSolution { get; set; }
        /// <summary>
        /// Was the curve solved.
        /// </summary>
        public bool Solved { get; set; }
        /// <summary>
        /// Solver message.
        /// </summary>
        public string Message { get; set; }
        /// <summary>
        /// Number of iterations before it was solved.
        /// </summary>
        public int Iterations { get; set; }
        /// <summary>
        /// The final cost of the solution.
        /// </summary>
        public double FinalCost { get; set; }

        /// <summary>
        /// Create a composte curve solver result.
        /// </summary>
        /// <param name="trials"></param>
        /// <param name="result"></param>
        /// <param name="solved"></param>
        /// <param name="message"></param>
        /// <param name="iterations"></param>
        /// <param name="finalCost"></param>
        public CompositeCurveSolverResult(CurveSpec[][] trials, CurveSpec[] result, bool solved, string message, int iterations, double finalCost)
        {
            TrialSolutions = trials;
            FinalSolution = result;
            Solved = solved;
            Message = message;
            Iterations = iterations;
            FinalCost = finalCost;
        }
    }

    /// <summary>
    /// Solve a CompositeCurveProblem using the 
    /// Levenberg–Marquardt algorithm
    /// </summary>
    public sealed class CompositeCurveSolver
    {
        /// <summary>
        /// Hard weight.
        /// </summary>
        public double HardWeight { get; set; } = 1e6;
        /// <summary>
        /// Used for finite difference approximation
        /// </summary>
        public double EpsFD { get; set; } = 1e-6;
        /// <summary>
        /// Max iterations of the solver
        /// </summary>
        public int MaxIter { get; set; } = 800;
        /// <summary>
        /// Solution is solved when the change in cost per step is lower than this value.
        /// </summary>
        public double RelTol { get; set; } = 1e-6;

        /// <summary>        /// 
        /// LMA: Pn+1 = Pn - ((hessian(Residual) + lambda(I))^-1)gradient(Residual)
        /// Where: Pn is the position vector at index n.
        /// Hessian(Residual) is the matrix of second partial derivates of the residual vector.
        /// Lambda is the weighting factor (chooses between Newton Raphson (small) and Gradient descent (large)).
        /// I is the identity matrix of the same dimension as the hessian matrix.
        /// Gradient(Residual) is the matrix of first order partial derivatives of the residual vector.
        /// </summary>
        /// <param name="problem"></param>
        /// <param name="optimizeP0">Include the initial position in the parameter vector.</param>
        /// <param name="optimizeLength">Include the length of each CurveSpec in the parameter vector.</param>
        /// <param name="optimizeR0">Include the R0 matrix in the parameter vector.</param>
        /// <returns></returns>
        public CompositeCurveSolverResult Solve(CompositeCurveProblem problem, bool optimizeP0 = true, bool optimizeLength = true, bool optimizeR0 = true)
        {
            List<CurveSpec[]> trialSolutions = new List<CurveSpec[]>();

            //Console.WriteLine("CompositeCurveProblem.Solve()");
            Pack pack = new Pack(problem.Seeds, optimizeP0, optimizeLength, optimizeR0);
            double[] currentParameters = pack.Pack2(problem.Seeds);
            double[] currentResiduals = BuildResidualOnly(problem.Seeds, problem.Constraints);
            double currentCost = Helpers.Dot(currentResiduals, currentResiduals);
            int _numParams = currentParameters.Length;

            trialSolutions.Add(problem.Seeds);

            double lambda = 1e-2;
            double lambdaFactor = 10;

            CurveSpec[] specs = problem.Seeds;
            for (int i = 0; i < MaxIter; i++)
            {
                specs = pack.Unpack(currentParameters);

                //Console.WriteLine();
                //Console.WriteLine($"==================== New Iteration: {i} ====================");
                // Calculate current residuals jacobian and hessian approximation

                double[,] Jacobian = BuildResidualsAndJacobian(specs, currentParameters, problem.Constraints, pack, out currentResiduals);
                //Helpers.PrintVector(currentResiduals);
                //Helpers.PrintMat(Jacobian);
                var JT = Helpers.Transpose(Jacobian);
                var hessian = Helpers.Multiply(JT, Jacobian);
                //Console.WriteLine($"Calculate current residuals jacobian and hessian approximation");

                // Calculate gradient
                var grad = Helpers.Multiply(JT, currentResiduals);
                //Console.WriteLine($"Calculate gradient");

                // (hess + lambda(I)) * dP = -grad | Ax=B
                var lambdaIdentity = Helpers.Multiply(lambda, RigidTransform.Identity(_numParams).R);
                var aMatrix = Helpers.Add(hessian, lambdaIdentity);
                var bMatrix = Helpers.Multiply(grad, -1);
                //Console.WriteLine($"(hess + lambda(I)) * dP = -grad");

                double[] dP;
                try
                {
                    // Solve for dp
                    var (LU, piv, pivSign) = LA.LUFactor(aMatrix);
                    dP = LA.LUSolve(LU, piv, bMatrix);
                }
                catch (Exception)
                {
                    return new CompositeCurveSolverResult(new CurveSpec[][] { specs }, specs, false, "Solution did not converge, matrix not factorizable.", i, currentCost);
                }
                //Console.WriteLine($"try solve for dP: ");
                //Helpers.PrintVector(dP);
                //Console.Write("Current parameters: ");
                //Helpers.PrintVector(currentParameters);
                //Console.WriteLine();

                // Calculate new parameters and new cost
                double[] newParameters = Helpers.Add(currentParameters, dP);
                //Console.WriteLine("New parameters: ");
                //Helpers.PrintVector(newParameters);
                //Console.WriteLine();
                specs = pack.Unpack(newParameters);

                /*Console.WriteLine($"New Specs");
                for (int j = 0; j < specs.Length; j++)
                {
                    specs[j].ShowInfo();
                    Console.WriteLine();
                }*/
                double[] newResiduals = BuildResidualOnly(specs, problem.Constraints);
                //Console.WriteLine("New residuals: ");
                //Helpers.PrintVector(newResiduals);
                //Console.WriteLine();
                double newCost = Helpers.Dot(newResiduals, newResiduals);
                //Console.WriteLine($"Calculate new parameters and new cost: {newCost}");

                // Decide whether to keep the new cost or not
                double costChange = currentCost - newCost;
                if (costChange > 0)
                {
                    trialSolutions.Add(specs);
                    currentParameters = newParameters;
                    //currentResiduals = newResiduals;
                    currentCost = newCost;
                    lambda /= lambdaFactor;

                    if (costChange < RelTol)
                    {
                        return new CompositeCurveSolverResult(trialSolutions.ToArray(), specs, true, "Solution converged, cost change below threshold.", i, currentCost);
                    }
                }
                else
                {
                    lambda *= lambdaFactor;
                }
                //Console.WriteLine($"Decide whether to keep the new cost or not");
            }
            //Console.Write($"Unpacking best theta: ");
            //Helpers.PrintVector(bestTheta);
            //Console.WriteLine();
            return new CompositeCurveSolverResult(trialSolutions.ToArray(), specs, false, "Solution did not converge, max iterations reached", MaxIter, currentCost);
        }

        private double[,] BuildResidualsAndJacobian(IReadOnlyList<CurveSpec> specs, double[] currenParams,
                                                                    List<ICompositeConstraint> constraints, Pack pack, out double[] residuals)
        {
            // Residual at the current constraints
            residuals = BuildResidualOnly(specs, constraints);
            int m = residuals.Length, p = currenParams.Length;
            var J = new double[m, p];
            // For each parameter add +h and calculate the FD derivative at i,j
            for (int j = 0; j < p; j++)
            {
                double h = EpsFD * (1.0 + System.Math.Abs(currenParams[j]));
                var params_h = (double[])currenParams.Clone();
                params_h[j] += h;
                var specs_h = pack.Unpack(params_h);
                // residual with the jth component of the parameter array increased by h
                var r_h = BuildResidualOnly(specs_h, constraints);
                int i = 0;
                try
                {
                    for (i = 0; i < m; i++) J[i, j] = (r_h[i] - residuals[i]) / h;
                } catch (Exception e)
                {
                    Console.WriteLine($"i, j => {i}, {j}");
                    Console.WriteLine($"J rows/cols => {m}/{p}");
                    Console.WriteLine($"len(r)/len(r_h) => {residuals.Length}/{r_h.Length}");
                    Console.WriteLine($"len(params)/len(params_h) => {p}/{params_h.Length}");
                    Console.WriteLine(e.Message);
                }
            }
            return J;
        }

        private double[] BuildResidualOnly(IReadOnlyList<CurveSpec> specs, List<ICompositeConstraint> constraints)
        {
            var all = new List<double>();
            foreach (ICompositeConstraint c in constraints)
            {
                //Console.WriteLine($"Building Residual for: ");
                //c.ShowInfo();
                var res = c.Residual(specs);
                //Console.WriteLine("Built residual only residual");
                double w = (c.Type == ConstraintType.Hard ? HardWeight : 1.0);
                if (c.Weight != 1.0) w *= c.Weight;
                if (w != 1.0) for (int i = 0; i < res.Length; i++) res[i] *= w;
                //Console.Write("residual: ");
                //Helpers.PrintVector(res);
                //Console.WriteLine();
                all.AddRange(res);
            }
            //Console.WriteLine("Returning residual...");
            return all.ToArray();
        }

        // ---------------- parameter pack across all segments ----------------
        private sealed class Pack
        {
            private readonly int _numSeg;
            private readonly int _N;
            private readonly bool _optP0, _optL, _optR0;
            private readonly CurveSpec[] _seed;
            private readonly IParamCurvatureLaw?[] _laws;
            private readonly int[] _pCounts;     // per-seg parameter counts
            private readonly int _totalP;
            private readonly int _totalL;

            public readonly int[] LengthIndices; // indices in the parameter vector that correspond with
                                                  // the length parameters of the individual curve specs
                                                  // keep this in case we want to limit length increases
                                                  // per step in the LMA loop.

            public Pack(CurveSpec[] seeds, bool optP0, bool optL, bool optR0)
            {
                _seed = seeds;
                _numSeg = seeds.Length;
                _N = seeds[0].N; // assume consistent N across segments
                _optP0 = optP0; 
                _optL = optL;
                _optR0 = optR0;
                // the segment laws
                _laws = new IParamCurvatureLaw?[_numSeg];
                // number of parameters in each segment
                _pCounts = new int[_numSeg];
                // total parameter count
                int total = 0;
                // total length parameter count
                int totalLIndices = 0;
                for (int s = 0; s < _numSeg; s++)
                {
                    var law = seeds[s].Kappa as IParamCurvatureLaw;
                    _laws[s] = law;

                    int cnt = 0;
                    if (_optP0) cnt += _N;
                    if (_optL)
                    {
                        cnt += 1;
                        totalLIndices++;
                    }
                    if (_optR0) cnt += SOParam.ParamCount(_N);
                    if (law != null) cnt += law.GetParams().Length;

                    _pCounts[s] = cnt;
                    total += cnt;
                }
                _totalP = total;
                _totalL = totalLIndices;

                LengthIndices = new int[_totalL];
            }

            /// <summary>
            /// Pack a CurveSpec list into a single dimensional array.
            /// Assuming we are optimizing P0, L and R0 it will look like this:
            /// [P00, P01, ... P0N, L, 
            /// </summary>
            /// <param name="specs"></param>
            /// <returns></returns>
            public double[] Pack2(CurveSpec[] specs)
            {
                var parameters = new double[_totalP];
                // current index in the parameter array
                // new parameters will go at this index.
                int k = 0;
                int l = 0;
                for (int s = 0; s < _numSeg; s++)
                {
                    var spec = specs[s];
                    if (_optP0) { Array.Copy(spec.P0, 0, parameters, k, _N); k += _N; }
                    if (_optL)
                    {
                        // After calling pack we can access the updated
                        // Length indices of the packed parameter vector
                        // with this property.
                        LengthIndices[l] = k;
                        l++;
                        parameters[k] = spec.Length; 
                        k++;
                    }
                    if (_optR0)
                    {
                        // start at zero offset; solver will find angles. Pack zeros.
                        int d = SOParam.ParamCount(_N);
                        for (int i = 0; i < d; i++) parameters[k++] = 0.0;
                    }
                    if (_laws[s] != null)
                    {
                        var p = _laws[s]!.GetParams();
                        Array.Copy(p, 0, parameters, k, p.Length);
                        k += p.Length;
                    }
                }
                return parameters;
            }

            /// <summary>
            /// Build CurveSpec list from packed parameter array.
            /// </summary>
            /// <param name="parameters"></param>
            /// <returns></returns>
            public CurveSpec[] Unpack(double[] parameters)
            {
                var specs = new CurveSpec[_numSeg];
                int k = 0; //keeps track of the indices of each parameter type
                int oldK = k;
                for (int s = 0; s < _numSeg; s++)
                {
                    var seed = _seed[s];

                    double[] P0 = _optP0 ? Slice(parameters, ref k, _N) : (double[])seed.P0.Clone();
                    //Console.WriteLine($"Unpack P0 | k ({oldK}=>{k}) | ({s} : {_numSeg-1})");
                    //Console.Write("P0 =");
                    //Helpers.PrintVector(P0);
                    //oldK = k;
                    double L = _optL ? Math.Abs(parameters[k++]) : seed.Length;
                    //Console.WriteLine($"Unpack L | k ({oldK}=>{k}) | ({s} : {_numSeg-1})");
                    //oldK = k;
                    double[,] R0 = (double[,])seed.R0.Clone();
                    if (_optR0)
                    {
                        int d = SOParam.ParamCount(_N);
                        var w = Slice(parameters, ref k, d);
                        var dR = SOParam.BuildRotation(_N, w);
                        R0 = Helpers.Multiply(dR, R0);
                    }
                    //Console.WriteLine($"Unpack R0 | k ({oldK}=>{k}) | ({s} : {_numSeg-1})");
                    //oldK = k;
                    ICurvatureLaw law;
                    if (_laws[s] != null)
                    {
                        //Console.Write("Is IParamCurvatureLaw => ");
                        int q = _laws[s]!.GetParams().Length;
                        var p = Slice(parameters, ref k, q);
                        law = _laws[s]!.CloneWithParams(p);
                    }
                    else law = seed.Kappa;
                    //Console.WriteLine($"Unpack kappa | k ({oldK}=>{k}) | ({s} : {_numSeg-1})");
                    //oldK = k;

                    specs[s] = new CurveSpec(_N, L, P0, R0, law, seed.Frame);
                    //string p0 = "";
                    //string r0 = "[";
                    //for (int i = 0; i < P0.Length; i++) p0 += P0[i].ToString();
                    //for (int i = 0; i < R0.GetLength(0); i++)
                    //{
                    //    r0 += "[";
                    //    for (int j = 0; j < R0.GetLength(1); j++)
                    //    {
                    //        r0 += R0[i, j].ToString();
                    //        r0 += j == R0.GetLength(1) ? "" : ", ";
                    //    }
                    //    r0 += "], ";
                    //}
                    //r0 += "]";
                    //Console.WriteLine($"Write spec number {s} | _N: {_N} | L: {L} | P0: [{p0}] | R0: {r0} | law: {law.ToString()} | frame: {seed.Frame.ToString()} | ({s} : {_numSeg-1})");
                    //Console.WriteLine();
                }
                //Console.WriteLine($"Returning specs (count: {specs.Length})");
                return specs;
            }


            private static double[] Slice(double[] v, ref int k, int len)
            {
                var a = new double[len];
                Array.Copy(v, k, a, 0, len);
                k += len;
                return a;
            }
        }
    }
}
