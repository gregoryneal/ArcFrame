using ArcFrame.Core.Constraints;
using ArcFrame.Core.Math;
using ArcFrame.Core.Params;

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
        public List<ICompositeConstraint> Constraints { get; } = new();
        public CompositeCurveProblem(CurveSpec[] seeds) { Seeds = seeds; }

        /// <summary>
        /// Wrap a single curve constraint in a SegmentConstraintWrapper for ICompositeConstraint
        /// </summary>
        /// <param name="single"></param>
        public void Add(ICurveConstraint single) => Constraints.Add(new SegmentConstraintWrapper(single));
        public void Add(ICompositeConstraint joint) => Constraints.Add(joint);
    }

    /// <summary>
    /// Solve a CompositeCurveProblem using the 
    /// Levenberg–Marquardt algorithm
    /// </summary>
    public sealed class CompositeCurveSolver
    {
        public double HardWeight { get; set; } = 1e6;
        public double EpsFD { get; set; } = 1e-6;
        public int MaxIter { get; set; } = 80;
        public double RelTol { get; set; } = 1e-6;

        /// <summary>
        /// TODO: Right now the algorithm is implemented incorrectly, it updates the cost
        /// no matter what and allows divergence....
        /// </summary>
        /// <param name="problem"></param>
        /// <param name="optimizeP0"></param>
        /// <param name="optimizeLength"></param>
        /// <param name="optimizeR0"></param>
        /// <returns></returns>
        public CurveSpec[] Solve(CompositeCurveProblem problem, bool optimizeP0 = true, bool optimizeLength = true, bool optimizeR0 = true)
        {
            //Console.WriteLine("CompositeCurveProblem.Solve()");
            var pack = new Pack(problem.Seeds, optimizeP0, optimizeLength, optimizeR0);
            var theta = pack.Pack2(problem.Seeds);

            double lambda = 1e-2;
            double bestCost = double.PositiveInfinity;
            var bestTheta = (double[])theta.Clone();

            for (int it = 0; it < MaxIter; it++)
            {
                //Console.WriteLine($"Solve ({it}, {MaxIter})");
                var specs = pack.Unpack(theta);
                //Console.WriteLine("Unpack");
                var (r, J) = BuildResidualsAndJacobian(specs, theta, problem.Constraints, pack);
                //Console.WriteLine("BuildResidualsAndJacobian");
                double cost = Helpers.Dot(r, r);
                if (cost < bestCost) { bestCost = cost; bestTheta = (double[])theta.Clone(); }

                var JT = Transpose(J);
                var JTJ = Helpers.Multiply(JT, J);
                AddLevenbergDamping(JTJ, lambda);
                //Console.WriteLine("AddLevenbergDamping");
                var g = Helpers.Multiply(JT, r);
                //Console.WriteLine("Multiply");
                for (int i = 0; i < g.Length; i++) g[i] = -g[i];
                //Console.WriteLine("Invert g");

                if (!Helpers.TryInvert(JTJ, out var JTJinv))
                {
                    //Console.WriteLine($"JTJ uninvertible, terminating");
                    //Helpers.PrintMat(JTJ);
                    break;
                }
                var delta = Helpers.Multiply(JTJinv!, g);

                if (Helpers.Len(delta) <= RelTol * (RelTol + Helpers.Len(theta)))
                {
                    //Console.WriteLine($"{Helpers.Len(delta)} <= {RelTol} * ({RelTol} + {Helpers.Len(theta)}");
                    break;
                }

                var thetaNew = Helpers.Add(theta, delta);
                var specsNew = pack.Unpack(thetaNew);
                //Console.WriteLine("Unpack thetaNew");
                var rNew = BuildResidualOnly(specsNew, problem.Constraints); // <-- error is here
                //Console.WriteLine("post BuildResidualOnly");
                double costNew = Helpers.Dot(rNew, rNew);

                double oldLambda = lambda;
                if (costNew < cost) { theta = thetaNew; lambda *= 0.5; bestTheta = (double[])thetaNew.Clone(); }
                else { lambda *= 4.0; }
                //Console.WriteLine($"Update lambda {oldLambda}=>{lambda}\n");
            }

            //Console.Write($"Unpacking best theta: ");
            //Helpers.PrintVector(bestTheta);
            //Console.WriteLine();
            return pack.Unpack(bestTheta);
        }

        private (double[] r, double[,] J) BuildResidualsAndJacobian(IReadOnlyList<CurveSpec> specs, double[] theta,
                                                                    List<ICompositeConstraint> constraints, Pack pack)
        {
            var r = BuildResidualOnly(specs, constraints);
            int m = r.Length, p = theta.Length;
            var J = new double[m, p];
            for (int j = 0; j < p; j++)
            {
                double h = EpsFD * (1.0 + System.Math.Abs(theta[j]));
                var theta_h = (double[])theta.Clone();
                theta_h[j] += h;
                var specs_h = pack.Unpack(theta_h);
                var r_h = BuildResidualOnly(specs_h, constraints);
                for (int i = 0; i < m; i++) J[i, j] = (r_h[i] - r[i]) / h;
            }
            return (r, J);
        }

        private double[] BuildResidualOnly(IReadOnlyList<CurveSpec> specs, List<ICompositeConstraint> constraints)
        {
            var all = new List<double>();
            foreach (ICompositeConstraint c in constraints)
            {
                //Console.WriteLine($"Building Residual for: ");
                //c.ShowInfo();
                var res = c.Residual(specs);
                //Console.WriteLine("Build residual only residual");
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

        private static double[,] Transpose(double[,] A)
        {
            int m = A.GetLength(0), n = A.GetLength(1);
            var AT = new double[n, m];
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    AT[j, i] = A[i, j];
            return AT;
        }

        private static void AddLevenbergDamping(double[,] H, double lambda)
        {
            int n = H.GetLength(0);
            for (int i = 0; i < n; i++) H[i, i] += lambda * H[i, i] + 1e-12;
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

            public Pack(CurveSpec[] seeds, bool optP0, bool optL, bool optR0)
            {
                _seed = seeds;
                _numSeg = seeds.Length;
                _N = seeds[0].N; // assume consistent N across segments
                _optP0 = optP0; _optL = optL; _optR0 = optR0;

                _laws = new IParamCurvatureLaw?[_numSeg];
                _pCounts = new int[_numSeg];
                int total = 0;
                for (int s = 0; s < _numSeg; s++)
                {
                    var law = seeds[s].Kappa as IParamCurvatureLaw;
                    _laws[s] = law;

                    int cnt = 0;
                    if (_optP0) cnt += _N;
                    if (_optL) cnt += 1;
                    if (_optR0) cnt += SOParam.ParamCount(_N);
                    if (law != null) cnt += law.GetParams().Length;

                    _pCounts[s] = cnt;
                    total += cnt;
                }
                _totalP = total;
            }

            public double[] Pack2(CurveSpec[] specs)
            {
                var theta = new double[_totalP];
                int k = 0;
                for (int s = 0; s < _numSeg; s++)
                {
                    var spec = specs[s];
                    if (_optP0) { Array.Copy(spec.P0, 0, theta, k, _N); k += _N; }
                    if (_optL) { theta[k++] = System.Math.Log(System.Math.Max(1e-9, spec.Length)); }
                    if (_optR0)
                    {
                        // start at zero offset; solver will find angles. Pack zeros.
                        int d = SOParam.ParamCount(_N);
                        for (int i = 0; i < d; i++) theta[k++] = 0.0;
                    }
                    if (_laws[s] != null)
                    {
                        var p = _laws[s]!.GetParams();
                        Array.Copy(p, 0, theta, k, p.Length);
                        k += p.Length;
                    }
                }
                return theta;
            }

            public CurveSpec[] Unpack(double[] theta)
            {
                var specs = new CurveSpec[_numSeg];
                int k = 0; //keeps track of the indices of each parameter type
                int oldK = k;
                for (int s = 0; s < _numSeg; s++)
                {
                    var seed = _seed[s];

                    oldK = k;
                    double[] P0 = _optP0 ? Slice(theta, ref k, _N) : (double[])seed.P0.Clone();
                    //Console.WriteLine($"Unpack P0 | k ({oldK}=>{k}) | ({s} : {_numSeg-1})");
                    //Console.Write("P0 =");
                    //Helpers.PrintVector(P0);
                    oldK = k;
                    double L = _optL ? System.Math.Exp(theta[k++]) : seed.Length;
                    //Console.WriteLine($"Unpack L | k ({oldK}=>{k}) | ({s} : {_numSeg-1})");

                    oldK = k;
                    double[,] R0 = (double[,])seed.R0.Clone();
                    if (_optR0)
                    {
                        int d = SOParam.ParamCount(_N);
                        var w = Slice(theta, ref k, d);
                        var dR = SOParam.BuildRotation(_N, w);
                        R0 = Helpers.Multiply(dR, R0);
                    }
                    //Console.WriteLine($"Unpack R0 | k ({oldK}=>{k}) | ({s} : {_numSeg-1})");
                    oldK = k;
                    ICurvatureLaw law;
                    if (_laws[s] != null)
                    {
                        //Console.Write("Is IParamCurvatureLaw => ");
                        int q = _laws[s]!.GetParams().Length;
                        var p = Slice(theta, ref k, q);
                        law = _laws[s]!.CloneWithParams(p);
                    }
                    else law = seed.Kappa;
                    //Console.WriteLine($"Unpack kappa | k ({oldK}=>{k}) | ({s} : {_numSeg-1})");

                    specs[s] = new CurveSpec(_N, L, P0, R0, law, seed.Frame);
                    string p0 = "";
                    string r0 = "[";
                    for (int i = 0; i < P0.Length; i++) p0 += P0[i].ToString();
                    for (int i = 0; i < R0.GetLength(0); i++)
                    {
                        r0 += "[";
                        for (int j = 0; j < R0.GetLength(1); j++)
                        {
                            r0 += R0[i, j].ToString();
                            r0 += j == R0.GetLength(1) ? "" : ", ";
                        }
                        r0 += "], ";
                    }
                    r0 += "]";
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
