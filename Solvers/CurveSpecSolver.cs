using ArcFrame.Core.Constraints;
using ArcFrame.Core.Math;
using ArcFrame.Core.Params;

namespace ArcFrame.Solvers
{
    /// <summary>
    /// Levenberg-Marquardt solver to fit a CurveSpec to a set of constraints.
    /// Optimizes P0, R0, Length, and curvature-law params if supported.
    /// </summary>
    public sealed class CurveSpecSolver
    {
        public double HardWeight { get; set; } = 1E6;
        public double EpsFD { get; set; } = 1E-6;
        public int MaxIter { get; set; } = 60;
        public double RelTol { get; set; } = 1E-6;
        public double SolveTol { get; set; } = 1E-8;

        public CurveSpec Solve(CurveSpec seed, IEnumerable<ICurveConstraint> constraints, bool optimizeP0 = true, bool optimizeLength = true, bool optimizeR0 = true)
        {
            var pack = new ParamPack(seed, optimizeP0, optimizeLength, optimizeR0);
            var theta = pack.Pack(seed);

            double lambda = 1e-2;
            double bestCost = double.PositiveInfinity;
            var bestTheta = (double[])theta.Clone();

            // Maybe remove the solve tolerance and just run all the iterations.. do some tests.
            if (SolveTol <= 0) SolveTol = 1E-12;

            for (int it = 0; it < MaxIter && bestCost > SolveTol; it++)
            {
                Console.WriteLine($"Solve iter: {it}/{MaxIter}");
                var spec = pack.Unpack(theta);
                var (r, J) = BuildResidualsAndJacobian(spec, theta, constraints, pack);
                double cost = Helpers.Dot(r, r);
                if (cost < bestCost)
                {
                    bestCost = cost;
                    bestTheta = (double[])theta.Clone();
                }
                Console.WriteLine("Solve 3 \n");

                // (J^T J + λ diag(J^T J)) Δ = -J^T r
                var JT = Transpose(J);
                var JTJ = Helpers.Multiply(JT, J);
                AddLevenbergDamping(JTJ, lambda);
                var g = Helpers.Multiply(Helpers.Multiply(JT, r), -1);
                //for (int i = 0; i < g.Length; i++) g[i] = -g[i];

                if (!Helpers.TryInvert(JTJ, out var JTJinv)) break;
                Console.WriteLine("Solve 4 \n");
                var delta = Helpers.Multiply(JTJinv!, g);

                if (Helpers.Len(delta) <= RelTol * (RelTol + Helpers.Len(theta))) break;

                var thetaNew = Helpers.Add(theta, delta);
                var specNew = pack.Unpack(thetaNew);
                Console.WriteLine("Solve 5 \n");
                var rNew = BuildResidualOnly(specNew, constraints);
                Console.WriteLine("Solve 6 \n");
                double costNew = Helpers.Dot(rNew, rNew);

                // adjust damping factor down as we approach the solution
                // this moves to Gauss-Newton method
                if (costNew < cost)
                {
                    theta = thetaNew;
                    lambda *= 0.5;

                    Console.WriteLine("Solve 7a \n");
                }
                // increase damping factor, this moves towards gradient descent
                else
                {
                    lambda *= 4.0;

                    Console.WriteLine("Solve 7b \n");
                }
                Console.WriteLine("Solve 8 \n");
            }
            Console.WriteLine("Solve 9 \n");

            return pack.Unpack(bestTheta);
        }

        private (double[] r, double[,] J) BuildResidualsAndJacobian(CurveSpec spec, double[] theta, IEnumerable<ICurveConstraint> constraints, ParamPack pack)
        {
            Console.WriteLine("BuildResidualsAndJacobian 1 \n");
            var r = BuildResidualOnly(spec, constraints);
            int m = r.Length;
            int p = theta.Length;
            var J = new double[m, p];

            for (int j = 0; j < p; j++)
            {
                double h = EpsFD * (1.0 + System.Math.Abs(theta[j]));
                var theta_h = (double[])theta.Clone();
                theta_h[j] += h;
                var spec_h = pack.Unpack(theta_h);
                var r_h = BuildResidualOnly(spec_h, constraints);
                for (int i = 0; i < m; i++) J[i, j] = (r_h[i] - r[i]) / h;
            }
            return (r, J);
        }

        private double[] BuildResidualOnly(CurveSpec spec, IEnumerable<ICurveConstraint> constraints)
        {
            Console.WriteLine("BuildResidualsOnly 1 \n");
            var all = new List<double>();
            foreach (var c in constraints)
            {
                if (c.SegmentIndex != 0) continue; // single-spec solver expects 0
                var res = c.Residual(spec);
                double w = (c.Type == ConstraintType.Hard ? HardWeight : 1.0);
                if (c.Weight != 1.0) w *= c.Weight;
                if (w != 1.0) for (int i = 0; i < res.Length; i++) res[i] *= w;
                all.AddRange(res);
            }
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

        private sealed class ParamPack
        {
            private readonly int _N;
            private readonly bool _optP0, _optL, _optR0;
            private readonly IParamCurvatureLaw? _kParam;
            private readonly CurveSpec _seed;

            public ParamPack(CurveSpec seed, bool optimizeP0, bool optimizeLength, bool optimizeR0)
            {
                _N = seed.N; _optP0 = optimizeP0; _optL = optimizeLength; _optR0 = optimizeR0;
                _kParam = seed.Kappa as IParamCurvatureLaw;
                _seed = seed;
            }

            /// <summary>
            /// Pack all of the curve parameters into a single dimensional array.
            /// </summary>
            /// <param name="spec"></param>
            /// <returns></returns>
            public double[] Pack(CurveSpec spec)
            {
                var list = new List<double>();
                if (_optP0) list.AddRange(spec.P0);
                if (_optL) list.Add(System.Math.Log(System.Math.Max(1e-9, spec.Length)));
                if (_optR0) list.AddRange(new double[SOParam.ParamCount(_N)]); // start at zero => identity offset
                if (_kParam != null) list.AddRange(_kParam.GetParams());
                return list.ToArray();
            }

            public CurveSpec Unpack(double[] theta)
            {
                int k = 0;
                double[] P0 = _optP0 ? Slice(theta, ref k, _N) : (double[])_seed.P0.Clone();
                double L = _optL ? System.Math.Exp(theta[k++]) : _seed.Length;

                double[,] R0 = (double[,])_seed.R0.Clone();
                if (_optR0)
                {
                    var w = Slice(theta, ref k, SOParam.ParamCount(_N));
                    var dR = SOParam.BuildRotation(_N, w);
                    R0 = Helpers.Multiply(dR, R0); // left-multiply offset onto seed R0
                }

                ICurvatureLaw law;
                if (_kParam != null)
                {
                    int q = _kParam.GetParams().Length;
                    var pv = Slice(theta, ref k, q);
                    law = _kParam.CloneWithParams(pv);
                }
                else law = _seed.Kappa;

                return new CurveSpec(_N, L, P0, R0, law, _seed.Frame);
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
