using ArcFrame.Core.Results;
using System;

namespace ArcFrame.Solvers.Core
{
    /// <summary>
    /// Very small, allocation-conscious Levenberg–Marquardt solver for
    /// low-dimensional nonlinear least-squares problems.
    /// 
    /// Designed for tiny parameter vectors (n <= ~16) and residuals,
    /// such as the 3D triple clothoid fitter parameter vector p.
    /// </summary>
    internal interface ISmallNonlinearLeastSquaresProblem
    {
        /// <summary>Number of parameters.</summary>
        int ParameterCount { get; }

        /// <summary>Number of residuals.</summary>
        int ResidualCount { get; }

        /// <summary>
        /// Evaluate residuals at the given parameter vector.
        /// The <paramref name="parameters"/> array is length ParameterCount.
        /// The <paramref name="residuals"/> array is length ResidualCount and
        /// must be fully overwritten by the implementation.
        /// 
        /// Implementations should NOT allocate per-call; write into the buffer.
        /// </summary>
        void Evaluate(double[] parameters, double[] residuals);
    }

    /// <summary>
    /// Tiny Levenberg–Marquardt implementation specialized for very small dense systems.
    /// 
    /// The goal is to avoid heap allocations inside the iteration loop and keep the
    /// implementation Unity/IL2CPP-friendly (no unsafe, no Span, no LINQ, no reflection).
    /// </summary>
    internal static class SmallLevenbergMarquardt
    {
        internal struct Options
        {
            /// <summary>Maximum LM iterations.</summary>
            public int MaxIterations;

            /// <summary>Initial damping parameter λ.</summary>
            public double InitialLambda;

            /// <summary>Multiplicative factor when a step is rejected (&gt; 1).</summary>
            public double LambdaIncrease;

            /// <summary>Multiplicative factor when a step is accepted (&lt; 1).</summary>
            public double LambdaDecrease;

            /// <summary>
            /// Finite-difference step factor for Jacobian:
            /// h_j = StepEps * (|p_j| + 1).
            /// </summary>
            public double StepEps;

            /// <summary>Relative cost tolerance.</summary>
            public double CostTolerance;

            /// <summary>Relative parameter-step tolerance.</summary>
            public double StepTolerance;

            /// <summary>Minimum allowed λ.</summary>
            public double MinLambda;

            /// <summary>Maximum allowed λ before we declare failure.</summary>
            public double MaxLambda;

            /// <summary>Return a set of reasonable defaults for tiny problems.</summary>
            public static Options Default
            {
                get
                {
                    return new Options
                    {
                        MaxIterations = 64,
                        InitialLambda = 1e-3,
                        LambdaIncrease = 10.0,
                        LambdaDecrease = 0.1,
                        StepEps = 1e-6,
                        CostTolerance = 1e-12,
                        StepTolerance = 1e-9,
                        MinLambda = 1e-12,
                        MaxLambda = 1e12
                    };
                }
            }
        }

        /// <summary>
        /// Solve a tiny nonlinear least-squares problem with LM.
        /// 
        /// The solver clones <paramref name="initialParameters"/> so the input is never modified.
        /// </summary>
        internal static Result Solve(
            ISmallNonlinearLeastSquaresProblem problem,
            double[] initialParameters,
            Options options)
        {
            if (problem == null) throw new ArgumentNullException(nameof(problem));
            if (initialParameters == null) throw new ArgumentNullException(nameof(initialParameters));

            int n = problem.ParameterCount;
            int m = problem.ResidualCount;

            if (n <= 0 || m <= 0)
                throw new ArgumentException("Problem must have positive parameter and residual counts.");

            if (initialParameters.Length != n)
                throw new ArgumentException("initialParameters.Length must equal problem.ParameterCount.");

            // Initialize options if caller passed a default-initialized struct.
            if (options.MaxIterations <= 0)
                options = Options.Default;

            // Current parameters and residuals
            double[] p = (double[])initialParameters.Clone();
            double[] r = new double[m];
            double[] rTrial = new double[m];
            double[] pTrial = new double[n];

            // Jacobian J (m x n) stored row-major in a flat array: J[i*n + j]
            double[] J = new double[m * n];

            // Normal equations: A = J^T J (n x n), g = J^T r (n)
            double[] A = new double[n * n];
            double[] g = new double[n];

            // LM step vector
            double[] delta = new double[n];

            // Evaluate initial residuals and cost
            problem.Evaluate(p, r);
            double cost = 0.5 * Dot(r, r);

            double lambda = options.InitialLambda > 0.0 ? options.InitialLambda : 1e-3;

            bool converged = false;
            int iter = 0;

            for (iter = 0; iter < options.MaxIterations; ++iter)
            {
                // Build Jacobian via forward finite differences (reuses rTrial as a scratch buffer)
                BuildJacobian(problem, p, r, J, options.StepEps, rTrial);

                // Build normal equations A = J^T J and g = J^T r
                BuildNormalEquations(J, r, m, n, A, g);

                // Apply LM damping: A_ii *= (1 + lambda)
                ApplyDamping(A, n, lambda);

                // Solve (A) * delta = -g
                for (int j = 0; j < n; ++j)
                    delta[j] = -g[j];

                bool ok = SolveSymmetricSystemInPlace(A, delta, n);
                if (!ok)
                {
                    // Matrix was singular or nearly so; treat as non-converged.
                    break;
                }

                // Compute trial parameters and norms
                double normDelta = 0.0;
                double normP = 0.0;
                for (int j = 0; j < n; ++j)
                {
                    double pj = p[j];
                    double dj = delta[j];
                    double trial = pj + dj;
                    pTrial[j] = trial;

                    normDelta += dj * dj;
                    normP += pj * pj;
                }
                normDelta = Math.Sqrt(normDelta);
                normP = Math.Sqrt(normP);

                // Early exit: tiny step
                if (normDelta <= options.StepTolerance * (normP + options.StepTolerance))
                {
                    converged = true;
                    break;
                }

                // Evaluate trial residuals and cost
                problem.Evaluate(pTrial, rTrial);
                double trialCost = 0.5 * Dot(rTrial, rTrial);
                double costChange = cost - trialCost;

                if (trialCost < cost)
                {
                    // Accept step
                    Array.Copy(pTrial, p, n);
                    Array.Copy(rTrial, r, m);
                    cost = trialCost;
                    lambda = Math.Max(lambda * options.LambdaDecrease, options.MinLambda);

                    // Check relative cost change
                    double relCostChange = costChange / (Math.Abs(cost) + 1e-12);
                    if (relCostChange < options.CostTolerance)
                    {
                        converged = true;
                        break;
                    }
                }
                else
                {
                    // Reject step, increase damping
                    lambda = Math.Min(lambda * options.LambdaIncrease, options.MaxLambda);
                    if (lambda >= options.MaxLambda)
                    {
                        // We're essentially stuck.
                        break;
                    }
                }
            }

            return new Result(converged, iter, cost, lambda, p);
        }

        private static double Dot(double[] v, double[] w)
        {
            double sum = 0.0;
            int len = v.Length;
            for (int i = 0; i < len; ++i)
                sum += v[i] * w[i];
            return sum;
        }

        private static void BuildJacobian(
            ISmallNonlinearLeastSquaresProblem problem,
            double[] p,
            double[] rAtP,
            double[] J,
            double stepEps,
            double[] workResidual)
        {
            int n = problem.ParameterCount;
            int m = problem.ResidualCount;

            for (int j = 0; j < n; ++j)
            {
                double pj = p[j];
                double h = stepEps * (Math.Abs(pj) + 1.0);
                if (h == 0.0) h = stepEps;

                p[j] = pj + h;
                problem.Evaluate(p, workResidual);
                p[j] = pj;

                double invH = 1.0 / h;
                int rowOffset = 0;
                for (int i = 0; i < m; ++i)
                {
                    double diff = (workResidual[i] - rAtP[i]) * invH;
                    J[rowOffset + j] = diff;
                    rowOffset += n;
                }
            }
        }

        private static void BuildNormalEquations(
            double[] J,
            double[] r,
            int m,
            int n,
            double[] A,
            double[] g)
        {
            // Zero A and g
            Array.Clear(A, 0, A.Length);
            Array.Clear(g, 0, g.Length);

            // A = J^T J, g = J^T r
            for (int i = 0; i < m; ++i)
            {
                int rowOffset = i * n;
                double ri = r[i];

                // g += J_i^T * r_i
                for (int j = 0; j < n; ++j)
                {
                    double Jij = J[rowOffset + j];
                    g[j] += Jij * ri;
                }

                // A += J_i^T * J_i  (rank-1 update)
                for (int j = 0; j < n; ++j)
                {
                    double Jij = J[rowOffset + j];
                    int aRow = j * n;
                    for (int k = j; k < n; ++k)
                    {
                        A[aRow + k] += Jij * J[rowOffset + k];
                    }
                }
            }

            // Fill lower triangle of A (since we only wrote upper)
            for (int j = 0; j < n; ++j)
            {
                int rowj = j * n;
                for (int k = j + 1; k < n; ++k)
                {
                    A[k * n + j] = A[rowj + k];
                }
            }
        }

        private static void ApplyDamping(double[] A, int n, double lambda)
        {
            for (int i = 0; i < n; ++i)
            {
                int idx = i * n + i;
                A[idx] = A[idx] * (1.0 + lambda);
            }
        }

        /// <summary>
        /// Solve A x = b in-place using Gaussian elimination with partial pivoting.
        /// 
        /// A is n x n in row-major form, b is length n. On success, b contains the solution x.
        /// On failure (singular or nearly singular matrix), returns false.
        /// </summary>
        private static bool SolveSymmetricSystemInPlace(double[] A, double[] b, int n)
        {
            const double pivotEps = 1e-15;

            // Forward elimination
            for (int k = 0; k < n; ++k)
            {
                // Find pivot row
                int pivotRow = k;
                double max = Math.Abs(A[k * n + k]);
                for (int i = k + 1; i < n; ++i)
                {
                    double val = Math.Abs(A[i * n + k]);
                    if (val > max)
                    {
                        max = val;
                        pivotRow = i;
                    }
                }

                if (max < pivotEps)
                {
                    return false; // singular
                }

                // Swap rows k and pivotRow in A and entries in b
                if (pivotRow != k)
                {
                    int rowK = k * n;
                    int rowP = pivotRow * n;
                    for (int j = 0; j < n; ++j)
                    {
                        double tmp = A[rowK + j];
                        A[rowK + j] = A[rowP + j];
                        A[rowP + j] = tmp;
                    }

                    double tb = b[k];
                    b[k] = b[pivotRow];
                    b[pivotRow] = tb;
                }

                // Eliminate below diag
                int rowKIndex = k * n;
                double diag = A[rowKIndex + k];
                double invDiag = 1.0 / diag;

                for (int i = k + 1; i < n; ++i)
                {
                    int rowI = i * n;
                    double factor = A[rowI + k] * invDiag;
                    A[rowI + k] = 0.0;

                    for (int j = k + 1; j < n; ++j)
                    {
                        A[rowI + j] -= factor * A[rowKIndex + j];
                    }

                    b[i] -= factor * b[k];
                }
            }

            // Back substitution
            for (int i = n - 1; i >= 0; --i)
            {
                int row = i * n;
                double sum = b[i];
                for (int j = i + 1; j < n; ++j)
                {
                    sum -= A[row + j] * b[j];
                }

                double diag = A[row + i];
                if (Math.Abs(diag) < pivotEps)
                    return false;

                b[i] = sum / diag;
            }

            return true;
        }
    }
}
