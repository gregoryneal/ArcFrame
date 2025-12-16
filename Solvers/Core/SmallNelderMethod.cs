using System;

namespace ArcFrame.Solvers.Core
{
    /// <summary>
    /// Very small, allocation-conscious Nelder–Mead solver for
    /// low-dimensional nonlinear least-squares problems.
    /// 
    /// Works with ISmallNonlinearLeastSquaresProblem by treating
    /// the objective as 0.5 * ||r(p)||^2.
    /// </summary>
    internal static class SmallNelderMead
    {
        internal struct Options
        {
            public int MaxIterations;

            /// <summary>Scale for building the initial simplex.</summary>
            public double InitialSimplexScale;

            /// <summary>Reflection coefficient (α).</summary>
            public double Reflection;

            /// <summary>Expansion coefficient (γ).</summary>
            public double Expansion;

            /// <summary>Contraction coefficient (ρ).</summary>
            public double Contraction;

            /// <summary>Shrink coefficient (σ).</summary>
            public double Shrink;

            /// <summary>Cost tolerance for convergence.</summary>
            public double CostTolerance;

            /// <summary>Simplex size tolerance (max distance from best).</summary>
            public double SimplexTolerance;

            public static Options Default
            {
                get
                {
                    return new Options
                    {
                        MaxIterations = 128,
                        InitialSimplexScale = 0.1,
                        Reflection = 1.0,
                        Expansion = 2.0,
                        Contraction = 0.5,
                        Shrink = 0.5,
                        CostTolerance = 1e-8,
                        SimplexTolerance = 1e-6
                    };
                }
            }
        }

        internal struct Result
        {
            public bool Converged;
            public int Iterations;
            public double Cost;
            public double[] Parameters;

            public Result(bool converged, int iterations, double cost, double[] parameters)
            {
                Converged = converged;
                Iterations = iterations;
                Cost = cost;
                Parameters = parameters;
            }
        }

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
                throw new ArgumentException("Invalid problem sizes.");

            if (initialParameters.Length != n)
                throw new ArgumentException("initialParameters.Length must equal problem.ParameterCount.");

            if (options.MaxIterations <= 0)
                options = Options.Default;

            int numVerts = n + 1;

            // Simplex vertices: flattened [vertexIndex * n + paramIndex]
            double[] simplex = new double[numVerts * n];
            double[] costs = new double[numVerts];

            // Work buffers
            double[] rWork = new double[m];
            double[] centroid = new double[n];
            double[] xr = new double[n];
            double[] xe = new double[n];
            double[] xc = new double[n];

            // Helper: evaluate cost = 0.5 * ||r(p)||^2
            double EvaluateCost(double[] p)
            {
                problem.Evaluate(p, rWork);
                double s = 0.0;
                for (int i = 0; i < m; ++i) s += rWork[i] * rWork[i];
                return 0.5 * s;
            }

            // Initialize simplex
            for (int j = 0; j < n; ++j)
                simplex[j] = initialParameters[j]; // vertex 0

            costs[0] = EvaluateCost(initialParameters);

            for (int v = 1; v < numVerts; ++v)
            {
                int offset = v * n;
                // base = p0
                for (int j = 0; j < n; ++j)
                    simplex[offset + j] = initialParameters[j];

                int jParam = v - 1;
                double step = options.InitialSimplexScale * (Math.Abs(initialParameters[jParam]) + 1.0);
                simplex[offset + jParam] += step;
                costs[v] = EvaluateCostVertex(simplex, offset, n, problem, rWork);
            }

            bool converged = false;
            int iter = 0;

            for (iter = 0; iter < options.MaxIterations; ++iter)
            {
                // Sort vertices by cost (ascending)
                SortSimplex(simplex, costs, numVerts, n);

                // Check cost convergence
                double bestCost = costs[0];
                double worstCost = costs[numVerts - 1];
                if (worstCost - bestCost < options.CostTolerance)
                {
                    converged = true;
                    break;
                }

                // Check simplex size convergence: max distance from best
                double maxDist = 0.0;
                for (int v = 1; v < numVerts; ++v)
                {
                    double d = 0.0;
                    int offV = v * n;
                    for (int j = 0; j < n; ++j)
                    {
                        double diff = simplex[offV + j] - simplex[j];
                        d += diff * diff;
                    }
                    if (d > maxDist) maxDist = d;
                }
                maxDist = Math.Sqrt(maxDist);
                if (maxDist < options.SimplexTolerance)
                {
                    converged = true;
                    break;
                }

                // Indices: best=0, worst=n, second worst=n-1
                int iBest = 0;
                int iWorst = numVerts - 1;
                int iSecondWorst = numVerts - 2;

                // Compute centroid of all vertices except worst
                for (int j = 0; j < n; ++j) centroid[j] = 0.0;
                for (int v = 0; v < numVerts - 1; ++v)
                {
                    int offV = v * n;
                    for (int j = 0; j < n; ++j)
                        centroid[j] += simplex[offV + j];
                }
                double invN = 1.0 / n;
                for (int j = 0; j < n; ++j)
                    centroid[j] *= invN;

                int offWorst = iWorst * n;

                // Reflection xr = c + α (c - x_worst)
                for (int j = 0; j < n; ++j)
                {
                    double cw = centroid[j] - simplex[offWorst + j];
                    xr[j] = centroid[j] + options.Reflection * cw;
                }
                double costR = EvaluateCost(xr);

                if (costR < costs[iBest])
                {
                    // Try expansion
                    for (int j = 0; j < n; ++j)
                        xe[j] = centroid[j] + options.Expansion * (xr[j] - centroid[j]);
                    double costE = EvaluateCost(xe);

                    if (costE < costR)
                    {
                        // Accept expansion
                        CopyVertex(xe, 0, simplex, offWorst, n);
                        costs[iWorst] = costE;
                    }
                    else
                    {
                        // Accept reflection
                        CopyVertex(xr, 0, simplex, offWorst, n);
                        costs[iWorst] = costR;
                    }
                }
                else if (costR < costs[iSecondWorst])
                {
                    // Accept reflection
                    CopyVertex(xr, 0, simplex, offWorst, n);
                    costs[iWorst] = costR;
                }
                else
                {
                    // Contraction
                    bool outside = costR < costs[iWorst];

                    if (outside)
                    {
                        // Outside contraction: xc = c + ρ (xr - c)
                        for (int j = 0; j < n; ++j)
                            xc[j] = centroid[j] + options.Contraction * (xr[j] - centroid[j]);
                    }
                    else
                    {
                        // Inside contraction: xc = c + ρ (x_worst - c)
                        for (int j = 0; j < n; ++j)
                            xc[j] = centroid[j] + options.Contraction * (simplex[offWorst + j] - centroid[j]);
                    }

                    double costC = EvaluateCost(xc);

                    if (costC < (outside ? costR : costs[iWorst]))
                    {
                        // Accept contraction
                        CopyVertex(xc, 0, simplex, offWorst, n);
                        costs[iWorst] = costC;
                    }
                    else
                    {
                        // Shrink towards best
                        int offBest = iBest * n;
                        for (int v = 1; v < numVerts; ++v)
                        {
                            int offV = v * n;
                            for (int j = 0; j < n; ++j)
                            {
                                simplex[offV + j] = simplex[offBest + j]
                                    + options.Shrink * (simplex[offV + j] - simplex[offBest + j]);
                            }
                            costs[v] = EvaluateCostVertex(simplex, offV, n, problem, rWork);
                        }
                    }
                }
            }

            // Best vertex is simplex[0..n-1]
            double[] bestP = new double[n];
            for (int j = 0; j < n; ++j)
                bestP[j] = simplex[j];

            // We already sorted in the last iteration or before
            SortSimplex(simplex, costs, numVerts, n);
            double bestCostFinal = costs[0];

            return new Result(converged, iter, bestCostFinal, bestP);
        }

        private static double EvaluateCostVertex(
            double[] simplex,
            int offset,
            int n,
            ISmallNonlinearLeastSquaresProblem problem,
            double[] rWork)
        {
            // copy vertex into temp param array
            // (we can alias directly since the simplex slice is contiguous)
            problem.Evaluate(slice: simplex, offset: offset, length: n, residuals: rWork);
            double s = 0.0;
            int m = problem.ResidualCount;
            for (int i = 0; i < m; ++i) s += rWork[i] * rWork[i];
            return 0.5 * s;
        }

        /// <summary>
        /// Extension on ISmallNonlinearLeastSquaresProblem to avoid allocations when
        /// evaluating a slice of a larger parameter array (simplex vertex).
        /// </summary>
        private static void Evaluate(
            this ISmallNonlinearLeastSquaresProblem problem,
            double[] slice,
            int offset,
            int length,
            double[] residuals)
        {
            // Copy into a small local parameter array if needed
            // (n is small, so this is cheap; but we can also allow
            // the ISmallNonlinearLeastSquaresProblem implementation
            // to manage its own scratch buffers if desired).
            int n = problem.ParameterCount;
            double[] p = new double[n];
            for (int j = 0; j < n; ++j)
                p[j] = slice[offset + j];

            problem.Evaluate(p, residuals);
        }

        private static void CopyVertex(double[] src, int srcOffset, double[] dst, int dstOffset, int n)
        {
            for (int j = 0; j < n; ++j)
                dst[dstOffset + j] = src[srcOffset + j];
        }

        private static void SortSimplex(double[] simplex, double[] costs, int numVerts, int n)
        {
            for (int i = 0; i < numVerts - 1; ++i)
            {
                int bestIdx = i;
                double bestCost = costs[i];
                for (int j = i + 1; j < numVerts; ++j)
                {
                    if (costs[j] < bestCost)
                    {
                        bestCost = costs[j];
                        bestIdx = j;
                    }
                }

                if (bestIdx != i)
                {
                    // swap costs
                    double tmpC = costs[i];
                    costs[i] = costs[bestIdx];
                    costs[bestIdx] = tmpC;

                    // swap vertices
                    int offI = i * n;
                    int offB = bestIdx * n;
                    for (int k = 0; k < n; ++k)
                    {
                        double tmp = simplex[offI + k];
                        simplex[offI + k] = simplex[offB + k];
                        simplex[offB + k] = tmp;
                    }
                }
            }
        }
    }
}
