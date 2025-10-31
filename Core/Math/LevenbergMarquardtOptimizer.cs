namespace ArcFrame.Core.Math
{
    /// <summary>
    /// The result of the optimization routine.
    /// </summary>
    public class OptimizationResult
    {
        /// <summary>
        /// Final parameters of the optimizer
        /// </summary>
        public double[] FinalParameters { get; set; } = { };
        /// <summary>
        /// Final cost of the optimization routine.
        /// </summary>
        public double FinalCost { get; set; } = 0;
        /// <summary>
        /// How many iterations the routine took.
        /// </summary>
        public int Iterations { get; set; } = 0;
        /// <summary>
        /// Was the final error within tolerance.
        /// </summary>
        public bool Success { get; set; } = false;
        /// <summary>
        /// The optimization message, this will tell you why a result failed.
        /// </summary>
        public string Message { get; set; } = "";
    }

    /// <summary>
    /// Create an optimizer that utilizes the LMA for least squares optimization.
    /// </summary>
    public class LevenbergMarquardtOptimizer
    {
        /// <summary>
        /// Delegate for the function that calculates the residuals.
        /// Takes the current parameters and returns an array of residuals.
        /// </summary>
        /// <param name="parameters"></param>
        /// <returns></returns>
        public delegate double[] ResidualFunction(double[] parameters);

        private readonly ResidualFunction _residualFunction;
        private readonly int _numParameters;
        private readonly int _numResiduals;

        /// <summary>
        /// Create new instance of the optimizer with the residual.
        /// </summary>
        /// <param name="residualFunction"></param>
        /// <param name="numParameters"></param>
        /// <param name="numResiduals"></param>
        public LevenbergMarquardtOptimizer(ResidualFunction residualFunction, int numParameters, int numResiduals)
        {
            _residualFunction = residualFunction;
            _numParameters = numParameters;
            _numResiduals = numResiduals;
        }

        /// <summary>
        /// Calculates the Jacobian matrix using central finite differences.
        /// J[i, j] = (r_i(p_j + h) - r_i(p_j - h)) / (2h)
        /// </summary>
        public double[,] CalculateJacobian(double[] parameters, double stepSize = 1e-6)
        {
            double[,] jacobian = new double[_numResiduals, _numParameters];
            double[] baseResiduals = _residualFunction(parameters);

            for (int j = 0; j < _numParameters; j++)
            {
                double[] paramsPlus = (double[])parameters.Clone();
                paramsPlus[j] += stepSize;
                double[] residualsPlus = _residualFunction(paramsPlus);

                double[] paramsMinus = (double[])parameters.Clone();
                paramsMinus[j] -= stepSize;
                double[] residualsMinus = _residualFunction(paramsMinus);

                for (int i = 0; i < _numResiduals; i++)
                {
                    jacobian[i, j] = (residualsPlus[i] - residualsMinus[i]) / (2 * stepSize);
                }
            }
            return jacobian;
        }

        /// <summary>
        /// Calculates the Gauss-Newton approximation of the Hessian matrix.
        /// This is simply the transpose of the Jacobian multiplied by the Jacobian itself (JᵀJ).
        /// </summary>
        public double[,] CalculateHessianApproximation(double[,] jacobian)
        {
            double[,] jacobianTranspose = Helpers.Transpose(jacobian);
            return Helpers.Multiply(jacobianTranspose, jacobian);
        }

        /// <summary>
        /// Calculates the sum of squares of the residuals (the cost function).
        /// </summary>
        private double CalculateCost(double[] residuals)
        {
            return Helpers.Dot(residuals, residuals);
        }

        /// <summary>
        /// Performs the Levenberg-Marquardt optimization.
        /// </summary>
        public OptimizationResult Optimize(double[] initialParameters)
        {
            double[] currentParameters = (double[])initialParameters.Clone();
            double[] currentResiduals = _residualFunction(currentParameters);
            double currentCost = CalculateCost(currentResiduals);

            double lambda = 1e-3; // Initial damping factor
            double lambdaFactor = 10.0;
            double tolerance = 1e-8;
            int maxIterations = 100;

            for (int iter = 0; iter < maxIterations; iter++)
            {
                // 1. Calculate Jacobian and Hessian approximation
                double[,] jacobian = CalculateJacobian(currentParameters);
                double[,] jacobianTranspose = Helpers.Transpose(jacobian);
                double[,] hessianApprox = Helpers.Multiply(jacobianTranspose, jacobian);

                // 2. Calculate the gradient (Jᵀr)
                double[] gradient = Helpers.Multiply(jacobianTranspose, currentResiduals);

                // 3. Form the augmented linear system: (H + λI) * delta_p = -gradient
                double[,] identity = RigidTransform.Identity(_numParameters).R;
                double[,] aMatrix = hessianApprox; // H
                for (int i = 0; i < _numParameters; i++)
                {
                    aMatrix[i, i] += lambda * identity[i, i]; // H + λI
                }

                double[] bVector = new double[_numParameters];
                for (int i = 0; i < _numParameters; i++)
                {
                    bVector[i] = -gradient[i]; // -gradient
                }

                // 4. Solve for the step delta_p
                double[] deltaP;
                try
                {
                    //solve linear equation
                    var (LU, piv, _) = LA.LUFactor(aMatrix);
                    deltaP = LA.LUSolve(LU, piv, bVector);
                    //deltaP = LA.S(aMatrix, bVector);
                }
                catch (System.Exception ex)
                {
                    return new OptimizationResult
                    {
                        FinalParameters = currentParameters,
                        FinalCost = currentCost,
                        Iterations = iter,
                        Success = false,
                        Message = $"Failed to solve linear system: {ex.Message}"
                    };
                }

                // 5. Calculate new parameters and new cost
                double[] newParameters = Helpers.Add(currentParameters, deltaP);
                double[] newResiduals = _residualFunction(newParameters);
                double newCost = CalculateCost(newResiduals);

                // 6. Decide whether to accept the step
                double costChange = currentCost - newCost;
                if (costChange > 0) // Good step, cost decreased
                {
                    currentParameters = newParameters;
                    currentResiduals = newResiduals;
                    currentCost = newCost;
                    lambda /= lambdaFactor; // Reduce damping, trust the model more

                    // Check for convergence
                    if (costChange < tolerance)
                    {
                        return new OptimizationResult
                        {
                            FinalParameters = currentParameters,
                            FinalCost = currentCost,
                            Iterations = iter + 1,
                            Success = true,
                            Message = "Converged based on cost change."
                        };
                    }
                }
                else // Bad step, cost increased
                {
                    lambda *= lambdaFactor; // Increase damping, act more like gradient descent
                }
            }

            return new OptimizationResult
            {
                FinalParameters = currentParameters,
                FinalCost = currentCost,
                Iterations = maxIterations,
                Success = false,
                Message = "Reached maximum iterations without converging."
            };
        }
    }
}
