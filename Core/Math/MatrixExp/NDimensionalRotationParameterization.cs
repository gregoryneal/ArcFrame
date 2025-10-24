namespace ArcFrame.Core.Math
{
    using System;

    /// <summary>
    /// Current constraint parameterization assumes a 2D rotation in ND space, this would 
    /// parameterize the entire SO(n) frame at once using exp(SO(n))
    /// </summary>
    public static class NDimensionalRotationParameterization
    {
        private static int _dimension;
        private static int _paramCount;
        private readonly static double _epsilon = 1e-8;
        private readonly static int _maxSeriesTerms = 50;

        /// <summary>
        /// Initializes the parameterization for a specific dimension.
        /// This must be called once before using Pack/Unpack.
        /// </summary>
        /// <param name="dimension">The dimension of the rotation space (N).</param>
        public static void Initialize(int dimension)
        {
            if (dimension < 2)
                throw new ArgumentException("Dimension must be 2 or greater.");
            _dimension = dimension;
            _paramCount = dimension * (dimension - 1) / 2;
        }

        /// <summary>
        /// Packs an N x N rotation matrix into a parameter vector of size N*(N-1)/2.
        /// This is the matrix logarithm, mapping from SO(N) to its tangent space.
        /// </summary>
        public static double[] Pack(double[,] rotationMatrix)
        {
            CheckInitialized();
            if (rotationMatrix.GetLength(0) != _dimension || rotationMatrix.GetLength(1) != _dimension)
                throw new ArgumentException($"Input matrix must be {_dimension}x{_dimension}.");

            // The matrix logarithm log(R) gives a skew-symmetric matrix A.
            // For small steps, log(I + E) ≈ E - E^2/2 + E^3/3 - ...
            // where E = R - I
            double[,] identity = MatrixIdentity(_dimension);
            double[,] E = MatrixSubtract(rotationMatrix, identity);
            double[,] A = MatrixZero(_dimension); // This will be our skew-symmetric matrix

            double[,] term = E.Clone() as double[,];
            double sign = 1.0;
            for (int k = 1; k <= _maxSeriesTerms; k++)
            {
                double factor = sign / k;
                A = MatrixAdd(A, MatrixScale(term, factor));

                // Check for convergence
                if (MatrixNorm(term) < _epsilon) break;

                term = MatrixMultiply(term, E);
                sign *= -1.0;
            }

            // Ensure the result is perfectly skew-symmetric
            for (int i = 0; i < _dimension; i++)
            {
                for (int j = 0; j < _dimension; j++)
                {
                    A[i, j] = 0.5 * (A[i, j] - A[j, i]);
                }
            }

            return SkewMatrixToVector(A);
        }

        /// <summary>
        /// Unpacks a parameter vector of size N*(N-1)/2 into an N x N rotation matrix.
        /// This is the matrix exponential, mapping from the tangent space to SO(N).
        /// </summary>
        public static double[,] Unpack(double[] parameters)
        {
            CheckInitialized();
            if (parameters.Length != _paramCount)
                throw new ArgumentException($"Parameter array must have length {_paramCount}.");

            // Convert the parameter vector to a skew-symmetric matrix A
            double[,] A = VectorToSkewMatrix(parameters);

            // The matrix exponential exp(A) gives the rotation matrix R.
            // R = I + A + A^2/2! + A^3/3! + ...
            double[,] R = MatrixIdentity(_dimension);
            double[,] term = MatrixIdentity(_dimension);
            double factorial = 1.0;

            for (int k = 1; k <= _maxSeriesTerms; k++)
            {
                factorial *= k;
                term = MatrixMultiply(term, A);
                R = MatrixAdd(R, MatrixScale(term, 1.0 / factorial));

                // Check for convergence
                if (MatrixNorm(term) / factorial < _epsilon) break;
            }

            return R;
        }

        #region Helper Methods

        private static void CheckInitialized()
        {
            if (_dimension == 0)
                throw new InvalidOperationException("NDimensionalRotationParameterization must be initialized first by calling Initialize(N).");
        }

        private static double[,] MatrixZero(int size)
        {
            return new double[size, size];
        }

        private static double[,] MatrixIdentity(int size)
        {
            double[,] I = new double[size, size];
            for (int i = 0; i < size; i++) I[i, i] = 1.0;
            return I;
        }

        private static double[,] MatrixAdd(double[,] A, double[,] B)
        {
            int size = A.GetLength(0);
            double[,] C = new double[size, size];
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    C[i, j] = A[i, j] + B[i, j];
            return C;
        }

        private static double[,] MatrixSubtract(double[,] A, double[,] B)
        {
            int size = A.GetLength(0);
            double[,] C = new double[size, size];
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    C[i, j] = A[i, j] - B[i, j];
            return C;
        }

        private static double[,] MatrixScale(double[,] M, double s)
        {
            int size = M.GetLength(0);
            double[,] C = new double[size, size];
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    C[i, j] = M[i, j] * s;
            return C;
        }

        private static double[,] MatrixMultiply(double[,] A, double[,] B)
        {
            int size = A.GetLength(0);
            double[,] C = new double[size, size];
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    for (int k = 0; k < size; k++)
                        C[i, j] += A[i, k] * B[k, j];
            return C;
        }

        private static double MatrixNorm(double[,] M)
        {
            // Frobenius norm
            double sum = 0;
            int size = M.GetLength(0);
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    sum += M[i, j] * M[i, j];
            return Math.Sqrt(sum);
        }

        private static double[,] VectorToSkewMatrix(double[] parameters)
        {
            double[,] A = new double[_dimension, _dimension];
            int paramIndex = 0;
            for (int i = 0; i < _dimension; i++)
            {
                for (int j = i + 1; j < _dimension; j++)
                {
                    A[i, j] = parameters[paramIndex];
                    A[j, i] = -parameters[paramIndex];
                    paramIndex++;
                }
            }
            return A;
        }

        private static double[] SkewMatrixToVector(double[,] A)
        {
            double[] parameters = new double[_paramCount];
            int paramIndex = 0;
            for (int i = 0; i < _dimension; i++)
            {
                for (int j = i + 1; j < _dimension; j++)
                {
                    parameters[paramIndex] = A[i, j];
                    paramIndex++;
                }
            }
            return parameters;
        }

        #endregion
    }
}
