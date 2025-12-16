
using ArcFrame.Core.Math;
using System;
using System.Runtime.CompilerServices;

namespace ArcFrame.Core.Kernels
{
    /// <summary>
    /// Doing math with unmanaged data and pointers. 
    /// </summary>
    public static unsafe class HpcMathNd
    {
        /// <summary>
        /// Nd dot product
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Dot(int n, double* a, double* b)
        {
            if (n <= 0) return 0;
            double s = 0.0;
            for (int i = 0; i < n; i++)
                s += a[i] * b[i];
            return s;
        }

        /// <summary>
        /// Nd unsafe addition
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="result"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void Add(int n, double* a, double* b, double* result)
        {
            if (n <= 0) return;
            for (int i = 0; i < n; i++)
                result[i] = a[i] + b[i];
        }

        /// <summary>
        /// Nd pointer subtraction
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="result"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void Sub(int n, double* a, double* b, double* result)
        {
            if (n <= 0) return;
            for (int i = 0; i < n; i++)
                result[i] = a[i] - b[i];
        }

        /// <summary>
        /// Nd vector scaling
        /// </summary>
        /// <param name="n"></param>
        /// <param name="s"></param>
        /// <param name="a"></param>
        /// <param name="result"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void Scale(int n, double s, double* a, double* result)
        {
            if (n <= 0) return;
            for (int i = 0; i < n; i++)
                result[i] = a[i] * s;
        }

        /// <summary>
        /// Multiply an n x m matrix with an m x o matrix
        /// Store result in n x o matrix
        /// </summary>
        /// <param name="n">num rows in M1</param>
        /// <param name="m">num cols in M1 and rows in M2</param>
        /// <param name="o">num cols in M2</param>
        /// <param name="M1"></param>
        /// <param name="M2"></param>
        /// <param name="result"></param>
        /// <exception cref="ArgumentException"></exception>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void MatMul(int n, int m, int o, double* M1, double* M2, double* result)
        {
            // M1: n x m, row-major
            // M2: m x o, row-major
            // result: n x o, row-major

            for (int i = 0; i < n; i++)
            {
                double* rowPtrM1 = M1 + i * m; // row i in M1

                for (int j = 0; j < o; j++)
                {
                    double sum = 0.0;

                    for (int k = 0; k < m; k++)
                    {
                        double a = rowPtrM1[k];          // M1[i, k]
                        double b = M2[k * o + j];        // M2[k, j]
                        sum += a * b;
                    }

                    result[i * o + j] = sum;             // result[i, j]
                }
            }
        }

        /// <summary>
        /// Row-major n x m matrix pointer, multiplied with m length vector pointer
        /// Note: always supply a fresh matrix with length numRows as the result, don't try
        /// to apply the result to the input vector. 
        /// </summary>
        /// <param name="numRows">number of rows in the matrix</param>
        /// <param name="numCols">number of columns in the matrix, and number of entries in the vector</param>
        /// <param name="M">pointer to first element of matrix</param>
        /// <param name="v">pointer to vector</param>
        /// <param name="result"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void MatVec(int numRows, int numCols, double* M, double* v, double* result)
        {
            // Bounds checking: Validate dimensions
            if (numRows <= 0 || numCols <= 0 || v == null || result == null)
            {
                throw new ArgumentException("Invalid matrix dimensions or null inputs.");
            }

            // Ensure the vector has the correct number of elements (must match numCols)
            if (v == null || numCols <= 0)
            {
                throw new ArgumentException("The vector must have a length equal to the number of columns in the matrix.");
            }

            // Ensure the result array has the correct length
            if (result == null)
            {
                throw new ArgumentNullException(nameof(result));
            }

            // Ensure the result array is large enough to hold the result (numRows elements)
            for (int i = 0; i < numRows; i++)
            {
                result[i] = 0;  // Initialize the result array to 0 (or any other default value)
            }

            // Pointer-based matrix-vector multiplication
            for (int row = 0; row < numRows; row++)
            {
                double sum = 0.0;
                double* rowPtr = M + (row * numCols);  // Pointer to the start of the current row in M

                // Vector v is accessed element-by-element using the index
                for (int col = 0; col < numCols; col++)
                {
                    sum += *(rowPtr + col) * *(v + col);  // Pointer arithmetic for matrix and vector
                }

                *(result + row) = sum;  // Store the result in the output vector
            }
        }

        /// <summary>
        /// Project vector v onto the vector u
        /// </summary>
        /// <param name="n"></param>
        /// <param name="v"></param>
        /// <param name="u"></param>
        /// <param name="result"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void Project(int n, double* v, double* u, double* result)
        {
            if (result == null || v == null || u == null || n <= 0) return;
            var denom = Dot(n, u, u); // len squared
            if (denom <= 0.0) return;
            var s = Dot(n, v, u) / denom;
            Scale(n, s, u, result);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void Lerp(int n, double t, double* a, double* b, double* result)
        {
            if (n <= 0) return;
            for (int i = 0; i < n; i++)
            {
                result[i] = a[i] + (b[i] - a[i]) * t;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double OneNorm(int rows, int cols, double* data, bool rowMajor)
        {
            if (rowMajor)
                return OneNormRowMajor(rows, cols, data);
            else
                return OneNormColMajor(rows, cols, data);
        }

        private static double OneNormRowMajor(int rows, int cols, double* data) // length = rows * cols, row-major
        {
            double maxColSum = 0.0;

            for (int col = 0; col < cols; col++)
            {
                double sum = 0.0;
                int idx = col; // starting at (row = 0, col)

                for (int row = 0; row < rows; row++)
                {
                    sum += System.Math.Abs(data[idx]);
                    idx += cols; // move to next row, same column: row*cols + col
                }

                if (sum > maxColSum)
                    maxColSum = sum;
            }

            return maxColSum;
        }

        private static double OneNormColMajor(int rows, int cols, double* data) // length = rows * cols, column-major
        {
            double maxColSum = 0.0;

            for (int col = 0; col < cols; col++)
            {
                double sum = 0.0;
                int idx = col * rows; // first element of the column

                for (int row = 0; row < rows; row++)
                {
                    sum += System.Math.Abs(data[idx + row]);
                }

                if (sum > maxColSum)
                    maxColSum = sum;
            }

            return maxColSum;
        }
    }
}
