using ArcFrame.Core.Math;
using ArcFrame.Core.Math.Geometry.Splines;
using System;

namespace ArcFrame.Core
{
    /// <summary>
    /// Con
    /// </summary>
    public static class Spline3dDescriptorFactory
    {
        /// <summary>
        /// Build a descriptor from a 3d spline. 
        /// </summary>
        /// <param name="spline"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static Spline3dDescriptor FromSpline(Spline spline)
        {
            if (spline.Dimension != 3)
                throw new ArgumentException("Spline3dDescriptor requires Dimension == 3");

            var desc = new Spline3dDescriptor
            {
                Degree = spline.Degree,
                ControlPointCount = spline.ControlPoints.Length,
                // ControlPointCount - Degree
                SegmentCount = spline.ComputeSegmentCount(),
                Length = spline.Length,
                Frame = spline.Frame,
                FastMode = spline.FastMode,
                // Each segment evaluation uses Degree+1 local control points
                Stride = spline.Degree + 1,
                // Flatten control points [x0,y0,z0, x1,y1,z1, ...]
                ControlPointsFlat = FlattenControlPoints(spline.ControlPoints),
                // Flatten BasisMatrix (row-major)
                BasisFlat = FlattenBasisMatrix(spline.BasisMatrix),
                // Use precomputed, cached GB blocks
                SegmentCoeff = FlattenGBCache(spline.GBCache, spline.Dimension, spline.Degree),
                // Optional constant R0 in FastMode
                R0Flat = Flatten3x3(spline.R0)
            };

            // Extract arc-length mapping table
            spline.Table.Export(out desc.ArcS, out desc.ArcT);

            return desc;
        }

        /// <summary>
        /// Create an array of flattened control points: [x0, y0, z0, ..., x1, y1, z1, ..., xn, yn, zn, ...]
        /// </summary>
        /// <returns></returns>
        public static double[] FlattenControlPoints(double[][] ControlPoints)
        {
            int n = ControlPoints.Length;
            if (n <= 0) return new double[] { };
            int m = ControlPoints[0].Length;
            if (m <= 0) return new double[] { };

            int len = n * m;
            double[] pts = new double[len];
            int ptsInd = 0;
            // iterate through the control points (length n)
            for (int i = 0; i < n; i++)
            {
                // find the start index of the control point in pts
                ptsInd = i * m;
                for (int j = 0; j < m; j++)
                {
                    pts[ptsInd + j] = ControlPoints[i][j];
                }
            }

            return pts;
        }

        /// <summary>
        /// Flatten our precomputed _gbCache into a single row-major array.
        ///
        /// Layout (for Dimension == 3):
        ///   SegmentCoeff[seg * (3*(deg+1)) + row*(deg+1) + col] = GB[seg][row, col]
        ///
        /// For general N:
        ///   len = SegmentCount * Dimension * (Degree+1)
        /// </summary>
        public static double[] FlattenGBCache(double[][,] _gbCache, int Dimension, int Degree)
        {
            if (_gbCache == null || _gbCache.Length == 0)
                return Array.Empty<double>();

            int segmentCount = _gbCache.Length;
            int dim = Dimension;
            int cols = Degree + 1;
            int segmentStride = dim * cols;

            var flat = new double[segmentCount * segmentStride];

            for (int seg = 0; seg < segmentCount; seg++)
            {
                var gb = _gbCache[seg];
                if (gb == null) continue; // should not happen, but be defensive

                int segOffset = seg * segmentStride;

                for (int r = 0; r < dim; r++)
                {
                    int rowOffset = segOffset + r * cols;
                    for (int c = 0; c < cols; c++)
                    {
                        flat[rowOffset + c] = gb[r, c];
                    }
                }
            }

            return flat;
        }

        private static double[] FlattenBasisMatrix(double[,] basis)
        {
            if (basis == null)
                return Array.Empty<double>();

            int rows = basis.GetLength(0);
            int cols = basis.GetLength(1);
            int len = rows * cols;

            var flat = new double[len];
            int idx = 0;

            for (int r = 0; r < rows; r++)
            {
                for (int c = 0; c < cols; c++)
                {
                    flat[idx++] = basis[r, c];
                }
            }

            return flat;
        }

        private static double[] Flatten3x3(double[,] R)
        {
            if (R == null)
                return Array.Empty<double>();

            if (R.GetLength(0) != 3 || R.GetLength(1) != 3)
                throw new ArgumentException("Flatten3x3 expects a 3x3 matrix.", nameof(R));

            var flat = new double[9];
            int idx = 0;

            for (int r = 0; r < 3; r++)
            {
                for (int c = 0; c < 3; c++)
                {
                    flat[idx++] = R[r, c];
                }
            }

            return flat;
        }
    }

    /// <summary>
    /// A description of a 3D spline with enough information to rebuild one in ArcFrame library.
    /// </summary>
    public sealed class Spline3dDescriptor
    {
        /// <summary>
        /// Spline degree
        /// </summary>
        public int Degree;
        /// <summary>
        /// The segment count, ControlPoint.Length - Degree
        /// </summary>
        public int SegmentCount;
        /// <summary>
        /// The computed length of the spline
        /// </summary>
        public double Length;
        /// <summary>
        /// Which frame type to use
        /// </summary>
        public FrameModel Frame;
        /// <summary>
        /// FastMode computes only the spline position and tangent at each point, no curvature data.
        /// </summary>
        public bool FastMode;

        /// <summary>
        /// How many control points are in our evaluation matrix
        /// </summary>
        public int Stride;

        /// <summary>
        /// Flattened control point array: [x0,y0,z0, x1,y1,z1, ...]
        /// Length: 3 * control point length
        /// </summary>
        public double[] ControlPointsFlat;

        /// <summary>
        /// The number of control points, depending on the type of spline, this might include points that the spline does not pass through. 
        /// </summary>
        public int ControlPointCount;

        /// <summary>
        /// Basis matrix: (deg+1)x(deg+1), row-major
        /// </summary>
        public double[] BasisFlat;

        /// <summary>
        /// Precomputed per-segment GB = G * BasisMatrix (3 x (deg+1))
        /// Flattened per segment: [seg0 3x(d+1), seg1 3x(d+1), ...] 
        /// len = SegmentCount * 3 * (Degree+1)
        /// </summary>
        public double[] SegmentCoeff;

        /// <summary>
        /// Arc-length mapping table: s_i ↔ t_i
        /// </summary>
        public double[] ArcS;                     // length = ArcSampleCount
        /// <summary>
        /// Arc-length mapping table: s_i ↔ t_i
        /// </summary>
        public double[] ArcT;                     // length = ArcSampleCount

        /// <summary>
        /// Constant frame in FastMode
        /// </summary>
        public double[]? R0Flat;                  // 3x3 row-major, or null if not used
    }
}
