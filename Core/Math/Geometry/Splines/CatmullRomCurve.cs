using ArcFrame.Core.Math;
using ArcFrame.Core.Math.Geometry.Splines;
using ArcFrame.Core.Results;

namespace ArcFrame.Core.Geometry.Splines
{
    /*
    /// <summary>
    /// Cubic Catmull-Rom curve as a Bezier-like with a 4x4 basis matrix.
    /// </summary>
    public sealed class CatmullRomCurve : BezierCurve
    {
        // Uniform Catmull–Rom in monomial basis:
        private static readonly double[,] M = new double[,]
        {
            {  0.0,   1.0,   0.0,  0.0 },
            { -0.5,   0.0,   0.5,  0.0 },
            {  1.0,  -2.5,   2.0, -0.5 },
            { -0.5,   1.5,  -1.5,  0.5 }
        };

        public CatmullRomCurve(IReadOnlyList<double[]> points, FrameModel frame = FrameModel.Bishop, double[]? upHint = null, int samples = 1024)
            : base(points, degree: 3, basisMatrix: M, frame: frame, upHint: upHint, samples: samples) { }
    }

    /// <summary>
    /// "cubic spline" placeholder — supply your preferred 4x4 basis matrix.
    /// </summary>
    public sealed class CubicSplineCurve : BezierCurve
    {
        public CubicSplineCurve(IReadOnlyList<double[]> points, double[,] splineBasis4x4,
                                FrameModel frame = FrameModel.Bishop, double[]? upHint = null, int samples = 1024)
            : base(points, degree: 3, basisMatrix: splineBasis4x4, frame: frame, upHint: upHint, samples: samples) { }
    }
    */
}