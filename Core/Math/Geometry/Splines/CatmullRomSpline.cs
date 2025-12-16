using ArcFrame.Core.Math;
using ArcFrame.Core.Math.Geometry.Splines;

namespace ArcFrame.Core.Geometry.Splines
{
    /// <summary>
    /// Class for a Catmull-Rom spline with adjustable tension.
    /// </summary>
    public class CatmullRomSpline : Spline
    {
        //public override int SegmentCount => ControlPoints.Length - Degree;
        /// <summary>
        /// Controls the tangent at the control points. 
        /// Higher values -> flatter tangents at the control points
        /// which might cause the curve to make large corrective bows.
        /// Lower values -> allows the spline to make sharper bends.
        /// </summary>
        public readonly double Tension;

        /// <summary>
        /// Create the Catmull-Rom spline with control points and tension.
        /// </summary>
        /// <param name="controlPoints"></param>
        /// <param name="tension"></param>
        /// <param name="frame"></param>
        /// <param name="fastMode"></param>
        /// <param name="cacheSamplesOverride"></param>
        public CatmullRomSpline(double[][] controlPoints, 
            double tension = 0.5, 
            FrameModel frame = FrameModel.Frenet,
            bool fastMode = false,
            int cacheSamplesOverride = 0) : base(controlPoints, ConstructBasis(tension), frame, fastMode, cacheSamplesOverride)
        { Tension = tension; }


        private static double[,] ConstructBasis(double alpha)
        {
            return new double[,]
            {
                { 0,        -alpha,     2*alpha,   -alpha },
                { 1,         0,         alpha-3,    2-alpha },
                { 0,         alpha,     3-2*alpha,  alpha-2 },
                { 0,         0,        -alpha,      alpha }
            };
        }
    }
}