using ArcFrame.Core.Math;
using ArcFrame.Core.Math.Geometry.Splines;

namespace ArcFrame.Core.Geometry.Splines
{
    public class CatmullRomSpline : Spline
    {
        public readonly double Tension;

        public CatmullRomSpline(double[][] controlPoints, double tension = 0.5, FrameModel frame = FrameModel.Frenet) : base(controlPoints, ConstructBasis(tension), frame)
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