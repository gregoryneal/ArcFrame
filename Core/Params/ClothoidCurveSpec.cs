using ArcFrame.Core.Math;

namespace ArcFrame.Core.Params
{
    /// <summary>
    /// Helper class for defining clothoid specific curve specs. Basically constrains ICurvatureLaw to be a LinearCurvatureLaw.
    /// </summary>
    public class ClothoidCurveSpec : CurveSpec
    {
        public ClothoidCurveSpec(int n, double length, double[] p0, double[,] r0, LinearCurvatureLaw kappa, FrameModel frame) : base(n, length, p0, r0, kappa, frame) { }
    }
}
