using ArcFrame.Core.Math.Geometry.Splines;

namespace ArcFrame.Core.Params.CurvatureLaws
{
    /*
    /// <summary>Base utility for Bezier curvature laws (value-only by default).</summary>
    public class BezierCurvatureLaw : ICurvatureLaw, ICurvatureJetProvider
    {
        protected readonly BezierCurve Curve;
        public int Order => System.Math.Max(0, Curve.Dimension - 1);
        public virtual int MaxSupportedJetOrder => Curve.Degree; // value-only
        public virtual bool IsLinear => false;
        public virtual bool IsConstant => false;

        public BezierCurvatureLaw(BezierCurve curve) => Curve = curve;

        public double[] Eval(double s)
        {
            s = System.Math.Clamp(s, 0.0, Curve.Length);
            double t = Curve.Length > 0 ? s / Curve.Length : 0.0;
            (int seg, double u) = Curve.Locate(t);
            // compute derivatives and build k vector [k1, k2, k3, ...] length = N-1
            // TODO: Test the N dimensional 
            Curve.ComputeSegmentJet(seg, u, MaxSupportedJetOrder, out _, out var jet);
            return Curve.GeneralizedCurvaturesFromJet(jet);
        }

        public double[][] EvalJet(double s, int jetOrder)
        {
            s = System.Math.Clamp(s, 0.0, Curve.Length);
            double t = Curve.Length > 0 ? s / Curve.Length : 0.0;
            (int seg, double u) = Curve.Locate(t);
            jetOrder = System.Math.Min(jetOrder, MaxSupportedJetOrder);
            Curve.ComputeSegmentJet(seg, u, jetOrder, out double[] P, out double[][] jets);
            var outJet = new double[jetOrder + 1][];
            outJet[0] = P;
            for (int i = 0; i < jets.Length; i++)
            {
                outJet[i + 1] = jets[i];
            }
            return outJet;
        }
    }
    */
}
