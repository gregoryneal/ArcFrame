using ArcFrame.Core.Geometry.Splines;

namespace ArcFrame.Core.Params.CurvatureLaws
{
    /// <summary>
    /// Create a curvature law to generate cubic splines.
    /// Given a list of knots and 
    /// </summary>
    public sealed class SplineCurvatureLaw : ICurvatureLaw, ICurvatureJetProvider
    {
        readonly CubicSpline1D[] comps; // one per curvature component
        readonly double sMin, sMax;
        public int Order => comps.Length;

        public int MaxSupportedJetOrder => 3;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="sKnots"></param>
        /// <param name="kSamples"></param>
        /// <exception cref="ArgumentException"></exception>
        public SplineCurvatureLaw(double[] sKnots, double[][] kSamples)
        {
            // kSamples[j][i] = component j at sKnots[i]
            int m = kSamples.Length;
            int n = sKnots.Length;
            if (m < 1 || n < 2) throw new ArgumentException();
            comps = new CubicSpline1D[m];
            for (int j = 0; j < m; j++) comps[j] = new CubicSpline1D(sKnots, kSamples[j]);
            sMin = sKnots[0]; 
            sMax = sKnots[^1];
        }
        public double[] Eval(double s)
        {
            double[] outK = new double[Order];
            s = System.Math.Max(sMin, System.Math.Min(sMax, s));
            for (int j = 0; j < comps.Length; j++) outK[j] = comps[j].Eval(s);
            return outK;
        }

        public double[][] EvalJet(double s, int order)
        {
            double[][] outJet = new double[1][];
            order = System.Math.Min(order, 3);
            for (int j = 0; j < comps.Length; j++)
            {
                comps[j].EvalJet(s, order, out double v0, out double v1, out double v2, out double v3);
                outJet[0][j] = v0;
                if (order >= 1) outJet[1][j] = v1;
                if (order >= 2) outJet[2][j] = v2;
                if (order >= 3) outJet[3][j] = v3;
            }
            return outJet;
        }
    }

}
