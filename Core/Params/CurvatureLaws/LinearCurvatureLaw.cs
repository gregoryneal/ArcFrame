using ArcFrame.Core.Params.CurvatureLaws;

namespace ArcFrame.Core.Params
{
    /// <summary>
    /// Curvature law where the curvature of the curve varies linearly with respect to arc length.
    /// k = k0 + s(dk)
    /// </summary>
    public sealed class LinearCurvatureLaw : ICurvatureLaw
    {
        private readonly double[] _k0, _dk; 
        public int Order => _k0.Length;
        public double[] Intercept => _k0;
        public double[] Slope => _dk;

        public LinearCurvatureLaw(double[] k0, double[] dk)
        { 
            _k0 = (double[])k0.Clone(); 
            _dk = (double[])dk.Clone(); 
        }
        public double[] Eval(double s)
        {
            double[] outK = new double[_k0.Length];
            for (int i = 0; i < _k0.Length; i++) outK[i] = _k0[i] + _dk[i] * s;
            return outK;
        }
    }
}
