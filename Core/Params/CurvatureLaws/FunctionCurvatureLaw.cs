namespace ArcFrame.Core.Params
{
    /// <summary>
    /// Curvature law where the generalized curvatures follow some function.
    /// Implent ICurvatureJetProvider in inherited classes to provide curvature jets. 
    /// </summary>
    public sealed class FunctionCurvatureLaw : ICurvatureLaw
    {
        readonly Func<double, double[]> _func;
        readonly int _o;
        public int Order => _o;

        public bool IsLinear => false;

        public bool IsConstant => false;

        public FunctionCurvatureLaw(int order, Func<double, double[]> kappa)
        {
            _o = order;
            _func = kappa;
        }

        public double[] Eval(double s)
        {
            return _func(s);
        }
    }
}
