namespace ArcFrame.Core.Params
{
    /// <summary>
    /// Curvature law where the generalized curvatures follow some function.
    /// </summary>
    public sealed class FunctionCurvatureLaw : ICurvatureLaw
    {
        readonly Func<double, double[]> _func;
        readonly int _o;
        public int Order => _o;

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
