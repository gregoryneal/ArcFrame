using System;

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
        /// <inheritdoc/>
        public int Order => _o;
        /// <inheritdoc/>
        public bool IsLinear => false;
        /// <inheritdoc/>
        public bool IsConstant => false;
        /// <summary>
        /// Create a curvature law with a predefined function for the curvatures.
        /// </summary>
        /// <param name="order"></param>
        /// <param name="kappa"></param>
        public FunctionCurvatureLaw(int order, Func<double, double[]> kappa)
        {
            _o = order;
            _func = kappa;
        }
        /// <inheritdoc/>
        public double[] Eval(double s)
        {
            return _func(s);
        }
    }
}
