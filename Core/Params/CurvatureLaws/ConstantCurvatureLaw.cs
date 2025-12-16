namespace ArcFrame.Core.Params
{
    /// <summary>
    /// Curvature law where the curvature does not change with respect to arc length.
    /// Use this for lines and arcs.
    /// </summary>
    public class ConstantCurvatureLaw : IParamCurvatureLaw
    {
        private readonly double[] _k0;
        /// <inheritdoc/>
        public int Order => _k0.Length;
        /// <inheritdoc/>
        public bool IsLinear => false;
        /// <inheritdoc/>
        public bool IsConstant => true;
        /// <summary>
        /// Create a constant curvature law
        /// </summary>
        /// <param name="k0"></param>
        public ConstantCurvatureLaw(double[] k0)
        {
            _k0 = (double[])k0.Clone();
        }
        /// <inheritdoc/>
        public double[] Eval(double s)
        {
            return (double[])_k0.Clone();
        }
        /// <inheritdoc/>
        public double[] GetParams()
        {
            return (double[])_k0.Clone();
        }
        /// <inheritdoc/>
        public void SetParams(double[] p)
        {
            for (int i = 0; i < p.Length; i++)
            {
                _k0[i] = p[i];
            }
        }
        /// <inheritdoc/>
        public ICurvatureLaw CloneWithParams(double[] p)
        {
            return new ConstantCurvatureLaw(p);
        }
    }
}
