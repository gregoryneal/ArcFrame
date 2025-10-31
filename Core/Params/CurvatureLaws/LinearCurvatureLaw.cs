namespace ArcFrame.Core.Params
{
    /// <summary>
    /// Curvature law where the curvature of the curve varies linearly with respect to arc length.
    /// k = k0 + s(dk)
    /// </summary>
    public sealed class LinearCurvatureLaw : IParamCurvatureLaw
    {
        private readonly double[] _k0, _dk;
        /// <inheritdoc/>
        public int Order => _k0.Length;
        /// <summary>
        /// The initial curvature vector.
        /// </summary>
        public double[] Intercept => _k0;
        /// <summary>
        /// The delta curvatures
        /// </summary>
        public double[] Slope => _dk;
        /// <inheritdoc/>
        public bool IsLinear => true;
        /// <inheritdoc/>
        public bool IsConstant => false;
        /// <summary>
        /// Create a curvature law where the curvature varies linearly with respect to arc length.
        /// </summary>
        /// <param name="k0"></param>
        /// <param name="dk"></param>
        public LinearCurvatureLaw(double[] k0, double[] dk)
        {
            _k0 = (double[])k0.Clone();
            _dk = (double[])dk.Clone();
        }
        /// <inheritdoc/>
        public double[] Eval(double s)
        {
            double[] outK = new double[_k0.Length];
            for (int i = 0; i < _k0.Length; i++) outK[i] = _k0[i] + _dk[i] * s;
            return outK;
        }
        /// <inheritdoc/>
        public double[] GetParams()
        {
            var p = new double[2 * Order];
            for (int i = 0; i < Order; i++) p[i] = Intercept[i];
            for (int i = 0; i < Order; i++) p[Order + i] = Slope[i];
            return p;
        }
        /// <inheritdoc/>
        public void SetParams(double[] p)
        {
            int m = Order;
            for (int i = 0; i < m; i++) _k0[i] = p[i];
            for (int i = 0; i < m; i++) _dk[i] = p[m + i];
        }
        /// <inheritdoc/>
        public ICurvatureLaw CloneWithParams(double[] p)
        {
            int m = Order;
            var k0 = new double[m];
            var dk = new double[m];
            for (int i = 0; i < m; i++) k0[i] = p[i];
            for (int i = 0; i < m; i++) dk[i] = p[m + i];
            return new LinearCurvatureLaw(k0, dk);
        }
    }
}
