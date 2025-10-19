namespace ArcFrame.Core.Params
{
    /// <summary>
    /// Curvature law where the curvature does not change with respect to arc length.
    /// Use this for lines and arcs.
    /// </summary>
    public class ConstantCurvatureLaw : IParamCurvatureLaw
    {
        private readonly double[] _k0;
        public int Order => _k0.Length;

        public bool IsLinear => false;

        public bool IsConstant => true;

        public ConstantCurvatureLaw(double[] k0)
        {
            _k0 = (double[])k0.Clone();
        }
        public double[] Eval(double s)
        {
            return (double[])_k0.Clone();
        }

        public double[] GetParams()
        {
            return (double[])_k0.Clone();
        }

        public void SetParams(double[] p)
        {
            for (int i = 0; i < p.Length; i++)
            {
                _k0[i] = p[i];
            }
        }

        public ICurvatureLaw CloneWithParams(double[] p)
        {
            return new ConstantCurvatureLaw((double[])p.Clone());
        }
    }
}
