namespace ArcFrame.Core.Params
{
    /// <summary>
    /// Curvature law where the curvature does not change with respect to arc length.
    /// Use this for lines and arcs.
    /// </summary>
    public class ConstantCurvatureLaw : ICurvatureLaw
    {
        private readonly double[] _k0;
        public int Order => _k0.Length;
        public ConstantCurvatureLaw(double[] k0)
        {
            _k0 = (double[])k0.Clone();
        }
        public double[] Eval(double s)
        {
            return (double[])_k0.Clone();
        }
    }
}
