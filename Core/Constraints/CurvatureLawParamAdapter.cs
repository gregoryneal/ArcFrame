using ArcFrame.Core.Params;

namespace ArcFrame.Core.Constraints
{
    /// <summary>
    /// Repackage the curvature law parameters into a flat array for optimization.
    /// </summary>
    public interface IParamCurvatureLaw : ICurvatureLaw
    {
        double[] GetParams();          // Flattened law params for optimization
        void SetParams(double[] p);    // Decode and apply
        ICurvatureLaw CloneWithParams(double[] p); // Return a new law instance with these params
    }

    /// <summary>
    /// Adapter for LinearCurvatureLaw (k(s) = k0 + s*dk).
    /// Pack as [k0..., dk...].
    /// </summary>
    public sealed class LinearCurvatureLawParamAdapter : IParamCurvatureLaw
    {
        private LinearCurvatureLaw _law;
        public LinearCurvatureLawParamAdapter(LinearCurvatureLaw law) { _law = law; }
        public int Order => _law.Order;
        public double[] Eval(double s) => _law.Eval(s);

        public double[] GetParams()
        {
            var p = new double[2 * _law.Order];
            for (int i = 0; i < _law.Order; i++) p[i] = _law.Intercept[i];
            for (int i = 0; i < _law.Order; i++) p[_law.Order + i] = _law.Slope[i];
            return p;
        }

        public void SetParams(double[] p)
        {
            int m = _law.Order;
            var k0 = new double[m];
            var dk = new double[m];
            for (int i = 0; i < m; i++) k0[i] = p[i];
            for (int i = 0; i < m; i++) dk[i] = p[m + i];
            _law = new LinearCurvatureLaw(k0, dk);
        }

        public ICurvatureLaw CloneWithParams(double[] p)
        {
            int m = _law.Order;
            var k0 = new double[m];
            var dk = new double[m];
            for (int i = 0; i < m; i++) k0[i] = p[i];
            for (int i = 0; i < m; i++) dk[i] = p[m + i];
            return new LinearCurvatureLaw(k0, dk);
        }
    }
}
