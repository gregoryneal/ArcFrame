using ArcFrame.Core.Params;

namespace ArcFrame.Solvers
{
    /// <summary>
    /// Solve a curve segments intrinsic parameters with specified endpoint parameters.
    /// </summary>
    public interface IFitterG2 : ISolverInfo
    {
        public CurveSpec Fit_G2(double x0, double y0, double t0, double k0, double x1, double y1, double t1, double k1);
    }
}
