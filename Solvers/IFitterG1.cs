using ArcFrame.Core.Params;

namespace ArcFrame.Solvers
{
    /// <summary>
    /// Solve for the intrinsic curve parameters given endpoint parameters.
    /// </summary>
    public interface IFitterG1 : ISolverInfo
    {
        public CurveSpec Fit_G1(double x0, double y0, double t0, double x1, double y1, double t1);
    }
}
