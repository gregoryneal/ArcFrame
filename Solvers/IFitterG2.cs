using ArcFrame.Core.Params;

namespace ArcFrame.Solvers
{
    /// <summary>
    /// Solve a curve segments intrinsic parameters with specified endpoint parameters.
    /// </summary>
    public interface IFitterG2 : ISolverInfo
    {
        /// <summary>
        /// Fit a 2D curve with G2 endpoints.
        /// </summary>
        /// <param name="x0"></param>
        /// <param name="y0"></param>
        /// <param name="t0"></param>
        /// <param name="k0"></param>
        /// <param name="x1"></param>
        /// <param name="y1"></param>
        /// <param name="t1"></param>
        /// <param name="k1"></param>
        /// <returns></returns>
        public CurveSpec Fit_G2(double x0, double y0, double t0, double k0, double x1, double y1, double t1, double k1);
    }
}
