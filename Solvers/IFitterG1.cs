using ArcFrame.Core.Params;

namespace ArcFrame.Solvers
{
    /// <summary>
    /// Solve for the intrinsic curve parameters given endpoint parameters.
    /// </summary>
    public interface IFitterG1 : ISolverInfo
    {
        /// <summary>
        /// Fit a 2D curve with G1 tangents.
        /// </summary>
        /// <param name="x0"></param>
        /// <param name="y0"></param>
        /// <param name="t0"></param>
        /// <param name="x1"></param>
        /// <param name="y1"></param>
        /// <param name="t1"></param>
        /// <returns></returns>
        public CurveSpec Fit_G1(double x0, double y0, double t0, double x1, double y1, double t1);
    }
}
