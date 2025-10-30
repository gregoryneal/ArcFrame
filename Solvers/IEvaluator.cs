using ArcFrame.Core.Params;
using ArcFrame.Core.Results;

namespace ArcFrame.Solvers
{
    /// <summary>
    /// Generate a Sample of a curve given the arc length and CurveSpec
    /// </summary>
    public interface IEvaluator : ISolverInfo
    {
        /// <summary>
        /// How many dimensions this evaluator is useful in
        /// </summary>
        public int TargetDimension { get; }
        /// <summary>
        /// Evaluate a curve described by intrinsic parameters at arc length s.
        /// </summary>
        /// <param name="p"></param>
        /// <param name="s"></param>
        /// <returns></returns>
        Sample Evaluate(CurveSpec p, double s);
    }
}
