using ArcFrame.Core.Geometry;
using ArcFrame.Core.Params;

namespace ArcFrame.Core.Results
{
    /// <summary>
    /// Class to transfer fitting results from an IFitter implementation.
    /// </summary>
    /// <typeparam name="IArcLengthCurve">An arc length parameterized curve that was fitted to the data</typeparam>
    public sealed class FitResult<IArcLengthCurve>
    {
        /// <summary>
        /// Did the FitResult find a good fit?
        /// </summary>
        public bool Ok { get; private set; } = false;
        /// <summary>
        /// Error message if applicable
        /// </summary>
        public string? Message { get; private set; } = "";
        /// <summary>
        /// The fit curve if one was found.
        /// </summary>
        public IArcLengthCurve Curve { get; private set; }
        /// <summary>
        /// Number of solver iterations, if applicable.
        /// </summary>
        public int Iterations { get; private set; } = 0;

        /// <summary>
        /// Successfully found a curve that fits.
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static FitResult<IArcLengthCurve> Success(IArcLengthCurve c)
        {
            var r = new FitResult<IArcLengthCurve>();
            r.Ok = true;
            r.Curve = c;
            r.Message = "Fit successful";
            return r;
        }
        /// <summary>
        /// We did not find a curve that fits.
        /// </summary>
        /// <param name="message"></param>
        /// <returns></returns>
        public static FitResult<IArcLengthCurve> Fail(string message)
        {
            var r = new FitResult<IArcLengthCurve>();
            r.Ok = false;
            r.Curve = default;
            r.Message = message;
            return r;
        }
    }

    /// <summary>
    /// Minimal solver result
    /// for curve spec parameter fitting.
    /// </summary>
    public struct Result
    {
        /// <summary>
        /// Did the result converge to a value.
        /// </summary>
        public bool Converged;
        /// <summary>
        /// Number of solver iterations to converge.
        /// </summary>
        public int Iterations;
        /// <summary>
        /// The final cost of the curve.
        /// </summary>
        public double FinalCost;
        /// <summary>
        /// The final lambda of the solver.
        /// </summary>
        public double Lambda;
        /// <summary>
        /// The final parameter set.
        /// </summary>
        public double[] Parameters;
        /// <summary>
        /// The final curve, if applicable. 
        /// </summary>
        public CompositeCurve? Curve;

        /// <summary>
        /// Build a Result object. Not sure what else I should put here.
        /// </summary>
        /// <param name="converged"></param>
        /// <param name="iterations"></param>
        /// <param name="finalCost"></param>
        /// <param name="lambda"></param>
        /// <param name="parameters"></param>
        /// <param name="curve"></param>
        public Result(bool converged, int iterations, double finalCost, double lambda, double[] parameters, CompositeCurve? curve = null)
        {
            Converged = converged;
            Iterations = iterations;
            FinalCost = finalCost;
            Lambda = lambda;
            Parameters = parameters;
            Curve = curve;
        }
    }
}
