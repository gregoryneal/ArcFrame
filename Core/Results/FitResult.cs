using ArcFrame.Core.Geometry;

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
        public string? Error { get; private set; } = "";
        /// <summary>
        /// The fit curve if one was found.
        /// </summary>
        public IArcLengthCurve Curve { get; private set; }
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
            r.Error = null;
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
            r.Error = message;
            return r;
        }
    }
}
