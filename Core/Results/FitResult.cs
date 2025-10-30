using ArcFrame.Core.Geometry;

namespace ArcFrame.Core.Results
{
    /// <summary>
    /// Class to transfer fitting results from an IFitter implementation.
    /// </summary>
    /// <typeparam name="IArcLengthCurve">An arc length parameterized curve that was fitted to the data</typeparam>
    public sealed class FitResult<IArcLengthCurve>
    {
        public bool Ok { get; private set; } = false;
        public string? Error { get; private set; } = "";
        public IArcLengthCurve Curve { get; private set; }

        public static FitResult<IArcLengthCurve> Success(IArcLengthCurve c)
        {
            var r = new FitResult<IArcLengthCurve>();
            r.Ok = true;
            r.Curve = c;
            r.Error = null;
            return r;
        }
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
