using ArcFrame.Core.Geometry;

namespace ArcFrame.Core.Results
{
    /// <summary>
    /// Class to transfer fitting results from an IFitter implementation.
    /// </summary>
    /// <typeparam name="TCurve">An arc length parameterized curve that was fitted to the data</typeparam>
    public sealed class FitResult<TCurve> where TCurve : IArcLengthCurve
    {
        public bool Ok { get; private set; }
        public string? Error { get; private set; }
        public TCurve? Curve { get; private set; }

        public static FitResult<TCurve> Success(TCurve c)
        {
            var r = new FitResult<TCurve>();
            r.Ok = true;
            r.Curve = c;
            r.Error = null;
            return r;
        }
        public static FitResult<TCurve> Fail(string message)
        {
            var r = new FitResult<TCurve>();
            r.Ok = false;
            r.Curve = default;
            r.Error = message;
            return r;
        }
    }
}
