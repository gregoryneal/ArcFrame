using ArcFrame.Core.Results;

namespace ArcFrame.Core.Geometry
{
    /// <summary>
    /// Interface for an ND arc length parameterized curve.
    /// </summary>
    public interface IArcLengthCurve
    {
        /// <summary>
        /// The spatial dimension of the curve.
        /// </summary>
        public int Dimension { get; }
        /// <summary>
        /// The total arc length of the curve.
        /// </summary>
        public double Length { get; }
        /// <summary>
        /// Clamp in [0, Length]
        /// </summary>
        public Sample Evaluate(double s);
        /// <summary>
        /// Shortcut to get the position. 
        /// Note: this does a full evaluation.
        /// If you need other information it 
        /// will be cheaper to call Evaluate(s)
        /// instead.
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        public double[] Position(double s)
        {
            return Evaluate(s).P;
        }
        /// <summary>
        /// Shortcut to grab the tangent of the Sample
        /// at s.
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        public double[] Tangent(double s)
        {
            return Evaluate(s).T;
        }

        /// <summary>
        /// Get a number of evenly spaced samples along the arc length
        /// </summary>
        /// <param name="count"></param>
        /// <returns></returns>
        public Sample[] GetSamples(int count)
        {
            Sample[] samples = new Sample[count];
            if (count <= 0) count = 1;
            double s;
            for (int i = 0; i < count; i++)
            {
                s = (i == count - 1) ? Length : i * (Length / System.Math.Max(1, count - 1));
                samples[i] = Evaluate(s);
            }
            return samples;
        }
    }
}
