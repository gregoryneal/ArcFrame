using ArcFrame.Core.Results;

namespace ArcFrame.Core.Geometry
{
    /// <summary>
    /// Interface for an ND arc length parameterized curve.
    /// </summary>
    public interface IArcLengthCurve
    {
        public int Dimension { get; }
        public double Length { get; }
        /// <summary>
        /// Clamp in [0, Length]
        /// </summary>
        public Sample Evaluate(double s);
        public double[] Position(double s)
        {
            return Evaluate(s).P;
        }

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
