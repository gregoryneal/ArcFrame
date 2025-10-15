using ArcFrame.Core.Math;
using ArcFrame.Core.Results;

namespace ArcFrame.Core.Geometry
{
    /// <summary>
    /// A class holding transformation and rotation data for an IArcLengthCurve
    /// These are applied on top of the generated curve.
    /// </summary>
    public sealed class TransformedCurve : IArcLengthCurve
    {
        private readonly IArcLengthCurve _inner;
        private readonly RigidTransform _transform;
        public double[] Position(double s) => Evaluate(s).P;
        public double[] Tangent(double s) => Evaluate(s).T;

        public TransformedCurve(IArcLengthCurve inner, RigidTransform transform)
        {
            if (inner.Dimension != transform.N) throw new ArgumentException($"Dimension mismatch between curve ({inner.Dimension}) and transform ({transform.N})");
            _inner = inner;
            _transform = transform;
        }

        /// <summary>
        /// Create transformed curve from the identity transform.
        /// </summary>
        /// <param name="inner"></param>
        /// <returns></returns>
        public TransformedCurve Identity(IArcLengthCurve inner)
        {
            return new TransformedCurve(inner, RigidTransform.Identity(inner.Dimension));
        }

        /// <summary>
        /// Create a transformation with an N dimensional rotation and translation matrices.
        /// </summary>
        /// <param name="inner">The inner curve</param>
        /// <param name="R">The NxN rotation matrix</param>
        /// <param name="T">The N dimensional translation vector</param>
        /// <returns></returns>
        public TransformedCurve FromRT(IArcLengthCurve inner, double[,] R, double[] T)
        {
            return new TransformedCurve(inner, new RigidTransform(R, T));
        }

        /// <summary>
        /// Create a transformation for rotating about the +Y axis and translating along the XZ plane.
        /// </summary>
        /// <param name="inner">Inner curve</param>
        /// <param name="angle">Rotation about the +Y axis</param>
        /// <param name="pivotX">Pivot for rotation X coordinate</param>
        /// <param name="pivotZ">Pivot for rotation Z coordinate</param>
        /// <param name="tx">Extra translation in X direction</param>
        /// <param name="tz">Extra translation in Z direction</param>
        /// <returns></returns>
        public TransformedCurve FromYawXZ(IArcLengthCurve inner, double angle, double pivotX, double pivotZ, double tx, double tz)
        {
            return new TransformedCurve(inner, RigidTransform.FromYawXZ(inner.Dimension, angle, pivotX, pivotZ, tx, tz));
        }

        public int Dimension => _inner.Dimension;

        public double Length => _inner.Length;

        public Sample Evaluate(double s)
        {
            var sp = _inner.Evaluate(s);
            var P = _transform.ApplyTransform(sp.P);
            var R = Helpers.Multiply(_transform.R, sp.R);
            //var T = _transform.ApplyRotation(sp.T);
            return new Sample(P, R, s, (double[])sp.k.Clone());
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
