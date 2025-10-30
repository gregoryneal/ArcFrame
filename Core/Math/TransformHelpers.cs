using ArcFrame.Core.Geometry;

using System;

namespace ArcFrame.Core.Math
{
    /// <summary>
    /// Helper methods to create transformed curves where the dimension of the input curve and the dimension of the rigid transform do not match.
    /// </summary>
    public static class TransformHelpers
    {
        /// <summary>
        /// Canonical promotion: embed curve in R^targetN along the first axes (pad zeros).
        /// </summary>
        /// <param name="curve"></param>
        /// <param name="xfHigher"></param>
        /// <param name="axisMap"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static TransformedCurve PromoteCurveAndTransform(IArcLengthCurve curve, RigidTransform xfHigher, int[]? axisMap = null)
        {
            if (xfHigher.N < curve.Dimension) throw new ArgumentException("Transform must have >= curve dimension for promotion.");
            var promoted = new PromotedCurve(curve, xfHigher.N, axisMap); // pads to N
            return new TransformedCurve(promoted, xfHigher);
        }

        /// <summary>
        /// Option: custom basis promotion using an embedding matrix E (N x d), columns orthonormal.
        /// </summary>
        /// <param name="curve">d dimensional curve</param>
        /// <param name="xfHigher">N dimensional transform</param>
        /// <param name="E">N x d embedding matrix</param>
        /// <param name="p0">Optional N dimensional offset</param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static TransformedCurve PromoteCurveWithBasisAndTransform(IArcLengthCurve curve, RigidTransform xfHigher, double[,] E /*N x d*/, double[]? p0 = null)
        {
            if (xfHigher.N < curve.Dimension) throw new ArgumentException("Transform too small.");
            var promoted = new BasisPromotedCurve(curve, E, p0); // embeds local d-dim curve into N-dim via E
            return new TransformedCurve(promoted, xfHigher);
        }
    }
}
