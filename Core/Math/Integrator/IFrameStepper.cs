using System;

namespace ArcFrame.Core.Math
{
    public interface IFrameStepper
    {
        /// <summary>
        /// Advance the frame pose (P, R) by a step size h at arc length s.
        /// </summary>
        /// <param name="P">Curve Position</param>
        /// <param name="R">Basis Vectors</param>
        /// <param name="kappa">Function to calculate intrinsic curvatures at s: K_s[] = kappa(s)</param>
        /// <param name="s">The current arc length position</param>
        /// <param name="h">The size of the step along the curve</param>
        /// <param name="frame">Frame model</param>
        public void Step(ref double[] P, ref double[,] R, Func<double, double[]> kappa, double s, double h, FrameModel frame);
    }
}
