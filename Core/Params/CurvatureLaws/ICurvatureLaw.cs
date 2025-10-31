namespace ArcFrame.Core.Params
{
    /// <summary>
    /// Generator for curvature(s) for intrinsic curves.
    /// -> Constant: Line/Arc   k = k0
    /// -> Linear: Clothoid     k = k0 + (k')s
    /// -> Polynomial: ?        k = k0 + (k')s + (k'')s^2 + ... (dk^n)s^n
    /// </summary>
    public interface ICurvatureLaw
    {
        /// <summary>
        /// Should be N-1 dimensions of target curve
        /// </summary>
        public int Order { get; }
        /// <summary>
        /// Calculate and return the generalized curvatures at the arc length s
        /// </summary>
        /// <param name="s">The arc length</param>
        /// <returns></returns>
        double[] Eval(double s);
        /// <summary>
        /// Does the curvature depend linearly on the arc length?
        /// </summary>
        public bool IsLinear { get; }
        /// <summary>
        /// Is the curvature constant along the arc length?
        /// </summary>
        public bool IsConstant { get; }

    }
}
