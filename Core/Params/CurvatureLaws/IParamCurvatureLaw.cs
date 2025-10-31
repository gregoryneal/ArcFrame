namespace ArcFrame.Core.Params
{
    /// <summary>
    /// Repackage the curvature law parameters into a flat array for optimization.
    /// </summary>
    public interface IParamCurvatureLaw : ICurvatureLaw
    {
        /// <summary>
        /// Get flattened parameters for optimization.
        /// </summary>
        /// <returns></returns>
        double[] GetParams();          // Flattened law params for optimization
        /// <summary>
        /// Decode the parameters and apply them.
        /// </summary>
        /// <param name="p"></param>
        void SetParams(double[] p);    // Decode and apply
        /// <summary>
        /// Return a new law instance with these parameters
        /// </summary>
        /// <param name="p"></param>
        /// <returns></returns>
        ICurvatureLaw CloneWithParams(double[] p); // Return a new law instance with these params
    }
}
