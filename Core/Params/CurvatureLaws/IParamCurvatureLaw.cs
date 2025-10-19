namespace ArcFrame.Core.Params
{
    /// <summary>
    /// Repackage the curvature law parameters into a flat array for optimization.
    /// </summary>
    public interface IParamCurvatureLaw : ICurvatureLaw
    {
        double[] GetParams();          // Flattened law params for optimization
        void SetParams(double[] p);    // Decode and apply
        ICurvatureLaw CloneWithParams(double[] p); // Return a new law instance with these params
    }
}
