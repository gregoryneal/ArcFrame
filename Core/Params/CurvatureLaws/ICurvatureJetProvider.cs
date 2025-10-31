namespace ArcFrame.Core.Params.CurvatureLaws
{
    /// <summary>
    /// Interface to provide curvature jets
    /// </summary>
    public interface ICurvatureJetProvider
    {
        // maximum derivative order supported exactly (use int.MaxValue if unbounded)
        /// <summary>
        /// Maximum derivative order supported exactly.
        /// </summary>
        int MaxSupportedJetOrder { get; }

        // outJet[o][j] = d^o/ds^o (kappa_j) at s, for o=0..jetOrder
        /// <summary>
        /// outJet[o][j] = d^o/ds^o (kappa_j) at s, for o=0..jetOrder
        /// </summary>
        /// <param name="s"></param>
        /// <param name="jetOrder"></param>
        /// <returns></returns>
        public double[][] EvalJet(double s, int jetOrder);
    }
}
