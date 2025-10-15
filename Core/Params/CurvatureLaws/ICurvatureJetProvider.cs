namespace ArcFrame.Core.Params.CurvatureLaws
{
    public interface ICurvatureJetProvider
    {
        // maximum derivative order supported exactly (use int.MaxValue if unbounded)
        int MaxSupportedJetOrder { get; }

        // outJet[o][j] = d^o/ds^o (kappa_j) at s, for o=0..jetOrder
        public double[][] EvalJet(double s, int jetOrder);
    }
}
