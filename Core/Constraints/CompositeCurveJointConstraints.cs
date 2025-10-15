using ArcFrame.Core.Geometry;
using ArcFrame.Core.Params;

namespace ArcFrame.Core.Constraints
{
    /// Enforce continuity between segment i (left) and i+1 (right).
    /// Choose any of C0 (pos), C1 (tangent), C2 (curvature).
    public sealed class CompositeCurveJointConstraint : ICompositeConstraint
    {
        public ConstraintType Type { get; init; } = ConstraintType.Hard;
        public double Weight { get; init; } = 1.0;

        public int LeftIndex { get; }     // segment i
        public int RightIndex { get; }    // segment i+1
        public bool C0 { get; init; } = true;
        public bool C1 { get; init; } = true;
        public bool C2 { get; init; } = false;
        public double wC0 { get; init; } = 1.0;
        public double wC1 { get; init; } = 1.0;
        public double wC2 { get; init; } = 1.0;

        public CompositeCurveJointConstraint(int leftIndex, int rightIndex, bool c0 = true, bool c1 = true, bool c2 = false, ConstraintType type = ConstraintType.Hard, double weight = 1.0, double wC0 = 1.0, double wC1 = 1.0, double wC2 = 1.0)
        {
            LeftIndex = leftIndex; 
            RightIndex = rightIndex;
            C0 = c0; 
            C1 = c1; 
            C2 = c2;
            Type = type; 
            Weight = weight; 
            this.wC0 = wC0; 
            this.wC1 = wC1; 
            this.wC2 = wC2;
        }

        public double[] Residual(IReadOnlyList<CurveSpec> specs)
        {
            var li = specs[LeftIndex];
            var ri = specs[RightIndex];

            var Lc = new CachedIntrinsicCurve(li);
            var Rc = new CachedIntrinsicCurve(ri);

            var le = Lc.Evaluate(li.Length);     // end of left
            var rs = Rc.Evaluate(0.0);           // start of right

            var list = new List<double>();

            if (C0)
            {
                for (int i = 0; i < li.N; i++) list.Add(wC0 * (le.P[i] - rs.P[i]));
            }
            if (C1)
            {
                for (int i = 0; i < li.N; i++) list.Add(wC1 * (le.T[i] - rs.T[i]));
            }
            if (C2)
            {
                var kl = li.Kappa.Eval(li.Length);
                var kr = ri.Kappa.Eval(0.0);
                int m = System.Math.Min(kl.Length, kr.Length);
                for (int i = 0; i < m; i++) list.Add(wC2 * (kl[i] - kr[i]));
            }

            var r = list.ToArray();
            if (Weight != 1.0) for (int i = 0; i < r.Length; i++) r[i] *= Weight;
            return r;
        }
    }
}
