using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;
using ArcFrame.Core.Params;

namespace ArcFrame.Core.Constraints
{
    public enum ConstraintType { Hard, Soft }

    /// <summary>
    /// Every constraint returns residuals; zeros mean "satisfied".
    /// </summary>
    public interface ICurveConstraint
    {
        ConstraintType Type { get; }
        /// <summary>Return residual vector for the given spec.</summary>
        double[] Residual(CurveSpec spec);
        /// <summary>Optional global weight to scale all residuals for this constraint.</summary>
        double Weight => 1.0;
        int SegmentIndex { get; }            // which segment (for composite problems)
    }

    /// Composite constraint (can see several CurveSpecs at once).
    public interface ICompositeConstraint
    {
        ConstraintType Type { get; }
        double Weight { get; }               // multiplies the whole residual block
        double[] Residual(IReadOnlyList<CurveSpec> specs);
    }

    /// Adapter to use any ICurveConstraint in a composite problem.
    public sealed class SegmentConstraintWrapper : ICompositeConstraint
    {
        private readonly ICurveConstraint _inner;
        public SegmentConstraintWrapper(ICurveConstraint inner) { _inner = inner; }
        public ConstraintType Type => _inner.Type;
        public double Weight => _inner.Weight;

        public double[] Residual(IReadOnlyList<CurveSpec> specs)
        {
            var res = _inner.Residual(specs[_inner.SegmentIndex]);
            if (Weight != 1.0)
                for (int i = 0; i < res.Length; i++) res[i] *= Weight;
            return res;
        }
    }

    /// <summary>
    /// "Hit this end pose" at s = L. Provide P (required) and optional T.
    /// </summary>
    public sealed class EndPoseConstraint : ICurveConstraint
    {
        public ConstraintType Type { get; init; } = ConstraintType.Hard;
        public double[] TargetP { get; init; }
        public double[]? TargetT { get; init; }
        public double wP { get; init; } = 1.0;
        public double wT { get; init; } = 1.0;
        public int SegmentIndex { get; init; }

        public EndPoseConstraint(int segmentIndex, double[] targetP, double[]? targetT = null, ConstraintType type = ConstraintType.Hard, double wp = 1.0, double wt = 1.0)
        {
            TargetP = (double[])targetP.Clone();
            TargetT = targetT == null ? null : (double[])targetT.Clone();
            Type = type;
            wP = wp; wT = wt;
            SegmentIndex = segmentIndex;
        }

        public double[] Residual(CurveSpec spec)
        {
            var curve = new CachedIntrinsicCurve(spec);
            var s = spec.Length;
            var smp = curve.Evaluate(s);
            var r = new List<double>(spec.N + (TargetT == null ? 0 : spec.N));
            for (int i = 0; i < spec.N; i++)
            {
                r.Add(wP * (smp.P[i] - TargetP[i]));
            }
            if (TargetT != null)
            {
                var T = smp.T;
                // Unit-length tangent; compare directly
                for (int i = 0; i < spec.N; i++)
                {
                    r.Add(wT * (T[i] - TargetT[i]));
                }
            }
            return r.ToArray();
        }
    }

    /// <summary>
    /// Pose constraint at an internal arclength s. Can include P and/or T.
    /// </summary>
    public sealed class PoseAtSConstraint : ICurveConstraint
    {
        public ConstraintType Type { get; init; } = ConstraintType.Soft;
        public double s { get; init; }
        public double[]? TargetP { get; init; }
        public double[]? TargetT { get; init; }
        public double wP { get; init; } = 1.0;
        public double wT { get; init; } = 1.0;
        public int SegmentIndex { get; init; }

        public PoseAtSConstraint(int segmentIndex,double s, double[]? targetP = null, double[]? targetT = null, ConstraintType type = ConstraintType.Soft, double wp = 1.0, double wt = 1.0)
        {
            SegmentIndex = segmentIndex;
            this.s = s; 
            TargetP = targetP == null ? null : (double[])targetP.Clone();
            TargetT = targetT == null ? null : (double[])targetT.Clone();
            Type = type; 
            wP = wp; 
            wT = wt;
        }

        public double[] Residual(CurveSpec spec)
        {
            var curve = new CachedIntrinsicCurve(spec);
            var smp = curve.Evaluate(System.Math.Clamp(s, 0, spec.Length));
            var r = new List<double>();
            if (TargetP != null)
            {
                for (int i = 0; i < spec.N; i++) r.Add(wP * (smp.P[i] - TargetP[i]));
            }
            if (TargetT != null)
            {
                var T = smp.T;
                for (int i = 0; i < spec.N; i++) r.Add(wT * (T[i] - TargetT[i]));
            }
            return r.ToArray();
        }
    }

    /// <summary>
    /// Keep the curve on (or on one side of) a plane: n · P(s) = d (equality) or ≥ d (one-sided).
    /// Samples M points in [s0, s1].
    /// </summary>
    public sealed class PlaneConstraint : ICurveConstraint
    {
        public ConstraintType Type { get; init; } = ConstraintType.Soft;
        public double[] n { get; init; }
        public double d { get; init; }
        public bool OneSided { get; init; } = false; // if true, enforce n·P(s) >= d via hinge
        public double s0 { get; init; } = 0;
        public double s1 { get; init; } = double.NaN;
        public int M { get; init; } = 9;
        public double Weight { get; init; } = 1.0;
        public int SegmentIndex { get; init; }

        public PlaneConstraint(int segmentIndex, double[] n, double d, bool oneSided = false, double s0 = 0, double s1 = double.NaN, int M = 9, ConstraintType type = ConstraintType.Soft, double weight = 1.0)
        {
            this.n = Helpers.Normalize(n);
            this.d = d; OneSided = oneSided; 
            this.s0 = s0; 
            this.s1 = s1; 
            this.M = System.Math.Max(2, M);
            Type = type; 
            Weight = weight;
            SegmentIndex = segmentIndex;
        }

        public double[] Residual(CurveSpec spec)
        {
            double L = spec.Length;
            double a = System.Math.Clamp(s0, 0, L);
            double b = double.IsNaN(s1) ? L : System.Math.Clamp(s1, 0, L);
            if (b < a) (a, b) = (b, a);

            var curve = new CachedIntrinsicCurve(spec);
            var res = new double[M];
            for (int i = 0; i < M; i++)
            {
                double s = a + (b - a) * i / (M - 1);
                var P = curve.Position(s);
                double t = Helpers.Dot(n, P) - d;
                res[i] = Weight * (OneSided ? System.Math.Min(0.0, t) : t);
            }
            return res;
        }
    }

    /// <summary>
    /// Soft cap on curvature norm: max(0, ||kappa|| - kMax) over samples in [s0,s1].
    /// </summary>
    public sealed class CurvatureBoundConstraint : ICurveConstraint
    {
        public ConstraintType Type { get; init; } = ConstraintType.Soft;
        public double s0 { get; init; } = 0;
        public double s1 { get; init; } = double.NaN;
        public int M { get; init; } = 9;
        public double kMax { get; init; } = 1.0;
        public double Weight { get; init; } = 1.0;
        public int SegmentIndex { get; init; }

        public CurvatureBoundConstraint(int segmentIndex, double kMax, double s0 = 0, double s1 = double.NaN, int M = 9, double weight = 1.0)
        {
            this.kMax = kMax; 
            this.s0 = s0; 
            this.s1 = s1; 
            this.M = System.Math.Max(2, M); 
            this.Weight = weight;
            SegmentIndex = segmentIndex;
        }

        public double[] Residual(CurveSpec spec)
        {
            double L = spec.Length;
            double a = System.Math.Clamp(s0, 0, L);
            double b = double.IsNaN(s1) ? L : System.Math.Clamp(s1, 0, L);
            if (b < a) (a, b) = (b, a);

            var res = new double[M];
            for (int i = 0; i < M; i++)
            {
                double s = a + (b - a) * i / (M - 1);
                var k = spec.Kappa.Eval(s);
                double kn = Helpers.Len(k);
                res[i] = Weight * System.Math.Max(0.0, kn - kMax);
            }
            return res;
        }
    }

    /// Position+tangent at s=0 for a segment (hard/soft).
    public sealed class StartPoseConstraint : ICurveConstraint
    {
        public ConstraintType Type { get; init; } = ConstraintType.Hard;
        public double Weight { get; init; } = 1.0;
        public int SegmentIndex { get; init; }

        public double[] TargetP { get; }
        public double[] TargetT { get; }
        public double wP { get; init; } = 1.0;
        public double wT { get; init; } = 1.0;

        public StartPoseConstraint(int segmentIndex, double[] targetP, double[] targetT,
                                   ConstraintType type = ConstraintType.Hard, double weight = 1.0,
                                   double wp = 1.0, double wt = 1.0)
        {
            SegmentIndex = segmentIndex; TargetP = (double[])targetP.Clone(); TargetT = (double[])targetT.Clone();
            Type = type; Weight = weight; wP = wp; wT = wt;
        }

        public double[] Residual(CurveSpec spec)
        {
            var cur = new CachedIntrinsicCurve(spec);
            var smp = cur.Evaluate(0.0);
            var r = new List<double>(spec.N * 2);
            for (int i = 0; i < spec.N; i++) r.Add(wP * (smp.P[i] - TargetP[i]));
            for (int i = 0; i < spec.N; i++) r.Add(wT * (smp.T[i] - TargetT[i]));
            if (Weight != 1.0) for (int i = 0; i < r.Count; i++) r[i] *= Weight;
            return r.ToArray();
        }
    }

    /// Curvature value at s=0 or s=L for a segment.
    public sealed class CurvatureBoundaryConstraint : ICurveConstraint
    {
        public enum Where { Start, End }
        public ConstraintType Type { get; init; } = ConstraintType.Hard;
        public double Weight { get; init; } = 1.0;
        public int SegmentIndex { get; init; }

        public Where Location { get; }
        public double TargetK { get; }

        public CurvatureBoundaryConstraint(int segmentIndex, double targetK, Where where,
                                           ConstraintType type = ConstraintType.Hard, double weight = 1.0)
        {
            SegmentIndex = segmentIndex; TargetK = targetK; Location = where; Type = type; Weight = weight;
        }

        public double[] Residual(CurveSpec spec)
        {
            // For planar curves the law returns a scalar curvature in [0].
            double k = Location == Where.Start
                ? spec.Kappa.Eval(0.0)[0]
                : spec.Kappa.Eval(spec.Length)[0];
            double r = (k - TargetK) * Weight;
            return new[] { r };
        }
    }
}
