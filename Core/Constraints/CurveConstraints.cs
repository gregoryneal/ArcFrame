using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using System.Collections.Generic;
using System;

namespace ArcFrame.Core.Constraints
{
    /// <summary>
    /// Constraint types control the overall weight of the constraint.
    /// </summary>
    public enum ConstraintType 
    {
        /// <summary>
        /// Hard constraints are constraints that cannot be altered, they are a non-negotiable part of the end system.
        /// </summary>
        Hard,
        /// <summary>
        /// Soft constraints can be changed to lower the overall cost of the resulting curve.
        /// </summary>
        Soft
    }

    /// <summary>
    /// Every constraint returns residuals; zeros mean "satisfied".
    /// </summary>
    public interface ICurveConstraint
    {
        /// <summary>
        /// A hard or soft target
        /// </summary>
        ConstraintType Type { get; }
        /// <summary>Return residual vector for the given spec.</summary>
        double[] Residual(CurveSpec spec);
        /// <summary>Optional global weight to scale all residuals for this constraint.</summary>
        double Weight => 1.0;
        /// <summary>
        /// Which segment this constraint acts on if it
        /// is part of a CompositeCurve, 0 otherwise.
        /// </summary>
        int SegmentIndex { get; }            // which segment this piece came from (for composite problems)
    }

    /// <summary>
    /// Composite constraint (can see several CurveSpecs at once).
    /// </summary>
    public interface ICompositeConstraint
    {
        /// <summary>
        /// A hard or soft target
        /// </summary>
        public ConstraintType Type { get; }
        /// <summary>Optional global weight to scale all residuals for this constraint.</summary>
        double Weight { get; }               // multiplies the whole residual block
        /// <summary>Return residual vector for the given spec.</summary>
        double[] Residual(IReadOnlyList<CurveSpec> specs);
        /// <summary>
        /// Display the information about this constraint in the console. For debug purposes.
        /// </summary>
        public void ShowInfo();
    }

    /// <summary>
    /// Represents a regularizer that applies a penalty to the slope of the curvature parameter in order to
    /// encourage smoother solutions.
    /// </summary>
    /// <remarks>This regularizer is used to penalize the magnitude of the change in curvature (dk)
    /// for a specific segment of the curve. It is typically applied in optimization problems to favor solutions
    /// with smaller slope variations, resulting in smoother curves.</remarks>
    public sealed class SlopeRegularizer : ICurveConstraint
    {
        /// <inheritdoc/>
        public ConstraintType Type { get; } = ConstraintType.Soft;
        /// <inheritdoc/>
        public double Weight { get; }
        /// <inheritdoc/>
        public int SegmentIndex { get; }
        /// <inheritdoc/>
        public SlopeRegularizer(int segIdx, double weight) { SegmentIndex = segIdx; Weight = weight; }
        /// <inheritdoc/>
        public double[] Residual(CurveSpec spec)
        {
            // Penalize dk magnitude: weight * dk
            var p = (spec.Kappa as IParamCurvatureLaw)?.GetParams();
            // For LinearCurvatureLawParamAdapter: p = [k0, dk]
            double dk = (p != null && p.Length >= 2) ? p[1] : 0.0;
            return new[] { Weight * dk };
        }
    }



    /// Position+tangent at s=0 for a segment (hard/soft).
    public sealed class StartPoseConstraint : ICurveConstraint
    {
        /// <inheritdoc/>
        public ConstraintType Type { get; } = ConstraintType.Hard;
        /// <inheritdoc/>
        public double Weight { get; } = 1.0;
        /// <inheritdoc/>
        public int SegmentIndex { get; }
        /// <summary>
        /// The target start position
        /// </summary>
        public double[] TargetP { get; }
        /// <summary>
        /// The target start tangent
        /// </summary>
        public double[] TargetT { get; }
        /// <summary>
        /// The weight of the position residual
        /// </summary>
        public double wP { get; } = 1.0;
        /// <summary>
        /// The weight of the tangent residual
        /// </summary>
        public double wT { get; } = 1.0;
        /// <inheritdoc/>
        public StartPoseConstraint(int segmentIndex, double[] targetP, double[] targetT,
                                   ConstraintType type = ConstraintType.Hard, double weight = 1.0,
                                   double wp = 1.0, double wt = 1.0)
        {
            SegmentIndex = segmentIndex; TargetP = (double[])targetP.Clone(); TargetT = (double[])targetT.Clone();
            Type = type; Weight = weight; wP = wp; wT = wt;
        }
        /// <inheritdoc/>
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

    /// <summary>
    /// "Hit this end pose" at s = L. Provide P and T.
    /// </summary>
    public sealed class EndPoseConstraint : ICurveConstraint
    {
        /// <inheritdoc/>
        public ConstraintType Type { get; } = ConstraintType.Hard;
        /// <summary>
        /// The target end position
        /// </summary>
        public double[] TargetP { get; }
        /// <summary>
        /// The target end tangent
        /// </summary>
        public double[] TargetT { get; }
        /// <summary>
        /// Weight of the position residual
        /// </summary>
        public double wP { get; } = 1.0;
        /// <summary>
        /// Weight of the tangent residual
        /// </summary>
        public double wT { get; } = 1.0;
        /// <inheritdoc/>
        public int SegmentIndex { get; }
        /// <inheritdoc/>
        public double Weight { get; }
        /// <inheritdoc/>
        public EndPoseConstraint(int segmentIndex, double[] targetP, double[] targetT,
                                   ConstraintType type = ConstraintType.Hard, double weight = 1.0,
                                   double wp = 1.0, double wt = 1.0)
        {
            TargetP = (double[])targetP.Clone();
            TargetT = (double[])targetT.Clone();
            Type = type;
            wP = wp; wT = wt;
            SegmentIndex = segmentIndex;
            Weight = weight;
        }
        /// <inheritdoc/>
        public double[] Residual(CurveSpec spec)
        {
            //Console.WriteLine("EndPoseConstraint");
            var cur = spec.GetOptimizedCurve();
            double l = spec.Length;
            //Console.WriteLine($"Length: {l} | curve type: {cur.GetType()}");
            var smp = cur.Evaluate(spec.Length);
            //Console.WriteLine(spec.N);
            //Helpers.PrintVector(smp.P);
            //Helpers.PrintVector(TargetP);
            var r = new List<double>(spec.N * 2);
            for (int i = 0; i < spec.N; i++) r.Add(wP * (smp.P[i] - TargetP[i]));
            for (int i = 0; i < spec.N; i++) r.Add(wT * (smp.T[i] - TargetT[i]));
            if (Weight != 1.0) for (int i = 0; i < r.Count; i++) r[i] *= Weight;
            return r.ToArray();
        }
    }

    /// <summary>
    /// Pose constraint at an internal arclength s. Can include P and/or T.
    /// </summary>
    public sealed class PoseAtSConstraint : ICurveConstraint
    {
        /// <inheritdoc/>
        public ConstraintType Type { get; } = ConstraintType.Soft;
        /// <summary>
        /// Target arc length
        /// </summary>
        public double s { get; }
        /// <summary>
        /// Target position at s
        /// </summary>
        public double[]? TargetP { get; }
        /// <summary>
        /// Target tangent at s
        /// </summary>
        public double[]? TargetT { get; }
        /// <summary>
        /// Weight of the position residual
        /// </summary>
        public double wP { get; } = 1.0;
        /// <summary>
        /// Weight of the tangent residual
        /// </summary>
        public double wT { get; } = 1.0;
        /// <inheritdoc/>
        public int SegmentIndex { get; }
        /// <inheritdoc/>
        public PoseAtSConstraint(int segmentIndex, double s, double[]? targetP = null, double[]? targetT = null, ConstraintType type = ConstraintType.Soft, double wp = 1.0, double wt = 1.0)
        {
            SegmentIndex = segmentIndex;
            this.s = s;
            TargetP = targetP == null ? null : (double[])targetP.Clone();
            TargetT = targetT == null ? null : (double[])targetT.Clone();
            Type = type;
            wP = wp;
            wT = wt;
        }
        /// <inheritdoc/>
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
        /// <inheritdoc/>
        public ConstraintType Type { get; } = ConstraintType.Soft;
        /// <summary>
        /// Normal vector of the plane in question.
        /// </summary>
        public double[] n { get; }
        /// <summary>
        /// Target distance from the plane.
        /// </summary>
        public double d { get; }
        /// <summary>
        /// If true, enforce n·P(s) >= d via hinge
        /// </summary>
        public bool OneSided { get; } = false; 
        /// <summary>
        /// Start arc length for sampling.
        /// </summary>
        public double s0 { get; } = 0;
        /// <summary>
        /// End arc length for sampling.
        /// </summary>
        public double s1 { get; } = double.NaN;
        /// <summary>
        /// Numer of samples
        /// </summary>
        public int M { get; } = 9;
        /// <inheritdoc/>
        public double Weight { get; } = 1.0;
        /// <inheritdoc/>
        public int SegmentIndex { get; }
        /// <inheritdoc/>
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
        /// <inheritdoc/>
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
        /// <inheritdoc/>
        public ConstraintType Type { get; } = ConstraintType.Soft;
        /// <summary>
        /// Start arc length for sampling
        /// </summary>
        public double s0 { get; } = 0;
        /// <summary>
        /// End arc length for sampling
        /// </summary>
        public double s1 { get; } = double.NaN;
        /// <summary>
        /// Number of samples
        /// </summary>
        public int M { get; } = 9;
        /// <summary>
        /// Target max curvature allowed
        /// </summary>
        public double kMax { get; } = 1.0;
        /// <inheritdoc/>
        public double Weight { get; } = 1.0;
        /// <inheritdoc/>
        public int SegmentIndex { get; }
        /// <inheritdoc/>
        public CurvatureBoundConstraint(int segmentIndex, double kMax, double s0 = 0, double s1 = double.NaN, int M = 9, double weight = 1.0)
        {
            this.kMax = kMax;
            this.s0 = s0;
            this.s1 = s1;
            this.M = System.Math.Max(2, M);
            this.Weight = weight;
            SegmentIndex = segmentIndex;
        }
        /// <inheritdoc/>
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

    /// Curvature value at s=0 or s=L for a segment.
    public sealed class CurvatureBoundaryConstraint : ICurveConstraint
    {
        /// <summary>
        /// Where on the curve should the bound be applied.
        /// </summary>
        public enum Where 
        {
            /// <summary>
            /// The start of the curve
            /// </summary>
            Start, 
            /// <summary>
            /// The end of the curve
            /// </summary>
            End 
        }
        /// <inheritdoc/>
        public ConstraintType Type { get; } = ConstraintType.Hard;
        /// <inheritdoc/>
        public double Weight { get; } = 1.0;
        /// <inheritdoc/>
        public int SegmentIndex { get; }
        /// <summary>
        /// Where to apply the constraint on the curve.
        /// </summary>
        public Where Location { get; }
        /// <summary>
        /// The target curvature
        /// </summary>
        public double TargetK { get; }
        /// <inheritdoc/>
        public CurvatureBoundaryConstraint(int segmentIndex, double targetK, Where where,
                                           ConstraintType type = ConstraintType.Hard, double weight = 1.0)
        {
            SegmentIndex = segmentIndex; TargetK = targetK; Location = where; Type = type; Weight = weight;
        }
        /// <inheritdoc/>
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

    /// <summary>
    /// tiny soft steering penalty (per segment): prefer small |dk|
    /// </summary>
    public sealed class LocalDkPenalty : ICurveConstraint
    {
        /// <inheritdoc/>
        public ConstraintType Type => ConstraintType.Soft;
        /// <inheritdoc/>
        public double Weight { get; }
        /// <inheritdoc/>
        public int SegmentIndex { get; }
        /// <inheritdoc/>
        public LocalDkPenalty(int segIdx, double weight) { SegmentIndex = segIdx; Weight = weight; }
        /// <inheritdoc/>
        public double[] Residual(CurveSpec spec)
        {
            var plc = spec.Kappa as IParamCurvatureLaw;
            if (plc == null) return new[] { 0.0 };
            var p = plc.GetParams(); // planar adapter packs [k0, dk]
            double dk = (p.Length >= 2) ? p[1] : 0.0;
            return new[] { Weight * dk };
        }
    }

    /// <summary>
    /// Force the length of the curve to be positive.
    /// </summary>
    public sealed class PositiveLengthPenalty : ICompositeConstraint
    {
        /// <inheritdoc/>
        public ConstraintType Type => ConstraintType.Hard;
        /// <inheritdoc/>
        public double Weight { get; } = 1.0;
        /// <inheritdoc/>
        public double[] Residual(IReadOnlyList<CurveSpec> specs)
        {
            double[] r = new double[specs.Count];
            for (int i = 0; i < r.Length; i++)
            {
                if (specs[i].Length > 0) r[i] = 0;
                else r[i] = Weight;
            }
            return r;
        }
        /// <inheritdoc/>
        public void ShowInfo()
        {
            Console.WriteLine("PositiveLengthPenalty");
        }
    }

    /// <summary>
    /// Enforce a small magnitude of the curvature vector [
    /// </summary>
    public sealed class CurvatureMagnitudeRegularizer : ICurveConstraint
    {
        /// <inheritdoc/>
        public ConstraintType Type => ConstraintType.Soft;
        /// <inheritdoc/>
        public double Weight { get; }
        /// <inheritdoc/>
        public int SegmentIndex { get; }
        /// <inheritdoc/>
        public CurvatureMagnitudeRegularizer(int segIdx, double weight) { SegmentIndex = segIdx; Weight = weight; }
        /// <inheritdoc/>
        public double[] Residual(CurveSpec spec)
        {
            var p = (spec.Kappa as IParamCurvatureLaw)?.GetParams(); // [k0, τ0, ..., dk, dτ, ...]
            if (p == null) return new[] { 0.0 };
            double mag = Helpers.Len(p);
            return new[] { Weight * mag };
        }
        /// <inheritdoc/>
        public void ShowInfo()
        {
            Console.WriteLine("CurvatureMagnitudeRegularizer");
        }
    }
}
