using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using ArcFrame.Core.Results;
using ArcFrame.Solvers;
using ArcFrame.Solvers.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.ConstrainedExecution;
using System.Threading;

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
        /// <summary>
        /// Build a residual from a given cached Sample object
        /// </summary>
        /// <param name="ctx"></param>
        /// <returns></returns>
        double[] Residual(IterationContext ctx);
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
        /// Build a residual from a given cached Sample object
        /// </summary>
        /// <param name="ctx"></param>
        /// <returns></returns>
        double[] Residual(IterationContext ctx);
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

        /// <inheritdoc/>
        public double[] Residual(IterationContext ctx)
        {
            var spec = ctx.Specs[0];
            return Residual(spec);
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

        private int? _sampleIndex;
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

        /// <inheritdoc/>
        public double[] Residual(IterationContext ctx)
        {
            int idx;
            if (_sampleIndex.HasValue) idx = _sampleIndex.Value;
            else
            {   
                // if we dont chain the curve specs we can optimize the position of each segment directly
                if (ctx.Mode == CompositeCurveSolver.CompositeSolverMode.Seperate)
                {
                    _sampleIndex = System.Math.Max(0, Array.IndexOf(ctx.SegIdx, SegmentIndex));
                }
                // otherwise we can only optimize the position of the first segment
                else
                {
                    _sampleIndex = 0;
                }

                idx = _sampleIndex.Value;
            }
            var smp = ctx.S[idx];
            var spec = ctx.Specs[SegmentIndex];
            var r = new List<double>(spec.N * 2);
            if (idx >= ctx.S.Length)
            {
                for (int i = 0; i < spec.N; i++)
                {
                    r.Add(wP * Weight);
                }
                for (int i = 0; i < spec.N; i++)
                {
                    r.Add(wT * Weight);
                }
                return r.ToArray();
            }
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

        private int? _sampleIndex;
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

        /// <inheritdoc/>
        public double[] Residual(IterationContext ctx)
        {
            int idx;
            var spec = ctx.Specs[SegmentIndex];
            var r = new List<double>(spec.N * 2);
            if (!_sampleIndex.HasValue)
            {
                if (SegmentIndex == ctx.Specs.Count - 1) _sampleIndex = ctx.S.Length - 1;
                // ensure we use IterationContext.Build_IncludeStartAndEnd
                else
                {
                    if (ctx.Mode == CompositeCurveSolver.CompositeSolverMode.Seperate)
                    {
                        // end of this segment is the final sample with this segment index
                        _sampleIndex = System.Math.Max(0, Array.LastIndexOf(ctx.SegIdx, SegmentIndex));
                    }
                    else
                    {
                        _sampleIndex = ctx.S.Length - 1;
                    }
                    if (_sampleIndex == 0)
                    {
                        // default to final sample
                        _sampleIndex = ctx.S.Length - 1;
                    }
                }
            }
            idx = _sampleIndex.Value;
            if (idx >= ctx.S.Length)
            {
                for (int i = 0; i < spec.N; i++)
                {
                    r.Add(wP * Weight);
                }
                for (int i = 0; i < spec.N; i++)
                {
                    r.Add(wT * Weight);
                }
                return r.ToArray();
            }

            //Console.WriteLine($"looking for index {idx} in array with max index {ctx.S.Length - 1}, seg count: {ctx.Specs.Count}");
            //Helpers.PrintVector(ctx.SegIdx.Select(i => (double)i).ToArray());
            var smp = ctx.S[idx];
            for (int i = 0; i < spec.N; i++)
            {
                r.Add(wP * (smp.P[i] - TargetP[i]));
            }
            for (int i = 0; i < spec.N; i++)
            {
                r.Add(wT * (smp.T[i] - TargetT[i]));
            }
            if (Weight != 1.0)
            {
                for (int i = 0; i < r.Count; i++) r[i] *= Weight;
            }
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

        /// <inheritdoc/>
        public double[] Residual(IterationContext ctx)
        {
            //int segIdx = System.Math.Clamp(SegmentIndex, 0, ctx.Specs.Count - 1);
            int N = ctx.Specs[0].N;
            (int sampleIdx, bool isExact) = Helpers.BinarySearch(ctx.SGlob, this.s);
            Sample smp;
            if (isExact)
            {
                smp = ctx.S[sampleIdx];
            }
            else
            {
                // if we didn't find the exact match the index tells us 
                // the index of the first item larger than value
                double interpDs = (s - ctx.SGlob[sampleIdx - 1]) / ctx.Ds;
                double sLoc = ctx.SGlob[sampleIdx - 1] + (interpDs * ctx.Ds);
                Sample s0 = ctx.S[sampleIdx - 1];
                Sample s1 = ctx.S[sampleIdx];
                smp = new Sample()
                {
                    k = Helpers.Add(s0.k, Helpers.Multiply(interpDs, Helpers.Subtract(s1.k, s0.k))),
                    P = Helpers.Add(s0.P, Helpers.Multiply(interpDs, Helpers.Subtract(s1.P, s0.P))),
                    R = Helpers.Add(s0.R, Helpers.Multiply(interpDs, Helpers.Subtract(s1.R, s0.R))),
                    s = sLoc,
                };
            }
            var r = new List<double>();

            if (TargetP != null)
            {
                for (int i = 0; i < N; i++) r.Add(wP * (smp.P[i] - TargetP[i]));
            }
            if (TargetT != null)
            {
                var T = smp.T;
                for (int i = 0; i < N; i++) r.Add(wT * (T[i] - TargetT[i]));
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
        /// <inheritdoc/>
        public double[] Residual(IterationContext ctx)
        {
            // 1. clamp the desired values
            // 2. find the index bounds that fully encompass desired arc length bounds
            // 3. recalculate arc length interval
            double s0_actual = System.Math.Clamp(s0, 0, ctx.Prefix[^1]);
            double s1_actual = System.Math.Clamp(s1, 0, ctx.Prefix[^1]);
            (int s0_idx, bool isExact0) = Helpers.BinarySearch(ctx.SGlob, s0_actual);
            // the start sample index
            s0_idx = System.Math.Clamp(s0_idx - 1, 0, ctx.SGlob.Length - 1);
            (int s1_idx, bool isExact1) = Helpers.BinarySearch(ctx.SGlob, s1_actual);
            // the end sample index
            s1_idx = System.Math.Clamp(s1_idx, 0, ctx.SGlob.Length - 1);

            if (s0_idx > s1_idx) (s0_idx, s1_idx) = (s1_idx, s0_idx);
            int M = s1_idx - s0_idx + 1;
            double[] res = new double[M];
            double s, t;
            int sIdx;
            Sample smp;
            double[] P;
            for (int i = 0; i < M; i++)
            {
                sIdx = s0_idx + i;
                smp = ctx.S[sIdx];
                s = smp.s;
                P = smp.P;
                t = Helpers.Dot(n, P) - d;
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

        /// <inheritdoc/>
        public double[] Residual(IterationContext ctx)
        {
            // 1. clamp the desired values
            // 2. find the index bounds that fully encompass desired arc length bounds
            // 3. recalculate arc length interval
            double s0_actual = System.Math.Clamp(s0, 0, ctx.Prefix[^1]);
            double s1_actual = System.Math.Clamp(s1, 0, ctx.Prefix[^1]);
            (int s0_idx, bool isExact0) = Helpers.BinarySearch(ctx.SGlob, s0_actual);
            // the start sample index
            s0_idx = System.Math.Clamp(s0_idx - 1, 0, ctx.SGlob.Length - 1);
            (int s1_idx, bool isExact1) = Helpers.BinarySearch(ctx.SGlob, s1_actual);
            // the end sample index
            s1_idx = System.Math.Clamp(s1_idx, 0, ctx.SGlob.Length - 1);

            if (s0_idx > s1_idx) (s0_idx, s1_idx) = (s1_idx, s0_idx);
            int M = s1_idx - s0_idx + 1;
            double[] res = new double[M];
            int sIdx;
            Sample smp;
            for (int i = 0; i < M; i++)
            {
                sIdx = s0_idx + i;
                smp = ctx.S[sIdx];
                double kn = Helpers.Len(smp.k);
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

        /// <inheritdoc/>
        public double[] Residual(IterationContext ctx)
        {
            double k = Location == Where.Start
                ? ctx.S[0].k[0]
                : ctx.S[^1].k[0];
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

        /// <inheritdoc/>
        public double[] Residual(IterationContext ctx)
        {
            int idx = System.Math.Clamp(SegmentIndex, 0, ctx.Specs.Count - 1);
            var plc = ctx.Specs[idx].Kappa as IParamCurvatureLaw;
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
        public double[] Residual(IterationContext ctx)
        {
            double[] r = new double[ctx.Specs.Count];
            for (int i = 0; i < r.Length; i++)
            {
                if (ctx.Specs[i].Length > 0) r[i] = 0;
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
        public double[] Residual(IterationContext ctx)
        {
            int idx = System.Math.Clamp(SegmentIndex, 0, ctx.Specs.Count - 1);
            var p = (ctx.Specs[idx].Kappa as IParamCurvatureLaw)?.GetParams(); // [k0, τ0, ..., dk, dτ, ...]
            if (p == null) return new[] { Weight };
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
