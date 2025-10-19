using ArcFrame.Core.Geometry;
using ArcFrame.Core.Params;
using static ArcFrame.Core.Constraints.LaneCorridorConstraint;

namespace ArcFrame.Core.Constraints
{
    /// <summary>
    /// Adapter to use any ICurveConstraint in a composite problem.
    /// </summary>
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

        public void ShowInfo()
        {
            Type t = GetType();
            Console.WriteLine($"SegmentConstraintWrapper -> {_inner.GetType().Name}");
        }
    }


    /// <summary>
    /// Represents a constraint that applies a regularization term based on the length of a specific curve segment.
    /// </summary>
    /// <remarks>This constraint is typically used in optimization scenarios to penalize or encourage certain
    /// segment lengths by applying a weighted regularization term. The regularization is applied to the segment
    /// specified by its index.</remarks>
    public sealed class LengthRegularizer : ICompositeConstraint
    {
        public ConstraintType Type { get; } = ConstraintType.Soft;
        public double Weight { get; }
        private readonly int _segIdx;

        public LengthRegularizer(int segIdx, double weight) { _segIdx = segIdx; Weight = weight; }

        public double[] Residual(IReadOnlyList<CurveSpec> specs)
        {
            return new[] { Weight * specs[_segIdx].Length };
        }

        public void ShowInfo()
        {
            Type t = GetType();
            Console.WriteLine($"LengthRegularizer");
        }
    }

    /// <summary>
    /// Enforce continuity between segment i (left) and i+1 (right).
    /// Choose any of C0 (pos), C1 (tangent), C2 (curvature).
    /// </summary>
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
            //Console.WriteLine($"LeftIndex: {LeftIndex} | RightIndex: {RightIndex}");
            //foreach (var spec in specs) spec.ShowInfo();
            var li = specs[LeftIndex];
            var ri = specs[RightIndex];

            IArcLengthCurve Lc = li.GetOptimizedCurve();
            IArcLengthCurve Rc = ri.GetOptimizedCurve();

            //Console.WriteLine($"CompositeCurveJointConstraint pre evaluate {li.Length} | {Lc.Length}");
            var le = Lc.Evaluate(li.Length);     // end of left
            //Console.WriteLine("CompositeCurveJointConstraint post evaluate left");
            var rs = Rc.Evaluate(0.0);           // start of right
            //Console.WriteLine("CompositeCurveJointConstraint post evaluate right");

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

        public void ShowInfo()
        {
            Type t = GetType();
            Console.WriteLine($"{t.Name}");
        }
    }

    /// <summary>
    /// dk jump smoother: encourage similar dk across composite joints
    /// </summary>
    public sealed class SlopeJumpRegularizer : ICompositeConstraint
    {
        public ConstraintType Type => ConstraintType.Soft;
        public double Weight { get; }
        public SlopeJumpRegularizer(double weight) { Weight = weight; }
        public double[] Residual(IReadOnlyList<CurveSpec> specs)
        {
            var r = new List<double>();
            for (int i = 0; i + 1 < specs.Count; i++)
            {
                double dki = GetDk(specs[i]);
                double dkj = GetDk(specs[i + 1]);
                r.Add(Weight * (dki - dkj));
            }
            return r.Count == 0 ? new[] { 0.0 } : r.ToArray();
        }
        private static double GetDk(CurveSpec s)
        {
            var plc = s.Kappa as IParamCurvatureLaw; if (plc == null) return 0.0;
            var p = plc.GetParams(); return (p.Length >= 2) ? p[1] : 0.0;
        }

        public void ShowInfo()
        {
            Type t = GetType();
            Console.WriteLine($"SlopeJumpRegularizer");
        }
    }

    /// <summary>
    /// Keep the path inside a varying "lane corridor":
    /// y in [centerY(u) - 0.5*width(u), centerY(u) + 0.5*width(u)], where u = s_global / L_total in [0,1].
    /// Hinge penalty for outside; zero residual if inside. Sampling is uniform in arclength.
    /// </summary>
    public sealed class LaneCorridorConstraint : ICompositeConstraint
    {
        public ConstraintType Type { get; } = ConstraintType.Soft;
        public double Weight { get; }
        private readonly double[] _centerY;   // samples over u in [0,1]
        private readonly double[] _width;     // samples over u in [0,1]
        private readonly int _samplesPerSegment;

        public void ShowInfo()
        {
            Type t = GetType();
            Console.WriteLine($"LaneCorridorConstraint");
        }

        public LaneCorridorConstraint(double[] laneCenterY, double[] laneWidth, double weight = 4.0, int samplesPerSegment = 21)
        {
            if (laneCenterY == null || laneCenterY.Length < 2) throw new ArgumentException("laneCenterY needs >=2 samples");
            if (laneWidth == null || laneWidth.Length < 2) throw new ArgumentException("laneWidth needs >=2 samples");
            _centerY = laneCenterY;
            _width = laneWidth;
            Weight = weight;
            _samplesPerSegment = System.Math.Max(4, samplesPerSegment);
        }

        public double[] Residual(IReadOnlyList<CurveSpec> specs)
        {
            // total length and per-segment prefix lengths
            double Ltot = 0.0;
            var prefix = new double[specs.Count + 1];
            prefix[0] = 0.0;
            for (int i = 0; i < specs.Count; i++) { Ltot += specs[i].Length; prefix[i + 1] = Ltot; }
            if (Ltot <= 0) return new double[] { 0.0 };

            var residuals = new List<double>();

            for (int seg = 0; seg < specs.Count; seg++)
            {
                var c = new CachedIntrinsicCurve(specs[seg]);
                for (int k = 0; k < _samplesPerSegment; k++)
                {
                    double sLocal = specs[seg].Length * k / (_samplesPerSegment - 1.0);
                    double sGlobal = prefix[seg] + sLocal;
                    double u = sGlobal / Ltot;

                    var P = c.Position(sLocal); // planar: [x, y]
                    double y = P[1];

                    // what the center should be
                    double center = PLInterp(u, _centerY);
                    double w = System.Math.Max(0.0, PLInterp(u, _width));
                    double lo = center - 0.5 * w;
                    double hi = center + 0.5 * w;

                    // hinge penalty: inside => 0; outside => distance to nearest bound
                    double pen = 0.0;
                    if (y < lo) pen = lo - y;
                    else if (y > hi) pen = y - hi;

                    if (pen != 0.0) residuals.Add(Weight * pen);
                }
            }

            return residuals.Count == 0 ? new double[] { 0.0 } : residuals.ToArray();
        }

        // piecewise-linear interpolation over u in [0,1]
        // given u in [0, 1] find what indices in a that falls
        // between. then return the interpolated value between
        // the values of a based on u.
        private static double PLInterp(double u, double[] a)
        {
            if (u <= 0) return a[0];
            if (u >= 1) return a[a.Length - 1];
            double t = u * (a.Length - 1);
            int i = (int)System.Math.Floor(t);
            double f = t - i;
            return a[i] * (1 - f) + a[i + 1] * f;
        }

        // ===============================
        // Polyline cache + projector (2D)
        // ===============================
        internal sealed class CurvePolylineCache
        {
            public readonly double Length;
            private readonly double[] _s;           // arclength samples (0..Length)
            private readonly double[][] _p;         // positions
            private readonly double[][] _t;         // tangents (unit)
            private readonly double[][] _n;         // normals (unit, +90° from T)

            // Build a uniform arclength polyline from a curve
            public CurvePolylineCache(IArcLengthCurve c, int samples = 256)
            {
                if (samples < 2) samples = 2;
                Length = c.Length;                      // assumes IArcLengthCurve exposes Length
                _s = new double[samples];
                _p = new double[samples][];
                _t = new double[samples][];
                _n = new double[samples][];

                for (int i = 0; i < samples; i++)
                {
                    double si = Length * i / (samples - 1.0);
                    var sm = c.Evaluate(si);
                    _s[i] = si;
                    _p[i] = sm.P;
                    _t[i] = sm.T;
                    _n[i] = new double[] { -_t[i][1], _t[i][0] }; // CCW normal in 2D
                }
            }

            // Progressive nearest-point projection: searches a small window around prevIdx
            // Returns: (idx, u, P(u), N(u)) where u is arclength along this cached curve
            public (int idx, double u, double[] pos, double[] nor) Project(double[] q, int prevIdx, int searchWindow = 8)
            {
                int n = _p.Length;
                if (n == 1) return (0, 0, _p[0], _n[0]);

                // Clamp prevIdx to [0, n-2] because we project to segments [i,i+1]
                int i0 = System.Math.Max(0, System.Math.Min(n - 2, prevIdx));
                int lo = System.Math.Max(0, i0 - searchWindow);
                int hi = System.Math.Min(n - 2, i0 + searchWindow);

                double bestD2 = double.PositiveInfinity;
                int bestI = i0;
                double bestAlpha = 0.0;

                for (int i = lo; i <= hi; i++)
                {
                    var a = _p[i];
                    var b = _p[i + 1];
                    double vx = b[0] - a[0], vy = b[1] - a[1];
                    double wx = q[0] - a[0], wy = q[1] - a[1];
                    double vv = vx * vx + vy * vy;
                    double t = vv > 0 ? (vx * wx + vy * wy) / vv : 0.0;
                    if (t < 0) t = 0; else if (t > 1) t = 1;

                    double px = a[0] + t * vx, py = a[1] + t * vy;
                    double dx = q[0] - px, dy = q[1] - py;
                    double d2 = dx * dx + dy * dy;
                    if (d2 < bestD2) { bestD2 = d2; bestI = i; bestAlpha = t; }
                }

                // Interpolate arclength, position, and normal at the footpoint
                double s0 = _s[bestI], s1 = _s[bestI + 1];
                double u = s0 + bestAlpha * (s1 - s0);

                var A = _p[bestI];
                var B = _p[bestI + 1];
                double[] pos = new double[] { A[0] + bestAlpha * (B[0] - A[0]), A[1] + bestAlpha * (B[1] - A[1]) };

                var n0 = _n[bestI];
                var n1 = _n[bestI + 1];
                double[] nor = Normalize(new double[] { n0[0] * (1 - bestAlpha) + n1[0] * bestAlpha,
                                                    n0[1] * (1 - bestAlpha) + n1[1] * bestAlpha });
                return (bestI, u, pos, nor);
            }

            private static double[] Normalize(double[] v)
            {
                double l = System.Math.Sqrt(v[0] * v[0] + v[1] * v[1]);
                return l > 0 ? new double[] { v[0] / l, v[1] / l } : new double[] { 1, 0 };
            }
        }
    }

    // ===================================================
    // CurveCorridorConstraint: centerline + width profile
    // ===================================================
    public sealed class CurveCorridorConstraint : ICompositeConstraint
    {
        public ConstraintType Type { get; } = ConstraintType.Soft;
        public double Weight { get; }

        public void ShowInfo()
        {
            Console.WriteLine($"CurveCorridorConstraint");
        }

        private readonly CurvePolylineCache _lane;
        private readonly double[] _widthProfile; // samples of w over normalized progress u/L ∈ [0,1]
        private readonly double _ds;             // sampling step along trial composite (meters)
        private readonly int _laneSearchWin;

        /// <param name="laneCenter">IArcLengthCurve centerline lane</param>
        /// <param name="widthProfile">piecewise-linear width samples over normalized lane progress</param>
        /// <param name="weight">penalty weight</param>
        /// <param name="ds">trial sampling step along composite (e.g., 0.5 m)</param>
        /// <param name="laneSamples">pre-sample count for the lane polyline cache</param>
        /// <param name="laneSearchWindow">±window size (in segments) for progressive projection</param>
        public CurveCorridorConstraint(IArcLengthCurve laneCenter,
                                       double[] widthProfile,
                                       double weight = 4.0,
                                       double ds = 0.5,
                                       int laneSamples = 512,
                                       int laneSearchWindow = 10)
        {
            if (widthProfile == null || widthProfile.Length < 2)
                throw new ArgumentException("widthProfile must have at least 2 samples.");
            _lane = new CurvePolylineCache(laneCenter, laneSamples);
            _widthProfile = widthProfile;
            Weight = weight;
            _ds = System.Math.Max(1e-3, ds);
            _laneSearchWin = System.Math.Max(2, laneSearchWindow);
        }

        public double[] Residual(IReadOnlyList<CurveSpec> specs)
        {
            // Prefix lengths for the trial composite
            var curves = new CachedIntrinsicCurve[specs.Count];
            var pre = new double[specs.Count + 1];
            pre[0] = 0;
            for (int i = 0; i < specs.Count; i++)
            {
                curves[i] = new CachedIntrinsicCurve(specs[i]);
                pre[i + 1] = pre[i] + specs[i].Length;
            }
            double Ltot = pre[specs.Count];
            if (Ltot <= 0) return new[] { 0.0 };

            var res = new List<double>();
            int lastLaneIdx = 0;

            // Sample uniformly in arclength along the entire composite
            int M = System.Math.Max(4, (int)System.Math.Ceiling(Ltot / _ds) + 1);
            for (int k = 0; k < M; k++)
            {
                double sGlob = Ltot * k / (M - 1.0);

                // Map to segment-local s
                int seg = UpperBound(pre, sGlob) - 1;                // seg with pre[seg] <= sGlob < pre[seg+1]
                double sLoc = sGlob - pre[seg];

                var P = curves[seg].Position(sLoc);
                // Progressive nearest point on lane center
                var pr = _lane.Project(P, lastLaneIdx, _laneSearchWin);
                lastLaneIdx = pr.idx;

                // Lateral offset is dot with lane normal at the station
                double dx = P[0] - pr.pos[0], dy = P[1] - pr.pos[1];
                double delta = dx * pr.nor[0] + dy * pr.nor[1];

                // Width at normalized lane progress
                double uNorm = (_lane.Length > 0) ? pr.u / _lane.Length : 0.0;
                double w = PLInterp(uNorm, _widthProfile);
                double half = 0.5 * System.Math.Max(0.0, w);

                double pen = 0.0;
                double absd = System.Math.Abs(delta);
                if (absd > half) pen = absd - half;

                if (pen > 0) res.Add(Weight * pen);
            }

            return res.Count == 0 ? new[] { 0.0 } : res.ToArray();
        }

        private static int UpperBound(double[] a, double x)
        {
            int lo = 0, hi = a.Length;
            while (lo < hi)
            {
                int mid = (lo + hi) >> 1;
                if (x < a[mid]) hi = mid; else lo = mid + 1;
            }
            return lo;
        }

        private static double PLInterp(double u, double[] v)
        {
            if (u <= 0) return v[0];
            if (u >= 1) return v[v.Length - 1];
            double t = u * (v.Length - 1);
            int i = (int)System.Math.Floor(t);
            double f = t - i;
            return v[i] * (1 - f) + v[i + 1] * f;
        }
    }

    // ==================================================
    // BoundedCurveConstraint: stay between two curves
    // (left and right boundaries; signed-distance hinge)
    // ==================================================
    public sealed class BoundedCurveConstraint : ICompositeConstraint
    {
        public ConstraintType Type { get; } = ConstraintType.Soft;
        public double Weight { get; }

        private readonly CurvePolylineCache _left;
        private readonly CurvePolylineCache _right;
        private readonly int _searchWin;
        private readonly double _ds;

        public void ShowInfo()
        {
            Type t = GetType();
            Console.WriteLine($"{t.Name}");
        }

        // inward normals (global flips applied at construction)
        private int _hintL = 0, _hintR = 0;
        private int _flipL = +1, _flipR = +1;

        /// <param name="leftBoundary">Left boundary (looking forward along the track)</param>
        /// <param name="rightBoundary">Right boundary</param>
        /// <param name="weight">penalty weight</param>
        /// <param name="ds">trial sampling step along composite</param>
        /// <param name="samplesPerBound">polyline samples per boundary curve</param>
        /// <param name="searchWindow">±window size for progressive projection</param>
        public BoundedCurveConstraint(IArcLengthCurve leftBoundary,
                                      IArcLengthCurve rightBoundary,
                                      double weight = 4.0,
                                      double ds = 0.5,
                                      int samplesPerBound = 512,
                                      int searchWindow = 10)
        {
            _left = new CurvePolylineCache(leftBoundary, samplesPerBound);
            _right = new CurvePolylineCache(rightBoundary, samplesPerBound);
            _searchWin = System.Math.Max(2, searchWindow);
            _ds = System.Math.Max(1e-3, ds);
            Weight = weight;

            // Calibrate inward normals (one-time): make left normals point toward right side and vice versa
            CalibrateNormals();
        }

        public double[] Residual(IReadOnlyList<CurveSpec> specs)
        {
            var curves = new CachedIntrinsicCurve[specs.Count];
            var pre = new double[specs.Count + 1];
            pre[0] = 0;
            for (int i = 0; i < specs.Count; i++)
            {
                curves[i] = new CachedIntrinsicCurve(specs[i]);
                pre[i + 1] = pre[i] + specs[i].Length;
            }
            double Ltot = pre[specs.Count];
            if (Ltot <= 0) return new[] { 0.0 };

            var res = new List<double>();
            for (int k = 0, M = System.Math.Max(4, (int)System.Math.Ceiling(Ltot / _ds) + 1); k < M; k++)
            {
                double sGlob = Ltot * k / (M - 1.0);
                int seg = UpperBound(pre, sGlob) - 1;
                double sLoc = sGlob - pre[seg];

                var P = curves[seg].Position(sLoc);

                // Signed distance to left (inward normal)
                var prL = _left.Project(P, _hintL, _searchWin);
                _hintL = prL.idx;
                double dxL = P[0] - prL.pos[0], dyL = P[1] - prL.pos[1];
                double dL = (dxL * prL.nor[0] + dyL * prL.nor[1]) * _flipL;

                // Signed distance to right (inward normal)
                var prR = _right.Project(P, _hintR, _searchWin);
                _hintR = prR.idx;
                double dxR = P[0] - prR.pos[0], dyR = P[1] - prR.pos[1];
                double dR = (dxR * prR.nor[0] + dyR * prR.nor[1]) * _flipR;

                // Inside if both >= 0; otherwise hinge back to 0
                double pen = 0.0;
                if (dL < 0) pen += -dL;
                if (dR < 0) pen += -dR;

                if (pen > 0) res.Add(Weight * pen);
            }

            return res.Count == 0 ? new[] { 0.0 } : res.ToArray();
        }

        private void CalibrateNormals()
        {
            // Take a few stations on left; ensure its normal points toward right boundary
            for (int i = 0; i < 5; i++)
            {
                double u = _left.Length * (i + 0.5) / 5.0;
                var smL = SampleOn(_left, u);
                var prR = _right.Project(smL.pos, 0, _searchWin);
                double vx = prR.pos[0] - smL.pos[0], vy = prR.pos[1] - smL.pos[1];
                if (vx * smL.nor[0] + vy * smL.nor[1] < 0) { _flipL = -1; break; }
            }
            // Do the symmetric check for right
            for (int i = 0; i < 5; i++)
            {
                double u = _right.Length * (i + 0.5) / 5.0;
                var smR = SampleOn(_right, u);
                var prL = _left.Project(smR.pos, 0, _searchWin);
                double vx = prL.pos[0] - smR.pos[0], vy = prL.pos[1] - smR.pos[1];
                if (vx * smR.nor[0] + vy * smR.nor[1] < 0) { _flipR = -1; break; }
            }
        }

        private static (double[] pos, double[] nor) SampleOn(CurvePolylineCache c, double u)
        {
            // simple nearest sample (good enough for calibration)
            int n = 256;
            int idx = (int)System.Math.Round((n - 1) * System.Math.Max(0, System.Math.Min(1, u / System.Math.Max(1e-9, c.Length))));
            idx = System.Math.Max(0, System.Math.Min(n - 1, idx));
            // Re-project onto cache for precise values
            var pr = c.Project(new double[] { 0, 0 }, idx, 0); // we don't use point, but want index clamped
                                                               // Use stored arrays
            return (pr.pos, pr.nor);
        }

        private static int UpperBound(double[] a, double x)
        {
            int lo = 0, hi = a.Length;
            while (lo < hi)
            {
                int mid = (lo + hi) >> 1;
                if (x < a[mid]) hi = mid; else lo = mid + 1;
            }
            return lo;
        }
    }
}
