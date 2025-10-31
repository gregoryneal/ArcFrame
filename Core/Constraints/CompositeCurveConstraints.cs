using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using System.Collections.Generic;
using System;

namespace ArcFrame.Core.Constraints
{
    internal sealed class CurvePolylineCache
    {
        public readonly double Length;
        public readonly double[] S;             // arclength samples (0..Length)
        public readonly double[][] P;           // positions
        public readonly double[][,] R;          // ON frames


        // Build a uniform arclength polyline from a curve
        public CurvePolylineCache(IArcLengthCurve c, int samples = 256)
        {
            if (samples < 2) samples = 2;
            Length = c.Length;
            S = new double[samples];
            P = new double[samples][];
            R = new double[samples][,];

            for (int i = 0; i < samples; i++)
            {
                double si = Length * i / (samples - 1.0);
                var sm = c.Evaluate(si);
                S[i] = si;
                P[i] = sm.P;
                R[i] = sm.R;
            }
        }

        // Progressive nearest-point projection: searches a small window around prevIdx
        /// <summary>
        /// Find the cached sample closest to the point q. 
        /// Only searches around the window around prevIdx. 
        /// </summary>
        /// <param name="q">The point to test</param>
        /// <param name="prevIdx">The index of this cached curve to search around</param>
        /// <param name="searchWindow"></param>
        /// <returns>(idx, u, P(u), N(u)) where u is arclength along this cached curve</returns>
        public (int idx, double u, double[] pos, double[,] frame) Project(double[] q, int prevIdx, int searchWindow = 8)
        {
            int n = P.Length;
            if (n == 1) return (0, 0, P[0], R[0]);

            // Clamp prevIdx to [0, n-2] because we project to segments [i,i+1]
            int i0 = System.Math.Max(0, System.Math.Min(n - 2, prevIdx));
            int lo = System.Math.Max(0, i0 - searchWindow);
            int hi = System.Math.Min(n - 2, i0 + searchWindow);

            double bestD2 = double.PositiveInfinity;
            int bestI = i0;
            double bestAlpha = 0.0;

            for (int i = lo; i <= hi; i++)
            {
                var a = P[i];
                var b = P[i + 1];
                
                double[] bSubA = Helpers.Subtract(b, a);
                //double vx = b[0] - a[0];
                //double vy = b[1] - a[1];
                double[] qSubA = Helpers.Subtract(q, a);
                //double wx = q[0] - a[0];
                //double wy = q[1] - a[1];
                double vv = Helpers.Dot(bSubA, bSubA);
                //double vv = vx * vx + vy * vy;
                double t = vv > 0 ? Helpers.Dot(bSubA, qSubA) / vv : 0;
                //double t = vv > 0 ? (vx * wx + vy * wy) / vv : 0.0;

                if (t < 0) t = 0;
                else if (t > 1) t = 1;

                // Interpolate between the two ends of the search window by t in [0, 1]
                // PP = a + t(b - a)
                double[] PP = Helpers.Add(a, Helpers.Multiply(t, bSubA));
                // DP = q - PP
                double[] DP = Helpers.Subtract(q, PP);
                double d2 = Helpers.Dot(DP, DP);
                //double px = a[0] + t * vx, py = a[1] + t * vy;
                //double dx = q[0] - px, dy = q[1] - py;
                //double d2 = dx * dx + dy * dy;
                if (d2 < bestD2) { bestD2 = d2; bestI = i; bestAlpha = t; }
            }

            // Interpolate arclength, position, and normal at the footpoint
            double s0 = S[bestI], s1 = S[bestI + 1];
            double u = s0 + bestAlpha * (s1 - s0);

            var A = P[bestI];
            var B = P[bestI + 1];
            double[] pos = Helpers.Add(A, Helpers.Multiply(bestAlpha, Helpers.Subtract(B, A)));
            //double[] pos = new double[] { A[0] + bestAlpha * (B[0] - A[0]), A[1] + bestAlpha * (B[1] - A[1]) };

            var R0 = R[bestI];
            var R1 = R[bestI + 1];
            double[,] nor = Helpers.Add(R0, Helpers.Multiply(bestAlpha, Helpers.Subtract(R1, R0)));
            //double[] nor = Normalize(new double[] { n0[0] * (1 - bestAlpha) + n1[0] * bestAlpha, n0[1] * (1 - bestAlpha) + n1[1] * bestAlpha });
            return (bestI, u, pos, nor);
        }
    }

    /// <summary>
    /// Adapter to use any ICurveConstraint in a composite problem.
    /// </summary>
    public sealed class SegmentConstraintWrapper : ICompositeConstraint
    {
        private readonly ICurveConstraint _inner;
        /// <inheritdoc/>
        public SegmentConstraintWrapper(ICurveConstraint inner) { _inner = inner; }
        /// <inheritdoc/>
        public ConstraintType Type => _inner.Type;
        /// <inheritdoc/>
        public double Weight => _inner.Weight;
        /// <inheritdoc/>
        public double[] Residual(IReadOnlyList<CurveSpec> specs)
        {
            var res = _inner.Residual(specs[_inner.SegmentIndex]);
            if (Weight != 1.0)
                for (int i = 0; i < res.Length; i++) res[i] *= Weight;
            return res;
        }
        /// <inheritdoc/>
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
        /// <inheritdoc/>
        public ConstraintType Type { get; } = ConstraintType.Soft;
        /// <inheritdoc/>
        public double Weight { get; }
        private readonly int _segIdx;
        /// <inheritdoc/>
        public LengthRegularizer(int segIdx, double weight) { _segIdx = segIdx; Weight = weight; }
        /// <inheritdoc/>
        public double[] Residual(IReadOnlyList<CurveSpec> specs)
        {
            return new[] { Weight * specs[_segIdx].Length };
        }
        /// <inheritdoc/>
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
        /// <inheritdoc/>
        public ConstraintType Type { get; } = ConstraintType.Hard;
        /// <inheritdoc/>
        public double Weight { get; } = 1.0;
        /// <summary>
        /// The first segment index
        /// </summary>
        public int LeftIndex { get; }     // segment i
        /// <summary>
        /// The second segment index
        /// </summary>
        public int RightIndex { get; }    // segment i+1
        /// <summary>
        /// Apply position constraint at joint
        /// </summary>
        public bool C0 { get; } = true;
        /// <summary>
        /// Apply tangent constraint at joint
        /// </summary>
        public bool C1 { get; } = true;
        /// <summary>
        /// Apply curvature constraint at joint
        /// </summary>
        public bool C2 { get; } = false;
        /// <summary>
        /// Weight of position residual
        /// </summary>
        public double wC0 { get; } = 1.0;
        /// <summary>
        /// Weight of tangent residual
        /// </summary>
        public double wC1 { get; } = 1.0;
        /// <summary>
        /// Weight of curvature residual
        /// </summary>
        public double wC2 { get; } = 1.0;
        /// <inheritdoc/>
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
        /// <inheritdoc/>
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
        /// <inheritdoc/>
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
        /// <inheritdoc/>
        public ConstraintType Type => ConstraintType.Soft;
        /// <inheritdoc/>
        public double Weight { get; }
        /// <inheritdoc/>
        public SlopeJumpRegularizer(double weight) { Weight = weight; }
        /// <inheritdoc/>
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
        /// <inheritdoc/>
        public void ShowInfo()
        {
            Type t = GetType();
            Console.WriteLine($"SlopeJumpRegularizer");
        }
    }

    // ==================================================
    // BoundedCurveConstraint: stay between two curves
    // (left and right boundaries; signed-distance hinge)
    // ==================================================
    /// <summary>
    /// Stay between two IArcLengthCurves. Curves and bounds
    /// must be in 2D XY plane, I think... TODO: Test bounded curves in 3D
    /// </summary>
    public sealed class BoundedCurveConstraint : ICompositeConstraint
    {
        /// <inheritdoc/>
        public ConstraintType Type { get; } = ConstraintType.Soft;
        /// <inheritdoc/>
        public double Weight { get; }

        private readonly CurvePolylineCache _left;
        private readonly CurvePolylineCache _right;
        private readonly int _searchWin;
        private readonly double _ds;
        private readonly double _buffer;
        /// <inheritdoc/>
        public void ShowInfo()
        {
            Type t = GetType();
            Console.WriteLine($"{t.Name}");
        }

        // inward normals (global flips applied at construction)
        private int _hintL = 0, _hintR = 0;
        private int _flipL = +1, _flipR = +1;

        /// <summary>
        /// Bound a curve between two other curves.
        /// </summary>
        /// <param name="leftBoundary">Left boundary (looking forward along the track)</param>
        /// <param name="rightBoundary">Right boundary</param>
        /// <param name="buffer">the minimum acceptable buffer inside of the boundary</param>
        /// <param name="weight">penalty weight</param>
        /// <param name="ds">trial sampling step along composite</param>
        /// <param name="samplesPerBound">polyline samples per boundary curve</param>
        /// <param name="searchWindow">±window size for progressive projection</param>
        public BoundedCurveConstraint(IArcLengthCurve leftBoundary,
                                      IArcLengthCurve rightBoundary,
                                      double buffer = 1,
                                      double weight = 4.0,
                                      double ds = 0.5,
                                      int samplesPerBound = 512,
                                      int searchWindow = 10)
        {
            _left = new CurvePolylineCache(leftBoundary, samplesPerBound);
            _right = new CurvePolylineCache(rightBoundary, samplesPerBound);
            _searchWin = System.Math.Max(2, searchWindow);
            _ds = System.Math.Max(1e-3, ds);
            _buffer = System.Math.Abs(buffer);
            Weight = weight;

            // Calibrate inward normals (one time): make left normals point toward right side and vice versa
            CalibrateNormals();

            // Test normals by exporting them
            /*
            double[][] normal_left;
            double[][] normal_right;
            int[] index = [0, _left.S.Length / 2, _left.S.Length - 1];
            
            for (int i = 0; i < index.Length; i++)
            {
                int j = index[i];
                double[] leftNormal = Helpers.Multiply(ONFrame.GetCol(_left.R[j], 1), _flipL);
                double[] rightNormal = Helpers.Multiply(ONFrame.GetCol(_right.R[j], 1), _flipR);

                double[] leftPos = _left.P[j];
                double[] rightPos = _right.P[j];
                double[] leftP2 = Helpers.Add(leftPos, leftNormal);
                double[] rightP2 = Helpers.Add(rightPos, rightNormal);

                Console.Write($"{j} => Left{i}0: ");
                Helpers.PrintVector(leftPos);

                Console.Write($"{j} => Left{i}1: ");
                Helpers.PrintVector(leftP2);

                Console.Write($"{j} => Right{i}0: ");
                Helpers.PrintVector(rightPos);

                Console.Write($"{j} => Right{i}1: ");
                Helpers.PrintVector(rightP2);
            }*/
        }
        /// <inheritdoc/>
        public double[] Residual(IReadOnlyList<CurveSpec> specs)
        {
            var curves = new CachedIntrinsicCurve[specs.Count];
            CompositeCurve curve = new CompositeCurve();
            // length prefixes (arc length at the start of each curve segment)
            // includes the total arc length at specs.Count
            var pre = new double[specs.Count + 1];
            pre[0] = 0;
            for (int i = 0; i < specs.Count; i++)
            {
                curves[i] = new CachedIntrinsicCurve(specs[i]);
                curve.AddG1(curves[i], out _);
                pre[i + 1] = pre[i] + specs[i].Length;
            }

            double Ltot = pre[specs.Count];
            // Break our trial curve up into M pieces (minimum 4) and iterate over them
            var M = System.Math.Max(4, (int)System.Math.Ceiling(Ltot / _ds));
            var res = new double[M+1];
            if (Ltot == 0) return res;
            if (Ltot < 0)
            {
                //discourage negative lengths
                for (int i = 0; i < M; i++) res[i] = -Ltot;
                return res;
            }

            _hintL = 0;
            _hintR = 0;
            
            for (int k = 0; k <= M; k++)
            {
                // Arc length at the start of the kth segment
                double sGlob = Ltot * k / M;
                // Map the global arc length at the kth segment 
                // to the index of the prefix array built earlier.
                //int seg = UpperBound(pre, sGlob) - 1;
                //double sLoc = sGlob - pre[seg];

                // Remap the index of the prefix array
                // to the index of the curve it belongs to
                // in the curve array
                //seg = System.Math.Max(0, seg - 1);

                // Test point for the residual
                double[] P;
                //Console.WriteLine($"seg: {seg} | pre length: {pre.Length} | sGlob: {sGlob} | Ltot: {Ltot} | k: {k} | M: {M} | Curves length: {curves.Length} | sLoc: {sLoc}");
                //Helpers.PrintVector(pre);

                P = curve.Position(sGlob);

                //if (seg < curves.Length) P = curve.Position(sGlob);
                //else P = curves[^1].Position(curves[^1].Length);


                // Find the closest point around from P to the search window in _left
                // In the window around index _hintL. Update _hintL so it doesn't fall
                // Behind out of a favorable search window.
                var (hint, _, pos_l, r_l) = _left.Project(P, _hintL, _searchWin);
                _hintL = hint;
                // (P - pos_l) • nor_l => >0 it is on the inside (same direction), <0 on the outside
                double[] n = Helpers.Multiply(ONFrame.GetCol(r_l, 1), _flipL);
                pos_l = Helpers.Add(pos_l, Helpers.Multiply(n, _buffer));
                double dL = Helpers.Dot(Helpers.Subtract(P, pos_l), n);
                //Console.WriteLine($"Bounded Curve left penalty: {dL}");
                //double dxL = P[0] - pos_l[0], dyL = P[1] - pos_l[1];
                //double dL = (dxL * nor_l[0] + dyL * nor_l[1]) * _flipL;

                // Signed distance to right (inward normal)
                var (hint2, _, pos_r, r_r) = _right.Project(P, _hintR, _searchWin);
                _hintR = hint2;
                n = Helpers.Multiply(ONFrame.GetCol(r_r, 1), _flipR);
                pos_r = Helpers.Add(pos_r, Helpers.Multiply(n, _buffer));
                double dR = Helpers.Dot(Helpers.Subtract(P, pos_r), n);
                //Console.WriteLine($"Bounded Curve right penalty: {dR}");
                //double dxR = P[0] - pos_r[0], dyR = P[1] - pos_r[1];
                //double dR = (dxR * nor_r[0] + dyR * nor_r[1]) * _flipR;

                // Inside if both >= 0; otherwise hinge back to 0
                double pen = 0.0;
                if (dL < 0) pen += -dL;
                if (dR < 0) pen += -dR;

                //Console.WriteLine($"Bounded Curve penalty at P: {pen}");
                //Helpers.PrintVector(P);

                if (pen > 0) res[k] = Weight * pen;
            }

            //Console.WriteLine("BoundedCurveConstraint residuals: ");
            //Helpers.PrintVector([.. res]);
            return res;
        }

        private void CalibrateNormals()
        {
            // Take a few stations on left; ensure its normal points toward right boundary
            int idx = 0;
            for (int i = 0; i < 5; i++)
            {
                double u = _left.Length * (i + 0.5) / 5.0;
                var (p_l, frame) = SampleOn(_left, u);
                double[] nor = ONFrame.GetCol(frame, 1);
                var (idx2, _, p_r, _) = _right.Project(p_l, idx, _searchWin);
                idx = idx2;
                double[] dv = Helpers.Subtract(p_r, p_l);
                double dot = Helpers.Dot(dv, nor);
                /*Console.WriteLine($"p_right - p_left: ");
                Helpers.PrintVector(p_r);
                Helpers.PrintVector(p_l);
                Helpers.PrintVector(dv);
                Helpers.PrintVector(nor);
                Console.WriteLine($"Left Normal dot: {dot}");*/
                if (dot < 0) { _flipL = -1; break; }
            }
            idx = 0;
            // Do the symmetric check for right
            for (int i = 0; i < 5; i++)
            {
                double u = _right.Length * (i + 0.5) / 5.0;
                var (pos, frame) = SampleOn(_right, u);
                double[] nor = ONFrame.GetCol(frame, 1);
                var (idx2, _, p, _) = _left.Project(pos, idx, _searchWin);
                idx = idx2;
                double[] dv = Helpers.Subtract(p, pos);
                if (Helpers.Dot(dv, nor) < 0) { _flipR = -1; break; }
            }

            //Console.WriteLine($"_flipL / _flipR: {_flipL} / {_flipR}");
        }

        // Find the sample with arc length closest to u in c.
        private static (double[] pos, double[,] frame) SampleOn(CurvePolylineCache c, double u)
        {
            int index = UpperBound(c.S, u) - 1;
            return (c.P[index], c.R[index]);
            /*
            // simple nearest sample (good enough for calibration)
            int n = 256;
            int idx = (int)System.Math.Round((n - 1) * System.Math.Max(0, System.Math.Min(1, u / System.Math.Max(1e-9, c.Length))));
            idx = System.Math.Max(0, System.Math.Min(n - 1, idx));
            // Re-project onto cache for precise values
            var pr = c.Project(new double[] { 0, 0 }, idx, 0); // we don't use point, but want index clamped
                                                               // Use stored arrays
            return (pr.pos, pr.frame);*/
        }

        /// <summary>
        /// Run a binary search on a.
        /// For any value x in [a[i], a[i+1])
        /// Return i+1.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        private static int UpperBound(double[] a, double x)
        {
            int lo = 0;
            int hi = a.Length;
            while (lo < hi)
            {
                int mid = (lo + hi) >> 1;
                if (x < a[mid]) hi = mid; 
                else lo = mid + 1;
            }
            return lo;
        }
    }
}
