using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using ArcFrame.Core.Results;
using System;
using System.Collections.Generic;

namespace ArcFrame.Core.Geometry
{
    /// <summary>
    /// Builder class for joining together IArcLengthCurve segments.
    /// It can handle automatic G1 continuity.
    /// </summary>
    public class CompositeCurve : IArcLengthCurve
    {
        /// <summary>
        /// Number of segments.
        /// </summary>
        public int Count => _segs.Count;
        private List<IArcLengthCurve> _segs = new List<IArcLengthCurve>();
        //the dimension of the curve, enforced across segments
        private int? _n;
        //total length of the curve
        private double _len;
        /// <summary>
        /// Arc length up to the start of the ith index.
        /// _prefix[0] = 0;
        /// </summary>
        private double[] _prefix = { };
        //call EnsurePrefix if true
        private bool _dirty = true;
        /// <inheritdoc/>
        public int Dimension
        {
            get
            {
                if (_n.HasValue) return _n.Value;
                if (_segs.Count == 0) return 0;
                _n = _segs[0].Dimension;
                return _n.Value;
            }
        }
        /// <inheritdoc/>
        public double Length
        {
            get
            {
                EnsurePrefix();
                return _len;
            }
        }

        /// <summary>
        /// Index the CompositeCurve, return the indexed segment.
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public IArcLengthCurve this[int index]
        {
            get
            {
                return _segs[index];
            }
        }
        /// <inheritdoc/>
        public double[] Position(double s) => Evaluate(s).P;
        /// <inheritdoc/>
        public double[] Tangent(double s) => Evaluate(s).T;

        /// <summary>
        /// Evaluate along the entire length of the curve.
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        public Sample Evaluate(double s)
        {
            //Console.WriteLine("Eval composite curve");
            EnsurePrefix();
            //Console.WriteLine("Ensure prefix composite curve");
            if (_segs.Count == 0)
            {
                //Console.WriteLine("Seg count is 0");
                double[,] R = new double[1, 1];
                return new Sample(new double[System.Math.Max(1, Dimension)], R, 0, new double[System.Math.Max(0, Dimension - 1)]);
            }
            if (s <= 0) s = 0;
            if (s >= _len) s = _len;

            int low = 0;
            int high = _segs.Count;
            int mid = 0;
            int u = 0;
            //binary search for segment
            while (low + 1 < high && ++u < _segs.Count)
            {
                mid = (low + high) >> 1;
                if (s >= _prefix[mid]) low = mid;
                else high = mid;
            }
            //Console.WriteLine("Binary search composite curve");

            /*if (u >= _segs.Count)
            {
                Console.WriteLine($"Binary search for composite segment failed! SegmentCount: {_segs.Count} | HighIndex: {high} | MidIndex: {mid} | LowIndex: {low} | s: {s} | Prefix: ");
                Helpers.PrintVector(_prefix);
            }*/

            double localS = s - _prefix[low];
            Sample sp = _segs[low].Evaluate(localS);
            //return sample with composite arc length not local
            return new Sample(sp.P, sp.R, s, sp.k);
        }

        /// <summary>
        /// Get a number of uniform samples along the curve.
        /// </summary>
        /// <param name="count"></param>
        /// <returns></returns>
        public Sample[] GetSamples(int count)
        {
            EnsurePrefix();
            if (count <= 0) count = 1;
            Sample[] pts = new Sample[count];
            double s;
            for (int i = 0; i < count; i++)
            {
                s = (i == count - 1) ? _len : i * (_len / System.Math.Max(1, count - 1));
                pts[i] = Evaluate(s);
            }

            return pts;
        }

        /// <summary>
        /// Add a segment with a full ON matching frame.
        /// </summary>
        /// <param name="segment"></param>
        /// <param name="xfUsed"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public CompositeCurve AddG1FullFrame(IArcLengthCurve segment, out RigidTransform xfUsed)
        {
            if (_segs.Count == 0)
            {
                xfUsed = RigidTransform.Identity(segment.Dimension);
                return Add(new TransformedCurve(segment, xfUsed));
            }
            if (segment.Dimension != Dimension)
                throw new ArgumentException("Dimension mismatch! New segment must match curve dimension.");

            EnsurePrefix();

            var prev = _segs[^1];
            var prevEnd = prev.Evaluate(prev.Length); // contains P, T, R, k, s
            var newStart = segment.Evaluate(0);

            // Full-frame alignment: find R so that R * R_new0 = R_prevEnd
            double[,] Ralign = Helpers.Multiply(prevEnd.R, Helpers.Transpose(newStart.R));

            // Translation so that the rotated start position sits on the previous end position
            double[] pRot = Helpers.Multiply(Ralign, newStart.P);
            double[] T = new double[Dimension];
            for (int i = 0; i < Dimension; i++) T[i] = prevEnd.P[i] - pRot[i];

            xfUsed = new RigidTransform(Ralign, T);
            return Add(new TransformedCurve(segment, xfUsed));
        }

        /// <summary>
        /// Takes a segment and creates a TransformedCurve that holds the input segment and a RigidTransform which ensures the new segment endpoint is aligned
        /// with the previous endpoint, and the new segment tangent is aligned with the previous tangent.
        /// </summary>
        /// <param name="segment"></param>
        /// <param name="xfUsed">The transform that was used to join the segment to the curve</param>
        /// <returns></returns>
        public CompositeCurve AddG1(IArcLengthCurve segment, out RigidTransform xfUsed)
        {
            if (_segs.Count == 0)
            {
                xfUsed = RigidTransform.Identity(segment.Dimension);
                return Add(new TransformedCurve(segment, xfUsed));
            }
            if (segment.Dimension != Dimension) throw new ArgumentException("Dimension mismatch! New segment must be the same dimensions as the input curve.");

            EnsurePrefix();
            IArcLengthCurve prev = _segs[^1];
            Sample prevEnd = prev.Evaluate(prev.Length);
            Sample newStart = segment.Evaluate(0);

            //R * new.T = prev.T
            double[,] TRotation = MinimumRotationAligning(prevEnd.T, newStart.T);
            //Rotate in place to find translation offset
            double[] rotatedPosition = Helpers.Multiply(TRotation, newStart.P);
            double[] newT = new double[Dimension];
            for (int i = 0; i < Dimension; i++)
            {
                newT[i] = prevEnd.P[i] - rotatedPosition[i];
            }

            xfUsed = new RigidTransform(TRotation, newT);
            return Add(new TransformedCurve(segment, xfUsed));
        }

        /// <summary>
        /// Adds a segment to the list.
        /// </summary>
        /// <param name="segment"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public CompositeCurve Add(IArcLengthCurve segment)
        {
            if (segment == null) throw new ArgumentException("Null segment");
            _segs.Add(segment);
            if (_n.HasValue && segment.Dimension != _n.Value)
            {
                throw new ArgumentException("Segment has different dimension than curve");
            }
            if (!_n.HasValue)
            {
                _n = segment.Dimension;
            }
            _dirty = true;
            return this;
        }

        /// <summary>
        /// Add a list of segments to this CompositeCurve
        /// </summary>
        /// <param name="segments"></param>
        /// <param name="g1">Enforce G1 joins (discard transformations)</param>
        /// <returns></returns>
        public CompositeCurve AddRange(IEnumerable<IArcLengthCurve> segments, bool g1 = true)
        {
            foreach (IArcLengthCurve segment in segments)
            {
                if (g1) AddG1(segment, out _);
                else Add(segment);
            }
            return this;
        }

        /// <summary>
        /// Count up the accumulated arc length at the start of each segment and store them in the _prefix array.
        /// </summary>
        private void EnsurePrefix()
        {
            if (!_dirty) return;
            int m = _segs.Count;
            _prefix = new double[m + 1];
            _prefix[0] = 0;
            for (int i = 1; i < _prefix.Length; i++)
            {
                _prefix[i] = _prefix[i - 1] + _segs[i - 1].Length;
            }
            _len = _prefix[m];
            _dirty = false;
        }

        /// <summary>
        /// Ensure R and t are compatible
        /// </summary>
        /// <param name="R"></param>
        /// <param name="t"></param>
        private void CheckRig(double[,] R, double[] t)
        {
            if (R == null || t == null) throw new ArgumentNullException();
            int n = t.Length;
            if (R.GetLength(0) != n || t.Length != n) throw new ArgumentException("R must be nxn and t must be length n");
            if (_n.HasValue && _n.Value != n)
            {
                throw new ArgumentException($"Dimension mismatch. Composite curve dimension: {_n.Value}, R/t dimension: {n}.");
            }
        }

        /// <summary>
        /// Find the minimum rotation matrix that aligns vector a with vector b.
        /// </summary>
        /// <param name="aRaw"></param>
        /// <param name="bRaw"></param>
        /// <returns></returns>
        private double[,] MinimumRotationAligning(double[] aRaw, double[] bRaw)
        {
            double tol = 1E-12;
            int n = aRaw.Length;
            var a = Helpers.Normalize((double[])aRaw.Clone());
            var b = Helpers.Normalize((double[])bRaw.Clone());
            double c = System.Math.Clamp(Helpers.Dot(a, b), -1.0, 1.0);       // cos(theta)
            if (c > 1 - tol) return RigidTransform.Identity(n).R;

            //build orthonormal basis vectors u1 and u2 such that span(u1, u2) = span(a, b)
            double[] u1 = a;
            double[] u2;
            double[] tmp = Helpers.Reject(b, a);
            double s = Helpers.Len(tmp);
            double dot;
            if (s <= tol)
            {
                //a and b are almost parallel
                u2 = Helpers.StdSeed(a);
                //make sure u2 is perpendicular to u1 (a)
                dot = Helpers.Dot(u2, a);
                for (int i = 0; i < n; i++)
                {
                    u2[i] -= dot * a[i];
                }
                Helpers.Normalize(u2);
                s = 0.0;
                c = -1.0;
            }
            else
            {
                u2 = Helpers.Multiply(1.0 / s, tmp);
            }

            //R = I + (cosθ - 1)(uu^T + vv^T) + sinθ(vu^T - uv^T)
            //rodriguez formula in N dimensions.
            //https://math.stackexchange.com/a/5024118
            double[,] R = RigidTransform.Identity(n).R;
            Helpers.AddOuterScaled(R, u1, u1, (c - 1));
            Helpers.AddOuterScaled(R, u2, u2, (c - 1));
            Helpers.AddOuterScaled(R, u2, u1, -s);
            Helpers.AddOuterScaled(R, u1, u2, s);

            return R;
        }

        /// <summary>
        /// Craft a CompositeCurve from a list of specs.
        /// </summary>
        /// <param name="specs"></param>
        /// <returns></returns>
        public static CompositeCurve FromSpecList(CurveSpec[] specs)
        {
            CompositeCurve c = new CompositeCurve();
            for (int i = 0; i < specs.Length; i++)
            {
                c.AddG1(specs[i].GetOptimizedCurve(), out _);
            }
            return c;
        }
    }
}
