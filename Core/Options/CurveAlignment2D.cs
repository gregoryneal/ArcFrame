using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;
using ArcFrame.Core.Results;

namespace ArcFrame.Core.Options
{
    /// <summary>
    /// Generate RigidTransform to align curves in 2D, rotate along axis and translate along 2D plane. Calculate center of mass and best rotation angle for curve alignment.
    /// </summary>
    public static class CurveAlignment2D
    {
        /// <summary>
        /// Align a curve to a polyline in the XZ plane.
        /// </summary>
        /// <param name="proto">The curve to be transformed, projected on to the XZ plane</param>
        /// <param name="polyline">The polyline to fit the curve to</param>
        /// <param name="samples">Number of samples to weight the pre transformed curve</param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static TransformedCurve AlignCOM_ByXZ(IArcLengthCurve proto, (double x, double z)[] polyline, int samples = 128)
        {
            if (proto.Dimension < 3) throw new ArgumentException("Curve must be at least 3D to rotate around the Y axis and translate along the XZ plane.");
            if (polyline == null || polyline.Length < 2) throw new ArgumentException("need at least 2 target points for alignment.");

            var A = new (double x, double z)[samples];
            double s;
            Sample sp;
            for (int i = 0; i < samples; i++)
            {
                s = (i == samples - 1) ? proto.Length : i * (proto.Length / (samples - 1));
                sp = proto.Evaluate(s);
                A[i] = (sp.P[0], sp.P[2]);
            }

            //Center of masses
            (double x, double z) comA = (0, 0), comB = (0, 0);
            foreach (var a in A)
            {
                comA.x += a.x;
                comA.z += a.z;
            }
            foreach (var b in polyline)
            {
                comB.x += b.x;
                comB.z += b.z;
            }
            comA.x /= samples;
            comA.z /= samples;
            comB.x /= polyline.Length;
            comB.z /= polyline.Length;

            // Centered sets (pair target to source by index, wrapping/clamping as needed)
            double Sxx = 0, Sxy = 0;
            int M = System.Math.Min(samples, polyline.Length);
            for (int i = 0; i < M; i++)
            {
                var a = (x: A[i].x - comA.x, z: A[i].z - comA.z);
                var b = (x: polyline[i].x - comB.x, z: polyline[i].z - comB.z);
                Sxx += a.x * b.x + a.z * b.z;         // dot
                Sxy += a.x * b.z - a.z * b.x;         // cross (signed area term)
            }
            double angle = System.Math.Atan2(Sxy, Sxx);

            // Build yaw transform about COM_A, then shift to COM_B
            var xfYaw = RigidTransform.FromYawXZ(proto.Dimension, angle, comA.x, comA.z);
            // after rotation about comA, we want comA to land on comB:
            // t_extra = comB - (R * comA + t) -> but FromYawXZ already rotates about comA,
            // so just add the delta in XZ.
            var t = (double[])xfYaw.T.Clone();
            t[0] += comB.x - (xfYaw.R[0, 0] * comA.x + xfYaw.R[0, 2] * comA.z + xfYaw.T[0]);
            t[2] += comB.z - (xfYaw.R[2, 0] * comA.x + xfYaw.R[2, 2] * comA.z + xfYaw.T[2]);
            var xf = new RigidTransform(xfYaw.R, t);

            return new TransformedCurve(proto, xf);
        }

        /// <summary>
        /// Same as AlignCOM_ByXZ(), but takes samples at equal arc lengths along the polyline and generated curve.
        /// </summary>
        /// <param name="proto"></param>
        /// <param name="target"></param>
        /// <param name="closed"></param>
        /// <param name="weightBySegments"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentNullException"></exception>
        /// <exception cref="ArgumentException"></exception>
        public static TransformedCurve AlignCOM_ByArcLength_XZ(IArcLengthCurve proto, (double x, double z)[] target, bool closed = false, bool weightBySegments = true)
        {
            if (proto == null) throw new ArgumentNullException(nameof(proto));
            if (proto.Dimension < 3) throw new ArgumentException("curve must be at least 3D (uses XZ plane).");
            if (target == null || target.Length < 2) throw new ArgumentException("need at least 2 target points.");

            // 1) arc-length of target polyline
            var L = CumulativeArcLengthsXZ(target, closed, out var cum, out var segLen);
            if (L <= 0) throw new ArgumentException("target polyline has zero length.");

            // 2) sample proto at the SAME arc-length fractions
            int M = target.Length;
            var A = new (double x, double z)[M];
            for (int i = 0; i < M; i++)
            {
                double f = cum[i] / L;                      // fraction along target
                double s = f * System.Math.Max(0.0, proto.Length); // corresponding arclength on proto
                var sp = proto.Evaluate(s);
                A[i] = (sp.P[0], sp.P[2]);
            }

            // 3) optional circular-shift search for CLOSED curves (phase alignment)
            int bestShift = 0;
            if (closed && M <= 4096) // cheap brute force; lift limit if needed
            {
                double bestScore = double.NegativeInfinity;
                for (int shift = 0; shift < M; shift++)
                {
                    var (Sxx, Sxy, wsum) = CrossStats(A, target, segLen, shift, weightBySegments);
                    double score = (wsum > 0) ? (Sxx * Sxx + Sxy * Sxy) : double.NegativeInfinity;
                    if (score > bestScore) { bestScore = score; bestShift = shift; }
                }
            }

            // 4) compute weighted Procrustes (2D rotation + translation) with the chosen pairing
            var (phi, tx, tz) = BestYawAndTranslation(A, target, segLen, bestShift, weightBySegments);

            // 5) build rigid transform (yaw in XZ, translation)
            var xf = RigidTransform.FromYawXZ(proto.Dimension, phi, 0, 0, tx, tz);

            return new TransformedCurve(proto, xf);
        }

        // --- helpers -------------------------------------------------------------

        // cumulative arc lengths; cum[i] is length from start up to point i (0..M-1)
        static double CumulativeArcLengthsXZ((double x, double z)[] pts, bool closed, out double[] cum, out double[] segLen)
        {
            int M = pts.Length;
            segLen = new double[M]; // segment i = dist(pts[i], pts[i+1]) (wrapping if closed)
            cum = new double[M];

            double total = 0.0;
            for (int i = 0; i < M - 1; i++)
            {
                double dx = pts[i + 1].x - pts[i].x;
                double dz = pts[i + 1].z - pts[i].z;
                segLen[i] = System.Math.Sqrt(dx * dx + dz * dz);
                total += segLen[i];
                cum[i + 1] = total;
            }
            if (closed)
            {
                double dx = pts[0].x - pts[M - 1].x;
                double dz = pts[0].z - pts[M - 1].z;
                segLen[M - 1] = System.Math.Sqrt(dx * dx + dz * dz);
                total += segLen[M - 1];
                // for closed curves we still leave cum[M-1] as accumulated (no extra element)
            }
            else segLen[M - 1] = 0.0;

            return total;
        }

        // Weighted cross stats for 2D Procrustes on XZ plane, with optional circular shift.
        static (double Sxx, double Sxy, double wsum) CrossStats((double x, double z)[] A, (double x, double z)[] B, double[] segLen, int shift, bool weightBySegments)
        {
            int M = A.Length;
            // centers (weighted if desired)
            double wsum = 0, Ax = 0, Az = 0, Bx = 0, Bz = 0;
            for (int i = 0; i < M; i++)
            {
                double w = weightBySegments ? SegWeight(segLen, i) : 1.0;
                int j = (i + shift) % M;
                Ax += w * A[i].x; Az += w * A[i].z;
                Bx += w * B[j].x; Bz += w * B[j].z;
                wsum += w;
            }
            if (wsum <= 0) return (0, 0, 0);
            Ax /= wsum; Az /= wsum; Bx /= wsum; Bz /= wsum;

            // cross-covariance terms
            double Sxx = 0, Sxy = 0;
            for (int i = 0; i < M; i++)
            {
                double w = weightBySegments ? SegWeight(segLen, i) : 1.0;
                int j = (i + shift) % M;
                double ax = A[i].x - Ax, az = A[i].z - Az;
                double bx = B[j].x - Bx, bz = B[j].z - Bz;
                Sxx += w * (ax * bx + az * bz);      // dot
                Sxy += w * (ax * bz - az * bx);      // cross
            }
            return (Sxx, Sxy, wsum);
        }

        /// <summary>
        /// Find best yaw rotation and XZ translation to fit data set A to B
        /// </summary>
        /// <param name="A">Data set to fit</param>
        /// <param name="B">Target data set</param>
        /// <param name="segLen"></param>
        /// <param name="shift"></param>
        /// <param name="weightBySegments"></param>
        /// <returns></returns>
        static (double phi, double tx, double tz) BestYawAndTranslation((double x, double z)[] A, (double x, double z)[] B, double[] segLen, int shift, bool weightBySegments)
        {
            // rotation
            var (Sxx, Sxy, wsum) = CrossStats(A, B, segLen, shift, weightBySegments);
            if (wsum <= 0) return (0, 0, 0);
            double phi = System.Math.Atan2(Sxy, Sxx);
            double c = System.Math.Cos(phi), s = System.Math.Sin(phi);

            // translation: t = COM_B - R * COM_A (same weights as above)
            int M = A.Length;
            double wtot = 0, Ax = 0, Az = 0, Bx = 0, Bz = 0;
            for (int i = 0; i < M; i++)
            {
                double w = weightBySegments ? SegWeight(segLen, i) : 1.0;
                int j = (i + shift) % M;
                Ax += w * A[i].x; Az += w * A[i].z;
                Bx += w * B[j].x; Bz += w * B[j].z;
                wtot += w;
            }
            Ax /= wtot; Az /= wtot; Bx /= wtot; Bz /= wtot;

            double rx = c * Ax - s * Az;
            double rz = s * Ax + c * Az;
            double tx = Bx - rx;
            double tz = Bz - rz;

            return (phi, tx, tz);
        }

        static double SegWeight(double[] segLen, int i)
        {
            // weight point i by the average length of adjacent segments (endpoints get one side)
            int M = segLen.Length;
            if (M == 1) return 1.0;
            double w = segLen[i];
            if (i > 0) w = 0.5 * (w + segLen[i - 1]);
            return System.Math.Max(w, 1e-12);
        }
    }
}
