using ArcFrame.Core.Math;
using System.Runtime.CompilerServices;

namespace ArcFrame.Core.Kernels
{
    /// <summary>
    /// Pointer version of <see cref="Spline3dDescriptor"/>
    /// </summary>
    public unsafe struct Spline3dRaw
    {
#pragma warning disable CS1591
        public int Degree;
        public int SegmentCount;
        public int ControlPointCount;
        public int ArcSampleCount;
        public double Length;
        public FrameModel Frame;
        public bool FastMode;

        // Flattened GB blocks: SegmentCount blocks of (3 x (Degree+1))
        public double* SegmentCoeff;

        // Arc-length mapping samples
        public double* ArcS;  // length = ArcSampleCount
        public double* ArcT;  // length = ArcSampleCount

        // Optional constant frame (FastMode path)
        public Mat3d R0;
#pragma warning restore
    }

    /// <summary>
    /// Kernel based spline evaluation methods
    /// </summary>
    public static unsafe class Spline3dKernels
    {
        /// <summary>
        /// Evaluate the 3d spline at its arc length s
        /// </summary>
        /// <param name="curve"></param>
        /// <param name="s"></param>
        /// <param name="p"></param>
        /// <param name="R"></param>
        /// <param name="kappa"></param>
        /// <param name="torsion"></param>
        public static void EvaluateAtS(
            in Spline3dRaw curve,
            double s,
            out Vec3d p,
            out Mat3d R,
            out double kappa,
            out double torsion)
        {
            // Clamp s
            if (s < 0) s = 0;
            if (s > curve.Length) s = curve.Length;

            // 1. s -> t using table
            double t = MapStoT(curve.ArcS, curve.ArcT, curve.ArcSampleCount, s);

            // 2. t -> segment index + local u
            int cpCount = curve.ControlPointCount;
            int degree = curve.Degree;
            int segCount = curve.SegmentCount; // same as your Locate() math

            double segReal = t * segCount;
            int seg = segReal < 0 ? 0 : (int)System.Math.Floor(segReal);
            if (seg >= segCount) seg = segCount - 1;
            double u = segReal - seg; // [0,1]

            // 3. Grab GB block for this segment (3 x (degree+1))
            int cols = degree + 1;
            int segStride = 3 * cols;
            double* GB = curve.SegmentCoeff + seg * segStride;

            // 4. Compute power monomials and their derivatives
            double* U0 = stackalloc double[cols]; // [1, u, u^2, ...]
            double* U1 = stackalloc double[cols]; // first derivative monomials
            double* U2 = stackalloc double[cols];
            double* U3 = stackalloc double[cols];

            PowerMonomialsAndDerivs(degree, u, U0, U1, U2, U3);

            // 5. Evaluate p(u), p'(u), p''(u), p'''(u)
            Vec3d r0 = Eval3(GB, degree, U0);
            Vec3d r1 = Eval3(GB, degree, U1);
            Vec3d r2 = Eval3(GB, degree, U2);
            Vec3d r3 = Eval3(GB, degree, U3);

            p = r0;
            // normalize tangent
            double tlen = r1.Length;
            Vec3d T = new Vec3d(r1.x, r1.y, r1.z) / tlen;

            if (curve.FastMode)
            {
                // Fast path: tangent only, zero curvature.
                R = curve.R0;
                R.m00 = T.x;
                R.m10 = T.y;
                R.m20 = T.z;
                kappa = 0.0;
                torsion = 0.0;
                return;
            }

            // 6. Frenet frame + curvature/torsion in 3D
            BuildFrenetFrameAndCurvature(in r1, in r2, in r3, out R, out kappa, out torsion);

            // TODO: switch on curve.Frame and
            // plug a 3D Bishop kernel in the same style.
        }

        // --- helpers used by both Burst and managed code --- //

        /// <summary>
        /// Map s to t uniform sampling.
        /// </summary>
        /// <param name="sSamples"></param>
        /// <param name="tSamples"></param>
        /// <param name="n"></param>
        /// <param name="s"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double MapStoT(double* sSamples, double* tSamples, int n, double s)
        {
            if (n <= 0) return 0.0;
            if (n == 1) return tSamples[0];

            // Clamp
            if (s <= sSamples[0]) return tSamples[0];
            if (s >= sSamples[n - 1]) return tSamples[n - 1];

            // Binary search for i such that sSamples[i] <= s < sSamples[i+1]
            int lo = 0;
            int hi = n - 1;
            while (hi - lo > 1)
            {
                int mid = (lo + hi) >> 1;
                if (sSamples[mid] <= s)
                    lo = mid;
                else
                    hi = mid;
            }

            double s0 = sSamples[lo];
            double s1 = sSamples[lo + 1];
            double t0 = tSamples[lo];
            double t1 = tSamples[lo + 1];

            double alpha = (s - s0) / (s1 - s0);
            return t0 + (t1 - t0) * alpha;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void PowerMonomialsAndDerivs(
            int degree,
            double u,
            double* U0,
            double* U1,
            double* U2,
            double* U3)
        {
            // U0[i] = u^i
            double p = 1.0;
            for (int i = 0; i <= degree; i++)
            {
                U0[i] = p;
                p *= u;
            }

            // U1: coefficients for first derivative: d/du of sum c_i u^i
            // So U1[i] corresponds to coefficient for c_i in r'(u): i * u^(i-1)
            U1[0] = 0.0;
            for (int i = 1; i <= degree; i++)
                U1[i] = i * U0[i - 1];

            // U2: second derivative: i*(i-1) u^(i-2)
            U2[0] = 0.0;
            if (degree >= 1) U2[1] = 0.0;
            for (int i = 2; i <= degree; i++)
                U2[i] = i * (i - 1) * U0[i - 2];

            // U3: third derivative
            U3[0] = 0.0;
            if (degree >= 1) U3[1] = 0.0;
            if (degree >= 2) U3[2] = 0.0;
            for (int i = 3; i <= degree; i++)
                U3[i] = i * (i - 1) * (i - 2) * U0[i - 3];
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static Vec3d Eval3(double* GB, int degree, double* U)
        {
            // GB is 3 x (degree+1), stored row-major in a contiguous block:
            // [x0,x1,...,xd, y0,y1,...,yd, z0,z1,...,zd]
            int cols = degree + 1;

            double px = 0, py = 0, pz = 0;

            int yOffset = cols;
            int zOffset = 2 * cols;

            for (int j = 0; j < cols; j++)
            {
                double uj = U[j];
                px += GB[j] * uj;
                py += GB[yOffset + j] * uj;
                pz += GB[zOffset + j] * uj;
            }

            return new Vec3d(px, py, pz);
        }

        private static void BuildFrenetFrameAndCurvature(
            in Vec3d r1, in Vec3d r2, in Vec3d r3,
            out Mat3d R,
            out double kappa,
            out double torsion)
        {
            // Classic 3D formulas with derivatives wrt parameter u:
            // κ = ||r' × r''|| / ||r'||^3
            // τ = det(r', r'', r''') / ||r' × r''||^2

            var Traw = r1;
            double lenR1 = HpcMath3d.Length(Traw);
            if (lenR1 <= 1e-12)
            {
                R = Mat3d.Identity;
                kappa = 0.0;
                torsion = 0.0;
                return;
            }

            Vec3d T = HpcMath3d.Normalize(Traw);
            Vec3d cross12 = HpcMath3d.Cross(r1, r2);
            double crossLen = HpcMath3d.Length(cross12);

            if (crossLen <= 1e-12)
            {
                // Nearly straight line: well-defined T, arbitrary N/B, zero curvature/torsion
                R = Mat3d.Identity;
                R.m00 = T.x; R.m10 = T.y; R.m20 = T.z;
                R.m11 = 1.0;
                R.m22 = 1.0;
                kappa = 0.0;
                torsion = 0.0;
                return;
            }

            kappa = crossLen / (lenR1 * lenR1 * lenR1);

            Vec3d B = new Vec3d(cross12.x / crossLen, cross12.y / crossLen, cross12.z / crossLen);
            Vec3d N = HpcMath3d.Cross(B, T);

            // Torsion via determinant
            double det = HpcMath3d.Det(r1, r2, r3);
            double denom = crossLen * crossLen;
            torsion = denom > 1e-12 ? det / denom : 0.0;

            // Build frame matrix with columns T, N, B
            R = default;
            R.m00 = T.x; R.m10 = T.y; R.m20 = T.z;
            R.m01 = N.x; R.m11 = N.y; R.m21 = N.z;
            R.m02 = B.x; R.m12 = B.y; R.m22 = B.z;
        }
    }
}
