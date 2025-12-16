using ArcFrame.Core.Geometry;
using ArcFrame.Core.Geometry.Splines;
using ArcFrame.Core.Kernels;
using ArcFrame.Core.Results;
using System;
using System.Collections.Generic;


namespace ArcFrame.Core.Math.Geometry.Splines
{
    /// <summary>
    /// General spline class using the matrix definition of a spline.
    /// Meant to be overridden with specific spline implementations.
    /// </summary>
    public class Spline : IArcLengthCurve
    {
        /// <summary>
        /// How many segments should Burst use?
        /// </summary>
        //public virtual int BurstSegmentCount => SegmentCount;

        /// <summary>
        /// Map segment index -> ControlPoint index
        /// </summary>
        /// <param name="segIndex"></param>
        /// <returns></returns>
        public virtual int GetControlPointStartForSegment(int segIndex) => segIndex;

        /// <summary>
        /// The number of spline segments. Default: sliding window 
        /// </summary>
        public virtual int ComputeSegmentCount() => System.Math.Max(1, ControlPoints.Length - Degree);

        /// <summary>
        /// The spline control points.
        /// </summary>
        public readonly double[][] ControlPoints;
        /// <summary>
        /// The arc length table which caches the curve on creation.
        /// </summary>
        protected readonly ArcLengthTable _table;
        /// <summary>
        /// Accessor for the internal ArcLengthTable
        /// </summary>
        public ArcLengthTable Table => _table;
        /// <summary>
        /// The frame model, frenet or RMF.
        /// </summary>
        public readonly FrameModel Frame;
        /// <summary>
        /// The start ON frame.
        /// </summary>
        public readonly double[,] R0;
        /// <summary>
        /// The cached ON frame along the curve.
        /// </summary>
        protected readonly double[][,] RCache;
        /// <summary>
        /// The cached tangents along the curve.
        /// </summary>
        protected readonly double[] T0;
        /// <summary>
        /// The cached normals along the curve.
        /// </summary>
        protected readonly double[]? N0; // used for RMF
        /// <summary>
        /// Number of cache samples to take.
        /// </summary>
        protected readonly int CacheSamples;
        /// <summary>
        /// Determines the complexity and smoothness of the curve.
        /// This is the highest power of the polynomial used to construct it.
        /// The degree must be less than or equal to the spatial Dimension
        /// parameter. This is because we build an initial SO(N) frame R0.
        /// </summary>
        public int Degree { get; }
        /// <summary>
        /// The basis matrix for the spline model.
        /// </summary>
        public double[,] BasisMatrix { get; }

        /// <summary>
        /// The spatial dimension the curve exists in.
        /// </summary>
        public int Dimension { get; }
        /// <summary>
        /// Cached arc length from the arc length table.
        /// </summary>
        public double Length => _table.Length;

        /// <summary>
        /// Fast mode disables frame calculations, ensuring only curve position is calculated.
        /// This is way faster because splines are kind of hacky in this library ngl.
        /// </summary>
        public bool FastMode { get; }

        /// <summary>
        /// Precomputed GB matrices per segment start index (cpi).
        /// Index = control point index returned by Locate(t).
        /// Each matrix is Dimension x (Degree+1).
        /// </summary>
        private readonly double[][,] _gbCache;

        /// <summary>
        /// Get the computed GB cache, G is the control point matrix, and B the basis matrix.
        /// </summary>
        public double[][,] GBCache => _gbCache;

        /// <summary>
        /// Scratch buffer for [1, u, u^2, ..., u^Degree].
        /// Reused across Position/Evaluate calls (not thread-safe).
        /// </summary>
        private readonly double[] _uScratch;

        /// <summary>
        /// Scratch buffer for RMF interpolation blending (Dimension x Dimension).
        /// Only used in non-fast Bishop mode.
        /// </summary>
        private readonly double[,] _rmfBlendScratch;

        /// <summary>
        /// Generate an arc length table
        /// </summary>
        public Spline(double[][] controlPoints, double[,] basisMatrix, FrameModel frame = FrameModel.Frenet, bool fastMode = false, int cacheSamplesOverride = 0)
        {
            //if (controlPoints.Length <= 0) return;
            BasisMatrix = basisMatrix;
            Degree = BasisMatrix.GetLength(0) - 1;

            ControlPoints = controlPoints;
            Dimension = ControlPoints[0].Length;
            Frame = frame;
            FastMode = fastMode;

            // ------------------------------------------------------------------
            // Precompute GB for all used segments
            // ------------------------------------------------------------------
            int _numSegments = ComputeSegmentCount();
            _gbCache = new double[_numSegments][,];
            for (int i = 0; i < _numSegments; i++)
            {
                int cpi = GetControlPointStartForSegment(i);
                double[,] Gseg = CreateG(cpi);
                _gbCache[i] = Helpers.Multiply(Gseg, BasisMatrix);
            }
            /*
            int maxSegStart = System.Math.Max(0, ControlPoints.Length - (Degree + 1));
            _gbCache = new double[maxSegStart + 1][,];

            for (int cpi = 0; cpi <= maxSegStart; cpi++)
            {
                // This uses the virtual CreateG, so Hermite/Bezier/etc. all
                // get the correct control point layout.
                double[,] Gseg = CreateG(cpi);
                _gbCache[cpi] = Helpers.Multiply(Gseg, BasisMatrix);
            }*/

            // Scratch for monomials and RMF interpolation
            _uScratch = new double[Degree + 1];
            _rmfBlendScratch = new double[Dimension, Dimension];

            // Scale arc length table number of samples by polyline length of control points
            double len = 0;
            for (int i = 0; i + 1 < controlPoints.Length; i++)
            {
                len += Helpers.Len(Helpers.Subtract(controlPoints[i + 1], controlPoints[i]));
            }

            // N samples per unit arc length
            if (cacheSamplesOverride > 0) CacheSamples = cacheSamplesOverride;
            else if (FastMode) CacheSamples = System.Math.Max(32, (int)System.Math.Round(len) * 64 + 1);
            else CacheSamples = System.Math.Max(64, (int)System.Math.Round(len) * 256 + 1);

            _table = new ArcLengthTable(Position, CacheSamples);

            if (!FastMode) {
                RCache = new double[CacheSamples][,];

                // FD derivative of the spline at t = 0
                double[][] derivatives = new double[Degree][];
                // Build ON frame at t = 0
                R0 = new double[Dimension, Degree];
                T0 = Helpers.Differentiate(Position, 0, 1, Dimension);
                ONFrame.SetCol(R0, 0, Helpers.Normalize(T0));
                derivatives[0] = T0;

                // Calculate derivatives beyong tangent
                for (int d = 1; d < Degree; d++)
                {
                    //d+1th derivative
                    derivatives[d] = Helpers.Differentiate(Position, 0, d + 1, Dimension);
                    ONFrame.SetCol(R0, d, Helpers.Normalize(derivatives[d]));
                }
                // Find the other ON bases if they exist, now R0 is square
                if (Dimension > Degree) R0 = ONFrame.R0_FromR_Complete(R0);

                // RMF frame is the same as frenet frame in Dimension <= 2
                if (Frame == FrameModel.Bishop && Dimension > 2)
                {
                    // In RMF we cache the frame along the curve
                    // and do a lookup + interpolation in Evaluate
                    // only if the curve is at least degree 2, a degree 1 curve 
                    // would just return the tangent between the correct control
                    // nodes in RMF and Frenet, since it would just be a polyline
                    if (Degree > 1)
                    {
                        N0 = ONFrame.GetCol(R0, 1);
                        double[] B0 = ONFrame.GetCol(R0, 2);
                        double[] previousP = Position(0);
                        double[] previousT = T0;
                        double[] previousN = N0;
                        double[] newT;
                        double[] calcT;
                        double[] newN;
                        double[] newB;
                        double[] v1;
                        double[] v2;

                        RCache[0] = R0;
                        for (int i = 1; i < CacheSamples; i++)
                        {
                            // Construct the sequential frame that minimizes rotation of N0 from the previous frame.
                            double t = (double)i / (CacheSamples - 1);

                            double[] p = Position(t);
                            v1 = Helpers.Subtract(p, previousP);
                            calcT = Helpers.Normalize(v1);
                            newT = Helpers.Subtract(previousT, Helpers.Multiply(2 * Helpers.Dot(previousT, v1) / Helpers.Dot(v1, v1), v1));
                            newN = Helpers.Subtract(previousN, Helpers.Multiply(2 * Helpers.Dot(previousN, v1) / Helpers.Dot(v1, v1), v1));

                            v2 = Helpers.Subtract(calcT, newT);
                            if (Helpers.Len(v2) > 1E-8)
                            {
                                newN = Helpers.Subtract(newN, Helpers.Multiply(2 * Helpers.Dot(newN, v2) / Helpers.Dot(v2, v2), v2));
                            }

                            newB = Helpers.Cross3(calcT, newN);
                            RCache[i] = ONFrame.R0_FromTNB_Complete(calcT, newN, newB);

                            previousN = newN;
                            previousT = calcT;
                            previousP = p;
                        }
                    }
                }
            }
            else
            {
                // minimal: just a fixed frame, no caches
                RCache = System.Array.Empty<double[,]>();
                //R0 = RigidTransform.Identity(Dimension).R; // or from first segment’s tangent
                T0 = Helpers.Differentiate(Position, 0, 1, Dimension);// new double[Dimension];       // optional
                R0 = ONFrame.R0_FromT_Complete(T0);
                N0 = null;
            }
        }

        

        /// <summary>
        /// Evaluating the Spline curve is done in a few steps.
        /// Goal: p(t) = GBt'
        /// p = position
        /// t = parameter in [0, 1]
        /// t' = power polynomial [1, t, t^2, t^3, ..., t^degree]^T
        /// G = degree + 1 control point column vectors whose first index is floor(t*(controlpoints.length - degree))
        /// B = basis matrix, controls the type of spline being generated.
        /// 
        /// Step 1: Get parameter t in [0, 1] from s in [0, Length] (s / Length)
        /// Step 2: Get degree+1 control points from t. Put this in Nx(d+1) matrix G (N is dimension).
        /// Step 3: Form t' with successive powers of t up to t^degree
        /// Step 4: Multiply the matrices to get the position p(t) = GBt'
        /// Step 5: Find all derivatives of the spline, and calculate the general curvatures from it.
        /// 5a) if our ON frame is Frenet and F = [T, N, B, ...] and our derivative frame is F' = [T', N', B', ...]
        ///     our general curvature matrix is K = [N•T', B•N', ...] or in general K[i] = F[i+1]•F'[i] for 0 <= i < N.
        /// 5b) if our ON frame is Bishop and F = [T, N, B, ... ] we must interpolate a cached array of rotation minimized frames.
        /// </summary>
        /// <param name="s">The arc length</param>
        /// <returns></returns>
        public Sample Evaluate(double s)
        {
            s = System.Math.Clamp(s, 0, Length);

            // Step 1.
            double t = _table.MapSToTUniform(s);//_table.MapStoT(s);
            (int seg, double u) = Locate(t);
            // Step 2.
            //double[,] G = CreateG(seg);
            double[,] GB = GetSegmentGB(seg);// Helpers.Multiply(G, BasisMatrix);
            // Step 4.
            double[] P = EvaluateSegmentPosition(seg, u);//Helpers.Multiply(GB, tp);
            // Step 5.
            double[,] R;
            double[] k;
            if (FastMode)
            {
                // We need T always for frame orientation
                var T = Helpers.Normalize(Helpers.Differentiate(Position, t, 1, 3));
                R = ONFrame.R0_FromT_Complete(T);
                k = new double[System.Math.Max(0, Dimension - 1)];
                return new Sample(P, R, s, k);
            }

            if (Frame == FrameModel.Bishop && RCache != null && RCache.Length == CacheSamples && RCache[0] != null)
            {
                // Interpolate cached RMF and re-orthonormalize (continuity-preserving)
                R = InterpolateCachedRMF(t);
            }
            else
            {
                // Build Frenet-like frame at s (any N) using successive u-derivatives
                R = BuildFrenetFrameAtS(GB, Degree, u, Dimension);
            }
            k = new double[System.Math.Max(1, Dimension - 1)]; // Curvatures are filled in different depending on the type of frame used, Frenet or RMF.

            if (Frame == FrameModel.Bishop)
            {
                double[] U1 = PowerMonomialsDerivN(Degree, 1, u);
                double[] U2 = PowerMonomialsDerivN(Degree, 2, u);

                // --- first two u-derivatives ---
                double[] v = Helpers.Multiply(GB, U1);   // dp/du
                double[] a = Helpers.Multiply(GB, U2);   // d^2p/du^2

                // --- Tangent and curvature vector ingredients (dimension-agnostic) ---
                double L = Helpers.Len(v);
                double invL = (L > 1e-16) ? 1.0 / L : 0.0;
                double invL2 = invL * invL;

                // Unit tangent T
                double[] T = (L > 1e-16) ? Helpers.Multiply(invL, v) : Helpers.Normalize(v);

                // dT/ds = (I - T T^T) * a / ||v||^2   (works in any N)
                double[] a_par_T = Helpers.Multiply(Helpers.Dot(a, T), T);
                double[] a_perp = Helpers.Subtract(a, a_par_T);
                double[] dTds = Helpers.Multiply(invL2, a_perp);

                // RMF: curvature components are simply dT/ds projected onto normal columns of R
                for (int i = 0; i < k.Length; i++)
                {
                    double[] Ni = ONFrame.GetCol(R, i + 1);
                    k[i] = Helpers.Dot(dTds, Ni);
                }
            }
            else
            {
                // Frenet (general N): extract κ_i from the connection Ω = R^T R'
                // Use small arc-length step h based on the arc-length table resolution.
                double h = System.Math.Max(Length / System.Math.Max(2, CacheSamples - 1), 1e-6 * System.Math.Max(1.0, Length));

                double sMinus = System.Math.Max(0.0, s - h);
                double sPlus = System.Math.Min(Length, s + h);
                if (sPlus == sMinus) sPlus = System.Math.Min(Length, sMinus + 1e-6);

                double[,] Rm = BuildFrameOnly(sMinus);  // no k; same method as above
                double[,] Rp = BuildFrameOnly(sPlus);

                // Align signs to avoid spurious flips (per-column)
                ProcrustesColumnSignFix(R, Rm);
                ProcrustesColumnSignFix(R, Rp);

                double invDelta = 1.0 / (sPlus - sMinus);
                double[,] Rprime = new double[Dimension, Dimension];
                for (int i = 0; i < Dimension; i++)
                    for (int j = 0; j < Dimension; j++)
                        Rprime[i, j] = (Rp[i, j] - Rm[i, j]) * invDelta;

                // Ω = R^T R'  (skew-symmetric for orthonormal R)
                double[,] Rt = Helpers.Transpose(R);
                double[,] Om = Helpers.Multiply(Rt, Rprime);

                // κ_i = |Ω_{i,i+1}|  (Cartan curvatures in general N)
                for (int i = 0; i < k.Length; i++)
                    k[i] = System.Math.Abs(Om[i, i + 1]);
            }


            return new Sample(P, R, s, k);
        }

        /// <summary>
        /// Test the 3d spline kernel
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        /// <exception cref="InvalidOperationException"></exception>
        public Frame3d Evaluate3dKernel(double s)
        {
            if (Dimension != 3)
                throw new InvalidOperationException("Evaluate3dKernel only valid for Dimension == 3.");

            // TODO: cache this
            var desc = Spline3dDescriptorFactory.FromSpline(this);

            unsafe
            {
                fixed (double* pSeg = desc.SegmentCoeff)
                fixed (double* pS = desc.ArcS)
                fixed (double* pT = desc.ArcT)
                {
                    var raw = new Spline3dRaw
                    {
                        Degree = desc.Degree,
                        SegmentCount = desc.SegmentCount,
                        ControlPointCount = desc.ControlPointCount,
                        ArcSampleCount = desc.ArcS.Length,
                        Length = desc.Length,
                        Frame = desc.Frame,
                        FastMode = desc.FastMode,
                        SegmentCoeff = pSeg,
                        ArcS = pS,
                        ArcT = pT,
                        R0 = desc.R0Flat != null ? Mat3d.FromRowMajor(desc.R0Flat) : Mat3d.Identity
                    };

                    Spline3dKernels.EvaluateAtS(in raw, s, out var p, out var R, out _, out _);

                    return new Frame3d { P = p, R = R };
                }
            }
        }

        /// <summary>
        /// Find a position along the spline given its parameter t in [0, 1]
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        public double[] Position(double t)
        {
            t = System.Math.Clamp(t, 0, 1);
            (int seg, double localS) = Locate(t);
            return EvaluateSegmentPosition(seg, localS);
        }

        /// <summary>
        /// Create a vector [1, u, u^2, ..., u^d]
        /// </summary>
        /// <param name="d"></param>
        /// <param name="u"></param>
        /// <returns></returns>
        protected static double[] PowerMonomials(int d, double u)
        {
            var U = new double[d + 1];
            double p = 1.0;
            for (int i = 0; i <= d; i++) { U[i] = p; p *= u; }
            return U;
        }

        /// <summary>
        /// Take the Nth derivative of the dth degree power monomial [1, u, u^2, ..., u^d]
        /// </summary>
        /// <param name="d"></param>
        /// <param name="N"></param>
        /// <param name="u"></param>
        /// <returns></returns>
        protected static double[] PowerMonomialsDerivN(int d, int N, double u)
        {
            double[] UN = new double[d + 1];
            if (N > d + 1) return UN; // all zeroes
            double p = 1.0;
            double coeff;
            // Everything before index N has derivative 0
            for (int i = N; i <= d; i++)
            {
                coeff = 1;
                // i * (i - 1) * (i - 2) * ... * (i - (N - 1))
                for (int j = 0; j < N; j++)
                {
                    coeff *= i - j;
                }
                UN[i] = coeff * p;
                p *= u;
            }
            return UN;
        }

        /// <summary>
        /// Find the starting control point given an interpolant along the curve.
        /// </summary>
        /// <param name="t"></param>
        /// <returns>Start control point index k, local interpolation along spline segment u.</returns>
        protected virtual (int k, double u) Locate(double t)
        {
            t = System.Math.Clamp(t, 0.0, 1.0);

            int segCount = ComputeSegmentCount();
            double seg = t * segCount;

            int i = System.Math.Min(segCount - 1, (int)System.Math.Floor(seg));
            double u = seg - i;

            // k is a segment index, not a CP index
            return (i, u);
            /*(
            if (t > 1) t = 1;
            if (t < 0) t = 0;
            double seg = t * (ControlPoints.Length - Degree);
            int k = System.Math.Min(ControlPoints.Length - (Degree + 1), (int)System.Math.Floor(seg));
            return (k, seg - k);*/
        }

        /// <summary>
        /// Create an Nx(degree+1) matrix of control points column vectors
        /// </summary>
        /// <param name="cpi"></param>
        /// <returns></returns>
        protected virtual double[,] CreateG(int cpi)
        {
            int n = Degree + 1;
            double[,] G = new double[Dimension, n];
            for (int col = 0; col < n; col++)
            {
                //Console.WriteLine($"Control point index: {cpi + col} | length of Control points: {ControlPoints.Length}");
                var P = ControlPoints[cpi + col];
                for (int row = 0; row < Dimension; row++) G[row, col] = P[row];
            }
            return G;
        }

        double[,] InterpolateCachedRMF(double t)
        {
            double pos = t * (CacheSamples - 1);
            int i0 = System.Math.Clamp((int)System.Math.Floor(pos), 0, CacheSamples - 1);
            int i1 = System.Math.Clamp(i0 + 1, 0, CacheSamples - 1);
            double a = System.Math.Clamp(pos - i0, 0, 1);

            double[,] R0c = RCache[i0];
            double[,] R1c = RCache[i1];

            double[,] Rlin = _rmfBlendScratch;
            for (int r = 0; r < Dimension; r++)
                for (int c = 0; c < Dimension; c++)
                    Rlin[r, c] = (1 - a) * R0c[r, c] + a * R1c[r, c];

            // Orthonormalize to fix interpolation drift
            return ONFrame.R0_FromR_Complete(Rlin);
        }
        double[,] BuildFrenetFrameAtS(double[,] GBm, int deg, double uu, int dim)
        {
            // Gather successive u-derivatives p^(k)(u) up to the dimension (or degree), then Gram–Schmidt.
            int maxK = System.Math.Min(deg, dim + 2); // a couple extra helps stability
            List<double[]> derivs = new List<double[]>(maxK);
            for (int ord = 1; ord <= maxK; ord++)
            {
                double[] Uk = PowerMonomialsDerivN(deg, ord, uu);
                derivs.Add(Helpers.Multiply(GBm, Uk));
            }

            // Modified Gram–Schmidt on {v, a, b, ...}
            double[,] Rloc = new double[dim, dim];
            int col = 0;
            for (int idx = 0; idx < derivs.Count && col < dim; idx++)
            {
                double[] w = (double[])derivs[idx].Clone();
                for (int j = 0; j < col; j++)
                {
                    double[] ej = ONFrame.GetCol(Rloc, j);
                    w = Helpers.Subtract(w, Helpers.Multiply(Helpers.Dot(w, ej), ej));
                }
                double nw = Helpers.Len(w);
                if (nw > 1e-12)
                {
                    w = Helpers.Multiply(1.0 / nw, w);
                    ONFrame.SetCol(Rloc, col, w);
                    col++;
                }
            }
            // Complete to SO(N) with a stable algorithm (uses R0 if available internally)
            return ONFrame.R0_FromR_Complete(Rloc);
        }

        double[,] BuildFrameOnly(double ss)
        {
            double tt = _table.MapSToTUniform(ss);//_table.MapStoT(ss);
            (int sseg, double suu) = Locate(tt);
            double[,] GBs = GetSegmentGB(sseg);
            return BuildFrenetFrameAtS(GBs, Degree, suu, Dimension);
        }

        void ProcrustesColumnSignFix(double[,] Rref, double[,] RtoFix)
        {
            int cols = Rref.GetLength(1);
            for (int c = 0; c < cols; c++)
            {
                double dp = 0.0;
                for (int r = 0; r < Rref.GetLength(0); r++) dp += Rref[r, c] * RtoFix[r, c];
                if (dp < 0)
                {
                    for (int r = 0; r < Rref.GetLength(0); r++) RtoFix[r, c] = -RtoFix[r, c];
                }
            }
        }

        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
        private double[] EvaluateSegmentPosition(int segStartIndex, double u)
        {
            // Clamp seg index defensively
            if (segStartIndex < 0) segStartIndex = 0;
            if (segStartIndex >= _gbCache.Length) segStartIndex = _gbCache.Length - 1;

            double[,] GB = _gbCache[segStartIndex];

            // Fill [1, u, u^2, ..., u^Degree] into scratch buffer
            double p = 1.0;
            for (int i = 0; i <= Degree; i++)
            {
                _uScratch[i] = p;
                p *= u;
            }

            // Multiply GB (Dimension x (Degree+1)) by _uScratch ((Degree+1) vector)
            double[] P = new double[Dimension];
            for (int r = 0; r < Dimension; r++)
            {
                double acc = 0.0;
                for (int c = 0; c <= Degree; c++)
                {
                    acc += GB[r, c] * _uScratch[c];
                }
                P[r] = acc;
            }

            return P;
        }

        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
        private double[,] GetSegmentGB(int segStartIndex)
        {
            if (segStartIndex < 0) segStartIndex = 0;
            if (segStartIndex >= _gbCache.Length) segStartIndex = _gbCache.Length - 1;
            return _gbCache[segStartIndex];
        }
    }
}
