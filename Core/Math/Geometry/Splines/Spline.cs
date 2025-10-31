using ArcFrame.Core.Geometry;
using ArcFrame.Core.Geometry.Splines;
using ArcFrame.Core.Results;
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
        /// The spline control points.
        /// </summary>
        protected readonly double[][] ControlPoints;
        /// <summary>
        /// The arc length table which caches the curve on creation.
        /// </summary>
        protected readonly ArcLengthTable _table;
        /// <summary>
        /// The frame model, frenet or RMF.
        /// </summary>
        protected readonly FrameModel Frame;
        /// <summary>
        /// The start ON frame.
        /// </summary>
        protected readonly double[,] R0;
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
        /// Generate an arc length table
        /// </summary>
        public Spline(double[][] controlPoints, double[,] basisMatrix, FrameModel frame = FrameModel.Frenet)
        {
            BasisMatrix = basisMatrix;
            Degree = BasisMatrix.GetLength(0) - 1;
            // Handle high degree splines in low dimensional spaces
            /*
            if (Degree > controlPoints[0].Length)
            {
                // Promote the control points into the new dimension, pad with 0s basically
                double[][] promotedControlPoints = new double[controlPoints.Length][];
                for (int i = 0; i < controlPoints.Length; i++)
                {
                    double[] v = new double[Degree];
                    for (int j = 0; j < controlPoints[i].Length; j++)
                    {
                        v[j] = controlPoints[i][j];
                    }
                    promotedControlPoints[i] = v;
                }

                controlPoints = promotedControlPoints;
            }*/

            ControlPoints = controlPoints;
            Dimension = ControlPoints[0].Length;
            Frame = frame;

            // Scale arc length table number of samples by polyline length of control points
            double len = 0;
            for (int i = 0; i + 1 < controlPoints.Length; i++)
            {
                len += Helpers.Len(Helpers.Subtract(controlPoints[i + 1], controlPoints[i]));
            }
            // N samples per unit arc length
            CacheSamples = System.Math.Max(64, (int)System.Math.Round(len) * 256 + 1);
            _table = new ArcLengthTable(Position, CacheSamples);
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

        /// <summary>
        /// Evaluating the Spline curve is done in a few steps.
        /// Goal: p(t) = GBt'
        /// p = position
        /// t = parameter in [0, 1]
        /// t' = power polynomial [1, t, t^2, t^3, ..., t^degree]^T
        /// G = degree + 1 control point column vectors whose first index is floor(t*(controlpoints.length - degree))
        /// B = basis matrix, controls the type of spline being generated.
        /// 
        /// Step 1: Get parameter t in [0, 1] from s in [0, Length] (s / Length)\
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
            double t = _table.MapStoT(s);
            (int seg, double u) = Locate(t);
            // Step 2.
            double[,] G = CreateG(seg);
            double[,] GB = Helpers.Multiply(G, BasisMatrix);
            // Step 3.
            double[] tp = PowerMonomials(Degree, u); //len degree + 1
            // Step 4.
            double[] P = Helpers.Multiply(GB, tp);
            // Step 5.
            double[,] R;

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
            double[] k = new double[System.Math.Max(1, Dimension - 1)]; // Curvatures are filled in different depending on the type of frame used, Frenet or RMF.

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
        /// Find a position along the spline given its parameter t in [0, 1]
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        public double[] Position(double t)
        {
            t = System.Math.Clamp(t, 0, 1);
            (int seg, double localS) = Locate(t);
            // Step 2.
            double[,] G = CreateG(seg);
            // Step 3.
            double[] tp = PowerMonomials(Degree, localS); //len degree + 1
            // Step 4.
            return Helpers.Multiply(Helpers.Multiply(G, BasisMatrix), tp);
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
        /// Fine the starting control point given an interpolant along the curve.
        /// </summary>
        /// <param name="t"></param>
        /// <returns>Start control point index k, local interpolation along spline segment u.</returns>
        protected virtual (int k, double u) Locate(double t)
        {
            if (t > 1) t = 1;
            if (t < 0) t = 0;
            double seg = t * (ControlPoints.Length - Degree);
            int k = System.Math.Min(ControlPoints.Length - (Degree + 1), (int)System.Math.Floor(seg));
            return (k, seg - k);
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

            double[,] Rlin = new double[Dimension, Dimension];
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
            double tt = _table.MapStoT(ss);
            (int sseg, double suu) = Locate(tt);
            double[,] Gs = CreateG(sseg);
            double[,] GBs = Helpers.Multiply(Gs, BasisMatrix);
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
    }
}
