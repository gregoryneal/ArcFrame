using ArcFrame.Core.Geometry;
using ArcFrame.Core.Geometry.Splines;
using ArcFrame.Core.Params;
using ArcFrame.Core.Params.CurvatureLaws;
using ArcFrame.Core.Results;
using System.ComponentModel.DataAnnotations;

namespace ArcFrame.Core.Math.Geometry.Splines
{
    public enum CubicBasis
    {
        BezierBernstein,
        CatmullRomUniform
    }
    /// <summary>
    /// Create a Bezier curve out of control points.
    /// This curve only supports a frenet frame at the moment.
    /// </summary>
    /*
    public class BezierCurve : IArcLengthCurve
    {
        public readonly double[][] ControlPoints;
        public readonly double[,] BasisMatrix;
        //TODO: Maybe cache the G matrix
        //public readonly double[,] G;
        public int Dimension => ControlPoints[0].Length;

        public double Length => _table.Length;

        private readonly ArcLengthTable _table;

        public BezierCurve(IReadOnlyList<double[]> points, double alpha = 0.5, int samples = 1024)
        {
            _table = new ArcLengthTable(Pos, samples);
            ControlPoints = new double[points.Count][];
            for (int i = 0; i < points.Count; i++)
            {
                ControlPoints[i] = points[i];
            }
        }

        /// <summary>
        /// Get a position along the curve given a parameter in [0,1]
        /// p(u) = GBu' 
        /// where G are the needed control points, B the basis matrix and u'
        /// the vector [1, u, u^2, u^3]
        /// </summary>
        /// <param name="t">Interpolation parameter</param>
        /// <returns></returns>
        public double[] Pos(double t)
        {
            (int seg, double localS)  = Locate(t);
            double s2 = localS * localS;
            double[] u = [1, localS, localS * localS, localS * s2];
            double[,] G = CreateG(seg);
            return Helpers.Multiply(Helpers.Multiply(G, BasisMatrix), u);
        }

        public Sample Evaluate(double s)
        {
            // Global t
            double t = _table.MapStoT(s);
            // Local t
            (int seg, double localT) = Locate(t);
            // Control matrix
            double[,] G = CreateG(seg);
            // Exact frenet frame
            double[] u = [1, localT, localT * localT, localT * localT * localT];
            double[] up = [0, 1, 2 * localT, 3 * localT * localT];
            double[] upp = [0, 0, 2, 6 * localT];
            double[,] GB = Helpers.Multiply(G, BasisMatrix);
            double[] P = Helpers.Multiply(GB, u);
            double[] GBup = Helpers.Multiply(GB, up);
            double[] T = Helpers.Normalize(GBup);
            double[] Tp = Helpers.Multiply(GB, upp);
            double[] B = Helpers.Normalize(Helpers.Cross3(T, Tp));
            double[] N = Helpers.Cross3(B, T);
            double[,] R = ONFrame.R0_FromTNB_Complete(T, N, B);
            double[] k = new double[System.Math.Max(0, Dimension - 1)];
            k[0] = Helpers.Len(GBup);
            k[1] = Helpers.Len(Tp);


            // Do graham schmidt orthonormalization here instead of R0_FromTNB_Complete
            // so we can use the intermediate vectors to calculate the curvatures.
            for (int i = 0; i < k.Length; i++)
            {
                // Preferably remove the static k[0] and k[1] reference and build everything dynamically
            }

            return new Sample(P, R, s, k);

            throw new NotImplementedException();
        }

        /// <summary>
        /// Fine the starting control point given an interpolant along the curve.
        /// </summary>
        /// <param name="t"></param>
        /// <returns>Start control point index k, local interpolation along spline segment u.</returns>
        (int k, double u) Locate(double t)
        {
            if (t > 1) t = 1;
            if (t < 0) t = 0;
            double seg = t * (ControlPoints.Length - 3); // segments = n-3
            int k = System.Math.Min(ControlPoints.Length - 4, (int)System.Math.Floor(seg));
            return (k, seg - k);
        }

        double[,] CreateG(int cpi)
        {
            int n = Dimension;
            double[,] G = new double[n, 4];
            for (int col = 0; col < 4; col++)
            {
                var P = ControlPoints[cpi + col];
                for (int r = 0; r < n; r++) G[r, col] = P[r];
            }
            return G;
        }
    }
    */

    /// <summary>
    /// Degree-d, N-D Bezier-like curve with pluggable basis per segment.
    /// Default basis is Bernstein for any degree. Matrix-defined bases (e.g., Catmull-Rom) can override EvaluateBasis().
    /// Provides Frenet or Bishop frames and exposes curvatures via ICurvatureLaw/ICurvatureJetProvider.
    /// </summary>
    /*
    public class BezierCurve : IArcLengthCurve
    {
        public readonly double[][] ControlPoints;
        public int Degree { get; }
        public int Dimension => ControlPoints[0].Length;

        // Framing
        private readonly FrameModel _frameModel;
        private readonly double[]? _upHint;

        // Basis configuration
        private readonly bool _useBernstein;          // default true
        protected readonly double[,]? BasisMatrix;    // optional (d+1)x(d+1) matrix for matrix-defined bases

        // Arc-length mapping
        private readonly ArcLengthTable _table;
        public double Length => _table.Length;

        // Curvature law surface for this extrinsic curve
        public ICurvatureLaw KappaLaw { get; }
        public ICurvatureJetProvider? KappaJet { get; }

        /// <summary>
        /// General constructor: Bernstein basis for any degree, Frenet frame by default.
        /// </summary>
        public BezierCurve(
            IReadOnlyList<double[]> points,
            int degree,
            FrameModel frame = FrameModel.Frenet,
            double[]? upHint = null,
            int samples = 1024)
        {
            if (points == null) throw new ArgumentNullException(nameof(points));
            if (points.Count < degree + 1) throw new ArgumentException("Not enough control points for the requested degree.");
            if (degree < 1) throw new ArgumentException("Degree must be >= 1.");

            Degree = degree;
            ControlPoints = new double[points.Count][];
            for (int i = 0; i < points.Count; i++) ControlPoints[i] = points[i];

            _frameModel = frame;
            _upHint = upHint;

            _useBernstein = true;
            BasisMatrix = null;                       // not used in Bernstein mode

            _table = new ArcLengthTable(Pos, samples); // s->t map using current Pos

            // Attach a curvature law implementation suitable for this degree.
            KappaLaw = new BezierCurvatureLaw(this);
            KappaJet = KappaLaw as ICurvatureJetProvider;
        }

        /// <summary>
        /// Alternate constructor: matrix-defined basis (e.g., Catmull–Rom for degree=3).
        /// </summary>
        public BezierCurve(
            IReadOnlyList<double[]> points,
            int degree,
            double[,] basisMatrix,
            FrameModel frame = FrameModel.Frenet,
            double[]? upHint = null,
            int samples = 1024)
        {
            if (points == null) throw new ArgumentNullException(nameof(points));
            if (points.Count < degree + 1) throw new ArgumentException("Not enough control points for the requested degree.");
            if (degree < 1) throw new ArgumentException("Degree must be >= 1.");
            if (basisMatrix.GetLength(0) != degree + 1 || basisMatrix.GetLength(1) != degree + 1)
                throw new ArgumentException("Basis matrix must be (degree+1)×(degree+1).");

            Degree = degree;
            ControlPoints = new double[points.Count][];
            for (int i = 0; i < points.Count; i++) ControlPoints[i] = points[i];

            _frameModel = frame;
            _upHint = upHint;

            _useBernstein = false;
            BasisMatrix = basisMatrix;

            _table = new ArcLengthTable(Pos, samples);

            KappaLaw = new BezierCurvatureLaw(this);
            KappaJet = KappaLaw as ICurvatureJetProvider;
        }

        // ---- IArcLengthCurve ----
        public Sample Evaluate(double s)
        {
            s = System.Math.Clamp(s, 0.0, Length);

            // 1) s -> (global) t in [0,1]
            double t = _table.MapStoT(s);

            // 2) Locate segment + local parameter u
            (int seg, double u) = Locate(t);

            // 3) Compute position and derivatives wrt local u (r', r'', r''')
            //    We compute up to 3rd derivative which supports κ1 (curvature) and κ2 (torsion).
            //    (Higher-degree curves can be extended to higher jets if/when needed.)
            ComputeSegmentJet(seg, u, 3, out var P, out var jets);

            // 4) Build frame (Frenet or Bishop)
            var T = Helpers.Normalize(jets[0]);
            double[,] R;
            if (_frameModel == FrameModel.Frenet)
            {
                // N1 ~ normalized projection of r'' orthogonal to T (falls back gracefully)
                var Nproj = Helpers.Reject(jets[1], T);
                double[]? N1 = Helpers.Len2(Nproj) > 1e-12 ? Helpers.Normalize(Nproj) : null;
                R = (N1 != null) ? ONFrame.R0_FromTN_Complete(T, N1) : ONFrame.R0_FromT_Complete(T, _upHint);
            }
            else
            {
                // Bishop (rotation-minimizing) – seeded by optional upHint; one-shot construction.
                R = ONFrame.R0_FromT_Complete(T, _upHint);
            }

            // 5) Curvatures via law (value-level)
            var k = KappaLaw.Eval(s);
            if (k == null || k.Length != System.Math.Max(0, Dimension - 1))
            {
                var kk = new double[System.Math.Max(0, Dimension - 1)];
                if (k != null) Array.Copy(k, kk, System.Math.Min(k.Length, kk.Length));
                k = kk;
            }

            return new Sample(P, R, s, k);
        }

        public double[] Position(double s) => Evaluate(s).P;

        // ---- Geometry core ----

        /// <summary>
        /// Get position P(t) in [0,1] using the configured basis.
        /// </summary>
        public double[] Pos(double t)
        {
            (int seg, double u) = Locate(t);
            ComputeSegmentJet(seg, u, 0, out var P, out _);
            return P;
        }

        /// <summary>
        /// Segment lookup: sliding window of (Degree+1) control points.
        /// </summary>
        public (int k, double u) Locate(double t)
        {
            if (t < 0) t = 0;
            if (t > 1) t = 1;

            int segments = ControlPoints.Length - Degree; // contiguous windows of size Degree+1
            if (segments <= 0) return (0, t);             // degenerate (single segment)

            double seg = t * segments;                    // [0,segments]
            int k = System.Math.Min(segments - 1, (int)System.Math.Floor(seg));
            return (k, seg - k);                          // local u in [0,1]
        }

        /// <summary>
        /// Compute P and derivatives up to 'maxOrder' for segment 'seg' at local parameter 'u'.
        /// Uses Bernstein for any degree by default; matrix-defined bases can override EvaluateBasis().
        /// </summary>
        public virtual void ComputeSegmentJet(int seg, double u, int maxOrder, out double[] P, out double[][] jets)
        {
            int n = Dimension;
            int d = Degree;

            // Collect degree+1 control points for this segment
            var Gcols = new double[d + 1][];
            for (int i = 0; i <= d; i++) Gcols[i] = ControlPoints[seg + i];

            // Get scalar basis arrays of length (d+1) for value and derivatives
            // 0 index is the function itself not the derivative
            double[][] bases = EvaluateBasis(d, u, maxOrder);
            Console.WriteLine("=====================");
            for (int i = 0; i < bases.Length; i++) Helpers.PrintVector(bases[i]);
            Console.WriteLine("=====================");
            P = new double[n];
            jets = new double[n-1][];
            for (int i = 0; i < n-1; i++)
            {
                jets[i] = new double[n];
            }

            // Accumulate
            for (int i = 0; i <= d; i++)
            {
                var Pi = Gcols[i];
                double[] b = new double[n];
                for (int j = 0; j < bases.Length; j++)
                {
                    //Console.WriteLine($"i: {i} | j: {j} | b len: {n} | bases[{j}] len: {bases[j].Length}");

                    b[j] = bases[j][i];
                }

                for (int a = 0; a < n; a++)
                {
                    double p = Pi[a];
                    P[a] += b[0] * p;
                    for (int j = 0; j < n-1; j++)
                    {
                        jets[j][a] += b[j+1] * p;
                    }
                }
            }
        }

        /// <summary>
        /// Return [B, dB, d2B, ... dNB]
        /// </summary>
        /// <param name="d"></param>
        /// <param name="u"></param>
        /// <param name="maxOrder"></param>
        /// <returns></returns>
        protected virtual double[][] EvaluateBasis(int d, double u, int maxOrder)
        {
            if (_useBernstein)
            {
                return BernsteinBasisAndDerivatives(d, u, maxOrder);
            }

            // p(u) = G * M * u 
            // G is the control point matrix
            // M the basis matrix
            // u the power monomial vector
            double[][] derivatives = new double[maxOrder + 1][];
            double[,] M = BasisMatrix!;
            double[] U = PowerMonomials(d, u);
            derivatives[0] = Helpers.Multiply(M, U);
            if (maxOrder > 0)
            {
                for (int i = 1; i > maxOrder; i++)
                {
                    //Console.WriteLine($"M size: [{M.GetLength(0)},{M.GetLength(1)}]");
                    //Console.WriteLine($"U size: [{U.Length}]");
                    //Console.WriteLine($"MaxOrder: {maxOrder} | i: {i}");
                    U = PowerMonomialsDerivN(d, i, u);
                    derivatives[i] = Helpers.Multiply(M, U);
                }
            }

            return derivatives;
        }

        // ---------- Basis helpers ----------

        /// <summary>
        /// Compute the Bernstein basis of degree d at u and its derivatives up to maxOrder.
        /// derivatives[0][i]   = B_i^d(u)
        /// derivatives[k][i]   = d^k/du^k B_i^d(u), for k = 1..maxOrder
        /// All rows have length (d+1). For k > d the row is all zeros.
        /// </summary>
        protected virtual double[][] BernsteinBasisAndDerivatives(int d, double u, int maxOrder)
        {
            if (d < 0) throw new ArgumentOutOfRangeException(nameof(d));
            if (maxOrder < 0) throw new ArgumentOutOfRangeException(nameof(maxOrder));

            int n = d;
            int orders = maxOrder + 1;

            // Allocate result: each derivative row has length (d+1)
            var derivatives = new double[orders][];
            for (int k = 0; k < orders; k++)
                derivatives[k] = new double[n + 1];

            // Early exit for n == 0
            if (n == 0)
            {
                derivatives[0][0] = 1.0;
                // higher derivatives remain zero
                return derivatives;
            }

            // Effective maximum derivative order that can be non-zero
            int effMax = System.Math.Min(maxOrder, n);

            // Precompute Bernstein rows Bdeg[k] = { B_0^{n-k}(u), ..., B_{n-k}^{n-k}(u) }
            // Start with degree n, then reduce degree step-by-step using exact degree-reduction.
            var Bdeg = new double[effMax + 1][];
            Bdeg[0] = ComputeBernsteinRow(n, u); // degree n
            for (int k = 1; k <= effMax; k++)
            {
                int mPrev = n - (k - 1);      // previous degree
                int m = mPrev - 1;            // new degree (n - k)
                var prev = Bdeg[k - 1];       // length mPrev + 1
                var curr = new double[m + 1]; // length m + 1

                // Exact relation: B_i^{m}(u) = ((m+1 - i)/(m+1)) B_i^{m+1}(u) + ((i+1)/(m+1)) B_{i+1}^{m+1}(u)
                double invDen = 1.0 / mPrev; // mPrev == (m + 1)
                for (int i = 0; i <= m; i++)
                {
                    curr[i] = ((mPrev - i) * invDen) * prev[i] + ((i + 1) * invDen) * prev[i + 1];
                }
                Bdeg[k] = curr;
            }

            // Fill the 0th derivative: just the degree-n Bernstein basis
            Array.Copy(Bdeg[0], derivatives[0], n + 1);

            // Precompute binomial coefficients up to effMax (tiny triangle)
            var binom = new double[effMax + 1][];
            for (int k = 0; k <= effMax; k++)
            {
                binom[k] = new double[k + 1];
                binom[k][0] = 1.0;
                binom[k][k] = 1.0;
                for (int j = 1; j < k; j++)
                    binom[k][j] = binom[k - 1][j - 1] + binom[k - 1][j];
            }

            // Falling factorials n^{\underline{k}} = n * (n-1) * ... * (n-k+1)
            var fall = new double[effMax + 1];
            fall[0] = 1.0;
            for (int k = 1; k <= effMax; k++)
                fall[k] = fall[k - 1] * (n - (k - 1));

            // Higher derivatives via the closed form:
            // d^k/du^k B_i^n(u) = n^{\underline{k}} * sum_{j=0..k} (-1)^{k-j} C(k,j) * B_{i-j}^{n-k}(u)
            for (int k = 1; k <= maxOrder; k++)
            {
                if (k > n)
                {
                    // row remains zeros
                    continue;
                }

                var rowNk = Bdeg[k];           // degree (n-k), length (n-k+1)
                int maxT = n - k;              // max valid index into rowNk

                for (int i = 0; i <= n; i++)
                {
                    // j must satisfy t = i - j in [0, n-k] and j in [0, k]
                    int minJ = System.Math.Max(0, i - maxT);
                    int maxJ = System.Math.Min(k, i);

                    double sum = 0.0;
                    for (int j = minJ; j <= maxJ; j++)
                    {
                        int t = i - j; // index into rowNk
                        double sign = (((k - j) & 1) == 1) ? -1.0 : 1.0; // (-1)^(k-j)
                        sum += sign * binom[k][j] * rowNk[t];
                    }

                    derivatives[k][i] = fall[k] * sum;
                }
            }

            return derivatives;

            // Compute {B_0^n(u), ..., B_n^n(u)} stably (O(n^2)) without explicit binomials.
            static double[] ComputeBernsteinRow(int nLocal, double uLocal)
            {
                var b = new double[nLocal + 1];
                b[0] = 1.0;
                double oneMinusU = 1.0 - uLocal;

                for (int j = 1; j <= nLocal; j++)
                {
                    double saved = 0.0;
                    for (int r = 0; r < j; r++)
                    {
                        double temp = b[r];
                        b[r] = oneMinusU * temp + saved;
                        saved = uLocal * temp;
                    }
                    b[j] = saved;
                }
                return b;
            }
        }


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

        // ---------- Curvature helpers (parameterization-invariant κ1, κ2) ----------
        /// <summary>
        /// Compute generalized curvatures [k1, k2, ...] from a curve jet at parameter u.
        /// jet[0] is optional (position); we only use jet[1]=r'(u), jet[2]=r''(u), ...
        /// Works in R^N and returns up to min(N-1, jet.Length-2) curvatures.
        /// - k1 is normal curvature (usual curvature magnitude).
        /// - k2 is torsion (signed in 3D; magnitude otherwise).
        /// - k3, k4, ... are higher Cartan curvatures (non-negative).
        /// </summary>
        /// <param name="jet">
        /// double[][] where jet[k] is the k-th derivative vector r^{(k)}(u).
        /// jet[1] must exist and all jet[k] have the same length N.
        /// </param>
        /// <returns>double[] curvatures of length R = min(N-1, availableDerivatives-1).</returns>
        internal double[] GeneralizedCurvaturesFromJet(double[][] jet)
        {
            if (jet == null || jet.Length < 2)
                throw new ArgumentException("jet must contain at least r'(u).", nameof(jet));

            // Identify derivative vectors: r'(u), r''(u), ...
            int firstDerIdx = 1;
            int numDerivs = jet.Length - firstDerIdx; // count of derivatives available
            if (numDerivs < 1)
                throw new ArgumentException("jet must contain at least r'(u).", nameof(jet));

            int N = jet[firstDerIdx].Length;
            for (int k = firstDerIdx; k < jet.Length; k++)
                if (jet[k] == null || jet[k].Length != N)
                    throw new ArgumentException("All derivative vectors must be non-null and have equal length.", nameof(jet));

            // We can only form Δ_k up to k = min(N, numDerivs).
            int K = System.Math.Min(N, numDerivs);
            if (K < 1)
                return Array.Empty<double>();

            // Build Δ_k using modified Gram-Schmidt on [r'(u), r''(u), ..., r^{(K)}(u)]
            // Property: if w_i are the orthogonalized vectors (not normalized),
            //           then Δ_k = Π_{i=1..k} ||w_i||^2.
            var q = new List<double[]>(capacity: K);     // orthonormal directions computed on the fly
            var wNorms = new double[K + 1];              // 1-based for readability; wNorms[i] = ||w_i||
            var Delta = new double[K + 1];               // Δ_0..Δ_K; Δ_0 := 1
            Delta[0] = 1.0;

            const double EPS = 1e-18;

            for (int i = 1; i <= K; i++)
            {
                double[] vi = (double[])jet[firstDerIdx + (i - 1)]!.Clone(); // r^{(i)}(u)

                // Subtract projections onto existing q_j (orthonormal)
                for (int j = 0; j < q.Count; j++)
                {
                    double dot = Dot(vi, q[j]);
                    Axpy(vi, q[j], -dot); // vi -= dot * q[j]
                }

                double norm = Norm(vi);
                // If degeneracy occurs, all higher Δ_k are zero too.
                if (norm <= EPS)
                {
                    for (int t = i; t <= K; t++)
                    {
                        wNorms[t] = 0.0;
                        Delta[t] = 0.0;
                    }
                    break;
                }

                wNorms[i] = norm;
                // Normalize to get q_i
                double inv = 1.0 / norm;
                for (int t = 0; t < N; t++) vi[t] *= inv;
                q.Add(vi);

                // Cumulative Δ_k
                Delta[i] = Delta[i - 1] * (norm * norm);
            }

            // How many curvatures can we return?
            // Need Δ_{r+1} to form k_r, so r_max = min(N-1, numDerivs-1).
            int rMax = System.Math.Min(N - 1, numDerivs - 1);
            if (rMax <= 0)
                return Array.Empty<double>();

            var kappa = new double[rMax];

            // Speed^2 = Δ_1
            double D1 = (K >= 1) ? Delta[1] : 0.0;
            if (D1 <= EPS)
                return kappa; // zero/degenerate velocity -> all curvatures treated as 0

            // k1 = sqrt(Δ2) / (Δ1)^(3/2)
            if (rMax >= 1 && K >= 2 && Delta[2] > EPS)
                kappa[0] = System.Math.Sqrt(Delta[2]) / System.Math.Pow(D1, 1.5);
            else if (rMax >= 1)
                kappa[0] = 0.0;

            // k_r = sqrt(Δ_{r+1} Δ_{r-1}) / (Δ_r * sqrt(Δ_1)), for r >= 2
            for (int r = 2; r <= rMax; r++)
            {
                if (r <= K - 1 && Delta[r - 1] > EPS && Delta[r] > EPS && Delta[r + 1] > EPS)
                {
                    kappa[r - 1] = System.Math.Sqrt(Delta[r + 1] * Delta[r - 1]) / (Delta[r] * System.Math.Sqrt(D1));
                }
                else
                {
                    kappa[r - 1] = 0.0;
                }
            }

            // Optional: sign the torsion in 3D (r=2) using the scalar triple product.
            if (N == 3 && rMax >= 2 && K >= 3 && Delta[2] > EPS)
            {
                double triple = ScalarTriple(jet[firstDerIdx], jet[firstDerIdx + 1], jet[firstDerIdx + 2]);
                double sign = (triple > 0) ? 1.0 : (triple < 0 ? -1.0 : 0.0);
                double torsionMag = System.Math.Sqrt(Delta[3]) / Delta[2]; // |det| / ||r'×r''||^2
                kappa[1] = sign * torsionMag;
            }

            return kappa;

            // ---- helpers ----
            static double Dot(double[] a, double[] b)
            {
                double s = 0.0;
                for (int i = 0; i < a.Length; i++) s += a[i] * b[i];
                return s;
            }

            static void Axpy(double[] y, double[] x, double alpha)
            {
                for (int i = 0; i < y.Length; i++) y[i] += alpha * x[i];
            }

            static double Norm(double[] a) => System.Math.Sqrt(Dot(a, a));

            static double ScalarTriple(double[] a, double[] b, double[] c)
            {
                // det([a b c]) = a · (b × c)
                double cx = b[1] * c[2] - b[2] * c[1];
                double cy = b[2] * c[0] - b[0] * c[2];
                double cz = b[0] * c[1] - b[1] * c[0];
                return a[0] * cx + a[1] * cy + a[2] * cz;
            }
        }


        internal void ComputeKappa12_Raw(double[] r1, double[] r2, double[] r3, out double k1, out double k2)
        {
            const double EPS = 1e-12;
            double r1r1 = Helpers.Dot(r1, r1);
            double r2r2 = Helpers.Dot(r2, r2);
            double r1r2 = Helpers.Dot(r1, r2);
            double cross2 = System.Math.Max(0.0, r1r1 * r2r2 - r1r2 * r1r2);   // ||r' ∧ r''||^2
            double denomK1 = r1r1 * System.Math.Sqrt(System.Math.Max(r1r1, EPS));     // ||r'||^3
            k1 = (denomK1 > EPS) ? System.Math.Sqrt(cross2) / denomK1 : 0.0;

            double r1r3 = Helpers.Dot(r1, r3);
            double r2r3 = Helpers.Dot(r2, r3);
            double[,] G3 = { { r1r1, r1r2, r1r3 }, { r1r2, r2r2, r2r3 }, { r1r3, r2r3, Helpers.Dot(r3, r3) } };
            double detG3 = Helpers.Det3(G3);
            k2 = (cross2 > EPS) ? System.Math.Sqrt(System.Math.Max(0.0, detG3)) / cross2 : 0.0;
        }
    }
    */
}
