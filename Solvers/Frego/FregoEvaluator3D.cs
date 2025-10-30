using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using ArcFrame.Core.Params.CurvatureLaws;
using ArcFrame.Core.Results;
using ArcFrame.Solvers.BertolazziFrego;
using System;

namespace ArcFrame.Solvers.Frego
{
    /// <summary>
    /// 3D evaluator from Marco Frego. I have not tried optimizing this
    /// but the unoptimized version is almost the same speed as the 
    /// built in IFrameStepper used on a normal IntrinsicCurve.
    /// In the canonical case the IFrameStepper is 2-3x faster.
    /// I probably need to dig in and optimize more.
    /// Also this clothoid evaluator does not support the frenet frame.
    /// If you supply that as input it will fall back to the default
    /// IFrameStepper.
    /// </summary>
    public class FregoEvaluator3D : IEvaluator
    {
        public string Key => "Clothoid3DEvaluator";

        public string Author => "Marco Frego";

        public string Reference => "https://doi.org/10.1016/j.amc.2021.126907";

        public int TargetDimension => 3;

        /// <summary>
        /// Frego luckily for me does a bunch of linear algebra, which I have already created an extensive
        /// library for. His method involves checking if the commutator condition holds: m = k't0 - t'k0 = 0
        /// there are closed form solutions. Otherwise I have to use the magnus expansion which involves
        /// creating matrices much larger than I would normally need to given the dimension of the space
        /// we are in. But Frego also has some methods to simplify these into just a few matrix operations.
        /// Something to be aware of is that in my skew symmetric matrices for the Frenet and bishop frames,
        /// I accidentally set everything up where they are transposed from how they are in common literature.
        /// I do not want to go through the trouble to fix everything just for that.
        /// This means I just need to take care when following matrix operations from the literature.
        /// 
        /// General process:
        /// If the clothoid does not satisfy the commutator condition, we must integrate the moving frame
        /// with a magnus expansion.
        /// 
        /// R'(s) = A(s)R(s) where R is the moving frame and A is the skew symmetric matrix with k, t on 
        /// sub and super diagonals. Some mathematician Magnus proposed a solution to this equation in the
        /// form:
        /// R(s) = exp(O(s))R(0) where
        /// O(s) is a sum from i = 1 to infinity of O_{i}(s)
        /// which is called the Magnus expansion.
        /// O_1 = integral from 0 to s of A(u)du
        /// O_2 = 0.5 * integral from 0 to s of [O_1(u), A(u)]du where
        /// [O_1(s), A(s)] =      ms^2  | 0  0 -1 |
        ///                       ____  | 0  0  0 |
        ///                         2   | 1  0  0 |
        /// m is defined above.
        /// We can derive 4th order magnus expansion from the first two terms.
        /// 
        /// Once we calculate O(s), we exponentiate it and then multiply with R(0) to get the frame at s.
        /// 
        /// Particulars:
        /// Solve above problem with X'(s) = B(s)X(s)
        ///     where X = | Tx  Nx  Bx  x(s) 0   0 |
        ///               | Ty  Ny  By   0  y(s) 0 |
        ///               | Tz  Nz  Bz   0   0 z(s)|
        ///     Take A(s) and embed it in a 4x4 matrix (final column all 0, final row = [1, 0, 0, 0]).
        ///     kronecker multiply that matrix by the 3x3 identity matrix to obtain a 12x12 matrix.
        ///     It's really freakin big. 
        ///     
        /// For our integration interval [0, L], we need to split it into n subintervals
        /// of length h = L/n and si = ih for i in [0, 1, 2, ... n]
        /// We calculate X_{i+1} = exp(E_i)X_i where E_i = O_1 + O_2.
        /// After all of the math simplification which I can not normally understand, we arrive
        /// at all of the little formulas we gotta consruct to build the final solution.
        /// 
        /// l1 = k'h^2(i + 0.5) + k0h
        /// l2 = t'h^2(i + 0.5) + t0h
        /// l3 = mh^3/12
        /// l4 = -k'h^3/12
        /// lambda = sqrt(l1^2 + l2^2 + l3^2)
        /// sl = sin(lambda)
        /// cl = cos(lambda)
        /// cl' = 1 - cos(lambda)
        /// q1 = [hl1^2 + l3(l2l4 + hl3)] / lambda
        /// q2 = [l4l1^2 + l2(l2l4 + hl3)] / lambda
        /// q3 = [-l1(hl2 - l3l4)] / lambda
        /// c41 = q1sl - l1l4cl' + l2(hl2 - l3l4)
        /// c42 = q2sl + hl1cl' - (hl2l3 - l4l3^2)
        /// c43 = q3sl + (l2l4 + hl3)cl' + (hl1l2 - l1l3l4)
        /// 
        /// now exp(E_i) = (1/lambda^2)C
        ///     | (l1^2 + l3^2)cl + l2^2  (lambda)l1sl - l2l3cl'  (lambda)l3sl + l1l2cl'  0 |
        /// C = | -(lambda)l1sl - l2l3cl' (l1^2 + l2^2)cl + l3^2  (lambda)l2sl - l1l3cl'  0 |
        ///     | -(lambda)l3sl - l1l2cl' -(lambda)l2sl - l1l3cl' (l2^2 + l3^2)cl + l1^2  0 |
        ///     |        c41                     c42                      c43       lambda^2|
        ///     
        /// do CX_i to find X_{i+1}
        /// 
        /// If the clothoid does satisfy the commutator condition we can use a different method:
        /// 
        /// Let G be a general clothoid that satisfies the commutator condition m = k't0 - t'k0 = 0
        /// Let C be the canonical clothoid (starts at the origin, k0 = t0 = 0, curve given by Fresnel integrals)
        /// Let N be a clothoid calculated with the heavy numerical recipe above: X_{i+1} = CX_i
        /// 
        /// [pCs, tCs, nCs, bCs] = C(s; p0, t0, n0, b0, k', t')
        /// [pGs, tGs, nGs, bGs] = G(s; p0, t0, n0, b0, k', t', k0, t0)
        /// [pNs, tNs, nNs, bNs] = N(s; p0, t0, n0, b0, k', t', k0, t0)
        /// 
        /// we can write 
        /// 
        /// G(s; pG0, tG0, nG0, bG0, k', t', k0, t0) = C(s - s'; pG0 - pCT, tC0, nC0, bC0, k', t')
        /// 
        /// where s' = -t0/t' = -k0/k'
        ///
        /// [pC0, tC0, nC0, bC0] = C(-s'; origin, tG0, nG0, bG0, -k', -t')
        /// [pCT, tCT, nCT, bCT] = C(-s'; origin, tC0, nC0, bC0,  k',  t') => only need pCT
        /// 
        /// The paper gives the canonical clothoid formula. It involves the Fresnel integrals.
        /// </summary>
        /// <param name="p"></param>
        /// <param name="s"></param>
        /// <returns></returns>
        /// <exception cref="System.NotImplementedException"></exception>
        public Sample Evaluate(CurveSpec p, double s)
        {
            // TODO: Add Bishop frame support, magnus expansion step can be done easily as we iterate, need to figure out something for the canonical method.
            if (p.N != 3) throw new ArgumentOutOfRangeException("Curve dimension should be 3 for 3D clothoid.");

            double[][] jet = p.Kappa.EvalJet(0); //defaults to first derivative only (always includes 0th order): { [k0, t0], [k', τ'] }
            //Helpers.PrintVector(jet[0]);
            //Helpers.PrintVector(jet[1]);
            double k0 = jet[0][0];
            double dk = jet[1][0];
            double t0 = jet[0][1];
            double dt = jet[1][1];

            //Console.WriteLine($"Input params: k0: {k0} | dk: {dk} | t0: {t0} | dt: {dt} | ");


            double[] P = (double[])p.P0.Clone();
            double[,] R = (double[,])p.R0.Clone();

            double[] origin = new double[] { 0, 0, 0 };

            // if the commutator condition m = 0 holds
            // we can do some fancy shit
            double m = (dk * t0) - (dt * k0);
            //Console.WriteLine($"m: {m}");
            double eps = 1E-12;
            if (Math.Abs(m) <= eps)
            {
                // pure canonical case
                if (Math.Abs(k0) <= eps && Math.Abs(t0) <= eps)
                {
                    var cc = CanonicalClothoid(s, R, dk, dt);
                    return new Sample(Helpers.Add(P, cc.P), cc.R, cc.s, p.Kappa.Eval(s));
                }
                // affine case
                double sbar = Math.Abs(dt) > Math.Abs(dk) ? -t0 / dt : -k0 / dk;
                double[,] C0 = CanonicalClothoid(-sbar, R, -dk, -dt).R;
                double[] CT = CanonicalClothoid(-sbar, C0, dk, dt).P;
                Sample S = CanonicalClothoid(s - sbar, C0, dk, dt);
                // final position
                double[] Pf = Helpers.Add(Helpers.Subtract(P, CT), S.P);
                return new Sample(Pf, S.R, s, p.Kappa.Eval(s));
            }

            // Integrate up to s using the truncated magnus expansion
            // state vector Y = [tx, ty, tz, nx, .., bz, x, y, z]
            // updated each step
            // Integration steps n >= 256 keep error: |t•n|+|t•b|+|n•b| around 1E-15 up to at least s = 20.
            double L = Math.Abs(s);
            int stepMultiplier = (int)Math.Ceiling(L / 20); // how many times we have a 20 arc length run in the evaluation rounded up
            int numSteps = Math.Max(1, 512 * stepMultiplier); // should test speed with 256, 512, 1024, etc up to 2^14

            // magnus guard, update per step
            // we choose h as the value that keeps integral from 0 to s of Q(e)de
            // which is a cubic in s that satisfies: as^3 + bs^2 + cs < pi
            // where
            // - a = (dk^2 + dt^2) / 3
            // - b = (kidk + tidt)
            // - c = ki^2 + ti^2
            double h;
            double si = 0;


            double[] Y = new double[] { R[0, 0], R[1, 0], R[2, 0], R[0, 1], R[1, 1], R[2, 1], R[0, 2], R[1, 2], R[2, 2], P[0], P[1], P[2] };
            //double h = L / numSteps; 
            double[,] C = new double[4, 4];
            double[,] I3 = RigidTransform.Identity(3).R;

            int i = 0; // counter for re-ON
            // Use a variable step size, different from the paper
            // I just learned that this is not good because she
            // converges for a rather large arc length, i need
            // to ensure we get a good sample depth.
            while (si < L)
            {
                h = Math.Min(1E-2, MagnusMaxStepAtS(si, k0, t0, dk, dt));

                double smid = si + (h / 2);

                //Console.WriteLine($"New h: {h}");
                double l1 = dk * h * (si + h / 2) + k0 * h;
                double l2 = dt * h * (si + h / 2) + t0 * h;
                double l3 = m * h * h * h / 12;
                double l4 = -dk * h * h * h / 12;
                double lambda = Math.Sqrt(l1 * l1 + l2 * l2 + l3 * l3);

                SafeTrig(lambda, out double sl, out double cl, out double clb); 
                double invLam = (Math.Abs(lambda) < 1e-12) ? 0.0 : 1.0 / lambda;
                double invLam2 = (Math.Abs(lambda) < 1e-12) ? 0.0 : 1.0 / (lambda * lambda);

                double q1 = (h * l1 * l1 + l3 * (l2 * l4 + h * l3)) * invLam;
                double q2 = (l4 * l1 * l1 + l2 * (l2 * l4 + h * l3)) * invLam;
                double q3 = (-l1 * (h * l2 - l3 * l4)) * invLam;
                double c41 = q1 * sl - l1 * l4 * clb + l2 * (h * l2 - l3 * l4);
                double c42 = q2 * sl + h * l1 * clb - (h * l2 * l3 - l4 * l3 * l3);
                double c43 = q3 * sl + (l2 * l4 + h * l3) * clb + (h * l1 * l2 - l1 * l3 * l4);

                C[0, 0] = (l1 * l1 + l3 * l3) * cl + l2 * l2;
                C[1, 0] = -lambda * l1 * sl - l2 * l3 * clb;
                C[2, 0] = -lambda * l3 * sl + l1 * l2 * clb;
                C[3, 0] = c41;
                C[0, 1] = lambda * l1 * sl - l2 * l3 * clb;
                C[1, 1] = (l1 * l1 + l2 * l2) * cl + l3 * l3;
                C[2, 1] = -lambda * l2 * sl - l1 * l3 * clb;
                C[3, 1] = c42;
                C[0, 2] = lambda * l3 * sl + l1 * l2 * clb;
                C[1, 2] = lambda * l2 * sl - l1 * l3 * clb;
                C[2, 2] = (l2 * l2 + l3 * l3) * cl + l1 * l1;
                C[3, 2] = c43;
                C[3, 3] = lambda * lambda;

                // maybe TODO: populate the 12x12 by hand?

                // exp(E_i) = 1 / lambda^2
                double[,] expE = Helpers.Multiply(invLam2, LA.KroneckerProduct(C, I3));
                Y = Helpers.Multiply(expE, Y);

                // re-orthonormalize t,n,b inside Y
                // every 32 steps
                if ((i % 31) == 0)
                {
                    double[] T = { Y[0], Y[1], Y[2] };
                    double[] N = { Y[3], Y[4], Y[5] };
                    double[] B;// = { Y[6], Y[7], Y[8] };

                    // if the normal is slightly off
                    // we find the part perpendicular
                    // to T and normalize it
                    N = Helpers.Normalize(Helpers.Reject(N, T));
                    B = Helpers.Cross3(T, N);

                    Y[0] = T[0];
                    Y[1] = T[1];
                    Y[2] = T[2];
                    Y[3] = N[0];
                    Y[4] = N[1];
                    Y[5] = N[2];
                    Y[6] = B[0];
                    Y[7] = B[1];
                    Y[8] = B[2];
                }

                i++;
                si += h;
            }

            // Extract final pose from Y
            R[0, 0] = Y[0];
            R[1, 0] = Y[1];
            R[2, 0] = Y[2];
            R[0, 1] = Y[3];
            R[1, 1] = Y[4];
            R[2, 1] = Y[5];
            R[0, 2] = Y[6];
            R[1, 2] = Y[7];
            R[2, 2] = Y[8];
            P[0] = Y[9];
            P[1] = Y[10];
            P[2] = Y[11];
            return new Sample(P, R, s, p.Kappa.Eval(s));
        }

        /// <summary>
        /// Evaluate the canonical 3d clothoid (P0 @ origin, k0 = t0 = 0)
        /// The paper only gives the parameterization in terms
        /// of the frenet serret frame, so i need to do some math
        /// when i implement this.
        /// </summary>
        /// <param name="s"></param>
        /// <param name="R0">the frame being integrated, this should be R at s0</param>
        /// <param name="dk"></param>
        /// <param name="dt"></param>
        /// <returns></returns>
        private static Sample CanonicalClothoid(double s, double[,] R0, double dk, double dt) 
        {
            double del = dk * dk + dt * dt;
            double rtdel = Math.Sqrt(del);
            double param = rtdel * s * s / 2;
            double cs = Math.Cos(param);
            double ss = Math.Sin(param);

            double[] t0 = ONFrame.GetCol(R0, 0);
            double[] n0 = ONFrame.GetCol(R0, 1);
            double[] b0 = ONFrame.GetCol(R0, 2);

            BFEvaluator.FresnelCS(rtdel / 2, out double C, out double S);

            // Frame at an arc length of s along the curve compared to R.
            double[] t = Helpers.Multiply(1 / del, Helpers.Add(Helpers.Add(Helpers.Multiply(dk * dk * cs + dt * dt, t0), Helpers.Multiply(rtdel * dk * ss, n0)), Helpers.Multiply(dk * dt * (1 - cs), b0)));
            double[] n = Helpers.Multiply(1 / rtdel, Helpers.Add(Helpers.Add(Helpers.Multiply(-dk * ss, t0), Helpers.Multiply(rtdel * cs, n0)), Helpers.Multiply(dt * ss, b0)));
            double[] b = Helpers.Multiply(1 / del, Helpers.Add(Helpers.Add(Helpers.Multiply(dk * dt * (1 - cs), t0), Helpers.Multiply(-rtdel * dt * ss, n0)), Helpers.Multiply(dt * dt * cs + dk * dk, b0)));
            double[] P = Helpers.Multiply(1 / del, Helpers.Add(Helpers.Add(Helpers.Multiply(dk * dk * C + dt * dt, t0), Helpers.Multiply(rtdel * dk * S, n0)), Helpers.Multiply(dk * dt * (s - C), b0)));

            LinearCurvatureLaw law = new LinearCurvatureLaw(new double[] { 0, 0 }, new double[] { dk, dt });
            return new Sample(P, ONFrame.R0_FromTNB(t, n, b), s, law.Eval(s));
        }

        /// <summary>
        /// Generalized Fresnel integrals for cos(a t^2), sin(a t^2)
        /// Uses the existing canonical FresnelCS (π/2 inside the trig) via scaling.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <param name="C"></param>
        /// <param name="S"></param>
        private static void FresnelCS(double a, double x, out double C, out double S)
        {
            // Handle a = 0 exactly: ∫0^x cos(0) dt = x, ∫0^x sin(0) dt = 0
            if (a == 0.0)
            {
                C = x;
                S = 0.0;
                return;
            }

            double b = Math.Abs(a);
            double signX = x < 0 ? -1.0 : 1.0;
            double signA = a < 0 ? -1.0 : 1.0;

            // Scale in/out
            double scaleIn = Math.Sqrt(2.0 * b / Math.PI);   // argument scale
            double scaleOut = Math.Sqrt(Math.PI / (2.0 * b)); // value scale

            // Call your existing canonical Fresnel on y ≥ 0
            double y = scaleIn * Math.Abs(x);
            double Cstd, Sstd;
            BFEvaluator.FresnelCS(y, out Cstd, out Sstd); 

            // Apply scales and signs:
            C = signX * (scaleOut * Cstd);            // C(a;x) is odd in x, even in a
            S = signX * signA * (scaleOut * Sstd);    // S(a;x) is odd in x and odd in a
        }

        static void SafeTrig(double lambda, out double sl, out double cl, out double clb)
        {
            if (Math.Abs(lambda) < 1e-6)
            {
                double l2 = lambda * lambda, l3 = l2 * lambda, l4 = l2 * l2, l6 = l4 * l2;
                sl = lambda - l3 / 6.0 + l6 / 120.0;
                cl = 1.0 - l2 / 2.0 + l4 / 24.0 - l6 / 720.0;
                clb = l2 / 2.0 - l4 / 24.0 + l6 / 720.0; // 1 - cos λ
            }
            else
            {
                sl = Math.Sin(lambda);
                cl = Math.Cos(lambda);
                clb = 1.0 - cl;
            }
        }

        /// <summary>
        /// Compute a safe local step h_max at s_i so that ∫_{s_i}^{s_i+h} (k^2 + t^2) ds ≤ π.
        /// Inputs: ki = κ(s_i), ti = τ(s_i), dk = κ', dt = τ'.
        /// Returns: safety * h_root (default safety=0.9).
        /// </summary>
        public static double MagnusMaxStepLocal(double ki, double ti, double dk, double dt, double safety = 0.9)
        {
            // Coefficients of I(h) = a h^3 + b h^2 + c h
            double a = (dk * dk + dt * dt) / 3.0;
            double b = ki * dk + ti * dt;
            double c = ki * ki + ti * ti;

            // We solve f(h) = a h^3 + b h^2 + c h - π = 0 for the largest positive root
            // with a safeguarded Newton method (bracket from the right).
            const double PI = Math.PI;
            const double eps = 1e-15;

            // Handle degeneracies cleanly
            a = Math.Max(a, 0.0);
            if (c < eps)
            {
                // c ≈ 0 ⇒ I(h) ≈ a h^3 + b h^2
                if (a < eps && b <= 0.0) return double.PositiveInfinity; // flat or decreasing ⇒ no finite bound
                if (a < eps)
                {                 // b h^2 = π
                    if (b <= 0.0) return double.PositiveInfinity;
                    return safety * Math.Sqrt(PI / b);
                }
                else if (b < eps)
                {          // a h^3 = π
                    return safety * Math.Cbrt(PI / a);
                }
                else
                {
                    // Use a simple bracket + Newton anyway
                }
            }

            double f(double h) => ((a * h + b) * h + c) * h - PI;
            double df(double h) => (3.0 * a * h + 2.0 * b) * h + c;

            // Initial guess: quadratic / linear fallbacks if useful
            double h0;
            if (a < 1e-12 && b > 0.0)
            {
                // Solve b h^2 + c h - π = 0 (positive root)
                double disc = c * c + 4.0 * b * PI;
                h0 = (-c + Math.Sqrt(Math.Max(disc, 0.0))) / (2.0 * b);
            }
            else if (c > eps)
            {
                h0 = PI / c;
            }
            else
            {
                h0 = 1.0;
            }
            if (!(h0 > 0.0) || double.IsInfinity(h0) || double.IsNaN(h0)) h0 = 1.0;

            // Find right bracket: f(0) = -π < 0. Increase hi until f(hi) > 0.
            double lo = 0.0, hi = h0;
            double fhi = f(hi);
            int it = 0;
            while (fhi <= 0.0 && hi < 1e16 && it++ < 60)
            {
                hi *= 2.0;
                fhi = f(hi);
            }
            if (fhi <= 0.0)
            {
                // As a last resort, pick a very large hi from h ~ cbrt(π/a)
                hi = (a > eps) ? 2.0 * Math.Cbrt(PI / a) : 1e8;
                fhi = f(hi);
                if (fhi <= 0.0) return safety * hi; // monotone but still under ⇒ accept huge (practically "no cap")
            }

            // Safeguarded Newton (fallback to bisection if Newton steps leave bracket)
            double h = Math.Min(h0, hi);
            if (!(h > 0.0) || f(h) <= 0.0) h = 0.5 * hi; // start inside (0, hi)
            for (int k = 0; k < 20; k++)
            {
                double fh = f(h);
                double dfh = df(h);
                // Converged?
                if (Math.Abs(fh) < 1e-12 * (1.0 + c * h)) break;

                // Newton trial
                double hN = h - fh / Math.Max(dfh, 1e-18);
                // Safeguard: keep in (lo, hi)
                if (!(hN > lo && hN < hi) || double.IsNaN(hN) || double.IsInfinity(hN))
                {
                    hN = 0.5 * (lo + hi); // bisection
                }

                // Update bracket
                double fN = f(hN);
                if (fN > 0.0) hi = hN; else lo = hN;
                h = hN;
            }

            return safety * h;
        }

        /// <summary>
        /// Convenience: compute local h_max at absolute arclength s_i using global (k0,t0,dk,dt).
        /// </summary>
        public static double MagnusMaxStepAtS(double s_i, double k0, double t0, double dk, double dt, double safety = 0.9)
        {
            double ki = k0 + dk * s_i;
            double ti = t0 + dt * s_i;
            return MagnusMaxStepLocal(ki, ti, dk, dt, safety);
        }

        /// <summary>
        /// Create the 12x12 matrix B evaluated at some s.
        /// B is defined as:
        /// 
        ///     | A(s)  0' |
        /// B = |  1    0  | X I_3
        /// 
        /// where
        /// - A is the 3x3 skew symmetric
        /// curvature matrix.
        /// - 0' is the 1x3 zero column vector
        /// - X is the kronecker product
        /// - I_3 is the 3x3 identity matrix
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        private double[,] ConstructB(double k, double t)
        {
            // construct B block by block
            double[,] B = new double[12, 12];
            B[3, 0] = -k;
            B[4, 1] = -k;
            B[5, 2] = -k;
            B[9, 0] = 1;
            B[10, 1] = 1;
            B[11, 2] = 1;
            B[0, 3] = k;
            B[1, 4] = k;
            B[2, 5] = k;
            B[6, 3] = -t;
            B[7, 4] = -t;
            B[8, 5] = -t;
            B[3, 6] = t;
            B[4, 7] = t;
            B[5, 8] = t;
            return B;
        }

        private double[,] ConstructE(double l1, double l2, double l3, double l4, double h)
        {
            double[,] E = new double[4, 4];
            E[1, 0] = -l1;
            E[2, 0] = -l3;
            E[3, 0] = h;
            E[0, 1] = l1;
            E[2, 1] = -l2;
            E[3, 1] = l4;
            E[0, 2] = l3;
            E[1, 2] = l2;
            return E;
        }
    }
}
