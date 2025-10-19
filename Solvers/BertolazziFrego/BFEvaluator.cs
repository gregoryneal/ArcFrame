using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using ArcFrame.Core.Params.CurvatureLaws;
using ArcFrame.Core.Results;

namespace ArcFrame.Solvers.BertolazziFrego
{
    public class BFEvaluator : IEvaluator
    {

        /// <summary>
        /// These are coefficients used in the initial guess of A in the root finding algorithm
        /// </summary>
        private static readonly double[] CF = { 2.989696028701907, 0.716228953608281, -0.458969738821509, -0.502821153340377, 0.261062141752652, -0.045854475238709 };
        static readonly double[] fn = {
            0.49999988085884732562,
            1.3511177791210715095,
            1.3175407836168659241,
            1.1861149300293854992,
            0.7709627298888346769,
            0.4173874338787963957,
            0.19044202705272903923,
            0.06655998896627697537,
            0.022789258616785717418,
            0.0040116689358507943804,
            0.0012192036851249883877
        };
        static readonly double[] fd = {
            1.0,
            2.7022305772400260215,
            4.2059268151438492767,
            4.5221882840107715516,
            3.7240352281630359588,
            2.4589286254678152943,
            1.3125491629443702962,
            0.5997685720120932908,
            0.20907680750378849485,
            0.07159621634657901433,
            0.012602969513793714191,
            0.0038302423512931250065
        };
        static readonly double[] gn = {
            0.50000014392706344801,
            0.032346434925349128728,
            0.17619325157863254363,
            0.038606273170706486252,
            0.023693692309257725361,
            0.007092018516845033662,
            0.0012492123212412087428,
            0.00044023040894778468486,
            -8.80266827476172521e-6,
            -1.4033554916580018648e-8,
            2.3509221782155474353e-10
        };
        static readonly double[] gd = {
            1.0,
            2.0646987497019598937,
            2.9109311766948031235,
            2.6561936751333032911,
            2.0195563983177268073,
            1.1167891129189363902,
            0.57267874755973172715,
            0.19408481169593070798,
            0.07634808341431248904,
            0.011573247407207865977,
            0.0044099273693067311209,
            -0.00009070958410429993314
        };
        /// <summary>
        /// Threshold for A to solve the fresnel equation with different solutions
        /// </summary>
        public static readonly double A_THRESHOLD = 1E-2;
        /// <summary>
        /// When A is small the momenta integrals are evaluated with an infinite series.
        /// Turns out this series converges very fast so we only need the first few terms.
        /// This is that number of terms.
        /// </summary>
        public static readonly int A_SMALL_SERIES_SIZE = 3;

        public string Key => throw new NotImplementedException();

        public string Author => throw new NotImplementedException();

        public string Reference => throw new NotImplementedException();

        /// <summary>
        /// The Generalized Fresnel integral in 2D as defined by Bertolazzi and Frego is given as:
        ///      1
        /// X = ∫cos(at²/2 + bt + c)
        ///      0
        ///      1
        /// Y = ∫sin(at²/2 + bt + c)
        ///      0
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        public static double[][] GeneralizedFresnelCS(int n, double a, double b, double c)
        {
            if (n < 1 || n > 3) throw new ArgumentOutOfRangeException($"Value of index: {n} | Expected value between 1 and 3 inclusive.");
            double cc = Math.Cos(c);
            double sc = Math.Sin(c);
            double[][] CS;

            if (Math.Abs(a) < A_THRESHOLD)
            {
                CS = EvalXYaSmall(n, a, b, A_SMALL_SERIES_SIZE);
            }
            else
            {
                CS = EvalXYaLarge(n, a, b);
            }

            double ci;
            double si;

            for (int i = 0; i < n; i++)
            {
                ci = CS[0][i];
                si = CS[1][i];
                CS[0][i] = ci * cc - si * sc;
                CS[1][i] = ci * sc + si * cc;
            }

            return CS;
        }

        /// <summary>
        /// Evaluate the GFI when a is small.
        /// </summary>
        /// <param name="k"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        private static double[][] EvalXYaSmall(int k, double a, double b, int p)
        {

            int nkk = k + 4 * p + 2;
            double[][] points = EvalXYaZero(b, nkk);
            double[] X0 = points[0];
            double[] Y0 = points[1];
            double[] X = new double[3];
            double[] Y = new double[3];

            for (int j = 0; j < k; j++)
            {
                X[j] = X0[j] - a / 2 * Y0[j + 2];
                Y[j] = Y0[j] + a / 2 * X0[j + 2];
            }

            double t = 1;
            double aa = -a * a / 4;
            double bf;
            int jj;
            for (int n = 1; n <= p; n++)
            {
                t *= aa / (2 * n * (2 * n - 1));
                bf = a / (4 * n + 2);

                for (int j = 0; j < k; j++)
                {
                    jj = 4 * n + j;
                    X[j] += t * (X0[jj] - bf * Y0[jj + 2]);
                    Y[j] += t * (Y0[jj] + bf * X0[jj + 2]);
                }
            }

            return [X, Y];
        }

        /// <summary>
        /// Evaluate the GFI when a is large.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        private static double[][] EvalXYaLarge(int n, double a, double b)
        {
            if (n > 3 || n < 0) throw new ArgumentOutOfRangeException();

            double s = a > 0 ? 1 : -1;
            double absa = Math.Abs(a);
            double z = Math.Sqrt(absa) / Math.Sqrt(Math.PI);
            double ell = s * b / Math.Sqrt(Math.PI * absa);
            double g = -s * b * b / (2 * absa);
            double cg = Math.Cos(g) / z;
            double sg = Math.Sin(g) / z;
            List<double[]> CS1 = FresnelCS(n, ell);
            List<double[]> CSz = FresnelCS(n, ell + z);
            double[] C1 = CS1[0];
            double[] S1 = CS1[1];
            double[] Cz = CSz[0];
            double[] Sz = CSz[1];

            double[] dC = new double[n];
            double[] dS = new double[n];

            for (int i = 0; i < n; i++)
            {
                dC[i] = Cz[i] - C1[i];
                dS[i] = Sz[i] - S1[i];
            }
            double[] X = new double[3];
            double[] Y = new double[3];

            X[0] = cg * dC[0] - s * sg * dS[0];
            Y[0] = sg * dC[0] + s * cg * dS[0];
            if (n > 1)
            {
                cg /= z;
                sg /= z;
                double DC = dC[1] - ell * dC[0];
                double DS = dS[1] - ell * dS[0];

                X[1] = cg * DC - s * sg * DS;
                Y[1] = sg * DC + s * cg * DS;

                if (n > 2)
                {
                    cg /= z;
                    sg /= z;
                    DC = dC[2] + ell * (ell * dC[0] - 2 * dC[1]);
                    DS = dS[2] + ell * (ell * dS[0] - 2 * dS[1]);
                    X[2] = cg * DC - s * sg * DS;
                    Y[2] = sg * DC + s * cg * DS;
                }
            }

            return [X, Y];
        }

        /// <summary>
        /// Evaluate the GFI when a is 0
        /// </summary>
        /// <param name="b"></param>
        /// <param name="k"></param>
        /// <param name="e">The allowable error</param>
        /// <returns></returns>
        private static double[][] EvalXYaZero(double b, int k)
        {
            double b2 = b * b;
            double sb = Math.Sin(b);
            double cb = Math.Cos(b);

            double[] X = new double[45];
            double[] Y = new double[45];

            if (Math.Abs(b) < 1E-3)
            {
                X[0] = 1 - b2 / 6 * (1 - b2 / 20 * (1 - b2 / 42));
                Y[0] = (b2 / 2) * (1 - b2 / 12 * (1 - b2 / 30));
            }
            else
            {
                X[0] = sb / b;
                Y[0] = (1 - cb) / b;
            }

            //int m = (int)Math.Min(Math.Max(1, Math.Floor(2 * b)), k);
            int m = (int)Math.Floor(2 * b);

            if (m >= k) m = k - 1;
            if (m < 1) m = 1;

            for (int i = 1; i < m; ++i)
            {
                X[i] = (sb - i * Y[i - 1]) / b;
                Y[i] = (i * X[i - 1] - cb) / b;
            }

            if (m < k)
            {
                double A = b * sb;
                double D = sb - b * cb;
                double B = b * D;
                double C = -b2 * sb;
                double rLa = RLommel(m + 0.5, 1.5, b);
                double rLd = RLommel(m + 0.5, 0.5, b);
                double rLb;
                double rLc;

                for (int i = m; i < k; ++i)
                {
                    rLb = RLommel(i + 1.5, 0.5, b);
                    rLc = RLommel(i + 1.5, 1.5, b);
                    X[i] = (i * A * rLa + B * rLb + cb) / (1 + i);
                    Y[i] = (C * rLc + sb) / (2 + i) + D * rLd;
                    rLa = rLc;
                    rLd = rLb;
                }

            }

            return [X, Y];
        }

        private static double RLommel(double mu, double v, double b)
        {

            double t = 1 / ((mu + v + 1) * (mu - v + 1));
            double r = t;
            double n = 1;
            double e = .00000000001;

            while (Math.Abs(t) > e * Math.Abs(r))
            {
                t *= -b / (2 * n + mu - v + 1) * (b / (2 * n + mu + v + 1));
                r += t;
                n++;
            }

            return r;
        }

        private static List<double[]> FresnelCS(int n, double t)
        {
            double[] C = new double[n];
            double[] S = new double[n];
            FresnelCS(t, out double c, out double s);
            C[0] = c;
            S[0] = s;

            if (n > 1)
            {
                double tt = Math.PI * t * t / 2;
                double stt = Math.Sin(tt);
                double ctt = Math.Cos(tt);
                C[1] = stt / Math.PI;
                S[1] = (1 - ctt) / Math.PI;
                if (n > 2)
                {
                    C[2] = (t * stt - S[0]) / Math.PI;
                    S[2] = (C[0] - t * ctt) / Math.PI;
                }
            }

            return new List<double[]>() { C, S };
        }


        /// <summary>
        /// Compute the Fresnel integral using a method defined by Venkata Sivakanth Telasula ~2005.
        /// </summary>
        /// <param name="arcLength"></param>
        /// <param name="C"></param>
        /// <param name="S"></param>
        /// <exception cref="Exception"></exception>
        private static void FresnelCS(double arcLength, out double C, out double S)
        {
            double x = Math.Abs(arcLength);
            double EPS = 1E-15;
            double PI_2 = Math.PI / 2;

            if (x < 1.0)
            {
                double term, sum;
                double s = PI_2 * (x * x);
                double t = -s * s;

                // Cosine integral series
                double twofn = 0.0;
                double fact = 1.0;
                double denterm = 1.0;
                double numterm = 1.0;
                sum = 1.0;
                do
                {
                    twofn += 2.0;
                    fact *= twofn * (twofn - 1.0);
                    denterm += 4.0;
                    numterm *= t;
                    term = numterm / (fact * denterm);
                    sum += term;
                } while (Math.Abs(term) > EPS * Math.Abs(sum));

                C = x * sum;

                // Sine integral series
                twofn = 1.0;
                fact = 1.0;
                denterm = 3.0;
                numterm = 1.0;
                sum = 1.0 / 3.0;
                do
                {
                    twofn += 2.0;
                    fact *= twofn * (twofn - 1.0);
                    denterm += 4.0;
                    numterm *= t;
                    term = numterm / (fact * denterm);
                    sum += term;
                } while (Math.Abs(term) > EPS * Math.Abs(sum));

                S = PI_2 * sum * (x * x * x);
            }
            else if (x < 6.0)
            {
                // Rational approximation for f
                double sumn = 0.0;
                double sumd = fd[11];
                for (int k = 10; k >= 0; --k)
                {
                    sumn = fn[k] + x * sumn;
                    sumd = fd[k] + x * sumd;
                }
                double f = sumn / sumd;

                // Rational approximation for g
                sumn = 0.0;
                sumd = gd[11];
                for (int k = 10; k >= 0; --k)
                {
                    sumn = gn[k] + x * sumn;
                    sumd = gd[k] + x * sumd;
                }
                double g = sumn / sumd;

                double U = PI_2 * (x * x);
                double SinU = Math.Sin(U);
                double CosU = Math.Cos(U);

                C = 0.5 + f * SinU - g * CosU;
                S = 0.5 - f * CosU - g * SinU;
            }
            else
            {
                // x >= 6.0; asymptotic expansions for f and g
                double s = Math.PI * x * x;
                double t = -1 / (s * s);

                // Expansion for f
                double numterm = -1.0;
                double term = 1.0;
                double sum = 1.0;
                double oldterm = 1.0;
                double eps10 = 0.1 * EPS;
                double absterm;

                do
                {
                    numterm += 4.0;
                    term *= numterm * (numterm - 2.0) * t;
                    sum += term;
                    absterm = Math.Abs(term);
                    if (oldterm < absterm)
                    {
                        throw new Exception($"In FresnelCS f not converged to eps, x = {x}, oldterm = {oldterm}, absterm = {absterm}");
                    }
                    oldterm = absterm;
                } while (absterm > eps10 * Math.Abs(sum));

                double f = sum / (Math.PI * x);

                // Expansion for g
                numterm = -1.0;
                term = 1.0;
                sum = 1.0;
                oldterm = 1.0;

                do
                {
                    numterm += 4.0;
                    term *= numterm * (numterm + 2.0) * t;
                    sum += term;
                    absterm = Math.Abs(term);
                    if (oldterm < absterm)
                    {
                        throw new Exception($"In FresnelCS g not converged to eps, x = {x}, oldterm = {oldterm}, absterm = {absterm}");
                    }
                    oldterm = absterm;
                } while (absterm > eps10 * Math.Abs(sum));

                double g = sum / (Math.PI * x * (Math.PI * x * x));

                double U = PI_2 * (x * x);
                double SinU = Math.Sin(U);
                double CosU = Math.Cos(U);

                C = 0.5 + f * SinU - g * CosU;
                S = 0.5 - f * CosU - g * SinU;
            }

            if (arcLength < 0)
            {
                C = -C;
                S = -S;
            }
        }

        /// <summary>
        /// Guess the initial value of A given the two values of phi. Phi is the angle difference between the tangent and the vector from start to end. 
        /// </summary>
        /// <param name="phi0"></param>
        /// <param name="phi1"></param>
        /// <returns></returns>
        public static double InitialGuessA(double phi0, double phi1)
        {
            double X = phi0 / Math.PI;
            double Y = phi1 / Math.PI;
            double xy = X * Y;
            X *= X;
            Y *= Y;
            return (phi0 + phi1) * (CF[0] + xy * (CF[1] + xy * CF[2]) + (CF[3] + xy * CF[4]) * (X + Y) + CF[5] * (X * X + Y * Y));
        }

        public static Sample Evaluate(double x0, double y0, double t0, double k0, double dk, double s, FrameModel frame)
        {
            double[][] XY = GeneralizedFresnelCS(1, dk * s * s, k0 * s, t0);
            double t1 = ((dk * s * s) / 2) + (k0 * s) + t0;

            double ct1 = Math.Cos(t1);
            double st1 = Math.Sin(t1);
            double[,] R;
            switch (frame)
            {
                case FrameModel.Frenet:
                    R = ONFrame.R0_FromTN([ct1, st1], [-st1, ct1]);
                    break;
                case FrameModel.Bishop:
                    //Rotate the initial frame by t1 in the XY plane
                    R = RigidTransform.Rotation2D(t1, 0, 0).R;
                    break;
                default:
                    R = ONFrame.R0_FromTN([ct1, st1], [-st1, ct1]);
                    break;
            }
            return new Sample([x0 + (s * XY[0][0]), y0 + (s * XY[1][0])], R, s, [k0 + (dk * s)]);
        }

        public Sample Evaluate(CurveSpec p, double s)
        {
            ICurvatureLaw l = p.Kappa;
            double[] k0 = l.Eval(0);
            //double[] k = l.Eval(s);
            //double endCurvature = k[0];
            //double sharpness = (endCurvature - k0[0]) / s;
            //or
            double sharpness;
            var jet = l.EvalJet(0);
            if (jet.Length > 1 && jet[1].Length > 0)
            {
                sharpness = jet[1][0];
            }
            else
            {
                sharpness = s != 0 ? (l.Eval(s)[0] - k0[0]) / s : 0;
            }
            double[] T0 = p.GetONAxis(0);
            double t0 = Math.Atan2(T0[1], T0[0]);
            double t1 = ((sharpness * s * s) / 2) + (k0[0] * s) + t0;
            double[][] XY = GeneralizedFresnelCS(1, sharpness * s * s, k0[0] * s, t0);

            double ct1 = Math.Cos(t1);
            double st1 = Math.Sin(t1);
            double[,] R;
            switch (p.Frame)
            {
                case FrameModel.Frenet:
                    R = ONFrame.R0_FromTN([ct1, st1], [-st1, ct1]);
                    break;
                case FrameModel.Bishop:
                    //Rotate the initial frame by t1 - t0 in the XY plane
                    R = Helpers.Multiply(RigidTransform.Rotation2D(t1 - t0, 0, 0).R, p.R0);
                    break;
                default:
                    R = ONFrame.R0_FromTN([ct1, st1], [-st1, ct1]);
                    break;
            }
            return new Sample([p.P0[0] + (XY[0][0] * s), p.P0[1] + (XY[1][0] * s)], R, s, l.Eval(s));
        }
    }
}
