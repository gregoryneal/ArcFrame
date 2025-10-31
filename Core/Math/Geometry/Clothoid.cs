using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using ArcFrame.Core.Results;
using ArcFrame.Solvers;

namespace ArcFrame.Core.Geometry
{
    /// <summary>
    /// N Dimensional clothoid curve handler class.
    /// There are no closed form solutions so to evaluate a Sample
    /// we must implement an IEvaluator that takes a CurveSpec and
    /// an arc length and calculates the Sample data.
    /// None of the papers i have read support bishop frames natively
    /// in the algorithm, so automatically fall back to a 
    /// LieGroupMidpointStepper in that case.
    /// </summary>
    public class Clothoid : IArcLengthCurve
    {
        private readonly CurveSpec _p;
        private readonly IEvaluator _eval;
        private readonly IFrameStepper? _stepper = null;
        private readonly IntegratorOptions? _opts = null;
        private readonly bool _useFrameStepper = false;

        /// <param name="P">Parameters for the clothoid</param>
        /// <param name="Eval">Evaluator for the clothoid, this generates Samples along the arc length s given P and s.</param>
        public Clothoid(CurveSpec P, IEvaluator Eval)
        {
            _p = P;
            _eval = Eval;
            _useFrameStepper = _p.Frame == FrameModel.Bishop || _eval.TargetDimension < Dimension;

            if (_useFrameStepper)
            {
                _stepper = new LieGroupMidpointStepper();
                _opts = IntegratorOptions.Default;
            }
        }
        /// <inheritdoc/>
        public int Dimension => _p.N;
        /// <inheritdoc/>
        public double Length => _p.Length;
        /// <inheritdoc/>
        public Sample Evaluate(double s)
        {
            s = System.Math.Clamp(s, 0.0, _p.Length);

            // Fall back to frame stepper for RMF or if our evaluator is not suitable for
            // solving in this spatial dimension
            if (_useFrameStepper)
            {
                var P = (double[])_p.P0.Clone();
                var R = (double[,])_p.R0.Clone();
                double si = 0;
                while (si < s) {
                    double h = _opts!.SuggestStep(_p.Kappa.Eval, si, s);
                    if (h <= 0) break;
                    _stepper!.Step(ref P, ref R, _p.Kappa.Eval, si, h, _p.Frame);
                    si += h;
                }
                double[] k = _p.Kappa.Eval(s);
                return new Sample(P, R, s, k);
            }

            return _eval.Evaluate(_p, s);
        }
        /// <inheritdoc/>
        public double[] Position(double s) => Evaluate(s).P;
        /// <inheritdoc/>
        public double[] Tangent(double s) => Evaluate(s).T;

        /// <summary>
        /// Get a number of evenly spaced samples along the arc length
        /// </summary>
        /// <param name="count"></param>
        /// <returns></returns>
        public Sample[] GetSamples(int count)
        {
            Sample[] samples = new Sample[count];
            if (count <= 0) count = 1;
            double s;
            for (int i = 0; i < count; i++)
            {
                s = (i == count - 1) ? Length : i * (Length / System.Math.Max(1, count - 1));
                samples[i] = Evaluate(s);
            }
            return samples;
        }

        /// <summary>
        /// Create from 2D clothoid case.
        /// </summary>
        /// <param name="x0"></param>
        /// <param name="y0"></param>
        /// <param name="t0"></param>
        /// <param name="k0"></param>
        /// <param name="dk"></param>
        /// <param name="L"></param>
        /// <param name="evaluator"></param>
        /// <param name="frame"></param>
        /// <returns></returns>
        public static Clothoid From2D(double x0, double y0, double t0, double k0, double dk, double L, IEvaluator evaluator, FrameModel frame = FrameModel.Frenet)
        {
            double[] p0 = { x0, y0 };
            double[,] R0 = new double[2, 2];
            double c = System.Math.Cos(t0);
            double s = System.Math.Sin(t0);

            //T [c, s]
            R0[0, 0] = c;
            R0[1, 0] = s;
            //N [-s, c]
            R0[0, 1] = -s;
            R0[1, 1] = c;

            double[] k0v = { k0 };
            double[] dkv = { dk };
            LinearCurvatureLaw law = new LinearCurvatureLaw(k0v, dkv);
            return new Clothoid(new ClothoidCurveSpec(2, L, p0, R0, law, frame), evaluator);
        }

        /// <summary>
        /// Create a 3D clothoid from a basis and curve parameters.
        /// </summary>
        /// <param name="P"></param>
        /// <param name="T"></param>
        /// <param name="N"></param>
        /// <param name="B"></param>
        /// <param name="k0"></param>
        /// <param name="dk"></param>
        /// <param name="tau0"></param>
        /// <param name="dtau"></param>
        /// <param name="L"></param>
        /// <param name="evaluator"></param>
        /// <param name="frame"></param>
        /// <returns></returns>
        public static Clothoid From3D(Vec3d P, Vec3d T, Vec3d N, Vec3d B, double k0, double dk, double tau0, double dtau, double L, IEvaluator evaluator, FrameModel frame = FrameModel.Frenet)
        {
            double[] p = new double[] { P.X, P.Y, P.Z };
            double[] t = new double[] { T.X, T.Y, T.Z };
            double[] n = new double[] { N.X, N.Y, N.Z };
            double[] b = new double[] { B.X, B.Y, B.Z }
            ;
            double[] k0v = { k0, tau0 };
            double[] dkv = { dk, dtau };
            LinearCurvatureLaw law = new LinearCurvatureLaw(k0v, dkv);
            return new Clothoid(new ClothoidCurveSpec(3, L, p, ONFrame.R0_FromTNB_Complete(t, n, b), law, frame), evaluator);
        }

        /// <summary>
        /// Create a 3D clothoid from a basis and clothoid parameters.
        /// </summary>
        /// <param name="P0"></param>
        /// <param name="R0"></param>
        /// <param name="k0"></param>
        /// <param name="dk"></param>
        /// <param name="tau0"></param>
        /// <param name="dtau"></param>
        /// <param name="L"></param>
        /// <param name="evaluator"></param>
        /// <param name="frame"></param>
        /// <returns></returns>
        public static Clothoid From3D(double[] P0, double[,] R0, double k0, double dk, double tau0, double dtau, double L, IEvaluator evaluator, FrameModel frame = FrameModel.Frenet)
        {
            double[] k0v = { k0, tau0 };
            double[] dkv = { dk, dtau };
            LinearCurvatureLaw law = new LinearCurvatureLaw(k0v, dkv);
            return new Clothoid(new ClothoidCurveSpec(3, L, P0, R0, law, frame), evaluator);
        }

//        public static Clothoid FitTo
    }
}
