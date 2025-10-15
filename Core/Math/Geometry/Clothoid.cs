using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using ArcFrame.Core.Results;
using ArcFrame.Solvers;
using System.Security.Cryptography;

namespace ArcFrame.Core.Geometry
{
    /// <summary>
    /// N Dimensional clothoid curve handler class.
    /// There are no closed form solutions so to evaluate a Sample
    /// we must implement an IEvaluator that takes a CurveSpec and
    /// an arc length and calculates the Sample data.
    /// </summary>
    /// <param name="P">Parameters for the clothoid</param>
    /// <param name="Eval">Evaluator for the clothoid, this generates Samples along the arc length s given P and s.</param>
    public class Clothoid(CurveSpec P, IEvaluator Eval) : IArcLengthCurve
    {
        private readonly CurveSpec _p = P;
        private readonly IEvaluator _eval = Eval;
        public int Dimension => _p.N;

        public double Length => _p.Length;

        public Sample Evaluate(double s) => _eval.Evaluate(_p, s);
        public double[] Position(double s) => Evaluate(s).P;
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

        public static Clothoid From3D(Vec3d P, Vec3d T, Vec3d N, Vec3d B, double k0, double dk, double tau0, double dtau, double L, IEvaluator evaluator, FrameModel frame = FrameModel.Frenet)
        {
            double[] p = [P.X, P.Y, P.Z];
            double[] t = [T.X, T.Y, T.Z];
            double[] n = [N.X, N.Y, N.Z];
            double[] b = [B.X, B.Y, B.Z];
            double[] k0v = { k0, tau0 };
            double[] dkv = { dk, dtau };
            LinearCurvatureLaw law = new LinearCurvatureLaw(k0v, dkv);
            return new Clothoid(new ClothoidCurveSpec(3, L, p, ONFrame.R0_FromTNB(t, n, b), law, frame), evaluator);
        }

        public static Clothoid From3D(double[] P0, double[,] R0, double k0, double dk, double tau0, double dtau, double L, IEvaluator evaluator, FrameModel frame = FrameModel.Frenet)
        {
            double[] k0v = { k0, tau0 };
            double[] dkv = { dk, dtau };
            LinearCurvatureLaw law = new LinearCurvatureLaw(k0v, dkv);
            return new Clothoid(new ClothoidCurveSpec(3, L, P0, R0, law, frame), evaluator);
        }
    }
}
