using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using ArcFrame.Solvers.BertolazziFrego;

namespace ArcFrame.Solvers.BertolazziFrego
{
    public class G1Fitter : IFitterG1
    {
        public string Key => "BF_G1";

        public string Author => "Bertolazzi";

        public string Reference => throw new NotImplementedException();
        /// <summary>
        /// Root finding tolerance
        /// </summary>
        private double ROOT_TOLERANCE = 1E-4;
        private int MAX_ITER = 30;
        /// <summary>
        /// The type of frame used.
        /// </summary>
        FrameModel Frame;

        public G1Fitter(FrameModel frame, double tolerance = 1E-4, int maxIterations = 30)
        {
            ROOT_TOLERANCE = tolerance;
            MAX_ITER = maxIterations;
            Frame = frame;
        }

        public CurveSpec Fit_G1(double x0, double y0, double t0, double x1, double y1, double t1)
        {
            double dx = x1 - x0;
            double dz = y1 - y0;
            double phi = Math.Atan2(dz, dx);
            double r = Math.Sqrt(dx * dx + dz * dz);
            double phi0 = NormalizeAngle(t0 - phi);
            double phi1 = NormalizeAngle(t1 - phi);

            double[] a = [x0, y0];
            double[] b = [x1, y1];
            double[] T;
            double[,] R0;

            if (Math.Abs(phi0) < ROOT_TOLERANCE && Math.Abs(phi1) == 0 || phi0 + Math.PI < ROOT_TOLERANCE && phi1 - Math.PI < ROOT_TOLERANCE || phi0 - Math.PI < ROOT_TOLERANCE && phi1 + Math.PI < ROOT_TOLERANCE)
            {

                T = Helpers.Normalize(Helpers.Subtract(b, a));
                R0 = ONFrame.R0_FromT_Complete(T);
                return new CurveSpec(2, r, a, R0, new ConstantCurvatureLaw([0]), Frame);
            }

            double d = phi1 - phi0;

            //Calculate the bounds of the solution A_max and T_max
            double absd = Math.Abs(d);
            double Tmax = Math.Max(0, Math.PI / 2 + Math.Sign(phi1) * phi0);
            double Amax = Tmax == 0 ? absd : absd + 2 * Tmax * (1 + Math.Sqrt(1 + absd / Tmax));

            double g;
            double dg;
            double[][] IntCS;
            int u = 0;
            double A = BFEvaluator.InitialGuessA(phi0, phi1);

            //print amax
            //Console.WriteLine($"{-Amax} <= {A} <= {Amax}");
            if (A > Amax || A < -Amax)
            {
                //print amax
                //Console.WriteLine($"{-Amax} <= {A} <= {Amax}\n");
                A = 0;
            }

            double change;
            double prevChange = double.MaxValue;

            do
            {
                IntCS = BFEvaluator.GeneralizedFresnelCS(3, 2 * A, d - A, phi0);
                g = IntCS[1][0];
                dg = IntCS[0][2] - IntCS[0][1];
                change = g / dg;

                //When the change is larger than the previous change, we are probably oscillating around the root, so cut the change in half
                //I noticed this happening when the initial guess was very close to the actual root. 
                if (Math.Abs(change) >= Math.Abs(prevChange))
                {
                    change = Math.Abs(prevChange) * Math.Sign(change) / 2;
                }
                prevChange = change;

                if (Math.Abs(g) > ROOT_TOLERANCE) A -= change;
                else break;
                //Console.WriteLine($"g: {g}, dg: {dg}, A: {A}, u: {u}");
            } while (++u < MAX_ITER && Math.Abs(g) > ROOT_TOLERANCE);

            double s = r / IntCS[0][0];
            double startCurvature = (d - A) / s;
            double sharpness = 2 * A / (s * s);
            double t = ((sharpness * s * s) / 2) + (startCurvature * s) + t0;
            T = [Math.Cos(t0), Math.Sin(t0)];
            R0 = ONFrame.R0_FromT_Complete(T);
            return new ClothoidCurveSpec(2, s, a, R0, new LinearCurvatureLaw([startCurvature], [sharpness]), Frame);
        }

        /// <summary>
        /// Normalize an angle in radians to be between -pi and pi.
        /// </summary>
        /// <param name="angle"></param>
        /// <returns></returns>
        private static double NormalizeAngle(double angle)
        {
            while (angle > Math.PI) angle -= 2 * Math.PI;
            while (angle < -Math.PI) angle += 2 * Math.PI;
            return angle;
        }       
    }
}
