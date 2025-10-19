using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;
using ArcFrame.Solvers.BertolazziFrego;

namespace ArcFrame.Core.Params
{
    /// <summary>
    /// Data class to hold intrisic curve data.
    /// </summary>
    public class CurveSpec
    {
        /// <summary>
        /// Dimension of the curve.
        /// </summary>
        public int N { get; }
        /// <summary>
        /// Total length of the curve.
        /// </summary>
        public double Length { get; }
        /// <summary>
        /// Initial position of the curve.
        /// </summary>
        public double[] P0 { get; }
        /// <summary>
        /// Initial orthonormal basis frame.
        /// </summary>
        public double[,] R0 { get; }
        /// <summary>
        /// The curvature model (how curvature will evolve with arc length)
        /// </summary>
        public ICurvatureLaw Kappa { get; }
        /// <summary>
        /// The type of frame (Frenet: curvature guided ONFrame | Bishop: tangent guided ONFrame)
        /// </summary>
        public FrameModel Frame { get; }

        public CurveSpec(int n, double length, double[] p0, double[,] r0, ICurvatureLaw kappa, FrameModel frame)
        {
            N = n;
            Length = length;
            P0 = (double[])p0.Clone();
            R0 = (double[,])r0.Clone();
            Kappa = kappa;
            Frame = frame;
        }

        /// <summary>
        /// Get a column from the R0 matrix.
        /// </summary>
        /// <param name="columnNumber"></param>
        /// <returns></returns>
        public double[] GetONAxis(int columnNumber)
        {
            int n = R0.GetLength(0);
            double[] v = new double[n];
            for (int i = 0; i < n; i++)
            {
                v[i] = R0[i, columnNumber];
            }
            return v;
        }

        /// <summary>
        /// Given a spec build a Line, Arc or 
        /// Build an optimized curve given the curve specs.
        /// Optimized curves are ones with significantly quicker
        /// evaluation methods than using the built in integrator.
        /// </summary>
        /// <returns></returns>
        public IArcLengthCurve GetOptimizedCurve()
        {
            double k = Kappa.Eval(0)[0];
            if (Kappa.IsLinear)
            {
                //Console.WriteLine("Optimizing for clothoid");
                var eval = new BFEvaluator();
                return new Clothoid(this, eval);
            }
            else if (Kappa.IsConstant)
            {
                if (k != 0)
                {
                    //Console.WriteLine("Optimizing for arc");
                    double[] t = GetONAxis(0);
                    double[] n = GetONAxis(1);
                    double startAngle = System.Math.Atan2(t[1], t[0]);
                    double deltaAngle = Length * k;
                    return new Arc(P0, t, n, 1 / k, startAngle, deltaAngle);
                }
                else
                {
                    //Console.WriteLine("Optimizing for line");
                    return new Line(P0, Helpers.Add(P0, Helpers.Multiply(Length, GetONAxis(0))));
                }
            }
            else
            {
                //Console.WriteLine("Curve optimization unavailable, fallback to IntrinsicCurve: ");
                //ShowInfo();
                return new IntrinsicCurve(this);
            }
        }

        public void ShowInfo()
        {
            Console.WriteLine();
            Console.WriteLine("============ Curve Info ============");
            Console.WriteLine($"Dimensions: {N}");
            Console.WriteLine($"Length: {Length}");
            Console.Write($"Start Position: ");
            Helpers.PrintVector(P0);
            Console.WriteLine();
            Console.WriteLine($"Start SO({N}): ");
            Helpers.PrintMat(R0);
            Console.WriteLine();
            Console.WriteLine($"FrameModel: {Frame}");
            Console.WriteLine($"CurvatureLaw: {Kappa.GetType().Name}");
            Console.WriteLine("============ Curve Info ============");
            Console.WriteLine();
        }
    }
}
