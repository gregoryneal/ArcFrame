using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using ArcFrame.Core.Results;

namespace ArcFrame.Examples
{
    /// <summary>
    /// Construct an IntrinsicCurve with a Bishop (rotation minimizing) moving frame.
    /// </summary>
    public class IntrinsicCurve_Bishop
    {
        public static void Main()
        {
            // 3d curve
            int dimension = 3;
            double arcLength = 5.2;
            double[] p0 = [0, 0, 0];
            // Build an orthonormal basis frame by rotating the XZ plane by 60deg about the Y axis.
            // does the same as rotating the curve by 60 deg along the Y axis.
            double[,] R0 = RigidTransform.FromYawXZ(dimension, Math.PI / 3, p0[0], p0[2]).R;

            FunctionCurvatureLaw law = new FunctionCurvatureLaw(dimension - 1, s => [s + 2, 2 / s]); //k = s + 2, tau = 2 / s
            CurveSpec spec = new CurveSpec(dimension, arcLength, p0, R0, law, FrameModel.Bishop);
            // Normally you might use a CachedIntrinsicCurve, as it does fewer computations when sampling.
            IntrinsicCurve curve = new IntrinsicCurve(spec);

            //Samples
            Sample[] samples = curve.GetSamples(100);
            double[][] positions = samples.Select(s => s.P).ToArray();
            double[][,] travelingONFrame = samples.Select(s => s.R).ToArray(); //orthonormal basis vectors, add them to positions at matching indices to view the frame
            double[][] curvatures = samples.Select(s => s.k).ToArray();
            int n1 = 0;
            int n2 = 50;
            int n3 = 100;
            //these are ON frames along the curve, add each column to the current position to generate the frame along the curve
            double[,] frame1 = travelingONFrame[n1];
            double[,] frame2 = travelingONFrame[n2];
            double[,] frame3 = travelingONFrame[n3];
            double[] position1 = positions[n1];
            double[] position2 = positions[n2];
            double[] position3 = positions[n3];
        }
    }
}
