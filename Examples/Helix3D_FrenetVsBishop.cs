using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using ArcFrame.Core.Results;
using System.Linq;

namespace ArcFrame.Examples
{
    /// <summary>
    /// Create a Helix in 3D, two different ways. The first way uses the Frenet Serret frame, which uses
    /// invariants of the curve like curvature and torsion, to evolve the frame at each step.
    /// The Bishop frame on the other hand is a rotation minimizing frame. The tangent goes with the curve
    /// like normal. A normal vector is chosen at P0. Then at each (small) step along the integrator path 
    /// we choose the frame that minimizes the rotation of this initial [T,N,..] frame
    /// </summary>
    public class Helix3D_FrenetVsBishop
    {
        /// <summary>
        /// Example
        /// </summary>
        public static void Main()
        {
            // Shared start pose
            var P0 = new[] { 0.0, 0.0, 0.0 };
            var R0 = RigidTransform.Identity(3).R; //identity matrix in 3D

            double radius = 1.0;
            double pitch = 0.5;
            double kappa = radius / (radius * radius + pitch * pitch);
            double tau = pitch / (radius * radius + pitch * pitch);
            double L = 40.0;

            ConstantCurvatureLaw fsLaw = new ConstantCurvatureLaw(new double[] { kappa, tau });
            CurveSpec specFS = new CurveSpec(3, L, P0, R0, fsLaw, FrameModel.Frenet);
            CachedIntrinsicCurve helixFS = new CachedIntrinsicCurve(specFS);

            //the normal vector doesn't spin around the tangent vector so the curvatures have to.
            FunctionCurvatureLaw bsLaw = new FunctionCurvatureLaw(2, s => new double[] { kappa* System.Math.Cos(tau * s), kappa* System.Math.Sin(tau * s)});
            CurveSpec specB = new CurveSpec(3, L, P0, R0, bsLaw, FrameModel.Bishop);
            CachedIntrinsicCurve helixBS = new CachedIntrinsicCurve(specB);

            //plot these samples
            double[][] samplesFS = helixFS.GetSamples(100).Select(s => s.P).ToArray(); // [ [p0x, p0y,..], [p1x, p1y, ...], ...]
            double[][] samplesBS = helixBS.GetSamples(100).Select(s => s.P).ToArray();
        }
    }
}
