using ArcFrame.Core.Geometry;
using ArcFrame.Core.Math;

namespace ArcFrame.Examples
{
    /// <summary>
    /// This example shows how to join segments, ensuring G1 continuity. using the composite curve AddG1 feature.
    /// </summary>
    public class CompositeCurve_G1Chain
    {
        public static void Main()
        {
            double[] center = [0, 0, 0];
            double[] endpoint1 = [2, 0, 0];
            double[] endpoint2 = [0, 2, 0];

            Line line1 = new Line(center, endpoint1);
            Line line2 = new Line(center, endpoint2); //lines are perpendicular
            Arc arc = new Arc(center, [1, 0, 0], [0, 0, 1], 3, Math.PI, Math.PI / 3); //60deg arc with radius 3 in the XZ plane centered on the origin

            CompositeCurve curve = new CompositeCurve();
            curve.Add(arc).AddG1(line1, out RigidTransform xfUsed1).AddG1(line2, out RigidTransform xfUsed2);

            //You can inspect the transforms that were done on the primitives
            //Pf = R*P0 + T
            xfUsed1.PrintR(); //print Rotation matrix to console 
            xfUsed1.PrintT(); //print Translation vector to console.
            xfUsed2.PrintR();
            xfUsed2.PrintT();

            //You can rebuild the same curve with the primitives like so
            CompositeCurve curve2 = new CompositeCurve();
            TransformedCurve tLine1 = new TransformedCurve(line1, xfUsed1);
            TransformedCurve tLine2 = new TransformedCurve(line2, xfUsed2);
            curve2.Add(arc).Add(tLine1).Add(tLine2); //add without G1, instead pre-transforming the lines first.
        }
    }
}
