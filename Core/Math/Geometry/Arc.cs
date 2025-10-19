using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using ArcFrame.Core.Results;

namespace ArcFrame.Core.Geometry
{
    /// <summary>
    /// Build an Arc in N dimensional space
    /// </summary>
    public sealed class Arc : IArcLengthCurve
    {
        //center and two basis vectors to establish the 2D subplane the circle inhabits
        public double[] _c, _e1, _e2;
        private double[,] _R0; // starting orthornormal frame, column vectors are R0 = [e1, e2, norm(e1xe2), ...]
        //radius, start angle and angle delta. If the arc is a cirle centered in the XZ plane, _a0 = 0 would correspond to the point (_r, 0), _ad > 0 corresponds to CCW rotations and _ad < 0 to CW rotations
        public double _r, _a0, _ad, _len, _sign;
        public int _n;
        public int Dimension => _n;
        public double Length => _len;

        /// <summary>
        /// An arc in N dimensional space
        /// </summary>
        /// <param name="center">Center of the arc in N dimensions</param>
        /// <param name="e1Unit">First unit basis vector of the local plane</param>
        /// <param name="e2Unit">Second unit basis vector of the local plane</param>
        /// <param name="radius">Arc radius</param>
        /// <param name="startAngle">Start angle, 0 corresponds to the point radius*e1</param>
        /// <param name="deltaAngle">The change in angle of the arc, positive -> CCW rotation when viewed from the e1Xe2 axis</param>
        public Arc(double[] center, double[] e1Unit, double[] e2Unit, double radius, double startAngle, double deltaAngle)
        {
            _n = center.Length;
            _c = (double[])center.Clone();
            _e1 = (double[])e1Unit.Clone();
            _e2 = (double[])e2Unit.Clone();
            _r = System.Math.Abs(radius);
            _a0 = startAngle;
            _ad = deltaAngle;
            _len = System.Math.Abs(_r * _ad);
            _sign = _ad >= 0 ? 1.0 : -1.0;

            double[] T0 = new double[_n];
            double[] N0 = new double[_n];
            double[] P0 = new double[_n];
            double sina = System.Math.Sin(_a0);
            double cosa = System.Math.Cos(_a0);
            for (int i = 0; i < _n; i++)
            {
                P0[i] = _c[i] + (_r * ((cosa * _e1[i]) + (sina * _e2[i])));
                T0[i] = _sign * ((-sina * _e1[i]) + (cosa * _e2[i]));
                N0[i] = _sign * ((-cosa * _e1[i]) + (-sina * _e2[i]));
            }
            _R0 = ONFrame.R0_FromT_Complete(T0, N0);
        }

        public Sample Evaluate(double s)
        {
            if (s <= 0) s = 0;
            else if (s >= _len) s = _len;
            double t = _len > 0 ? s / _len : 0;
            double a = _a0 + (_ad * t);
            double cosa = System.Math.Cos(a);
            double sina = System.Math.Sin(a);

            double[] P = new double[_n];
            double[] T = new double[_n];
            double[] N = new double[_n];
            double[,] R = (double[,])_R0.Clone();


            for (int i = 0; i < _n; i++)
            {
                P[i] = _c[i] + (_r * ((cosa * _e1[i]) + (sina * _e2[i])));
                // only tangent and normals change with an arc, other normals do not
                R[i, 0] = _sign * ((-sina * _e1[i]) + (cosa * _e2[i])); // T
                R[i, 1] = _sign * ((-cosa * _e1[i]) + (-sina * _e2[i]));// N
            }

            double[] K = new double[System.Math.Max(0, _n - 1)];
            if (K.Length > 0) K[0] = 1 / _r;
            return new Sample(P, R, s, K);
        }

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
        /// Build an arc from XZ plane-parameterized arc.
        /// </summary>
        /// <param name="cx"></param>
        /// <param name="cz"></param>
        /// <param name="radius"></param>
        /// <param name="startAngle"></param>
        /// <param name="deltaAngle"></param>
        /// <returns></returns>
        public static Arc From3D_XZ(double cx, double cz, double radius, double startAngle, double deltaAngle)
        {
            double[] e1 = { 1, 0, 0 };
            double[] e2 = { 0, 0, 1 };
            double[] c = { cx, 0, cz };

            return new Arc(c, e1, e2, radius, startAngle, deltaAngle);
        }

        /// <summary>
        /// Build an arc from XZ plane-parameterized arc.
        /// </summary>
        /// <param name="cx"></param>
        /// <param name="cz"></param>
        /// <param name="radius"></param>
        /// <param name="startAngle"></param>
        /// <param name="deltaAngle"></param>
        /// <returns></returns>
        public static Arc From2D(double cx, double cy, double radius, double startAngle, double deltaAngle)
        {
            double[] e1 = { 1, 0 };
            double[] e2 = { 0, 1 };
            double[] c = { cx, cy };

            return new Arc(c, e1, e2, radius, startAngle, deltaAngle);
        }

        /// <summary>
        /// Build from an arc in 3D space.
        /// </summary>
        /// <param name="center"></param>
        /// <param name="planeNormal"></param>
        /// <param name="radius"></param>
        /// <param name="startAngle"></param>
        /// <param name="deltaAngle"></param>
        /// <returns></returns>
        public static Arc From3D(double[] center, double[] planeNormal3, double radius, double startAngle, double deltaAngle)
        {
            double[] normalizedPlaneNormal = Helpers.Normalize(planeNormal3);
            double[] tmp = System.Math.Abs(normalizedPlaneNormal[0]) < .9 ? new double[] { 1, 0, 0 } : new double[] { 0, 1, 0 };
            double[] e1 = Helpers.Cross3(normalizedPlaneNormal, tmp);
            double[] e2 = Helpers.Normalize(Helpers.Cross3(normalizedPlaneNormal, e1));

            return new Arc(center, e1, e2, radius, startAngle, deltaAngle);
        }

        /// <summary>
        /// Build an arc from a CurveSpec. The angle delta will always be positive.
        /// </summary>
        /// <param name="spec"></param>
        /// <param name="arc"></param>
        /// <returns></returns>
        public static bool TryFromSpec(CurveSpec spec, out Arc? arc)
        {
            if (spec.N != 2)
            {
                arc = null;
                return false;
            }
            double k = spec.Kappa.Eval(0)[0];
            if (k != 0)
            {
                double[] t = spec.GetONAxis(0);
                double[] n = spec.GetONAxis(1);
                double startAngle = System.Math.Atan2(t[1], t[0]);
                double deltaAngle = spec.Length * k;
                arc = new Arc(spec.P0, t, n, 1 / k, startAngle, deltaAngle);
                return true;
            }
            arc = null;
            return false;
        }

        public double[] Position(double s) => Evaluate(s).P;
        public double[] Tangent(double s) => Evaluate(s).T;
    }
}
