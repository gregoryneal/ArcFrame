using System;

namespace ArcFrame.Core.Math
{
    /// <summary>
    /// Holds rigid transformation data to be applied on a vector.
    /// </summary>
    public sealed class RigidTransform
    {
        /// <summary>
        /// The dimension of the transformation.
        /// </summary>
        public int N { get; }
        //rotation in N dimensions
        /// <summary>
        /// The rotation matrix of the transformation.
        /// </summary>
        public double[,] R { get; }
        //translation in N dimensions
        /// <summary>
        /// The translation vector of the transformation.
        /// </summary>
        public double[] T { get; }

        /// <summary>
        /// Create a new transform given a rotation matrix and translation.
        /// </summary>
        /// <param name="R"></param>
        /// <param name="T"></param>
        /// <exception cref="ArgumentException"></exception>
        public RigidTransform(double[,] R, double[] T)
        {
            int n0 = R.GetLength(0);
            int n1 = R.GetLength(1);
            if (n0 != n1) throw new ArgumentException("R must be square");
            if (T.Length != n0) throw new ArgumentException("T length must match R size");
            this.R = (double[,])R.Clone();
            this.T = (double[])T.Clone();
            this.N = n0;
        }

        /// <summary>
        /// Compose: apply this transform, then apply other.
        /// f ∘ g = f(g(x)) where this transform is g and other is f.
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public RigidTransform Compose(RigidTransform other)
        {
            if (other.N != N) throw new ArgumentException("dimension mismatch!");
            var RR = new double[N, N];
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    //dot the ith row of other with the jth column of this rotation matrix
                    double dot = 0;
                    for (int k = 0; k < N; k++)
                    {
                        dot += other.R[i, k] * R[k, j];
                    }
                    RR[i, j] = dot;
                }
            }

            //R * t => apply rotation then translation to T
            var TT = new double[N];
            for (int i = 0; i < N; i++)
            {
                //dot product of row i with T
                double dot = 0;
                for (int j = 0; j < N; j++)
                {
                    dot += other.R[i, j] * T[j];
                }

                TT[i] = dot + other.T[i];
            }

            return new RigidTransform(RR, TT);
        }

        /// <summary>
        /// Apply rotation and translation to point
        /// 
        /// p1 = (R * p0) + T
        /// </summary>
        /// <param name="p">N dimensional point</param>
        /// <returns></returns>
        public double[] ApplyTransform(double[] p)
        {
            if (p.Length != N) throw new ArgumentException("dimension mismatch!");
            var r = new double[N];
            double dot;
            for (int i = 0; i < N; i++)
            {
                dot = 0;
                for (int j = 0; j < N; j++)
                {
                    dot += R[i, j] * p[j];
                }
                r[i] = dot + T[i];
            }
            return r;
        }

        /// <summary>
        /// Apply rotation only, useful for tangents.
        /// </summary>
        /// <param name="p"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public double[] ApplyRotation(double[] p)
        {
            if (p.Length != N) throw new ArgumentException("dimension mismatch!");
            var r = new double[N];
            double dot;
            for (int i = 0; i < N; i++)
            {
                dot = 0;
                for (int j = 0; j < N; j++)
                {
                    dot += R[i, j] * p[j];
                }
                r[i] = dot;
            }
            return r;
        }

        /// <summary>
        /// Print the rotation matrix to the console.
        /// </summary>
        public void PrintR()
        {
            Console.WriteLine("[");
            for (int i = 0; i < N; i++)
            {
                Console.Write("[");
                for (int j = 0; j < N; j++)
                {
                    Console.Write(R[i, j]);
                    if (j < N - 1) Console.Write(", ");
                }
                if (i < N - 1) Console.WriteLine("],");
                else Console.WriteLine("]");
            }
            Console.WriteLine("]");
        }

        /// <summary>
        /// Helper to print the translation vector to the console.
        /// </summary>
        public void PrintT()
        {
            Console.Write("[");
            for (int i = 0; i < N; i++)
            {
                Console.Write(T[i]);
                if (i < N - 1) Console.Write(", ");
            }
            Console.WriteLine("]");
        }

        /// <summary>
        /// Identity transformation
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static RigidTransform Identity(int n)
        {
            var R = new double[n, n];
            for (int i = 0; i < n; i++) R[i, i] = 1;
            var T = new double[n];
            return new RigidTransform(R, T);
        }

        /// <summary>
        /// Build an XZ yaw rotation (angle about +Y) embedded in N>=3, with pivot and extra translation.
        /// </summary>
        /// <param name="n">n dimensions</param>
        /// <param name="angle">angle about y (dimension 2)</param>
        /// <param name="pivotX">center of rotation X (dimension 1)</param>
        /// <param name="pivotZ">center of rotation Z (dimension 3)</param>
        /// <param name="tx">X translation</param>
        /// <param name="tz">Z translation</param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static RigidTransform FromYawXZ(int n, double angle, double pivotX, double pivotZ, double tx = 0, double tz = 0)
        {
            if (n < 3) throw new ArgumentException("N must be >= 3 for XZ yaw");
            var R = Identity(n).R; // start as identity
            double c = System.Math.Cos(angle), s = System.Math.Sin(angle);
            // 2x2 rotation in (x,z) i.e., indices (0,2)
            R[0, 0] = c; R[0, 2] = -s;
            R[2, 0] = s; R[2, 2] = c;

            // translation: t = (I - R)*pivot + extra
            var t = new double[n];
            t[0] = (1 - R[0, 0]) * pivotX - R[0, 2] * pivotZ + tx;   // (I - R)*p
            t[2] = -R[2, 0] * pivotX + (1 - R[2, 2]) * pivotZ + tz;
            // others zero
            return new RigidTransform(R, t);
        }

        /// <summary>
        /// Create a roll transformation in 3D in the XY plane.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="angle"></param>
        /// <param name="pivotX"></param>
        /// <param name="pivotY"></param>
        /// <param name="tx"></param>
        /// <param name="ty"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static RigidTransform FromRollXY(int n, double angle, double pivotX, double pivotY, double tx = 0, double ty = 0)
        {
            if (n < 3) throw new ArgumentException("N must be >= 3 for XY roll");
            var R = Identity(n).R; // start as identity
            double c = System.Math.Cos(angle), s = System.Math.Sin(angle);
            // 2x2 rotation in (x,y) i.e., indices (0,1)
            R[0, 0] = c; R[0, 1] = -s;
            R[1, 0] = s; R[1, 1] = c;

            // translation: t = (I - R)*pivot + extra
            var t = new double[n];
            t[0] = (1 - R[0, 0]) * pivotX - R[0, 1] * pivotY + tx;   // (I - R)*p
            t[1] = -R[1, 0] * pivotX + (1 - R[1, 1]) * pivotY + ty;
            // others zero
            return new RigidTransform(R, t);
        }

        /// <summary>
        /// Create a pitch transformation in 3D in the YZ plane.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="angle"></param>
        /// <param name="pivotY"></param>
        /// <param name="pivotZ"></param>
        /// <param name="ty"></param>
        /// <param name="tz"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static RigidTransform FromPitchYZ(int n, double angle, double pivotY, double pivotZ, double ty = 0, double tz = 0)
        {
            if (n < 3) throw new ArgumentException("N must be >= 3 for XZ yaw");
            var R = Identity(n).R; // start as identity
            double c = System.Math.Cos(angle), s = System.Math.Sin(angle);
            // 2x2 rotation in (y,z) i.e., indices (1,2)
            R[1, 1] = c; R[1, 2] = -s;
            R[2, 1] = s; R[2, 2] = c;

            // translation: t = (I - R)*pivot + extra
            var t = new double[n];
            t[1] = (1 - R[1, 1]) * pivotY - R[1, 2] * pivotZ + ty;   // (I - R)*p
            t[2] = -R[2, 1] * pivotY + (1 - R[2, 2]) * pivotZ + tz;
            // others zero
            return new RigidTransform(R, t);
        }

        /// <summary>
        /// Create a 2D rotation transorm.
        /// </summary>
        /// <param name="angle"></param>
        /// <param name="pivotX"></param>
        /// <param name="pivotY"></param>
        /// <param name="tx"></param>
        /// <param name="ty"></param>
        /// <returns></returns>
        public static RigidTransform Rotation2D(double angle, double pivotX, double pivotY, double tx = 0, double ty = 0)
        {
            var R = Identity(2).R;
            double c = System.Math.Cos(angle);
            double s = System.Math.Sin(angle);
            R[0, 0] = c;
            R[1, 0] = s;
            R[0, 1] = -s;
            R[1, 1] = c;

            // translation: t = (I - R)*pivot + extra
            var t = new double[2];
            t[0] = (1 - R[0, 0]) * pivotX - R[0, 1] * pivotY + tx;   // (I - R)*p
            t[1] = -R[1, 0] * pivotX + (1 - R[1, 1]) * pivotY + ty;
            return new RigidTransform(R, t);
        }
    }

    /// <summary>
    /// Extension methods for the rigid transforms.
    /// </summary>
    public static class RigitTransformExtensions
    {
        // Expand a smaller transform to 'targetN' by block-diagonal identity and zero padding.
        // axisMap selects which curve axes the small transform acts on (default: first M axes).

        /// <summary>
        /// Expand a smaller transformation into a target dimension N, with optional axis remapping.
        /// </summary>
        /// <param name="xfSmall"></param>
        /// <param name="targetN"></param>
        /// <param name="axisMap"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static RigidTransform ExpandedTo(this RigidTransform xfSmall, int targetN, int[]? axisMap = null)
        {
            int M = xfSmall.N;
            if (targetN < M) throw new ArgumentException("targetN must be >= transform dimension.");
            var map = axisMap ?? DefaultMap(M, targetN);

            var R = new double[targetN, targetN];
            for (int i = 0; i < targetN; i++) R[i, i] = 1.0; // identity
                                                             // insert small rotation block into selected axes
            for (int i = 0; i < M; i++)
                for (int j = 0; j < M; j++)
                    R[map[i], map[j]] = xfSmall.R[i, j];

            var t = new double[targetN];
            for (int i = 0; i < M; i++) t[map[i]] = xfSmall.T[i];

            return new RigidTransform(R, t);

            static int[] DefaultMap(int m, int N)
            {
                var a = new int[m];
                for (int i = 0; i < m; i++) a[i] = i;
                return a;
            }
        }
    }
}
