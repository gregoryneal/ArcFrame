using System.Numerics;
using System;

namespace ArcFrame.Core.Math
{
    /// <summary>
    /// 3D vector class using doubles.
    /// </summary>
    public readonly struct Vec3d
    {
        /// <summary>
        /// Zero vector
        /// </summary>
        public static readonly Vec3d Zero = new Vec3d(0, 0, 0);
        /// <summary>
        /// Unity vector
        /// </summary>
        public static readonly Vec3d One = new Vec3d(1, 1, 1);
        /// <summary>
        /// Unit Y vector
        /// </summary>
        public static readonly Vec3d UnitY = new Vec3d(0, 1, 0);
        /// <summary>
        /// Unit X vector
        /// </summary>
        public static readonly Vec3d UnitX = new Vec3d(1, 0, 0);
        /// <summary>
        /// Unit Z vector
        /// </summary>
        public static readonly Vec3d UnitZ = new Vec3d(0, 0, 1);
        /// <summary>
        /// X coordinate.
        /// </summary>
        public double X { get => _vals[0]; }
        /// <summary>
        /// Y coordinate.
        /// </summary>
        public double Y { get => _vals[1]; }
        /// <summary>
        /// Z coordinate.
        /// </summary>
        public double Z { get => _vals[2]; }
        private readonly double[] _vals;
        /// <summary>
        /// Create a new vector at the given position.
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <param name="Z"></param>
        public Vec3d(double X, double Y, double Z)
        {
            _vals = new double[3] { X, Y, Z };
        }
        /// <summary>
        /// Create a new vector from another vector.
        /// </summary>
        /// <param name="vector"></param>
        public Vec3d(Vector3 vector)
        {
            _vals = new double[3] { vector.X, vector.Y, vector.Z };
        }
        /// <summary>
        /// Length of the vector
        /// </summary>
        public double Length
        {
            get
            {
                return System.Math.Sqrt((X * X) + (Y * Y) + (Z * Z));
            }
        }
        /// <summary>
        /// Convert to System.Numerics.Vector3
        /// </summary>
        /// <returns></returns>
        public Vector3 ToV3()
        {
            return new Vector3((float)X, (float)Y, (float)Z);
        }
        /// <summary>
        /// Display the properties of the matrix to 6 decimals.
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return $"<{System.Math.Round(X, 6)}, {System.Math.Round(Y, 6)}, {System.Math.Round(Z, 6)}>";
        }
        /// <summary>
        /// Index the coordinate values.
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        /// <exception cref="IndexOutOfRangeException"></exception>
        public double this[int index]
        {
            get
            {
                if (index < 0 || index >= 3) throw new IndexOutOfRangeException("Index must be 0, 1, or 2.");
                return _vals[index];
            }
        }
        /// <summary>
        /// Rotate a vector around a vector axis.
        /// </summary>
        /// <param name="v"></param>
        /// <param name="axis"></param>
        /// <param name="angle"></param>
        /// <returns></returns>
        public static Vec3d RotateAboutAxis(Vec3d v, Vec3d axis, double angle)
        {
            // Rodrigues' rotation formula
            double cosTheta = System.Math.Cos(angle);
            double sinTheta = System.Math.Sin(angle);
            double dot = Dot(v, axis);
            return v * cosTheta + (Cross(axis, v) * sinTheta) + (axis * dot * (1 - cosTheta));
        }
        /// <summary>
        /// Normalize a vector (set its length to 1).
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        /// <exception cref="DivideByZeroException"></exception>
        public static Vec3d Normalize(Vec3d v)
        {
            double length = v.Length;
            if (length == 0) throw new DivideByZeroException("Cannot normalize a zero vector.");
            return v / length;
        }

        /// <summary>
        /// Add two vectors.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Vec3d operator +(Vec3d a, Vec3d b)
        {
            return new Vec3d(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
        }
        /// <summary>
        /// Subtract two vectors.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Vec3d operator -(Vec3d a, Vec3d b)
        {
            return new Vec3d(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
        }
        /// <summary>
        /// Multiply a vector by a scalar.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="scalar"></param>
        /// <returns></returns>
        public static Vec3d operator *(Vec3d a, double scalar)
        {
            return new Vec3d(a.X * scalar, a.Y * scalar, a.Z * scalar);
        }
        /// <summary>
        /// Multiply a vector by a scalar.
        /// </summary>
        /// <param name="scalar"></param>
        /// <param name="a"></param>
        /// <returns></returns>
        public static Vec3d operator *(double scalar, Vec3d a)
        {
            return a * scalar;
        }
        /// <summary>
        /// Divide a vector by a scalar.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="scalar"></param>
        /// <returns></returns>
        /// <exception cref="DivideByZeroException"></exception>
        public static Vec3d operator /(Vec3d a, double scalar)
        {
            if (scalar == 0) throw new DivideByZeroException("Cannot divide by zero.");
            return new Vec3d(a.X / scalar, a.Y / scalar, a.Z / scalar);
        }
        /// <summary>
        /// Negate a vector.
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static Vec3d operator -(Vec3d a)
        {
            return new Vec3d(-a.X, -a.Y, -a.Z);
        }
        /// <summary>
        /// Dot product of two vectors a.x * b.x + a.y * b.y + a.z * b.z
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static double Dot(Vec3d a, Vec3d b)
        {
            return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
        }
        /// <summary>
        /// Cross product of two vectors.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Vec3d Cross(Vec3d a, Vec3d b)
        {
            return new Vec3d(
                a.Y * b.Z - a.Z * b.Y,
                a.Z * b.X - a.X * b.Z,
                a.X * b.Y - a.Y * b.X
            );
        }
        /// <summary>
        /// Lerp between two vectors.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public static Vec3d Lerp(Vec3d a, Vec3d b, double t)
        {
            return a + (b - a) * t;
        }
        /// <summary>
        /// Distance between two vectors.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static double Distance(Vec3d a, Vec3d b)
        {
            if (Vector.IsHardwareAccelerated)
            {
                Vec3d diff = a - b;
                return System.Math.Sqrt(Dot(diff, diff));
            }
            else
            {
                double dx = a.X - b.X;
                double dy = a.Y - b.Y;
                double dz = a.Z - b.Z;
                double ls = dx * dx + dy * dy + dz * dz;
                return System.Math.Sqrt(ls);
            }
        }
        /// <summary>
        /// Implicitly convert to a System.Numerics.Vector3
        /// </summary>
        /// <param name="v"></param>
        public static implicit operator Vector3(Vec3d v)
        {
            return v.ToV3();
        }

    }

}
