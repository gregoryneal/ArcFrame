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
        public readonly double x;
        /// <summary>
        /// Y coordinate.
        /// </summary>
        public readonly double y;
        /// <summary>
        /// Z coordinate.
        /// </summary>
        public readonly double z;
        /// <summary>
        /// Create a new vector at the given position.
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <param name="Z"></param>
        public Vec3d(double X, double Y, double Z)
        {
            x = X;
            y = Y;
            z = Z;
        }

        /// <summary>
        /// Create a new vector from another vector.
        /// </summary>
        /// <param name="vector"></param>
        public Vec3d(Vector3 vector)
        {
            x = vector.X;
            y = vector.Y;
            z = vector.Z;
        }
        /// <summary>
        /// Length of the vector
        /// </summary>
        public double Length => System.Math.Sqrt((x * x) + (y * y) + (z * z));

        /// <summary>
        /// Display the properties of the matrix to 6 decimals.
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return $"<{System.Math.Round(x, 6)}, {System.Math.Round(y, 6)}, {System.Math.Round(z, 6)}>";
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
                if (index == 0) return x;
                if (index == 1) return y;
                else return z;
            }
        }

        /// <summary>
        /// Add two vectors.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Vec3d operator +(Vec3d a, Vec3d b)
        {
            return new Vec3d(a.x + b.x, a.y + b.y, a.z + b.z);
        }
        /// <summary>
        /// Subtract two vectors.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Vec3d operator -(Vec3d a, Vec3d b)
        {
            return new Vec3d(a.x - b.x, a.y - b.y, a.z - b.z);
        }
        /// <summary>
        /// Multiply a vector by a scalar.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="scalar"></param>
        /// <returns></returns>
        public static Vec3d operator *(Vec3d a, double scalar)
        {
            return new Vec3d(a.x * scalar, a.y * scalar, a.z * scalar);
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
            return new Vec3d(a.x / scalar, a.y / scalar, a.z / scalar);
        }
        /// <summary>
        /// Negate a vector.
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static Vec3d operator -(Vec3d a)
        {
            return new Vec3d(-a.x, -a.y, -a.z);
        }
    }
}
