using System.Numerics;
using System;

namespace ArcFrame.Core.Math
{
    /// <summary>
    /// 3D vector class using doubles.
    /// </summary>
    public readonly struct Vec3d
    {
        public static readonly Vec3d Zero = new Vec3d(0, 0, 0);
        public static readonly Vec3d One = new Vec3d(1, 1, 1);
        public static readonly Vec3d UnitY = new Vec3d(0, 1, 0);
        public static readonly Vec3d UnitX = new Vec3d(1, 0, 0);
        public static readonly Vec3d UnitZ = new Vec3d(0, 0, 1);

        public double X { get => _vals[0]; }
        public double Y { get => _vals[1]; }
        public double Z { get => _vals[2]; }
        private readonly double[] _vals;
        public Vec3d(double X, double Y, double Z)
        {
            _vals = new double[3] { X, Y, Z };
        }

        public Vec3d(Vector3 vector)
        {
            _vals = new double[3] { vector.X, vector.Y, vector.Z };
        }

        public double Length
        {
            get
            {
                return System.Math.Sqrt((X * X) + (Y * Y) + (Z * Z));
            }
        }

        public Vector3 ToV3()
        {
            return new Vector3((float)X, (float)Y, (float)Z);
        }

        public override string ToString()
        {
            return $"<{System.Math.Round(X, 6)}, {System.Math.Round(Y, 6)}, {System.Math.Round(Z, 6)}>";
        }

        public double this[int index]
        {
            get
            {
                if (index < 0 || index >= 3) throw new IndexOutOfRangeException("Index must be 0, 1, or 2.");
                return _vals[index];
            }
        }

        public static Vec3d RotateAboutAxis(Vec3d v, Vec3d axis, double angle)
        {
            // Rodrigues' rotation formula
            double cosTheta = System.Math.Cos(angle);
            double sinTheta = System.Math.Sin(angle);
            double dot = Dot(v, axis);
            return v * cosTheta + (Cross(axis, v) * sinTheta) + (axis * dot * (1 - cosTheta));
        }

        public static Vec3d Normalize(Vec3d v)
        {
            double length = v.Length;
            if (length == 0) throw new DivideByZeroException("Cannot normalize a zero vector.");
            return v / length;
        }

        public static Vec3d operator +(Vec3d a, Vec3d b)
        {
            return new Vec3d(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
        }
        public static Vec3d operator -(Vec3d a, Vec3d b)
        {
            return new Vec3d(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
        }
        public static Vec3d operator *(Vec3d a, double scalar)
        {
            return new Vec3d(a.X * scalar, a.Y * scalar, a.Z * scalar);
        }
        public static Vec3d operator *(double scalar, Vec3d a)
        {
            return a * scalar;
        }
        public static Vec3d operator /(Vec3d a, double scalar)
        {
            if (scalar == 0) throw new DivideByZeroException("Cannot divide by zero.");
            return new Vec3d(a.X / scalar, a.Y / scalar, a.Z / scalar);
        }
        public static Vec3d operator -(Vec3d a)
        {
            return new Vec3d(-a.X, -a.Y, -a.Z);
        }

        public static double Dot(Vec3d a, Vec3d b)
        {
            return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
        }
        public static Vec3d Cross(Vec3d a, Vec3d b)
        {
            return new Vec3d(
                a.Y * b.Z - a.Z * b.Y,
                a.Z * b.X - a.X * b.Z,
                a.X * b.Y - a.Y * b.X
            );
        }

        public static Vec3d Lerp(Vec3d a, Vec3d b, double t)
        {
            return a + (b - a) * t;
        }

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

        public static implicit operator Vector3(Vec3d v)
        {
            return v.ToV3();
        }

    }

}
