using ArcFrame.Core.Math;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace ArcFrame.Core.Kernels
{
    

    

    /// <summary>
    /// High performance computing 3d math for blittable types
    /// </summary>
    public static class HpcMath3d
    {
        /// <summary>
        /// Dot product in 3d
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Dot(in Vec3d a, in Vec3d b) => a.x * b.x + a.y * b.y + a.z * b.z;

        /// <summary>
        /// Cross product in 3d
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3d Cross(in Vec3d a, in Vec3d b) => new Vec3d(a.y * b.z - a.z * b.y,
                                                                       a.z * b.x - a.x * b.z,
                                                                       a.x * b.y - a.y * b.x);

        /// <summary>
        /// Rotate a vector around a vector axis.
        /// </summary>
        /// <param name="v"></param>
        /// <param name="axis"></param>
        /// <param name="angle"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3d RotateAboutAxis(Vec3d v, Vec3d axis, double angle)
        {
            // Rodrigues' rotation formula
            double cosTheta = System.Math.Cos(angle);
            double sinTheta = System.Math.Sin(angle);
            double dot = Dot(v, axis);
            return v * cosTheta + (Cross(axis, v) * sinTheta) + (axis * dot * (1 - cosTheta));
        }

        /// <summary>
        /// Distance between two vectors.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Distance(Vec3d a, Vec3d b)
        {
            if (Vector.IsHardwareAccelerated)
            {
                Vec3d diff = a - b;
                return System.Math.Sqrt(Dot(diff, diff));
            }
            else
            {
                double dx = a.x - b.x;
                double dy = a.y - b.y;
                double dz = a.z - b.z;
                double ls = dx * dx + dy * dy + dz * dz;
                return System.Math.Sqrt(ls);
            }
        }

        /// <summary>
        /// Squared length of <see cref="Vec3d"/>
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double LengthSq(in Vec3d v) => Dot(v, v);

        /// <summary>
        /// Length of <see cref="Vec3d"/>
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Length(in Vec3d v) => System.Math.Sqrt(LengthSq(v));

        /// <summary>
        /// Normalize a <see cref="Vec3d"/>
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3d Normalize(in Vec3d v)
        {
            var len = Length(v);
            if (len <= 0.0) return v;
            var inv = 1.0 / len;
            return new Vec3d(v.x * inv, v.y * inv, v.z * inv);
        }

        /// <summary>
        /// Determinant of 3D column vectors a, b, c
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Det(in Vec3d a, in Vec3d b, in Vec3d c)
        {
            // determinant of [a b c] columns (or rows, just be consistent)
            return
                a.x * (b.y * c.z - b.z * c.y) -
                a.y * (b.x * c.z - b.z * c.x) +
                a.z * (b.x * c.y - b.y * c.x);
        }

        /// <summary>
        /// Lerp between two <see cref="Vec3d"/>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="u"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3d Lerp(in Vec3d a, in Vec3d b, double u)
            => new Vec3d(
                a.x + (b.x - a.x) * u,
                a.y + (b.y - a.y) * u,
                a.z + (b.z - a.z) * u
            );

        /// <summary>
        /// Find the projection of v onto the vector onto
        /// </summary>
        /// <param name="v"></param>
        /// <param name="onto"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3d Project(in Vec3d v, in Vec3d onto)
        {
            var denom = LengthSq(onto);
            if (denom <= 0.0) return new Vec3d(0, 0, 0);
            var s = Dot(v, onto) / denom;
            return onto * s;
        }

        /// <summary>
        /// Find the vector that is perpendicular to u whose endpoint ends at v's endpoint.
        /// 
        ///       v   ↗^
        ///          / |
        ///         /  |
        ///        /   |
        ///       /    |
        ///      /     | Reject_u(v)
        ///     /      |
        ///    /       |
        ///   /        |
        ///  /         |
        ///  --------->----------------> onto
        ///  Project_onto(v)
        /// 
        /// Reject_u(onto) = v - Dot(onto, v)onto
        /// </summary>
        /// <param name="v">Normalized ND vector</param>
        /// <param name="onto">Normalized ND vector, this will be perpendicular to the return vector.</param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3d Reject(in Vec3d v, in Vec3d onto)
        {
            var proj = Project(v, onto);
            return new Vec3d(v.x - proj.x, v.y - proj.y, v.z - proj.z);
        }

        /// <summary>
        /// Multiply a <see cref="Mat3d"/> and a <see cref="Vec3d"/>
        /// </summary>
        /// <param name="m"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3d Mul(in Mat3d m, in Vec3d v)
            => new Vec3d(
                m.m00 * v.x + m.m01 * v.y + m.m02 * v.z,
                m.m10 * v.x + m.m11 * v.y + m.m12 * v.z,
                m.m20 * v.x + m.m21 * v.y + m.m22 * v.z
            );

        /// <summary>
        /// Multiply two <see cref="Mat3d"/>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Mat3d Mul(in Mat3d a, in Mat3d b)
        {
            Mat3d r;

            r.m00 = a.m00 * b.m00 + a.m01 * b.m10 + a.m02 * b.m20;
            r.m01 = a.m00 * b.m01 + a.m01 * b.m11 + a.m02 * b.m21;
            r.m02 = a.m00 * b.m02 + a.m01 * b.m12 + a.m02 * b.m22;

            r.m10 = a.m10 * b.m00 + a.m11 * b.m10 + a.m12 * b.m20;
            r.m11 = a.m10 * b.m01 + a.m11 * b.m11 + a.m12 * b.m21;
            r.m12 = a.m10 * b.m02 + a.m11 * b.m12 + a.m12 * b.m22;

            r.m20 = a.m20 * b.m00 + a.m21 * b.m10 + a.m22 * b.m20;
            r.m21 = a.m20 * b.m01 + a.m21 * b.m11 + a.m22 * b.m21;
            r.m22 = a.m20 * b.m02 + a.m21 * b.m12 + a.m22 * b.m22;

            return r;
        }

        /// <summary>
        /// Outer product of two <see cref="Vec3d"/>
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Mat3d Outer(in Vec3d u, in Vec3d v)
        {
            Mat3d m;
            m.m00 = u.x * v.x; m.m01 = u.x * v.y; m.m02 = u.x * v.z;
            m.m10 = u.y * v.x; m.m11 = u.y * v.y; m.m12 = u.y * v.z;
            m.m20 = u.z * v.x; m.m21 = u.z * v.y; m.m22 = u.z * v.z;
            return m;
        }
    }

}
