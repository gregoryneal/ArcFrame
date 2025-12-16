namespace ArcFrame.Core.Math
{
    /// <summary>
    /// 3x3 Matrix type
    /// </summary>
    public struct Mat3d
    {
#pragma warning disable
        // Row-major 3x3 matrix
        public double m00, m01, m02;
        public double m10, m11, m12;
        public double m20, m21, m22;
#pragma warning restore

        /// <summary>
        /// Identity 3x3
        /// </summary>
        public static Mat3d Identity => new Mat3d
        {
            m00 = 1.0,
            m11 = 1.0,
            m22 = 1.0
        };

        /// <summary>
        /// Extract 3x3 matrix from length 9 array
        /// </summary>
        /// <returns></returns>
        public static Mat3d FromRowMajor(double[] rowMajorVector)
        {
            if (rowMajorVector.Length < 9) throw new System.ArgumentOutOfRangeException("Row major vector needs 9 values to form 3x3 matrix.");

            return new Mat3d
            {
                m00 = rowMajorVector[0],
                m01 = rowMajorVector[1],
                m02 = rowMajorVector[2],
                m10 = rowMajorVector[3],
                m11 = rowMajorVector[4],
                m12 = rowMajorVector[5],
                m20 = rowMajorVector[6],
                m21 = rowMajorVector[7],
                m22 = rowMajorVector[8],
            };
        }
    }
}
