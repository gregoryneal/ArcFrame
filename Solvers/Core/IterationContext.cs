using ArcFrame.Core.Geometry;
using ArcFrame.Core.Params;
using ArcFrame.Core.Results;
using System;
using System.Collections.Generic;

namespace ArcFrame.Solvers.Core
{
    /// <summary>
    /// Per-iteration shared data: prefix lengths, total length, uniform global samples,
    /// and the single-pass P/T evaluations for all constraints to reuse.
    /// 
    /// The idea is that at each iteration we build one of these and pass it to each
    /// Constraint who can evaluate it at the cached location instead of building its
    /// own IntrinsicCurve each time it needs to caculate the residual. At the end of
    /// the loop we update each spec and recalculate the samples again.
    /// It isn't implemented yet but one day I will try to use it for optimization.
    /// </summary>
    public sealed class IterationContext : IDisposable
    {
        /// <summary>
        /// Which mode is used by the parent composite curve solver. This is mainly for residuals to see.
        /// </summary>
        public CompositeCurveSolver.CompositeSolverMode Mode = CompositeCurveSolver.CompositeSolverMode.Seperate;
        /// <summary>
        /// These CurveSpecs represent a contiguous concatenated curve,
        /// not multiple curves starting from the same position
        /// </summary>
        public IReadOnlyList<CurveSpec> Specs { get; }
        /// <summary>
        /// Curves
        /// </summary>
        public IArcLengthCurve[] Curves { get; }

        /// <summary>
        /// The total arc length
        /// at the start of each segment.
        /// Also includes the total arc length
        /// at the end.
        /// </summary>
        public double[] Prefix { get; }    // prefix[i] = sum_{k< i} Length_k
        /// <summary>
        /// Total curve length
        /// </summary>
        public double TotalLength { get; }
        /// <summary>
        /// Sampling step
        /// </summary>
        public double Ds { get; }        // sampling step used
        /// <summary>
        /// Samples to take
        /// </summary>
        public int SampleCount { get; }
        /// <summary>
        /// use a large ds vs a small one.
        /// </summary>
        public bool FastMode { get; }  // handy toggle for coarse vs fine

        // Global samples in arclength over the concatenated curve
        /// <summary>
        /// Global sample arc length
        /// </summary>
        public double[] SGlob { get; }     // size M
        /// <summary>
        /// Segment index
        /// </summary>
        public int[] SegIdx { get; }    // size M
        /// <summary>
        /// Local arc lengths
        /// </summary>
        public double[] SLocal { get; }    // size M

        // Shared evaluations (size M)
        /// <summary>
        /// Sample cache.
        /// </summary>
        public Sample[] S { get; }       // samples

        /// <summary>
        /// Build an IterationContext with a cached set of
        /// precomputed samples along a CompositeCurve
        /// comprised of the input CurveSpecs
        /// </summary>
        /// <param name="specs"></param>
        /// <param name="curves"></param>
        /// <param name="prefix"></param>
        /// <param name="Ltot"></param>
        /// <param name="ds"></param>
        /// <param name="fast"></param>
        /// <param name="sGlob"></param>
        /// <param name="segIdx"></param>
        /// <param name="sLocal"></param>
        /// <param name="s"></param>
        private IterationContext(IReadOnlyList<CurveSpec> specs,
                                 IArcLengthCurve[] curves,
                                 double[] prefix, double Ltot,
                                 double ds, bool fast,
                                 double[] sGlob, int[] segIdx, double[] sLocal,
                                 Sample[] s)
        {
            Specs = specs;
            Curves = curves;
            Prefix = prefix;
            TotalLength = Ltot;
            Ds = ds;
            FastMode = fast;
            SGlob = sGlob;
            SegIdx = segIdx;
            SLocal = sLocal;
            S = s;
            SampleCount = sGlob.Length;
        }

        /// <summary>
        /// Build the context from specs and a sample step size.
        /// Does not guarantee endpoint samples.
        /// </summary>
        /// <param name="specs"></param>
        /// <param name="N">number of samples to take, including the start and endpoint</param>
        /// <param name="fastMode"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static IterationContext Build(IReadOnlyList<CurveSpec> specs, int N, bool fastMode)
        {
            if (specs == null || specs.Count == 0)
                throw new ArgumentException("IterationContext.Build: specs empty.");

            // Prefix & curves
            var curves = new IArcLengthCurve[specs.Count];
            var prefix = new double[specs.Count + 1];
            prefix[0] = 0.0;
            for (int i = 0; i < specs.Count; i++)
            {
                curves[i] = new CachedIntrinsicCurve(specs[i], specs[i].Length / 20);
                prefix[i + 1] = prefix[i] + specs[i].Length;
            }
            double Ltot = prefix[specs.Count];
            // Num sample points, ensure at least start and endpoint
            int M = Math.Max(2, N);

            var sGlob = new double[M];
            var segIdx = new int[M];
            var sLocal = new double[M];
            var S = new Sample[M];

            // Walk segments once; avoid per-sample binary search
            int seg = 0;
            for (int m = 0; m < M; m++)
            {
                // global s
                double sg = (Ltot * m) / (M - 1.0);
                // seg => which segment we are in by index
                while (seg + 1 < prefix.Length && sg >= prefix[seg + 1]) seg++;
                seg = System.Math.Clamp(seg, 0, curves.Length - 1); // hack
                // local s
                double sl = sg - prefix[seg];

                //Console.WriteLine($"m: {m} | M: {M} | sg: {sg} | sl: {sl} | prefix.Len: {prefix.Length} | curves len: {curves.Length} | index: {seg}");
                var smp = curves[seg].Evaluate(sl); // single expensive call per sample
                // cache the global and local s, the segment index and the Sample itself
                sGlob[m] = sg;
                segIdx[m] = seg;
                sLocal[m] = sl;
                S[m] = smp;
            }

            return new IterationContext(specs, curves, prefix, Ltot, Ltot / (M - 1), fastMode, sGlob, segIdx, sLocal, S)
            {
                Mode = CompositeCurveSolver.CompositeSolverMode.Seperate
            };
        }

        /// <summary>
        /// Build the context from specs and a sample step size.
        /// This one is different in that it ensures that the
        /// start and end of each segment is sampled. 
        /// The solver will include the position and frame of all segments as optimization parameters.
        /// </summary>
        /// <param name="specs"></param>
        /// <param name="N">Desired number of samples. This is not explicitly respected in this method. 
        /// If you want an exact number of samples use Build()</param>
        /// <param name="fastMode"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static IterationContext Build_IncludeStartAndEnd(IReadOnlyList<CurveSpec> specs, int N, bool fastMode)
        {
            if (specs == null || specs.Count == 0)
                throw new ArgumentException("IterationContext.Build: specs empty.");

            // Prefix & curves
            var curves = new IArcLengthCurve[specs.Count];
            var prefix = new double[specs.Count + 1];
            prefix[0] = 0.0;
            for (int i = 0; i < specs.Count; i++)
            {
                curves[i] = specs[i].GetOptimizedCurve();//new CachedIntrinsicCurve(specs[i], specs[i].Length / 20);
                prefix[i + 1] = prefix[i] + specs[i].Length;
            }
            double Ltot = prefix[specs.Count];
            // Num sample points, ensure the start and end of each segment
            // we need at least 2*spec.Count samples
            int M = Math.Max(specs.Count * 2, N);
            double ds = Ltot / M;

            // Now we start at the first segment
            // and take samples at increments of ds
            // until we surpass the segment length.
            // then we take an endpoint sample and
            // move to the next segment. 
            var sGlob = new List<double>();
            var segIdx = new List<int>();
            var sLocal = new List<double>();
            var S = new List<Sample>();
            double sglob;
            for (int segCount = 0; segCount < curves.Length; segCount++)
            {
                var seg = curves[segCount];
                var s = 0.0;
                sglob = prefix[segCount];
                var segL = seg.Length;
                while (s < segL)
                {
                    sGlob.Add(sglob);
                    sLocal.Add(s);
                    segIdx.Add(segCount);
                    S.Add(seg.Evaluate(s));

                    s += ds;
                    sglob += ds;
                }

                // include segment endpoint
                sglob = prefix[segCount + 1];
                sGlob.Add(sglob);
                sLocal.Add(segL);
                segIdx.Add(segCount);
                S.Add(seg.Evaluate(segL));
            }

            return new IterationContext(specs, curves, prefix, Ltot, Ltot / (M - 1), fastMode, sGlob.ToArray(), segIdx.ToArray(), sLocal.ToArray(), S.ToArray())
            { 
                Mode = CompositeCurveSolver.CompositeSolverMode.Seperate
            };
        }

        /// <summary>
        /// Builds a G1FullFrame <see cref="CompositeCurve"/> and captures at least N samples.
        /// A sample at the start of each segment and at the final endpoint is guaranteed.
        /// The solver will include the position and frame of only the first segment.
        /// </summary>
        /// <param name="specs"></param>
        /// <param name="N"></param>
        /// <param name="fastMode"></param>
        /// <returns></returns>
        public static IterationContext Build_Chained_IncludeStartAndEnd(IReadOnlyList<CurveSpec> specs, int N, bool fastMode)
        {
            if (specs == null || specs.Count == 0)
                throw new ArgumentException("IterationContext.Build: specs empty.");

            // Prefix & curves
            var curves = new IArcLengthCurve[specs.Count];
            var prefix = new double[specs.Count + 1];
            prefix[0] = 0.0;
            CompositeCurve c = new CompositeCurve();
            //Console.WriteLine("Building Chained Context");
            for (int i = 0; i < specs.Count; i++)
            {
                curves[i] = specs[i].GetOptimizedCurve();// new CachedIntrinsicCurve(specs[i], specs[i].Length / 20);
                prefix[i + 1] = prefix[i] + specs[i].Length;
                c.AddG1FullFrame(curves[i], out _);
                //specs[i].ShowInfo();
            }
            double Ltot = prefix[specs.Count];
            // Num sample points, ensure the endpoint of each segment
            // we need at least spec.Count + 1 samples
            int M = Math.Max(specs.Count + 1, N);
            double ds = Ltot / M;
            double ds2 = ds / 2;
            var sGlob = new List<double>();
            var segIdx = new List<int>();
            var sLocal = new List<double>();
            var S = new List<Sample>();
            double sglob;
            for (int segCount = 0; segCount < curves.Length; segCount++)
            {
                var seg = curves[segCount];
                var s = 0.0;
                sglob = prefix[segCount];
                var segL = seg.Length;
                // if the distance from our sample to the endpoint of this segment
                // is greater than ds / 2
                while (segL - s > ds2)
                {
                    sGlob.Add(sglob);
                    sLocal.Add(s);
                    segIdx.Add(segCount);
                    S.Add(c.Evaluate(sglob));

                    s += ds;
                    sglob += ds;
                }
            }

            // include final sample
            var sl = curves[^1].Length;
            sGlob.Add(Ltot);
            sLocal.Add(sl);
            segIdx.Add(curves.Length - 1);
            S.Add(c.Evaluate(Ltot));

            return new IterationContext(specs, curves, prefix, Ltot, Ltot / (M - 1), fastMode, sGlob.ToArray(), segIdx.ToArray(), sLocal.ToArray(), S.ToArray()) 
            { 
                Mode = CompositeCurveSolver.CompositeSolverMode.Chained 
            };
        }

        /// <summary>
        /// Dispose of any disposables here
        /// </summary>
        public void Dispose()
        {
            // if we use Array.Pool or something to optimize even more we might need to dispose of them.
        }
    }
}
