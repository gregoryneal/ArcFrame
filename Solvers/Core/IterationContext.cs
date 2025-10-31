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
        /// Specs
        /// </summary>
        public IReadOnlyList<CurveSpec> Specs { get; }
        /// <summary>
        /// Curves
        /// </summary>
        public IArcLengthCurve[] Curves { get; }

        /// <summary>
        /// The prefix of the curve
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
        /// </summary>
        /// <param name="specs"></param>
        /// <param name="ds"></param>
        /// <param name="fastMode"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static IterationContext Build(IReadOnlyList<CurveSpec> specs, double ds, bool fastMode)
        {
            if (specs == null || specs.Count == 0)
                throw new ArgumentException("IterationContext.Build: specs empty.");

            // Prefix & curves
            var curves = new IArcLengthCurve[specs.Count];
            var prefix = new double[specs.Count + 1];
            prefix[0] = 0.0;
            for (int i = 0; i < specs.Count; i++)
            {
                curves[i] = specs[i].GetOptimizedCurve();
                prefix[i + 1] = prefix[i] + specs[i].Length;
            }
            double Ltot = prefix[specs.Count];
            // Num sample points
            int M = Math.Max(2, (int)Math.Ceiling(Ltot / Math.Max(1e-6, ds)) + 1);

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
                // local s
                double sl = sg - prefix[seg];

                var smp = curves[seg].Evaluate(sl); // single expensive call per sample
                // cache the global and local s, the segment index and the Sample itself
                sGlob[m] = sg;
                segIdx[m] = seg;
                sLocal[m] = sl;
                S[m] = smp;
            }

            return new IterationContext(specs, curves, prefix, Ltot, ds, fastMode, sGlob, segIdx, sLocal, S);
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
