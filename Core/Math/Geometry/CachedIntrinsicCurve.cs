using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using ArcFrame.Core.Results;
using System;
using System.Collections.Generic;

namespace ArcFrame.Core.Geometry
{
    /// <summary>
    /// A version of IntrinsicCurve that caches computed points along the curve.
    /// </summary>
    public sealed class CachedIntrinsicCurve : IArcLengthCurve
    {
        private readonly IntrinsicCurve _inner;
        private readonly IFrameStepper _stepper;
        private readonly IntegratorOptions _opt;
        private readonly CurveSpec _spec;
        private readonly double _ds; // target checkpoint spacing

        // checkpoint arrays
        private readonly List<double> _S = new List<double>();
        private readonly List<double[]> _P = new List<double[]>();
        private readonly List<double[,]> _R = new List<double[,]>();

        /// <summary>
        /// Create a cached intrinsic curve with optional integrator.
        /// </summary>
        /// <param name="spec"></param>
        /// <param name="checkpointSpacing"></param>
        /// <param name="stepper"></param>
        /// <param name="opt"></param>
        public CachedIntrinsicCurve(CurveSpec spec, double checkpointSpacing = 0.25, IFrameStepper? stepper = null, IntegratorOptions? opt = null)
        {
            _spec = spec;
            _inner = new IntrinsicCurve(spec, stepper, opt);
            _stepper = stepper ?? new LieGroupMidpointStepper();
            _opt = opt ?? IntegratorOptions.Default;
            _ds = System.Math.Max(1e-6, checkpointSpacing);

            // seed s=0
            _S.Add(0.0);
            _P.Add((double[])spec.P0.Clone());
            _R.Add((double[,])spec.R0.Clone());
        }
        /// <inheritdoc/>
        public int Dimension => _inner.Dimension;
        /// <inheritdoc/>
        public double Length => _inner.Length;
        /// <inheritdoc/>
        public Sample Evaluate(double s)
        {
            s = System.Math.Clamp(s, 0.0, Length);

            // ensure checkpoints up to s
            EnsureCheckpointsUpTo(s);

            // find nearest previous checkpoint
            int idx = UpperBound(_S, s) - 1; // _S[idx] <= s < _S[idx+1]
            var P = (double[])_P[idx].Clone();
            var R = (double[,])_R[idx].Clone();

            // short hop from _S[idx] to s
            double s0 = _S[idx];
            while (s0 < s - 0)
            {
                double h = _opt.SuggestStep(_spec.Kappa.Eval, s0, s);
                if (s0 + h > s) h = s - s0;
                _stepper.Step(ref P, ref R, _spec.Kappa.Eval, s0, h, _spec.Frame);
                s0 += h;
            }

            var k = _spec.Kappa.Eval(s);
            return new Sample(P, R, s, k);
        }

        /// <summary>
        /// Check if we have cached a Sample beyond arc length s. If not we sample from the last cached point up to s.
        /// </summary>
        /// <param name="s"></param>
        private void EnsureCheckpointsUpTo(double s)
        {
            double sLast = _S[^1];
            if (s <= sLast) return;

            var P = (double[])_P[^1].Clone();
            var R = (double[,])_R[^1].Clone();

            while (sLast < s - 0)
            {
                double sNext = System.Math.Min(Length, sLast + _ds);
                // integrate from sLast to sNext with adaptive steps
                double s0 = sLast;
                while (s0 < sNext)
                {
                    double h = _opt.SuggestStep(_spec.Kappa.Eval, s0, sNext);
                    if (s0 + h > sNext) h = sNext - s0;
                    _stepper.Step(ref P, ref R, _spec.Kappa.Eval, s0, h, _spec.Frame);
                    s0 += h;
                }

                _S.Add(sNext);
                _P.Add((double[])P.Clone());
                _R.Add((double[,])R.Clone());
                sLast = sNext;
            }
        }

        static int UpperBound(List<double> a, double x)
        {
            int lo = 0, hi = a.Count;
            while (lo < hi) 
            { 
                int mid = (lo + hi) >> 1;
                if (a[mid] <= x) lo = mid + 1; 
                else hi = mid; 
            }
            return lo;
        }
        /// <inheritdoc/>
        public double[] Position(double s) => Evaluate(s).P;
        /// <inheritdoc/>
        public double[] Tangent(double s) => Evaluate(s).T;

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
    }

}
