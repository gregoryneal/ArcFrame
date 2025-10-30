using ArcFrame.Core.Math;
using ArcFrame.Core.Params;
using ArcFrame.Core.Results;
using System;

namespace ArcFrame.Core.Geometry
{
    /// <summary>
    /// A curve described using its intrinsic parameters. These are generally integrated via frame-stepping to find points along its arc length.
    /// </summary>
    public class IntrinsicCurve : IArcLengthCurve
    {
        private readonly CurveSpec _spec;
        private readonly IFrameStepper _stepper;
        private readonly IntegratorOptions _opt;

        public IntrinsicCurve(CurveSpec spec, IFrameStepper? stepper = null, IntegratorOptions? opt = null)
        {
            _spec = spec;
            _opt = opt ?? IntegratorOptions.Default;
            _stepper = stepper ?? new LieGroupMidpointStepper();
        }

        public int Dimension => _spec.N;

        public double Length => _spec.Length;

        public Sample Evaluate(double s)
        {
            s = System.Math.Clamp(s, 0.0, _spec.Length);

            var P = (double[])_spec.P0.Clone();
            var R = (double[,])_spec.R0.Clone();

            double s0 = 0.0;
            while (s0 < s)
            {
                double h = _opt.SuggestStep(_spec.Kappa.Eval, s0, s);
                //Console.WriteLine($"Suggested step size: {h} ({s0},{s})");
                if (h <= 0) break;
                _stepper.Step(ref P, ref R, _spec.Kappa.Eval, s0, h, _spec.Frame);
                s0 += h;
            }

            // final sample
            var k = _spec.Kappa.Eval(s);
            return new Sample(P, R, s, k);
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
                //Console.WriteLine("New position sample:");
                //Helpers.PrintVector(samples[i].P);
            }
            return samples;
        }

        public double[] Position(double s) => Evaluate(s).P;

        public double[] Tangent(double s) => Evaluate(s).T;
    }
}
