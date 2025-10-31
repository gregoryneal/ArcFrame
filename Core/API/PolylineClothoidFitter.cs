﻿using ArcFrame.Core.Geometry;
using ArcFrame.Core.Params;
using ArcFrame.Core.Results;
using ArcFrame.Solvers;
using System;
using System.Collections.Generic;

namespace ArcFrame.Core.API
{
    /// <summary>
    /// A class that attempts to fit a clothoid or clothoid segments to a set of points.
    /// Note: This is a redundant special case of the curve constraint solvers.
    /// </summary>
    public static class PolylineClothoidFitter
    {
        /// <summary>
        /// Fit a polyline (x,z, theta)[] with N points into N-1 G1 continuous clothoid segments (planar, τ≡0),
        /// lifted to R^n (n>=2). Uses a paper-specific 2D fitter & evaluator that you pass in.
        /// </summary>
        /// <param name="points">Polyline in XZ (Y assumed 0).</param>
        /// <param name="fit2D">
        ///   Delegate: given two endpoint states (x,z,theta)_0 and (x, z, theta)_1 return a CurveSpec 
        /// </param>
        /// <param name="closed">If true, connects last→first (N segments).</param>
        public static FitResult<CompositeCurve> FitPolyline_G1_XY(IReadOnlyList<(double x, double y, double t)> points, Func<double, double, double, double, double, double, CurveSpec> fit2D, IEvaluator evaluator, bool closed = false)
        {
            if (points == null || points.Count < 2) return FitResult<CompositeCurve>.Fail("Not enough points in polyline");

            try
            {
                CompositeCurve c = new CompositeCurve();

                var comp = new CompositeCurve();

                int segCount = closed ? points.Count : points.Count - 1;

                for (int i = 0; i < segCount; i++)
                {
                    int j = (i + 1) % points.Count;
                    var A = (points[i].x, points[i].y, points[i].t);
                    var B = (points[j].x, points[j].y, points[j].t);

                    // 3) Fit a 2D clothoid between endpoints A→B (paper specific)
                    var cparams = fit2D(A.x, A.y, A.t, B.x, B.y, B.t);

                    // Validate (light)
                    if (double.IsNaN(cparams.Length) || cparams.Length <= 0)
                    {
                        // Fallback: straight line (ND)
                        var line = Line.From2D_XZ(A.x, A.y, B.x, B.y);
                        comp.Add(line);
                        continue;
                    }

                    // 5) Wrap with 2D evaluator adapter
                    var seg = new Clothoid(cparams, evaluator);

                    // 6) Add to composite 
                    comp.Add(seg);
                }

                return FitResult<CompositeCurve>.Success(comp);
            }
            catch (Exception ex)
            {
                return FitResult<CompositeCurve>.Fail("Polyline fit failed: " + ex.Message);
            }
        }

        /// <summary>
        /// Build a CompositeCurve out of a set of clothoid with associated endpoints, tangend and curvatures.
        /// </summary>
        /// <param name="points"></param>
        /// <param name="fit2D"></param>
        /// <param name="evaluator"></param>
        /// <param name="closed"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static FitResult<CompositeCurve> FitPolyline_G2_XY(IReadOnlyList<(double x, double y, double t, double k)> points, Func<double, double, double, double, double, double, double, double, ClothoidCurveSpec> fit2D, IEvaluator evaluator, bool closed, int n)
        {
            if (points == null || points.Count < 2) return FitResult<CompositeCurve>.Fail("Not enough points in polyline");
            if (n < 2) n = 2;

            try
            {
                CompositeCurve c = new CompositeCurve();

                var comp = new CompositeCurve();

                int segCount = closed ? points.Count : points.Count - 1;

                for (int i = 0; i < segCount; i++)
                {
                    int j = (i + 1) % points.Count;
                    var A = (points[i].x, points[i].y, points[i].t, points[i].k);
                    var B = (points[j].x, points[j].y, points[j].t, points[j].k);

                    // 3) Fit a 2D clothoid between endpoints A→B (paper-specific)
                    var cparams = fit2D(A.x, A.y, A.t, A.k, B.x, B.y, B.t, B.k);

                    // Validate (light)
                    if (double.IsNaN(cparams.Length) || cparams.Length <= 0)
                    {
                        // Fallback: straight line (ND)
                        var line = Line.From2D_XZ(A.x, A.y, B.x, B.y);
                        comp.Add(line);
                        continue;
                    }

                    // 5) Wrap with 2D evaluator adapter → ND curve
                    var seg = new Clothoid(cparams, evaluator);

                    comp.Add(seg);
                }

                return FitResult<CompositeCurve>.Success(comp);
            }
            catch (Exception ex)
            {
                return FitResult<CompositeCurve>.Fail("Polyline fit failed: " + ex.Message);
            }
        }
    }
}
