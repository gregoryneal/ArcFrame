# ArcFrame

ArcFrame is a C# library for generating, evaluating, and composing arc-length-parameterized curves in 2D, 3D, and higher dimensions. It centers on intrinsic descriptions of curves (curvature as a function of arc length) and integrates moving frames on SO(N) to recover position and orientation with stable, adaptive numerics. The library includes robust clothoid utilities based on the Bertolazzi–Frego formulations for G1 interpolation and extensions toward G2 spline construction.

## What it provides

### Unified curve model
IArcLengthCurve with Evaluate(s) → (P, R, k) where P is position, R is an orthonormal frame [T, N, …], and k are principal curvatures. This contract ensures the curve can be
embedded in any spatial dimension.

### Intrinsic curves
IntrinsicCurve builds curves from a CurveSpec that couples initial pose (ND), frame model (Frenet or Bishop), and a curvature law (N-1 components).
<img width="1432" height="683" alt="IntrinsicCurveDemo" src="https://github.com/user-attachments/assets/656fccb1-cc99-485f-9d67-71782e9cf95a" />

### Lie-group frame stepping
A midpoint integrator on SO(N) with adaptive step control for stable propagation of pose over arc length.
TODO: implement RK4 frame stepping.

### Constraint Solver
Uses the Levenberg-Marquardt algorithm to solve a CompositeCurveProblem given a set of constraints.
<img width="1396" height="613" alt="Example_BuildBoundedClothoid" src="https://github.com/user-attachments/assets/cc1f612a-b41e-4c0e-8a20-dca5f23ac01a" />
<img width="640" height="700" alt="3ArcClothoidSolver3D" src="https://github.com/user-attachments/assets/5be5c4ea-1f93-478c-9211-c8ab66c4c9d0" />

### Clothoids and curvature laws
Linear curvature laws and generalized Fresnel integral evaluators support accurate clothoid construction and sampling.

### Native Spline model
A generalized matrix based Spline object, works as a native IArcLengthCurve by caching an ArcLengthTable on creation. Works with any matrix based spline model.

Position found via right-multiplication: P(t) = GBt' where G is the control point window, B is the basis matrix and t' is the power monomial vector [1, t, t^2, ...]^T

Built in cubic B-Splines include: CatmullRomSpline, BezierSpline, BSpline, HermiteSpline. Native support for any degree: linear, quadratic, cubic, quartic, and beyond.
<img width="1455" height="766" alt="CatmullRomSpline" src="https://github.com/user-attachments/assets/ce37cbc6-f717-47b0-926f-ba606cbb9cc6" />

### Composition and transforms
CompositeCurve concatenates segments with automatic G1 alignment. TransformedCurve, RigidTransform, and PromotedCurve embed and position curves across dimensions.

### Performance helpers
CachedIntrinsicCurve checkpoints along arc length for fast re-evaluation and uniform sampling.

## Use cases

Path planning and motion primitives for robotics and vehicles

Track, road, and rail geometry design

Game and simulation tooling that needs smooth, controllable splines

CAD/CAM workflows that benefit from curvature-aware construction
