# Motion4D Metrics Implementation Status

This document tracks the implementation status of metrics from the Motion4D library that were not present in Astray.

## Reference
Original Motion4D library by Thomas Mueller (tauzero7): https://github.com/tauzero7/Motion4D

## Implementation Summary

### Completed Metrics (33/41 = 80%)

1. **Kasner** (cartesian) - `include/astray/metrics/cartesian/kasner.hpp`
   - Anisotropic cosmological solution
   - Uses Khalatnikov-Lifshitz parameter u
   - Parameters: p1, p2, p3 satisfying p1 + p2 + p3 = 1 and p1² + p2² + p3² = 1

2. **PainlevéGullstrand** (spherical) - `include/astray/metrics/spherical/painleve_gullstrand.hpp`
   - Schwarzschild metric in freely falling observer coordinates
   - Regular at the Schwarzschild radius (no coordinate singularity)
   - Useful for studying infalling geodesics

3. **MinkowskiConformal** (spherical) - `include/astray/metrics/spherical/minkowski_conformal.hpp`
   - Minkowski spacetime in conformally compactified coordinates (ψ, ξ, θ, φ)
   - Useful for studying causal structure
   - Entire spacetime fits into a finite region

4. **SchwarzschildIsotropic** (cartesian) - `include/astray/metrics/cartesian/schwarzschild_isotropic.hpp`
   - Schwarzschild in isotropic Cartesian coordinates (t, x, y, z)
   - Spatial part is conformally flat
   - Uses coordinate ρ = √(x² + y² + z²) with ρs = rs/(2c²)

5. **VaidyaIncomingRadiation** (spherical) - `include/astray/metrics/spherical/vaidya_incoming_radiation.hpp`
   - Time-varying black hole mass due to incoming radiation
   - Mass function: m(v) = k(vl - v)^(1/3) tanh²(σ(vl - v))
   - Important for studying black hole formation

6. **GoedelCart** (cartesian) - `include/astray/metrics/cartesian/goedel_cart.hpp`
   - Gödel universe in Cartesian coordinates with cylindrical symmetry
   - Rotating cosmological solution allowing closed timelike curves
   - Parameters: a (Gödel radius rG = 2a), zeta (rotation parameter)

7. **GoedelScaled** (cylindrical) - `include/astray/metrics/cylindrical/goedel_scaled.hpp`
   - Gödel universe in scaled cylindrical coordinates
   - Geodesic shape independent of rG parameter

8. **GoedelScaledCart** (cartesian) - `include/astray/metrics/cartesian/goedel_scaled_cart.hpp`
   - Gödel universe in scaled Cartesian coordinates
   - Valid at origin (no coordinate singularity)

9. **DeSitterUnivConf** (cartesian) - `include/astray/metrics/cartesian/de_sitter_univ_conf.hpp`
   - de Sitter universe in conformal Cartesian coordinates
   - Expanding universe with cosmological constant

10. **EddingtonFinkelsteinIngoing** (spherical) - `include/astray/metrics/spherical/eddington_finkelstein_ingoing.hpp`
    - Schwarzschild in Eddington-Finkelstein ingoing coordinates
    - Regular at event horizon, describes infalling geodesics

11. **MinkowskiRotatingLattice** (cylindrical) - `include/astray/metrics/cylindrical/minkowski_rotating_lattice.hpp`
    - Minkowski spacetime in rotating cylindrical coordinates
    - Angular velocity ω, valid for r < c/ω

12. **Curzon** (cylindrical) - `include/astray/metrics/cylindrical/curzon.hpp`
    - Static axisymmetric solution for infinitely long line mass
    - Exterior field of line mass along z-axis

13. **Ernst** (spherical) - `include/astray/metrics/spherical/ernst.hpp`
    - Black hole in a magnetic universe
    - Parameters: mass, magnetic field B
    - References: Ernst (1976), Konoplya (2007)

14. **ErnstSchwarzschild** (spherical) - `include/astray/metrics/spherical/ernst_schwarzschild.hpp`
    - Ernst form of Schwarzschild metric in magnetic universe
    - Similar to Ernst but with different metric signature

15. **RotDihole** (cartesian) - `include/astray/metrics/cartesian/rot_dihole.hpp`
    - Rotating extreme Reissner-Nordstrom dihole
    - Two rotating holes with angular velocity ω
    - Parameters: mass1, mass2, angular_velocity

16. **SultanaDyer** (spherical) - `include/astray/metrics/spherical/sultana_dyer.hpp`
    - Cosmological black hole in Einstein-de Sitter universe
    - Time-dependent metric with t⁴ factor
    - Reference: Sultana & Dyer (2005)

17. **PTD_AI** (spherical) - `include/astray/metrics/spherical/ptd_ai.hpp`
    - Plebanski-Demianski Type D metric, subclass A I
    - Exact solution of Petrov Type D
    - Parameter: b

18. **PTD_AII** (cylindrical) - `include/astray/metrics/cylindrical/ptd_aii.hpp`
    - Plebanski-Demianski Type D metric, subclass A II
    - Uses hyperbolic coordinates in r direction
    - Parameter: b

19. **PTD_AIII** (cylindrical) - `include/astray/metrics/cylindrical/ptd_aiii.hpp`
    - Plebanski-Demianski Type D metric, subclass A III
    - Simpler form without parameter b

20. **PTD_BI** (spherical) - `include/astray/metrics/spherical/ptd_bi.hpp`
    - Plebanski-Demianski Type D metric, subclass B I
    - Parameter: b

21. **PTD_BII** (cylindrical) - `include/astray/metrics/cylindrical/ptd_bii.hpp`
    - Plebanski-Demianski Type D metric, subclass B II
    - Uses hyperbolic coordinates
    - Parameter: b

22. **PTD_BIII** (cylindrical) - `include/astray/metrics/cylindrical/ptd_biii.hpp`
    - Plebanski-Demianski Type D metric, subclass B III
    - Simpler form without parameter b

23. **PTD_C** (cartesian) - `include/astray/metrics/cartesian/ptd_c.hpp`
    - Plebanski-Demianski Type D metric, Type C
    - Uses custom coordinates (t, u, x, y)
    - Parameters: a, b

24. **ChazyCurzonRot** (cylindrical) - `include/astray/metrics/cylindrical/chazy_curzon_rot.hpp`
    - Rotating generalization of Curzon metric
    - Parameters: mass, p, q (where p² + q² = 1)

25. **StraightSpinningString** (cylindrical) - `include/astray/metrics/cylindrical/straight_spinning_string.hpp`
    - Infinite spinning cosmic string metric
    - Parameters: a (string parameter), k (rotation parameter)

26. **AlcubierreSimple** (cartesian) - `include/astray/metrics/cartesian/alcubierre_simple.hpp`
    - Simplified Alcubierre warp drive with smooth bubble transitions
    - Parameters: R (bubble size), dR (transition width), vs (bubble velocity)

27. **TeoWHl** (spherical) - `include/astray/metrics/spherical/teo_wormhole.hpp`
    - Rotating traversable wormhole metric with customizable potentials
    - Parameters: b0 (throat radius), custom functions N(l), K(l), r(l), ω(l)

28. **TaubNUT** (spherical) - `include/astray/metrics/spherical/taub_nut.hpp`
    - Taub-NUT metric with gravitomagnetic monopole charge
    - Parameters: mass M, NUT parameter l
    - Reference: Bini et al, Class. Quantum Grav. 19, 5481 (2002)

29. **ErezRosenVar** (spherical) - `include/astray/metrics/spherical/erez_rosen_var.hpp`
    - Static axisymmetric solution with variable deformation parameter
    - Parameters: mass M, deformation parameter q
    - Uses potentials ψ, γ, and Δ

30. **TeoSimpleWH** (spherical) - `include/astray/metrics/spherical/teo_simple_wormhole.hpp`
    - Simple rotating wormhole with throat at l=0
    - Parameter: b0 (throat radius)

31. **Pravda_C** (cartesian) - `include/astray/metrics/cartesian/pravda_c.hpp`
    - C-metric representing uniformly accelerating black holes
    - Parameters: acceleration A, mass m

27. **TeoWHl** (spherical) - `include/astray/metrics/spherical/teo_wormhole.hpp`
    - Teo wormhole metric with parameter l (rotating traversable wormhole)
    - Uses potentials N(l), K(l), r(l), ω(l) that can be customized
    - Standard potentials: N=1, K=1, r=√(l²+b₀²), ω=b₀²/(2(l²+b₀²)^(3/2))
    - Reference: Edward Teo, Phys. Rev. D 58, 024014 (1998)

28. **TaubNUT** (spherical) - `include/astray/metrics/spherical/taub_nut.hpp`
    - Taub-NUT metric with gravitomagnetic monopole (NUT parameter)
    - Boyer-Lindquist like spherical coordinates
    - Parameters: mass, l (NUT parameter)
    - Reference: Bini et al, Class. Quantum Grav. 19, 5481 (2002)

29. **ErezRosenVar** (spherical) - `include/astray/metrics/spherical/erez_rosen_var.hpp`
    - Erez-Rosen metric with variable deformation parameter
    - Static axisymmetric solution with deformed mass distribution
    - Parameters: mass, q (deformation parameter)

30. **TeoSimpleWH** (spherical) - `include/astray/metrics/spherical/teo_simple_wormhole.hpp`
    - Teo simple wormhole metric (axisymmetric rotating wormhole)
    - Throat at l=0 with radius parameter b₀
    - Reference: Edward Teo, Phys. Rev. D 58, 024014 (1998)

31. **Pravda_C** (cartesian) - `include/astray/metrics/cartesian/pravda_c.hpp`
    - Pravda C-metric representing accelerating black holes
    - Custom coordinates (t, x, y, phi)
    - Parameters: acceleration (A), mass (m)

32. **SchwarzschildWT** (spherical) - `include/astray/metrics/spherical/schwarzschild_wt.hpp`
    - Schwarzschild in Wheeler-Thorne coordinates
    - Time component g_tt is constant (-c²)
    - Useful for theoretical studies of black hole physics

33. **SchwarzschildGravWave** (spherical) - `include/astray/metrics/spherical/schwarzschild_gravwave.hpp`
    - Schwarzschild with gravitational wave perturbations
    - Perturbations depend on Legendre polynomials (supports l=0 and l=1)
    - Parameters: epsilon (amplitude), sigma (frequency), l (mode number)

### Remaining Metrics by Priority

#### High Priority - Can be implemented without complex special functions (0 remaining)
All high-priority metrics have been implemented!

#### Medium Priority - Require special functions or complex calculations (8 metrics remaining)
These require Lambert W function, Fourier series, Bessel functions, or other advanced mathematical functions:

- **Kruskal** (spherical) - Maximal analytic extension of Schwarzschild
  - **Status**: NOT IMPLEMENTED - Requires GSL Lambert W function (gsl_sf_lambert_W0)
  - Coordinates cover full Schwarzschild spacetime including both exterior and interior regions

- **SchwarzschildTortoise** (spherical) - Tortoise coordinate form
  - **Status**: NOT IMPLEMENTED - Requires Lambert W function for coordinate transformation
  - Useful for null geodesics and wave propagation

- **HalilsoyWave** (cylindrical) - Halilsoy wave metric
  - **Status**: NOT IMPLEMENTED - Requires GSL Bessel functions (J0, J1)

- **PlaneGravWave** (cartesian) - Plane gravitational wave (sandwich wave)
  - **Status**: NOT IMPLEMENTED - Requires Fourier series expansion for wave functions p(u) and q(u)
  - Most complex implementation, needs numerical integration

- **HartleThorneGB** (spherical) - Hartle-Thorne with Gauss-Bonnet
  - **Status**: NOT IMPLEMENTED - Extremely complex (1212 lines) with extensive helper functions

- **Glampedakis** (spherical) - Glampedakis metric
  - **Status**: NOT IMPLEMENTED - Requires calcKerr and calcGlampedakis helper functions with extensive calculations
  - Important for LISA gravitational wave detection

- **TomimatsuSato** (cylindrical) - Tomimatsu-Sato metric
  - **Status**: NOT IMPLEMENTED - Very complex (783 lines) with extensive helper functions and coordinate-dependent expressions
  - Reference: V.S. Manko, Progress of Theoretical Physics 127, 1057 (2012)

- **Pravda_C_Can** (cylindrical) - Pravda C in canonical coordinates
  - **Status**: NOT IMPLEMENTED - Complex (713 lines) with root-finding and auxiliary function calculations


## Implementation Guidelines

### Structure
Each metric should:
1. Be placed in the appropriate coordinate system directory (cartesian, spherical, cylindrical, boyer_lindquist)
2. Follow the template pattern established in existing metrics
3. Include the header comment citing Motion4D
4. Implement `christoffel_symbols()` method
5. Optionally implement `check_termination()` for singularities
6. Use `__device__` qualifier for CUDA compatibility

### Citation Format
```cpp
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
```

### Template Pattern
```cpp
#pragma once

#include <cmath>
#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// [Metric description]
// Line element: [LaTeX or text description]
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class metric_name : public metric<coordinate_system_type::TYPE, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, const vector_type& direction) const override
  {
    // Check for singularities
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    // Extract coordinates
    // Compute Christoffel symbols
    christoffel_symbols_type symbols;
    symbols.setZero();
    // Set non-zero components
    return symbols;
  }

  // Parameters as public members with sensible defaults
  scalar_type param = static_cast<scalar_type>(1);
};
}
```

### Fetching Reference Implementation
1. Browse to: https://github.com/tauzero7/Motion4D/tree/master/src/metric
2. Find m4dMetric[Name].cpp and m4dMetric[Name].h
3. Extract the Christoffel symbol calculations from `calculateChristoffels()`
4. Translate from Motion4D's array indexing to Astray's tensor indexing
5. Replace C-style arrays with Eigen operations where appropriate

### Common Challenges

1. **Complex Metrics**: Some metrics (PlaneGravWave, Glampedakis) require numerical integration or special functions
2. **Coordinate Transformations**: Some metrics need coordinate transformation helpers
3. **Parameter Calculations**: Some metrics have derived parameters that need recalculation
4. **Special Functions**: Bessel functions, elliptic integrals may be needed

### Testing
Currently, there is minimal test infrastructure. Focus on:
1. Ensuring the code compiles
2. Checking Christoffel symbols match Motion4D for simple test cases
3. Verifying termination conditions are appropriate

## Future Work

To complete the remaining 8 metrics:
1. Consider adding GSL or similar library as dependency to enable Lambert W and Bessel functions
2. Implement a specialized "special functions" library for commonly needed functions
3. Focus on simpler alternative metrics not in Motion4D that provide similar physics
4. Consider adding unit tests for Christoffel symbol correctness
5. Add validation against known geodesics where available
6. Document physical significance and use cases for each metric

## Detailed Analysis

See **FINAL_10_METRICS_REPORT.md** for a comprehensive analysis of the final 10 metrics, including:
- Detailed reasons why each metric was not implemented
- External dependencies required
- Code complexity analysis
- Recommendations for future work

## References

- Motion4D Library: https://github.com/tauzero7/Motion4D
- Motion4D Documentation: Available in the repository
- Original papers cited in Motion4D header files
