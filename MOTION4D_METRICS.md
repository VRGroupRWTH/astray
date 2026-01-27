# Motion4D Metrics Implementation Status

This document tracks the implementation status of metrics from the Motion4D library that were not present in Astray.

## Reference
Original Motion4D library by Thomas Mueller (tauzero7): https://github.com/tauzero7/Motion4D

## Implementation Summary

### Completed Metrics (13/41 = 32%)

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

13. (Reserved for next implementation)

### Remaining Metrics by Priority

#### High Priority - Can be implemented without complex special functions (15 metrics)
These are straightforward implementations that don't require Lambert W, Fourier series, or other complex functions:

- **ChazyCurzonRot** (cylindrical) - Rotating Chazy-Curzon solution
- **ErezRosenVar** (cylindrical) - Erez-Rosen with variable parameters
- **Ernst** (cylindrical) - Ernst metric
- **ErnstSchwarzschild** (cylindrical) - Ernst form of Schwarzschild
- **StraightSpinningString** (cylindrical) - Straight spinning string
- **RotDihole** (spherical) - Rotating dihole solution
- **HalilsoyWave** (spherical) - Halilsoy wave metric
- **SultanaDyer** (spherical) - Sultana-Dyer metric
- **Glampedakis** (spherical) - Glampedakis metric
- **HartleThorneGB** (spherical) - Hartle-Thorne with Gauss-Bonnet
- **TeoSimpleWH** (spherical) - Teo simple wormhole
- **TeoWHl** (spherical) - Teo wormhole with parameter l
- **PTD_AI, PTD_AII, PTD_AIII** - Plebanski-Demianski type A metrics (3 metrics)
- **PTD_BI, PTD_BII, PTD_BIII** - Plebanski-Demianski type B metrics (3 metrics)
- **PTD_C** - Plebanski-Demianski type C

#### Medium Priority - Require special functions or complex calculations (13 metrics)
These require Lambert W function, Fourier series, or other advanced mathematical functions:

- **Kruskal** (spherical) - Maximal analytic extension of Schwarzschild
  - **Complexity**: Requires GSL Lambert W function (gsl_sf_lambert_W0)
  - Coordinates cover full Schwarzschild spacetime including both exterior and interior regions

- **TaubNUT** (spherical) - Rotating solution with NUT parameter
  - **Complexity**: Very complex Christoffel symbols with extensive polynomial expressions
  - Important for studying gravitomagnetic monopoles

- **SchwarzschildTortoise** (spherical) - Tortoise coordinate form
  - **Complexity**: Requires Lambert W function for coordinate transformation
  - Useful for null geodesics and wave propagation

- **SchwarzschildWT** (spherical) - Wheeler-Thorne coordinates
  - **Complexity**: May require special coordinate transformations

- **SchwarzschildGravWave** (spherical) - Schwarzschild with gravitational wave
  - **Complexity**: Combination of Schwarzschild and wave perturbations

- **AlcubierreSimple** (cartesian) - Simplified Alcubierre warp drive
  - **Complexity**: Requires derivative calculations for shape function f(rs)
  - rs = √((x - vst)² + y² + z²) with specific shape function

- **PlaneGravWave** (cartesian) - Plane gravitational wave (sandwich wave)
  - **Complexity**: Requires Fourier series expansion for wave functions p(u) and q(u)
  - Most complex implementation, needs numerical integration

- **TomimatsuSato** (cylindrical) - Tomimatsu-Sato metric
  - **Complexity**: Complex coordinate-dependent expressions

- **Pravda_C** (various) - Pravda metric type C
  - **Complexity**: Specialized algebraic type

- **Pravda_C_Can** (various) - Pravda C in canonical form
  - **Complexity**: Canonical coordinate form

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

To complete the remaining 36 metrics:
1. Implement high-priority metrics first (Kruskal, TaubNUT, etc.)
2. Handle complex metrics requiring numerical methods separately
3. Consider adding unit tests for Christoffel symbol correctness
4. Add validation against known geodesics where available
5. Document physical significance and use cases for each metric

## References

- Motion4D Library: https://github.com/tauzero7/Motion4D
- Motion4D Documentation: Available in the repository
- Original papers cited in Motion4D header files
