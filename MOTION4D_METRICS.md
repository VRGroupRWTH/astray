# Motion4D Metrics Implementation Status

This document tracks the implementation status of metrics from the Motion4D library that were not present in Astray.

## Reference
Original Motion4D library by Thomas Mueller (tauzero7): https://github.com/tauzero7/Motion4D

## Implementation Summary

### Completed Metrics (5/41 = 12%)

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

### Remaining Metrics by Priority

#### High Priority (6 metrics)
These are well-known, important metrics that should be implemented next:

- **Kruskal** (spherical) - Maximal analytic extension of Schwarzschild
- **TaubNUT** (spherical) - Rotating solution with NUT parameter
- **DeSitterUnivConf** (cartesian) - Conformal de Sitter universe
- **GoedelCart** (cartesian) - Gödel universe in Cartesian coordinates
- **GoedelScaled** (cartesian) - Scaled Gödel universe
- **GoedelScaledCart** (cartesian) - Scaled Gödel in Cartesian coordinates

#### Medium Priority (15 metrics)
More specialized but still useful:

- **AlcubierreSimple** (cartesian) - Simplified Alcubierre warp drive
- **PlaneGravWave** (cartesian) - Plane gravitational wave (complex: needs Fourier series)
- **Curzon** (cylindrical) - Static axisymmetric solution
- **ChazyCurzonRot** (cylindrical) - Rotating Chazy-Curzon
- **ErezRosenVar** (cylindrical) - Erez-Rosen with variable parameters
- **Ernst** (cylindrical) - Ernst metric
- **ErnstSchwarzschild** (cylindrical) - Ernst form of Schwarzschild
- **SchwarzschildTortoise** (spherical) - Tortoise coordinate form
- **SchwarzschildWT** (spherical) - Wheeler-Thorne coordinates
- **SchwarzschildGravWave** (spherical) - Schwarzschild with gravitational wave
- **RotDihole** (spherical) - Rotating dihole solution
- **HalilsoyWave** (spherical) - Halilsoy wave metric
- **SultanaDyer** (spherical) - Sultana-Dyer metric

#### Lower Priority (20 metrics)
Highly specialized metrics:

- **EddFinkIn** (spherical) - Eddington-Finkelstein ingoing
- **Glampedakis** (spherical) - Glampedakis metric
- **HartleThorneGB** (spherical) - Hartle-Thorne with Gauss-Bonnet
- **TeoSimpleWH** (spherical) - Teo simple wormhole
- **TeoWHl** (spherical) - Teo wormhole with parameter l
- **TomimatsuSato** (cylindrical) - Tomimatsu-Sato metric
- **MinkRotLattice** (cartesian) - Minkowski with rotating lattice
- **StraightSpinningString** (cylindrical) - Straight spinning string
- **PTD_AI, PTD_AII, PTD_AIII** - Plebanski-Demianski type A metrics
- **PTD_BI, PTD_BII, PTD_BIII** - Plebanski-Demianski type B metrics
- **PTD_C** - Plebanski-Demianski type C
- **Pravda_C**, **Pravda_C_Can** - Pravda metrics

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
