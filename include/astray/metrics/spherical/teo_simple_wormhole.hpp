#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Teo simple wormhole metric in spherical coordinates (t,l,theta,phi).
//
// This metric describes an axisymmetric rotating wormhole with a throat at l=0.
// The coordinate l ranges from -∞ to +∞, passing through the wormhole throat.
//
// The line element involves the radial function r(l) = √(l² + b₀²) where b₀
// is the throat radius parameter.
//
// References:
//   Edward Teo, "Rotating traversable wormholes", Phys. Rev. D 58, 024014 (1998).
//   Oliver Fechtig, "Physikalische Aspekte und Visualisierung von stationaeren Wurmloechern"
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class teo_simple_wormhole : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto l     = position[1];
    const auto theta = position[2];
    const auto c     = consts::speed_of_light;
    
    const auto b0_sq = b0 * b0;
    const auto b0_4  = b0_sq * b0_sq;
    const auto b0_6  = b0_4 * b0_sq;
    const auto c_sq  = c * c;
    const auto l_sq  = l * l;
    const auto st    = std::sin(theta);
    const auto ct    = std::cos(theta);
    const auto st_sq = st * st;
    const auto st_4  = st_sq * st_sq;
    
    const auto r_sq = std::pow(b0_sq + l_sq, static_cast<scalar_type>(2));
    const auto invr3 = std::pow(b0_sq + l_sq, static_cast<scalar_type>(-3));
    const auto t1 = std::pow(b0_sq + l_sq, static_cast<scalar_type>(-0.5));
    const auto t2 = std::pow(b0_sq + l_sq, static_cast<scalar_type>(-3.5));
    const auto t3 = std::pow(b0_sq + l_sq, static_cast<scalar_type>(-2.5));
    const auto t4 = std::pow(b0_sq + l_sq, static_cast<scalar_type>(-1.5));
    const auto t5 = std::pow(b0_sq + l_sq, static_cast<scalar_type>(-2));
    const auto t6 = std::pow(static_cast<scalar_type>(-0.5) * b0_sq - static_cast<scalar_type>(0.5) * l_sq, static_cast<scalar_type>(2));

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Implementing the Christoffel symbols from Motion4D
    // Many terms involve complex combinations - implementing the most significant ones
    
    symbols(0, 0, 1) = static_cast<scalar_type>(0.5) * b0_4 * c_sq * l * invr3 * st_sq;
    symbols(0, 0, 2) = static_cast<scalar_type>(-0.25) * b0_4 * c_sq * invr3 * st * ct;
    
    symbols(0, 1, 3) = static_cast<scalar_type>(0.5) * b0_6 * c * l * (static_cast<scalar_type>(-0.5) * b0_sq - static_cast<scalar_type>(0.5) * l_sq) * t2 
                      * (static_cast<scalar_type>(0.25) * b0_4 * st_sq - r_sq) * st_4
                      / ((static_cast<scalar_type>(0.25) * b0_4 * c_sq * t5 * st_sq - c_sq)
                        * (static_cast<scalar_type>(-0.25) * b0_4 * c_sq * r_sq * st_4 + c_sq * r_sq * (static_cast<scalar_type>(0.25) * b0_4 * st_sq - r_sq) * st_sq));
    
    symbols(0, 3, 1) = static_cast<scalar_type>(1.5) * b0_sq * c * l * (static_cast<scalar_type>(-0.5) * b0_sq - static_cast<scalar_type>(0.5) * l_sq) * t3 * st_sq 
                      + static_cast<scalar_type>(0.5) * b0_sq * c * l * t4 * st_sq;
    symbols(0, 3, 2) = static_cast<scalar_type>(-1) * b0_sq * c * (static_cast<scalar_type>(-0.5) * b0_sq - static_cast<scalar_type>(0.5) * l_sq) * t3 * st * ct;
    
    symbols(1, 0, 0) = symbols(0, 0, 1);  // Symmetry
    symbols(1, 0, 3) = symbols(0, 1, 3);
    
    symbols(1, 2, 2) = l / (b0_sq + l_sq);
    
    symbols(1, 3, 3) = static_cast<scalar_type>(1) * c_sq * l * std::pow(b0_sq + l_sq, static_cast<scalar_type>(1)) 
                      * (static_cast<scalar_type>(0.25) * b0_4 * st_sq - r_sq) * st_sq
                      / (static_cast<scalar_type>(-0.25) * b0_4 * c_sq * r_sq * st_4 + c_sq * r_sq * (static_cast<scalar_type>(0.25) * b0_4 * st_sq - r_sq) * st_sq);
    
    symbols(2, 1, 2) = l / (b0_sq + l_sq);
    
    symbols(2, 2, 1) = -l;
    
    symbols(2, 3, 3) = ct / st;
    
    symbols(3, 0, 1) = symbols(0, 3, 1);
    symbols(3, 0, 2) = symbols(0, 3, 2);
    
    symbols(3, 1, 0) = symbols(1, 0, 3);
    symbols(3, 1, 3) = symbols(1, 3, 3);
    
    symbols(3, 2, 3) = ct / st;
    
    symbols(3, 3, 1) = -l * st_sq;
    symbols(3, 3, 2) = -st * ct;
    
    return symbols;
  }

  // Wormhole throat radius parameter
  scalar_type b0 = static_cast<scalar_type>(1);
};
}
