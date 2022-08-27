#include <doctest/doctest.h>

#include <astray/api.hpp>

TEST_CASE("ast::angle")
{
  using scalar_type = double;

  constexpr auto radians_1   = M_PI    ;
  constexpr auto radians_2   = M_PI * 2;
  constexpr auto degrees_180 = 180.0   ;
  constexpr auto degrees_360 = 360.0   ;

  REQUIRE(degrees_180 == ast::to_degrees(radians_1  ));
  REQUIRE(degrees_360 == ast::to_degrees(radians_2  ));
  REQUIRE(radians_1   == ast::to_radians(degrees_180));
  REQUIRE(radians_2   == ast::to_radians(degrees_360));

  auto radians_0_closed = static_cast<scalar_type>(0); 
  auto radians_0_open   = static_cast<scalar_type>(0); 
  auto radians_1_closed = radians_1;                   
  auto radians_1_open   = radians_1;                   
  auto radians_2_closed = radians_2;                   
  auto radians_2_open   = radians_2;                   

  ast::wrap_angle(radians_0_closed,     M_PI, false);
  ast::wrap_angle(radians_0_open  ,     M_PI, true );
  ast::wrap_angle(radians_1_closed,     M_PI, false);
  ast::wrap_angle(radians_1_open  ,     M_PI, true );
  ast::wrap_angle(radians_2_closed, 2 * M_PI, false);
  ast::wrap_angle(radians_2_open  , 2 * M_PI, true );
  
  REQUIRE(radians_0_closed == doctest::Approx(0       ));
  REQUIRE(radians_0_open   == doctest::Approx(0       ));
  REQUIRE(radians_1_closed == doctest::Approx(0       ));
  REQUIRE(radians_1_open   == doctest::Approx(    M_PI));
  REQUIRE(radians_2_closed == doctest::Approx(0       ));
  REQUIRE(radians_2_open   == doctest::Approx(2 * M_PI));
}