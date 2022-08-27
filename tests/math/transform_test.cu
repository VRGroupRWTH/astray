#include <doctest/doctest.h>

#include <astray/api.hpp>

TEST_CASE("ast::transform")
{
  using transform_type = ast::transform<double>;
  using scalar_type    = transform_type::scalar_type;
  using vector_type    = transform_type::vector_type;
  
  transform_type transform;

  REQUIRE(transform.right  () == vector_type::UnitX());
  REQUIRE(transform.up     () == vector_type::UnitY());
  REQUIRE(transform.forward() == vector_type::UnitZ());

  transform.rotation_from_euler(vector_type(ast::to_radians(static_cast<scalar_type>(270)), 0, 0));
  REQUIRE((transform.forward()[0] == doctest::Approx(vector_type::UnitY()[0]) &&
           transform.forward()[1] == doctest::Approx(vector_type::UnitY()[1]) &&
           transform.forward()[2] == doctest::Approx(vector_type::UnitY()[2])));

  transform.rotation_from_euler(vector_type(0, ast::to_radians(90.0), 0));
  REQUIRE((transform.forward()[0] == doctest::Approx(vector_type::UnitX()[0]) &&
           transform.forward()[1] == doctest::Approx(vector_type::UnitX()[1]) &&
           transform.forward()[2] == doctest::Approx(vector_type::UnitX()[2])));

  transform.rotation_from_euler(vector_type(0, 0, ast::to_radians(90.0)));
  REQUIRE((transform.forward()[0] == doctest::Approx(vector_type::UnitZ()[0]) &&
           transform.forward()[1] == doctest::Approx(vector_type::UnitZ()[1]) &&
           transform.forward()[2] == doctest::Approx(vector_type::UnitZ()[2])));

  transform.translation = -10 * vector_type::UnitX();
  transform.look_at(vector_type::Zero());
  REQUIRE((transform.forward()[0] == doctest::Approx(vector_type::UnitX()[0]) &&
           transform.forward()[1] == doctest::Approx(vector_type::UnitX()[1]) &&
           transform.forward()[2] == doctest::Approx(vector_type::UnitX()[2])));

  transform.translation = -10 * vector_type::UnitY();
  transform.look_at(vector_type::Zero(), -vector_type::UnitZ());
  REQUIRE((transform.forward()[0] == doctest::Approx(vector_type::UnitY()[0]) &&
           transform.forward()[1] == doctest::Approx(vector_type::UnitY()[1]) &&
           transform.forward()[2] == doctest::Approx(vector_type::UnitY()[2])));

  transform.translation = -10 * vector_type::UnitZ();
  transform.look_at(vector_type::Zero());
  REQUIRE((transform.forward()[0] == doctest::Approx(vector_type::UnitZ()[0]) &&
           transform.forward()[1] == doctest::Approx(vector_type::UnitZ()[1]) &&
           transform.forward()[2] == doctest::Approx(vector_type::UnitZ()[2])));
}