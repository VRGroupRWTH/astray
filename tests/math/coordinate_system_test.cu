#include <doctest/doctest.h>

#include <astray/api.hpp>

TEST_CASE("ast::coordinate_system")
{
  using scalar_type = double;
  using vector_type = ast::vector4<scalar_type>;

  auto device                 = std::random_device                         ();
  auto mersenne               = std::mt19937                               (device());
  auto cartesian_distribution = std::uniform_real_distribution<scalar_type>(1  , 10);
  auto radius_distribution    = std::uniform_real_distribution<scalar_type>(0.1, 10);
  auto pi_distribution        = std::uniform_real_distribution<scalar_type>(0.1, M_PI     - 0.1);
  auto two_pi_distribution    = std::uniform_real_distribution<scalar_type>(0.1, M_PI * 2 - 0.1);
  auto sigma_distribution     = std::uniform_real_distribution<scalar_type>( 1 , 10);
  auto tau_distribution       = std::uniform_real_distribution<scalar_type>(-1 , 1 );
  auto alpha_distribution     = std::uniform_real_distribution<scalar_type>(0.1, 10);

  for (auto i = 0; i < 10000; ++i)
  {
    auto cartesian = ast::ray<vector_type>
    {
      vector_type(0, cartesian_distribution(mersenne), cartesian_distribution(mersenne), cartesian_distribution(mersenne)),
      vector_type(0, cartesian_distribution(mersenne), cartesian_distribution(mersenne), cartesian_distribution(mersenne))
    };

    auto alpha                           = alpha_distribution(mersenne);
    auto cartesian_to_spherical          = cartesian                      ; ast::convert_ray<ast::coordinate_system_type::cartesian         , ast::coordinate_system_type::spherical>         (cartesian_to_spherical  );
    auto spherical_to_cartesian          = cartesian_to_spherical         ; ast::convert_ray<ast::coordinate_system_type::spherical         , ast::coordinate_system_type::cartesian>         (spherical_to_cartesian  );
    auto cartesian_to_cylindrical        = cartesian                      ; ast::convert_ray<ast::coordinate_system_type::cartesian         , ast::coordinate_system_type::cylindrical>       (cartesian_to_cylindrical);
    auto cylindrical_to_cartesian        = cartesian_to_cylindrical       ; ast::convert_ray<ast::coordinate_system_type::cylindrical       , ast::coordinate_system_type::cartesian>         (cylindrical_to_cartesian);
    auto cartesian_to_prolate_spheroidal = cartesian                      ; ast::convert_ray<ast::coordinate_system_type::cartesian         , ast::coordinate_system_type::prolate_spheroidal>(cartesian_to_prolate_spheroidal, alpha);
    auto prolate_spheroidal_to_cartesian = cartesian_to_prolate_spheroidal; ast::convert_ray<ast::coordinate_system_type::prolate_spheroidal, ast::coordinate_system_type::cartesian>         (prolate_spheroidal_to_cartesian, alpha);
    auto cartesian_to_boyer_lindquist    = cartesian                      ; ast::convert_ray<ast::coordinate_system_type::cartesian         , ast::coordinate_system_type::boyer_lindquist>   (cartesian_to_boyer_lindquist   , alpha);
    auto boyer_lindquist_to_cartesian    = cartesian_to_boyer_lindquist   ; ast::convert_ray<ast::coordinate_system_type::boyer_lindquist   , ast::coordinate_system_type::cartesian>         (boyer_lindquist_to_cartesian   , alpha);
    
    REQUIRE((cartesian.position [0] == doctest::Approx(spherical_to_cartesian         .position [0]) &&
             cartesian.position [1] == doctest::Approx(spherical_to_cartesian         .position [1]) &&
             cartesian.position [2] == doctest::Approx(spherical_to_cartesian         .position [2]) &&
             cartesian.position [3] == doctest::Approx(spherical_to_cartesian         .position [3])));
    REQUIRE((cartesian.position [0] == doctest::Approx(cylindrical_to_cartesian       .position [0]) &&
             cartesian.position [1] == doctest::Approx(cylindrical_to_cartesian       .position [1]) &&
             cartesian.position [2] == doctest::Approx(cylindrical_to_cartesian       .position [2]) &&
             cartesian.position [3] == doctest::Approx(cylindrical_to_cartesian       .position [3])));
    REQUIRE((cartesian.position [0] == doctest::Approx(prolate_spheroidal_to_cartesian.position [0]) &&
             cartesian.position [1] == doctest::Approx(prolate_spheroidal_to_cartesian.position [1]) &&
             cartesian.position [2] == doctest::Approx(prolate_spheroidal_to_cartesian.position [2]) &&
             cartesian.position [3] == doctest::Approx(prolate_spheroidal_to_cartesian.position [3])));
    REQUIRE((cartesian.position [0] == doctest::Approx(boyer_lindquist_to_cartesian   .position [0]) &&
             cartesian.position [1] == doctest::Approx(boyer_lindquist_to_cartesian   .position [1]) &&
             cartesian.position [2] == doctest::Approx(boyer_lindquist_to_cartesian   .position [2]) &&
             cartesian.position [3] == doctest::Approx(boyer_lindquist_to_cartesian   .position [3])));
    REQUIRE((cartesian.direction[0] == doctest::Approx(spherical_to_cartesian         .direction[0]) &&
             cartesian.direction[1] == doctest::Approx(spherical_to_cartesian         .direction[1]) &&
             cartesian.direction[2] == doctest::Approx(spherical_to_cartesian         .direction[2]) &&
             cartesian.direction[3] == doctest::Approx(spherical_to_cartesian         .direction[3])));
    REQUIRE((cartesian.direction[0] == doctest::Approx(cylindrical_to_cartesian       .direction[0]) &&
             cartesian.direction[1] == doctest::Approx(cylindrical_to_cartesian       .direction[1]) &&
             cartesian.direction[2] == doctest::Approx(cylindrical_to_cartesian       .direction[2]) &&
             cartesian.direction[3] == doctest::Approx(cylindrical_to_cartesian       .direction[3])));
    REQUIRE((cartesian.direction[0] == doctest::Approx(prolate_spheroidal_to_cartesian.direction[0]) &&
             cartesian.direction[1] == doctest::Approx(prolate_spheroidal_to_cartesian.direction[1]) &&
             cartesian.direction[2] == doctest::Approx(prolate_spheroidal_to_cartesian.direction[2]) &&
             cartesian.direction[3] == doctest::Approx(prolate_spheroidal_to_cartesian.direction[3])));
    REQUIRE((cartesian.direction[0] == doctest::Approx(boyer_lindquist_to_cartesian   .direction[0]) &&
             cartesian.direction[1] == doctest::Approx(boyer_lindquist_to_cartesian   .direction[1]) &&
             cartesian.direction[2] == doctest::Approx(boyer_lindquist_to_cartesian   .direction[2]) &&
             cartesian.direction[3] == doctest::Approx(boyer_lindquist_to_cartesian   .direction[3])));
  }
  
  for (auto i = 0; i < 10000; ++i)
  {
    auto spherical                       = vector_type(0, radius_distribution(mersenne), pi_distribution(mersenne), two_pi_distribution(mersenne));
    auto alpha                           = alpha_distribution(mersenne);
    auto spherical_to_cartesian          = spherical                      ; ast::convert<ast::coordinate_system_type::spherical         , ast::coordinate_system_type::cartesian>         (spherical_to_cartesian);
    auto cartesian_to_spherical          = spherical_to_cartesian         ; ast::convert<ast::coordinate_system_type::cartesian         , ast::coordinate_system_type::spherical>         (cartesian_to_spherical);
    auto spherical_to_cylindrical        = spherical                      ; ast::convert<ast::coordinate_system_type::spherical         , ast::coordinate_system_type::cylindrical>       (spherical_to_cylindrical);
    auto cylindrical_to_spherical        = spherical_to_cylindrical       ; ast::convert<ast::coordinate_system_type::cylindrical       , ast::coordinate_system_type::spherical>         (cylindrical_to_spherical);
    auto spherical_to_prolate_spheroidal = spherical                      ; ast::convert<ast::coordinate_system_type::spherical         , ast::coordinate_system_type::prolate_spheroidal>(spherical_to_prolate_spheroidal, alpha);
    auto prolate_spheroidal_to_spherical = spherical_to_prolate_spheroidal; ast::convert<ast::coordinate_system_type::prolate_spheroidal, ast::coordinate_system_type::spherical>         (prolate_spheroidal_to_spherical, alpha);
    auto spherical_to_boyer_lindquist    = spherical                      ; ast::convert<ast::coordinate_system_type::spherical         , ast::coordinate_system_type::boyer_lindquist>   (spherical_to_boyer_lindquist   , alpha);
    auto boyer_lindquist_to_spherical    = spherical_to_boyer_lindquist   ; ast::convert<ast::coordinate_system_type::boyer_lindquist   , ast::coordinate_system_type::spherical>         (boyer_lindquist_to_spherical   , alpha);
    REQUIRE((spherical[0] == doctest::Approx(cartesian_to_spherical         [0]) &&
             spherical[1] == doctest::Approx(cartesian_to_spherical         [1]) &&
             spherical[2] == doctest::Approx(cartesian_to_spherical         [2]) &&
             spherical[3] == doctest::Approx(cartesian_to_spherical         [3])));
    REQUIRE((spherical[0] == doctest::Approx(cylindrical_to_spherical       [0]) &&
             spherical[1] == doctest::Approx(cylindrical_to_spherical       [1]) &&
             spherical[2] == doctest::Approx(cylindrical_to_spherical       [2]) &&
             spherical[3] == doctest::Approx(cylindrical_to_spherical       [3])));
    REQUIRE((spherical[0] == doctest::Approx(prolate_spheroidal_to_spherical[0]) &&
             spherical[1] == doctest::Approx(prolate_spheroidal_to_spherical[1]) &&
             spherical[2] == doctest::Approx(prolate_spheroidal_to_spherical[2]) &&
             spherical[3] == doctest::Approx(prolate_spheroidal_to_spherical[3])));
    REQUIRE((spherical[0] == doctest::Approx(boyer_lindquist_to_spherical   [0]) &&
             spherical[1] == doctest::Approx(boyer_lindquist_to_spherical   [1]) &&
             spherical[2] == doctest::Approx(boyer_lindquist_to_spherical   [2]) &&
             spherical[3] == doctest::Approx(boyer_lindquist_to_spherical   [3])));
  }

  for (auto i = 0; i < 10000; ++i)
  {
    auto cylindrical                       = vector_type(0, radius_distribution(mersenne), two_pi_distribution(mersenne), cartesian_distribution(mersenne));
    auto alpha                             = alpha_distribution(mersenne);
    auto cylindrical_to_cartesian          = cylindrical                      ; ast::convert<ast::coordinate_system_type::cylindrical       , ast::coordinate_system_type::cartesian>         (cylindrical_to_cartesian);
    auto cartesian_to_cylindrical          = cylindrical_to_cartesian         ; ast::convert<ast::coordinate_system_type::cartesian         , ast::coordinate_system_type::cylindrical>       (cartesian_to_cylindrical);
    auto cylindrical_to_spherical          = cylindrical                      ; ast::convert<ast::coordinate_system_type::cylindrical       , ast::coordinate_system_type::spherical>         (cylindrical_to_spherical);
    auto spherical_to_cylindrical          = cylindrical_to_spherical         ; ast::convert<ast::coordinate_system_type::spherical         , ast::coordinate_system_type::cylindrical>       (spherical_to_cylindrical);
    auto cylindrical_to_prolate_spheroidal = cylindrical                      ; ast::convert<ast::coordinate_system_type::cylindrical       , ast::coordinate_system_type::prolate_spheroidal>(cylindrical_to_prolate_spheroidal, alpha);
    auto prolate_spheroidal_to_cylindrical = cylindrical_to_prolate_spheroidal; ast::convert<ast::coordinate_system_type::prolate_spheroidal, ast::coordinate_system_type::cylindrical>       (prolate_spheroidal_to_cylindrical, alpha);
    auto cylindrical_to_boyer_lindquist    = cylindrical                      ; ast::convert<ast::coordinate_system_type::cylindrical       , ast::coordinate_system_type::boyer_lindquist>   (cylindrical_to_boyer_lindquist   , alpha);
    auto boyer_lindquist_to_cylindrical    = cylindrical_to_prolate_spheroidal; ast::convert<ast::coordinate_system_type::boyer_lindquist   , ast::coordinate_system_type::cylindrical>       (boyer_lindquist_to_cylindrical   , alpha);
    REQUIRE((cylindrical[0] == doctest::Approx(cartesian_to_cylindrical         [0]) &&
             cylindrical[1] == doctest::Approx(cartesian_to_cylindrical         [1]) &&
             cylindrical[2] == doctest::Approx(cartesian_to_cylindrical         [2]) &&
             cylindrical[3] == doctest::Approx(cartesian_to_cylindrical         [3])));
    REQUIRE((cylindrical[0] == doctest::Approx(spherical_to_cylindrical         [0]) &&
             cylindrical[1] == doctest::Approx(spherical_to_cylindrical         [1]) &&
             cylindrical[2] == doctest::Approx(spherical_to_cylindrical         [2]) &&
             cylindrical[3] == doctest::Approx(spherical_to_cylindrical         [3])));
    REQUIRE((cylindrical[0] == doctest::Approx(prolate_spheroidal_to_cylindrical[0]) &&
             cylindrical[1] == doctest::Approx(prolate_spheroidal_to_cylindrical[1]) &&
             cylindrical[2] == doctest::Approx(prolate_spheroidal_to_cylindrical[2]) &&
             cylindrical[3] == doctest::Approx(prolate_spheroidal_to_cylindrical[3])));
    //REQUIRE((cylindrical[0] == doctest::Approx(boyer_lindquist_to_cylindrical   [0]) &&
    //         cylindrical[1] == doctest::Approx(boyer_lindquist_to_cylindrical   [1]) &&
    //         cylindrical[2] == doctest::Approx(boyer_lindquist_to_cylindrical   [2]) &&
    //         cylindrical[3] == doctest::Approx(boyer_lindquist_to_cylindrical   [3])));
  }
  
  for (auto i = 0; i < 10000; ++i)
  {
    auto prolate_spheroidal                = vector_type(0, sigma_distribution(mersenne), tau_distribution(mersenne), two_pi_distribution(mersenne));
    auto alpha                             = alpha_distribution(mersenne);
    auto prolate_spheroidal_to_cartesian   = prolate_spheroidal               ; ast::convert<ast::coordinate_system_type::prolate_spheroidal, ast::coordinate_system_type::cartesian>         (prolate_spheroidal_to_cartesian  , alpha);
    auto cartesian_to_prolate_spheroidal   = prolate_spheroidal_to_cartesian  ; ast::convert<ast::coordinate_system_type::cartesian         , ast::coordinate_system_type::prolate_spheroidal>(cartesian_to_prolate_spheroidal  , alpha);
    auto prolate_spheroidal_to_cylindrical = prolate_spheroidal               ; ast::convert<ast::coordinate_system_type::prolate_spheroidal, ast::coordinate_system_type::cylindrical>       (prolate_spheroidal_to_cylindrical, alpha);
    auto cylindrical_to_prolate_spheroidal = prolate_spheroidal_to_cylindrical; ast::convert<ast::coordinate_system_type::cylindrical       , ast::coordinate_system_type::prolate_spheroidal>(cylindrical_to_prolate_spheroidal, alpha);
    auto prolate_spheroidal_to_spherical   = prolate_spheroidal               ; ast::convert<ast::coordinate_system_type::prolate_spheroidal, ast::coordinate_system_type::spherical>         (prolate_spheroidal_to_spherical  , alpha);
    auto spherical_to_prolate_spheroidal   = prolate_spheroidal_to_spherical  ; ast::convert<ast::coordinate_system_type::spherical         , ast::coordinate_system_type::prolate_spheroidal>(spherical_to_prolate_spheroidal  , alpha);
    REQUIRE((prolate_spheroidal[0] == doctest::Approx(cartesian_to_prolate_spheroidal  [0]) &&
             prolate_spheroidal[1] == doctest::Approx(cartesian_to_prolate_spheroidal  [1]) &&
             prolate_spheroidal[2] == doctest::Approx(cartesian_to_prolate_spheroidal  [2]) &&
             prolate_spheroidal[3] == doctest::Approx(cartesian_to_prolate_spheroidal  [3])));
    REQUIRE((prolate_spheroidal[0] == doctest::Approx(cylindrical_to_prolate_spheroidal[0]) &&
             prolate_spheroidal[1] == doctest::Approx(cylindrical_to_prolate_spheroidal[1]) &&
             prolate_spheroidal[2] == doctest::Approx(cylindrical_to_prolate_spheroidal[2]) &&
             prolate_spheroidal[3] == doctest::Approx(cylindrical_to_prolate_spheroidal[3])));
    REQUIRE((prolate_spheroidal[0] == doctest::Approx(spherical_to_prolate_spheroidal  [0]) &&
             prolate_spheroidal[1] == doctest::Approx(spherical_to_prolate_spheroidal  [1]) &&
             prolate_spheroidal[2] == doctest::Approx(spherical_to_prolate_spheroidal  [2]) &&
             prolate_spheroidal[3] == doctest::Approx(spherical_to_prolate_spheroidal  [3])));
  }
  
  for (auto i = 0; i < 10000; ++i)
  {
    auto boyer_lindquist                = vector_type(0, radius_distribution(mersenne), pi_distribution(mersenne), two_pi_distribution(mersenne));
    auto alpha                          = alpha_distribution(mersenne);
    auto boyer_lindquist_to_cartesian   = boyer_lindquist               ; ast::convert<ast::coordinate_system_type::boyer_lindquist, ast::coordinate_system_type::cartesian>      (boyer_lindquist_to_cartesian  , alpha);
    auto cartesian_to_boyer_lindquist   = boyer_lindquist_to_cartesian  ; ast::convert<ast::coordinate_system_type::cartesian      , ast::coordinate_system_type::boyer_lindquist>(cartesian_to_boyer_lindquist  , alpha);
    auto boyer_lindquist_to_cylindrical = boyer_lindquist               ; ast::convert<ast::coordinate_system_type::boyer_lindquist, ast::coordinate_system_type::cylindrical>    (boyer_lindquist_to_cylindrical, alpha);
    auto cylindrical_to_boyer_lindquist = boyer_lindquist_to_cylindrical; ast::convert<ast::coordinate_system_type::cylindrical    , ast::coordinate_system_type::boyer_lindquist>(cylindrical_to_boyer_lindquist, alpha);
    auto boyer_lindquist_to_spherical   = boyer_lindquist               ; ast::convert<ast::coordinate_system_type::boyer_lindquist, ast::coordinate_system_type::spherical>      (boyer_lindquist_to_spherical  , alpha);
    auto spherical_to_boyer_lindquist   = boyer_lindquist_to_spherical  ; ast::convert<ast::coordinate_system_type::spherical      , ast::coordinate_system_type::boyer_lindquist>(spherical_to_boyer_lindquist  , alpha);
    REQUIRE((boyer_lindquist[0] == doctest::Approx(cartesian_to_boyer_lindquist  [0]) &&
             boyer_lindquist[1] == doctest::Approx(cartesian_to_boyer_lindquist  [1]) &&
             boyer_lindquist[2] == doctest::Approx(cartesian_to_boyer_lindquist  [2]) &&
             boyer_lindquist[3] == doctest::Approx(cartesian_to_boyer_lindquist  [3])));
    REQUIRE((boyer_lindquist[0] == doctest::Approx(cylindrical_to_boyer_lindquist[0]) &&
             boyer_lindquist[1] == doctest::Approx(cylindrical_to_boyer_lindquist[1]) &&
             boyer_lindquist[2] == doctest::Approx(cylindrical_to_boyer_lindquist[2]) &&
             boyer_lindquist[3] == doctest::Approx(cylindrical_to_boyer_lindquist[3])));
    REQUIRE((boyer_lindquist[0] == doctest::Approx(spherical_to_boyer_lindquist  [0]) &&
             boyer_lindquist[1] == doctest::Approx(spherical_to_boyer_lindquist  [1]) &&
             boyer_lindquist[2] == doctest::Approx(spherical_to_boyer_lindquist  [2]) &&
             boyer_lindquist[3] == doctest::Approx(spherical_to_boyer_lindquist  [3])));
  }
}
