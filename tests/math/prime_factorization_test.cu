#include <doctest/doctest.h>

#include <astray/api.hpp>

TEST_CASE("ast::prime_factorization")
{
  {
    const auto result   = ast::prime_factorize<std::size_t>(2);
    const auto expected = std::vector<std::size_t>({2});
    REQUIRE(result == expected);
  }
  {
    const auto result   = ast::prime_factorize<std::size_t>(3);
    const auto expected = std::vector<std::size_t>({3});
    REQUIRE(result == expected);
  }
  {
    const auto result   = ast::prime_factorize<std::size_t>(4);
    const auto expected = std::vector<std::size_t>({2, 2});
    REQUIRE(result == expected);
  }
  {
    const auto result   = ast::prime_factorize<std::size_t>(5);
    const auto expected = std::vector<std::size_t>({5});
    REQUIRE(result == expected);
  }
  {
    const auto result   = ast::prime_factorize<std::size_t>(6);
    const auto expected = std::vector<std::size_t>({2, 3});
    REQUIRE(result == expected);
  }
  {
    const auto result   = ast::prime_factorize<std::size_t>(7);
    const auto expected = std::vector<std::size_t>({7});
    REQUIRE(result == expected);
  }
  {
    const auto result   = ast::prime_factorize<std::size_t>(8);
    const auto expected = std::vector<std::size_t>({2, 2, 2});
    REQUIRE(result == expected);
  }
  {
    const auto result   = ast::prime_factorize<std::size_t>(9);
    const auto expected = std::vector<std::size_t>({3, 3});
    REQUIRE(result == expected);
  }
  {
    const auto result   = ast::prime_factorize<std::size_t>(10);
    const auto expected = std::vector<std::size_t>({2, 5});
    REQUIRE(result == expected);
  }
}