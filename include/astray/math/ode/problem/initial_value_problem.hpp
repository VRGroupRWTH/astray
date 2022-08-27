#pragma once

#include <astray/parallel/thrust.hpp>

namespace ast
{
template <typename time_type_, typename value_type_, typename function_type_ = device_function<value_type_(time_type_, const value_type_&)>>
struct initial_value_problem
{
  using time_type     = time_type_    ;
  using value_type    = value_type_   ;
  using function_type = function_type_;

  time_type     time    = time_type(0);
  value_type    value   = {};
  function_type function;
};
}