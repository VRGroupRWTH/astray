#pragma once

namespace ast
{
enum class termination_reason
{
  none                ,
  constraint_violation,
  numeric_error       ,
  out_of_bounds       ,
  spacetime_breakdown ,
};
}