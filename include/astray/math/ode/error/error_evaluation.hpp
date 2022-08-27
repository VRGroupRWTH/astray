#pragma once

namespace ast
{
template <typename type>
struct error_evaluation
{
  bool accept;
  type next_step_size;
};
}