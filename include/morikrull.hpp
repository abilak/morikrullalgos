#pragma once

#include <string>
#include <vector>

namespace morikrull {

struct Result {
  bool holds = false;
  std::string property;
  std::string reason;
  std::vector<std::string> gens;
  std::vector<std::string> hilbert_basis;
  std::string hole;
};

// one pos algebraic number, given by its minimal poly over Q with the constant
// c first, either integers or rational strings
//
// root picks which root is meant when the poly has two pos ones, as a decimal
// or a fraction near the intended val; leave it empty when only one root is pos
struct Generator {
  std::vector<std::string> poly;
  std::string root;
};

// one gen is handled in any deg, several only when they generate a quadratic field
Result mori(const std::vector<Generator>& gens);
Result krull(const std::vector<Generator>& gens);

} // namespace morikrull
