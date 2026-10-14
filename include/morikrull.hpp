#pragma once

#include <string>
#include <vector>

namespace morikrull {

struct Result {
  bool holds = false;
  std::string property;
  std::string reason;
  std::vector<std::string> generators;
  std::vector<std::string> hilbert_basis;
  std::string hole;
};

// one positive algebraic number, given by its minimal polynomial over Q with the constant
// coefficient first, either integers or rational strings
//
// root picks which root is meant when the polynomial has two positive ones, as a decimal
// or a fraction near the intended value; leave it empty when only one root is positive
struct Generator {
  std::vector<std::string> polynomial;
  std::string root;
};

// one generator is handled in any degree, several only when they generate a quadratic field
Result mori(const std::vector<Generator>& generators);
Result krull(const std::vector<Generator>& generators);

} // namespace morikrull
