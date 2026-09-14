#pragma once

#include <string>
#include <vector>

namespace morikrull {

struct Answer {
  bool value = false;
  std::string property;
  std::string reason;
  std::vector<std::string> generators;
  std::vector<std::string> hilbert_generators;
  std::string hole;
};

// A quadratic generator is written a + b*sqrt(D).  Both a and b are rational
// strings, for example {"1/2", "1/2"} for the golden ratio in Q(sqrt(5)).
struct QuadraticGenerator {
  std::string a;
  std::string b;
};

Answer quadratic_mori(long long D, const std::vector<QuadraticGenerator>& generators);
Answer quadratic_krull(long long D, const std::vector<QuadraticGenerator>& generators);

// Coefficients are constant first.  They may be integers or rational strings.
// The polynomial is assumed to be the minimal polynomial over Q.
Answer cyclic_mori(const std::vector<std::string>& polynomial);
Answer cyclic_krull(const std::vector<std::string>& polynomial);

} // namespace morikrull

