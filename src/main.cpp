#include "morikrull.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using morikrull::Answer;

namespace {

std::vector<std::string> split(const std::string& s, char delimiter)
{
  std::vector<std::string> out;
  std::stringstream in(s);
  std::string part;
  while (std::getline(in, part, delimiter)) out.push_back(part);
  return out;
}

void print(const Answer& answer)
{
  std::cout << answer.property << ": " << (answer.value ? "true" : "false") << '\n';
  std::cout << "reason: " << answer.reason << '\n';

  if (!answer.generators.empty()) {
    std::cout << "additive generators (" << answer.generators.size() << "):\n";
    for (const auto& g : answer.generators) std::cout << "  " << g << '\n';
  }
  if (!answer.hilbert_generators.empty()) {
    std::cout << "cone lattice generators (" << answer.hilbert_generators.size() << "):\n";
    for (const auto& h : answer.hilbert_generators) std::cout << "  " << h << '\n';
  }
  if (!answer.hole.empty()) std::cout << "lattice hole: " << answer.hole << '\n';
}

[[noreturn]] void usage()
{
  std::cerr
    << "usage:\n"
    << "  morikrull quadratic-mori  D a,b [a,b ...]\n"
    << "  morikrull quadratic-krull D a,b [a,b ...]\n"
    << "  morikrull cyclic-mori  a0,a1,...,ad\n"
    << "  morikrull cyclic-krull a0,a1,...,ad\n\n"
    << "Polynomial coefficients are constant first. Rational values use p/q.\n"
    << "Example: morikrull quadratic-krull 5 0,1 1/2,1/2\n"
    << "Example: morikrull cyclic-krull -2,0,1\n";
  std::exit(2);
}

} // namespace

int main(int argc, char** argv)
{
  try {
    if (argc < 3) usage();
    std::string command = argv[1];

    if (command == "cyclic-mori" || command == "cyclic-krull") {
      if (argc != 3) usage();
      auto polynomial = split(argv[2], ',');
      print(command == "cyclic-mori"
        ? morikrull::cyclic_mori(polynomial)
        : morikrull::cyclic_krull(polynomial));
      return 0;
    }

    if (command == "quadratic-mori" || command == "quadratic-krull") {
      if (argc < 4) usage();
      long long D = std::stoll(argv[2]);
      std::vector<morikrull::QuadraticGenerator> generators;
      for (int i = 3; i < argc; ++i) {
        auto pair = split(argv[i], ',');
        if (pair.size() != 2) usage();
        generators.push_back({pair[0], pair[1]});
      }
      print(command == "quadratic-mori"
        ? morikrull::quadratic_mori(D, generators)
        : morikrull::quadratic_krull(D, generators));
      return 0;
    }

    usage();
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return 2;
  }
}
