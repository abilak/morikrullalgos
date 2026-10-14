#include "morikrull.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using morikrull::Generator;
using morikrull::Result;

namespace {

std::vector<std::string> split(const std::string& text, char separator)
{
  std::vector<std::string> pieces;
  std::stringstream line(text);
  std::string piece;
  while (std::getline(line, piece, separator)) pieces.push_back(piece);
  return pieces;
}

Generator parse(const std::string& argument)
{
  size_t at = argument.find('@');
  if (at == std::string::npos) return {split(argument, ','), ""};
  return {split(argument.substr(0, at), ','), argument.substr(at + 1)};
}

void print_result(const Result& result)
{
  std::cout << result.property << ": " << (result.holds ? "true" : "false") << '\n';
  std::cout << "reason: " << result.reason << '\n';

  if (!result.generators.empty()) {
    std::cout << "additive generators (" << result.generators.size() << "):\n";
    for (const auto& generator : result.generators) std::cout << "  " << generator << '\n';
  }
  if (!result.hilbert_basis.empty()) {
    std::cout << "Hilbert basis (" << result.hilbert_basis.size() << "):\n";
    for (const auto& point : result.hilbert_basis) std::cout << "  " << point << '\n';
  }
  if (!result.hole.empty()) std::cout << "lattice hole: " << result.hole << '\n';
}

[[noreturn]] void usage()
{
  std::cerr
    << "usage:\n"
    << "  morikrull mori  <generator> [<generator> ...]\n"
    << "  morikrull krull <generator> [<generator> ...]\n\n"
    << "A generator is the minimal polynomial of a positive algebraic number, constant\n"
    << "coefficient first, p/q for rationals. Append @value to pick a root when the\n"
    << "polynomial has two positive ones.\n\n"
    << "One generator works in any degree. Several have to generate a quadratic field.\n\n"
    << "  morikrull krull -2,0,1              sqrt 2\n"
    << "  morikrull krull -1,-1,1            the golden ratio\n"
    << "  morikrull krull -5,0,1 -1,-1,1     sqrt 5 and the golden ratio\n"
    << "  morikrull krull 5,-5,1@3.6 -5,0,1  (5 + sqrt 5)/2 and sqrt 5\n";
  std::exit(2);
}

int run(int argc, char** argv)
{
  if (argc < 3) usage();
  std::string mode = argv[1];
  if (mode != "mori" && mode != "krull") usage();

  std::vector<Generator> generators;
  for (int i = 2; i < argc; ++i) generators.push_back(parse(argv[i]));

  print_result(mode == "mori" ? morikrull::mori(generators) : morikrull::krull(generators));
  return 0;
}

} // namespace

int main(int argc, char** argv)
{
  try {
    return run(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
