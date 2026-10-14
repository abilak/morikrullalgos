#include "morikrull.hpp"
#include <gmpxx.h> 
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/qqbar.h>
#include <flint/ulong_extras.h>
#include <libnormaliz/cone.h>
#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace morikrull {
namespace {

using Point = std::vector<mpz_class>;
using Points = std::vector<Point>;

mpq_class rational(const std::string& text) {
  mpq_class value;
  if (value.set_str(text, 10) != 0) throw std::invalid_argument("bad coefficient " + text);
  if (value.get_den() == 0) throw std::invalid_argument("zero denominator in " + text);
  value.canonicalize();
  return value;
}

mpq_class rational(const mpz_class& numerator, const mpz_class& denominator) {
  mpq_class value(numerator, denominator);
  value.canonicalize();
  return value;
}

void lcm_denominator(mpz_class& common, const mpq_class& value) {
  mpz_lcm(common.get_mpz_t(), common.get_mpz_t(), value.get_den().get_mpz_t());
}

std::string point_string(const Point& point) {
  std::ostringstream out;
  out << '[';
  for (size_t i = 0; i < point.size(); ++i) out << (i ? ", " : "") << point[i];
  out << ']';
  return out.str();
}

bool contains(Points generators, const Point& target) {
  auto vanishes = [](const Point& p) {
    return std::all_of(p.begin(), p.end(), [](const mpz_class& x) { return x == 0; });
  };
  generators.erase(std::remove_if(generators.begin(), generators.end(), vanishes), generators.end());
  if (vanishes(target)) return true;
  if (generators.empty()) return false;

  std::map<libnormaliz::Type::InputType, Points> input;
  for (size_t coordinate = 0; coordinate < target.size(); ++coordinate) {
    Point equation;
    for (const auto& generator : generators) equation.push_back(generator[coordinate]);
    equation.push_back(-target[coordinate]);
    input[libnormaliz::Type::inhom_equations].push_back(std::move(equation));
  }
  input[libnormaliz::Type::signs] = {Point(generators.size(), 1)};
  return libnormaliz::Cone<mpz_class>(input).getNrModuleGenerators() != 0;
}

// prop 5.10 and cor 5.11 -- cone_and_lattice makes the ambient lattice q(H), where closedness is tested
template <typename Show>
void krull_step(Result& result, const Points& lattice, Show show)
{
  libnormaliz::Cone<mpz_class> cone(libnormaliz::Type::cone_and_lattice, lattice);
  for (const auto& point : cone.getHilbertBasis()) result.hilbert_basis.push_back(show(point));

  if (!cone.isIntegrallyClosed()) {
    result.reason = "not integer divisor closed in q(H)";
    result.hole = show(cone.getWitnessNotIntegrallyClosed());
    return;
  }
  result.holds = true;
  result.reason = "finitely generated and integer divisor closed in q(H)";
}

struct Quad {
  mpq_class a, b; // a + b*sqrt(D)
  bool operator<(const Quad& y) const { return std::tie(a, b) < std::tie(y.a, y.b); }
};

struct Field {
  mpq_class D;

  Quad mul(const Quad& x, const Quad& y) const {
    return {x.a * y.a + D * x.b * y.b, x.a * y.b + x.b * y.a};
  }

  Quad conjugate(const Quad& x) const { return {x.a, -x.b}; }

  int sign(const Quad& x) const {
    int sa = sgn(x.a), sb = sgn(x.b);
    if (sa == 0 || sb == 0 || sa == sb) return sa != 0 ? sa : sb;
    int larger = cmp(x.a * x.a, D * x.b * x.b);
    return larger == 0 ? 0 : (larger > 0 ? sa : sb);
  }

  bool integral(const Quad& x) const {
    mpq_class trace = 2 * x.a, norm = x.a * x.a - D * x.b * x.b;
    return trace.get_den() == 1 && norm.get_den() == 1;
  }

  // |a - b*sqrt(D)| <= a + b*sqrt(D) for positive x
  bool weak_perron(const Quad& x) const { return x.a >= 0 && x.b >= 0; }

  std::string text(const Quad& x) const {
    std::ostringstream out;
    mpq_class magnitude = abs(x.b);
    out << x.a.get_str() << (sgn(x.b) < 0 ? " - " : " + ") << magnitude.get_str()
        << "*sqrt(" << D.get_str() << ')';
    return out.str();
  }
};

std::pair<mpz_class, Points> to_lattice(const std::vector<Quad>& values) {
  mpz_class common = 1;
  for (const auto& x : values) {
    lcm_denominator(common, x.a);
    lcm_denominator(common, x.b);
  }

  Points points;
  for (const auto& x : values)
    points.push_back({x.a.get_num() * (common / x.a.get_den()), x.b.get_num() * (common / x.b.get_den())});
  return {common, points};
}

bool contains(const std::vector<Quad>& generators, const Quad& target)
{
  std::vector<Quad> together = generators;
  together.push_back(target);
  Points points = to_lattice(together).second;
  Point goal = std::move(points.back());
  points.pop_back();
  return contains(std::move(points), goal);
}

// search at end of section 4: G holds the monomials of degree at most m, and m grows until G is closed under multiplication by every alpha_i
std::vector<Quad> generating_set(const Field& field, const std::vector<Quad>& alphas)
{
  const Quad one{1, 0};
  std::set<Quad> G{one}, degree{one}; // G is G_m, degree holds the monomials of degree m

  for (;;) {
    std::vector<Quad> monomials(G.begin(), G.end());
    bool closed = true;
    for (const auto& g : monomials)
      for (const auto& alpha : alphas)
        closed = closed && contains(monomials, field.mul(g, alpha));
    if (closed) return monomials;

    std::set<Quad> next;
    for (const auto& x : degree)
      for (const auto& alpha : alphas) next.insert(field.mul(x, alpha));
    degree = std::move(next);
    G.insert(degree.begin(), degree.end());
  }
}

struct QuadraticCheck {
  bool finite = false;
  std::string reason;
  std::vector<Quad> alphas;
};

// finite generation classification, thm 4.2
QuadraticCheck quadratic_check(const Field& field, std::vector<Quad> alphas)
{
  QuadraticCheck check;
  bool negative_conjugate = false;

  for (const auto& alpha : alphas) {
    if (!field.integral(alpha)) {
      check.reason = field.text(alpha) + " is not an algebraic integer";
      return check;
    }
    if (!field.weak_perron(alpha)) {
      check.reason = field.text(alpha) + " is not weak Perron";
      return check;
    }
    negative_conjugate = negative_conjugate || field.sign(field.conjugate(alpha)) < 0;
  }

  if (!negative_conjugate) {
    check.reason = "no generator has a negative conjugate";
    return check;
  }

  check.alphas = std::move(alphas);
  check.finite = true;
  check.reason = "every generator is a weak Perron algebraic integer and one has a negative conjugate";
  return check;
}

struct Roots {
  int positive = 0;
  bool bigger_conjugate = false;
};

Roots inspect(const std::vector<mpz_class>& coefficients) {
  fmpz_poly_t polynomial;
  fmpz_poly_init(polynomial);
  fmpz_t value;
  fmpz_init(value);
  for (size_t i = 0; i < coefficients.size(); ++i) {
    fmpz_set_mpz(value, coefficients[i].get_mpz_t());
    fmpz_poly_set_coeff_fmpz(polynomial, i, value);
  }
  fmpz_clear(value);

  fmpz_poly_factor_t factors;
  fmpz_poly_factor_init(factors);
  fmpz_poly_factor(factors, polynomial);
  bool irreducible = factors->num == 1 && factors->exp[0] == 1;
  fmpz_poly_factor_clear(factors);
  if (!irreducible) {
    fmpz_poly_clear(polynomial);
    throw std::invalid_argument("polynomial is reducible");
  }

  slong degree = fmpz_poly_degree(polynomial);
  qqbar_ptr roots = _qqbar_vec_init(degree);
  qqbar_roots_fmpz_poly(roots, polynomial, QQBAR_ROOTS_IRREDUCIBLE);

  Roots result;
  qqbar_ptr alpha = nullptr;
  for (slong i = 0; i < degree; ++i)
    if (qqbar_is_real(roots + i) && qqbar_sgn_re(roots + i) > 0) {
      alpha = roots + i;
      ++result.positive;
    }
  if (result.positive == 1)
    for (slong i = 0; i < degree; ++i) result.bigger_conjugate |= qqbar_cmpabs(roots + i, alpha) > 0;

  _qqbar_vec_clear(roots, degree);
  fmpz_poly_clear(polynomial);
  return result;
}

struct CyclicCheck {
  bool finite = false;
  std::string reason;
  std::vector<mpz_class> monic; 
};

// cor 3.13
CyclicCheck cyclic_check(const std::vector<std::string>& coefficients) {
  std::vector<mpq_class> polynomial;
  for (const auto& text : coefficients) polynomial.push_back(rational(text));
  while (polynomial.size() > 1 && polynomial.back() == 0) polynomial.pop_back();
  if (polynomial.size() < 2) throw std::invalid_argument("need degree at least one");

  mpz_class common = 1;
  for (const auto& coefficient : polynomial) lcm_denominator(common, coefficient);
  std::vector<mpz_class> cleared;
  for (const auto& coefficient : polynomial)
    cleared.push_back(coefficient.get_num() * (common / coefficient.get_den()));
  Roots roots = inspect(cleared);

  CyclicCheck check;
  for (const auto& coefficient : polynomial) {
    mpq_class scaled = coefficient / polynomial.back();
    if (scaled.get_den() != 1) {
      check.reason = "alpha is not an algebraic integer";
      return check;
    }
    check.monic.push_back(scaled.get_num());
  }

  if (roots.positive != 1) {
    check.reason = "alpha is not the only positive root";
    return check;
  }
  if (roots.bigger_conjugate) {
    check.reason = "a conjugate has modulus larger than alpha";
    return check;
  }
  check.finite = true;
  check.reason = "alpha is a weak Perron algebraic integer and the only positive root";
  return check;
}

// section 5.2, in the basis 1, alpha, ..., alpha^(d-1): v_(n+d) = -sum a_i v_(n+i), and k is the first power lying in the span of the earlier ones
Points cyclic_generating_set(const std::vector<mpz_class>& monic) {
  size_t d = monic.size() - 1;
  Points powers;
  for (size_t i = 0; i < d; ++i) {
    Point unit(d);
    unit[i] = 1;
    powers.push_back(std::move(unit));
  }

  for (size_t s = d;; ++s) {
    Point next(d);
    for (size_t i = 0; i < d; ++i)
      for (size_t j = 0; j < d; ++j) next[j] -= monic[i] * powers[s - d + i][j];

    if (contains(powers, next)) return powers;
    powers.push_back(std::move(next));
  }
}

mpz_class squarefree_part(const mpz_class& n) {
  fmpz_t value;
  fmpz_init(value);
  fmpz_set_mpz(value, n.get_mpz_t());
  fmpz_factor_t factors;
  fmpz_factor_init(factors);
  fmpz_factor(factors, value);

  mpz_class part = 1;
  for (slong i = 0; i < factors->num; ++i)
    if (factors->exp[i] % 2) {
      mpz_class prime;
      fmpz_get_mpz(prime.get_mpz_t(), factors->p + i);
      part *= prime;
    }

  fmpz_factor_clear(factors);
  fmpz_clear(value);
  return part;
}

mpq_class decimal(const std::string& text) {
  size_t dot = text.find('.');
  if (dot == std::string::npos) return rational(text);

  mpz_class scale = 1;
  for (size_t i = dot + 1; i < text.size(); ++i) scale *= 10;
  return rational(mpz_class(text.substr(0, dot) + text.substr(dot + 1)), scale);
}

// a generator with its root picked and its field known
struct Root {
  Quad value; 
  mpz_class D;
};

Root parse_root(const Generator& given) {
  std::vector<mpq_class> p;
  for (const auto& text : given.polynomial) p.push_back(rational(text));
  while (p.size() > 1 && p.back() == 0) p.pop_back();
  if (p.size() < 2 || p.size() > 3)
    throw std::invalid_argument("with several generators every polynomial must have degree 1 or 2");

  if (p.size() == 2) {
    mpq_class value = -p[0] / p[1];
    if (value <= 0) throw std::invalid_argument("generators must be positive");
    return {{value, 0}, 0};
  }

  mpz_class common = 1;
  for (const auto& c : p) lcm_denominator(common, c);
  if (p[2] < 0) common = -common; // keep the leading coefficient positive
  mpz_class a0 = p[0].get_num() * (common / p[0].get_den());
  mpz_class a1 = p[1].get_num() * (common / p[1].get_den());
  mpz_class a2 = p[2].get_num() * (common / p[2].get_den());

  mpz_class disc = a1 * a1 - 4 * a0 * a2;
  if (disc < 0) throw std::invalid_argument("root is not real");
  mpz_class D = disc == 0 ? 1 : squarefree_part(disc);
  if (D == 1) throw std::invalid_argument("polynomial has rational roots");

  mpq_class centre = rational(-a1, 2 * a2), spread = rational(sqrt(disc / D), 2 * a2);
  Field field{mpq_class(D)};
  Quad low{centre, -spread}, high{centre, spread};

  if (!given.root.empty()) {
    Quad picked = decimal(given.root) > centre ? high : low;
    if (field.sign(picked) <= 0) throw std::invalid_argument("generators must be positive");
    return {picked, D};
  }

  if (field.sign(low) > 0) throw std::invalid_argument("two positive roots, pick one with @");
  if (field.sign(high) <= 0) throw std::invalid_argument("generators must be positive");
  return {high, D};
}

Result decide_quadratic(std::vector<Quad> alphas, const mpz_class& D, bool krull) {
  Field field{mpq_class(D)};
  QuadraticCheck check = quadratic_check(field, std::move(alphas));

  Result result;
  result.property = krull ? "Krull" : "Mori";
  if (!check.finite) {
    result.reason = check.reason + "; not finitely generated";
    return result;
  }

  std::vector<Quad> generators = generating_set(field, check.alphas);
  for (const auto& x : generators) result.generators.push_back(field.text(x));
  if (!krull) {
    result.holds = true;
    result.reason = check.reason + "; Mori is the same as finite generation here";
    return result;
  }

  auto [common, lattice] = to_lattice(generators);
  krull_step(result, lattice, [&](const Point& point) {
    return field.text({rational(point[0], common), rational(point[1], common)});
  });
  return result;
}

Result decide_cyclic(const std::vector<std::string>& polynomial, bool krull) {
  CyclicCheck check = cyclic_check(polynomial);

  Result result;
  result.property = krull ? "Krull" : "Mori";
  if (!check.finite) {
    result.reason = check.reason + "; not finitely generated";
    return result;
  }

  Points generators = cyclic_generating_set(check.monic);
  for (size_t i = 0; i < generators.size(); ++i)
    result.generators.push_back("alpha^" + std::to_string(i) + " = " + point_string(generators[i]));
  if (!krull) {
    result.holds = true;
    result.reason = check.reason + "; Mori is the same as finite generation";
    return result;
  }

  krull_step(result, generators, point_string);
  return result;
}

// one generator is the cyclic case and works in any degree
// several only work when they all live in one quadratic field, which is as far as thm 4.2 reaches
Result decide(const std::vector<Generator>& given, bool krull) {
  if (given.empty()) throw std::invalid_argument("need at least one generator");
  if (given.size() == 1) {
    if (!given.front().root.empty())
      throw std::invalid_argument("@ does nothing with one generator");
    return decide_cyclic(given.front().polynomial, krull);
  }

  std::vector<Quad> alphas;
  mpz_class D = 0;
  for (const auto& generator : given) {
    Root root = parse_root(generator);
    if (root.D != 0) {
      if (D != 0 && D != root.D)
        throw std::invalid_argument("generators lie in different quadratic fields");
      D = root.D;
    }
    alphas.push_back(std::move(root.value));
  }
  if (D == 0) throw std::invalid_argument("generators only span Q");

  return decide_quadratic(std::move(alphas), D, krull);
}

} // namespace

Result mori(const std::vector<Generator>& generators) { return decide(generators, false); }

Result krull(const std::vector<Generator>& generators) { return decide(generators, true); }

} // namespace morikrull
