#include "morikrull.hpp"
#include <gmpxx.h> 
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/qqbar.h>
#include <libnormaliz/cone.h>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace morikrull {
namespace {

using std::invalid_argument;
using std::string;
using std::vector;

using Point = vector<mpz_class>;
using Points = vector<Point>;

mpq_class rational(const string& inp) {
  mpq_class val;
  if(val.set_str(inp, 10) != 0) { 
    throw invalid_argument("bad coeff " + inp);
  }
  if(val.get_den() == 0) { 
    throw invalid_argument("zero denom in " + inp);
  }
  val.canonicalize();
  return val;
}

// pair<vector<mpz_class>, mpz_class> clear_denoms(const vector<mpq_class>& vals);
vector<mpz_class> clear_denoms(const vector<mpq_class>& vals, mpz_class& lcm) {
  lcm = 1;
  for(const auto& v : vals) {
    mpz_lcm(lcm.get_mpz_t(), lcm.get_mpz_t(), v.get_den().get_mpz_t()); 
  } 

  vector<mpz_class> new_nums;
  for(const auto& v : vals) {
    new_nums.push_back(v.get_num() * (lcm / v.get_den()));
  } 
  return new_nums;
}

vector<mpq_class> read_poly(const vector<string>& coeffs) {
  vector<mpq_class> p;
  for(const auto& text : coeffs) {
    p.push_back(rational(text));
  }
  while(p.size() > 1 && p.back() == 0) {
    p.pop_back();
  }
  return p;
}

bool contains(const Points& gens, const Point& target) {
  std::map<libnormaliz::Type::InputType, Points> input;
  for(size_t coord = 0; coord < target.size(); ++coord) {
    Point eq;
    for(const auto& gen : gens) {
      eq.push_back(gen[coord]);
    }
    eq.push_back(-target[coord]);
    input[libnormaliz::Type::inhom_equations].push_back(std::move(eq));
  }
  input[libnormaliz::Type::signs] = {Point(gens.size(), 1)};

  return libnormaliz::Cone<mpz_class>(input).getNrModuleGenerators() != 0;
}

// prop 5.10 and cor 5.11
template <typename Show>
void krull_step(Result& res, const Points& lattice, Show show) {
  libnormaliz::Cone<mpz_class> cone(libnormaliz::Type::cone_and_lattice, lattice);
  for(const auto& point : cone.getHilbertBasis()) {
    res.hilbert_basis.push_back(show(point));
  }

  if(!cone.isIntegrallyClosed()) {
    res.reason = "not idv closed in q(H)";
    res.hole = show(cone.getWitnessNotIntegrallyClosed());
    return;
  }
  res.holds = true;
  res.reason = "fin gen and idv closed in q(H)";
}


struct Quad {
  mpq_class a, b; // a + b*sqrt(D)
  bool operator<(const Quad& y) const {
    return std::tie(a, b) < std::tie(y.a, y.b);
  }
};

struct Field {
  mpq_class D;

  Quad mul(const Quad& x, const Quad& y) const {
    return {x.a*y.a + D*x.b*y.b, x.a*y.b + x.b*y.a};
  }

  Quad conjugate(const Quad& x) const {
    return {x.a, -x.b};
  }

  int sign(const Quad& x) const {
    int sa = sgn(x.a), sb = sgn(x.b);
    if(sa == 0 || sb == 0 || sa == sb) {
      return sa != 0 ? sa : sb;
    }
    int larger = cmp(x.a*x.a, D*x.b*x.b);
    return larger == 0 ? 0 : (larger > 0 ? sa : sb);
  }

  bool integral(const Quad& x) const {
    mpq_class trace = 2*x.a, norm = x.a*x.a - D*x.b*x.b;
    return trace.get_den() == 1 && norm.get_den() == 1;
  }

  // |a - b*sqrt(D)| <= a + b*sqrt(D) for pos x
  bool weak_perron(const Quad& x) const {
    return x.a >= 0 && x.b >= 0;
  }

  string text(const Quad& x) const {
    std::ostringstream out;
    mpq_class mag = abs(x.b);
    out << x.a.get_str() << (sgn(x.b) < 0 ? " - " : " + ") << mag.get_str()
        << "*sqrt(" << D.get_str() << ')';
    return out.str();
  }
};

mpz_class sqfree(const mpz_class& n) {
  fmpz_t val;
  fmpz_init(val);
  fmpz_set_mpz(val, n.get_mpz_t());
  fmpz_factor_t facs;
  fmpz_factor_init(facs);
  fmpz_factor(facs, val);

  mpz_class part = 1;
  for(slong i = 0; i < facs->num; i++) {
    if(facs->exp[i] % 2) {
      mpz_class p;
      fmpz_get_mpz(p.get_mpz_t(), facs->p + i);
      part *= p;
    }
  }

  fmpz_factor_clear(facs);
  fmpz_clear(val);
  return part;
}
mpq_class decimal(const string& text) {
  size_t dot = text.find('.');
  if(dot == string::npos) {
    return rational(text);
  }

  mpz_class scale = 1;
  for(size_t i = dot + 1; i < text.size(); i++) {
    scale *= 10;
  }
  return rational(text.substr(0, dot) + text.substr(dot + 1)) / scale;
}

Quad parse_root(const Generator& given, mpz_class& D) {
  vector<mpq_class> p = read_poly(given.poly);
  if(p.size() < 2 || p.size() > 3) {
    throw invalid_argument("with several generators every polynomial must have degree 1 or 2");
  }

  if(p.size() == 2) {
    mpq_class val = -p[0] / p[1];
    if(val <= 0) {
      throw invalid_argument("generators must be positive");
    }
    D = 0;
    return {val, 0};
  }

  mpz_class lcm;
  vector<mpz_class> c = clear_denoms(p, lcm);
  if(c[2] < 0) {
    for(auto& x : c) {
      x = -x;
    }
  }

  mpz_class disc = c[1]*c[1] - 4*c[0]*c[2];
  if(disc < 0) {
    throw invalid_argument("root is not real");
  }
  D = disc == 0 ? 1 : sqfree(disc);
  if(D == 1) {
    throw invalid_argument("polynomial has rational roots");
  }

  mpq_class centre = mpq_class(-c[1]) / (2*c[2]), spread = mpq_class(sqrt(disc / D)) / (2*c[2]);
  Field field{mpq_class(D)};
  Quad low{centre, -spread}, high{centre, spread};

  if(!given.root.empty()) {
    Quad picked = decimal(given.root) > centre ? high : low;
    if(field.sign(picked) <= 0) {
      throw invalid_argument("generators must be positive");
    }
    return picked;
  }

  if(field.sign(low) > 0) {
    throw invalid_argument("two positive roots, pick one with @");
  }
  if(field.sign(high) <= 0) {
    throw invalid_argument("generators must be positive");
  }
  return high;
}

Points to_lattice(const vector<Quad>& vals, mpz_class& lcm) {
  vector<mpq_class> flat;
  for(const auto& x : vals) {
    flat.push_back(x.a);
    flat.push_back(x.b);
  }
  vector<mpz_class> cleared = clear_denoms(flat, lcm);

  Points pts;
  for(size_t i = 0; i < vals.size(); ++i) {
    pts.push_back({cleared[2*i], cleared[2*i + 1]});
  }
  return pts;
}

bool contains(const vector<Quad>& gens, const Quad& target) {
  vector<Quad> together = gens;
  together.push_back(target);
  mpz_class lcm;
  Points pts = to_lattice(together, lcm);
  Point goal = std::move(pts.back());
  pts.pop_back();
  return contains(pts, goal);
}

// search at end of section 4, which was G holds the monos of deg at most m, and m grows until G is closed under multiplication by every alpha_i
vector<Quad> gen_set(const Field& field, const vector<Quad>& alphas) {
  const Quad one{1, 0}; // thm 4.2 is what makes this terminate
  std::set<Quad> G{one}, deg{one}; 

  for(;;) {
    vector<Quad> monos(G.begin(), G.end());
    bool closed = true;
    for(const auto& g : monos) {
      for(const auto& alpha : alphas) {
        closed = closed && contains(monos, field.mul(g, alpha));
      }
    }
    if(closed) {
      return monos;
    }

    std::set<Quad> next;
    for(const auto& x : deg) {
      for(const auto& alpha : alphas) {
        next.insert(field.mul(x, alpha));
      }
    }
    deg = std::move(next);
    G.insert(deg.begin(), deg.end());
  }
}

// fin gen classification, thm 4.2
bool quad_check(const Field& field, const vector<Quad>& alphas, string& reason) {
  bool neg_conj = false;

  for(const auto& alpha : alphas) {
    if(!field.integral(alpha)) {
      reason = field.text(alpha) + " is not an algebraic integer";
      return false;
    }
    if(!field.weak_perron(alpha)) {
      reason = field.text(alpha) + " is not weak Perron";
      return false;
    }
    neg_conj = neg_conj || field.sign(field.conjugate(alpha)) < 0;
  }

  if(!neg_conj) {
    reason = "no generator has a negative conjugate";
    return false;
  }

  reason = "every generator is a weak Perron algebraic integer and one has a negative conjugate";
  return true;
}


// how many roots are real and positive; bigger says some conjugate beats alpha in modulus
// Roots inspect(const vector<mpz_class>& coeffs);
int inspect(const vector<mpz_class>& coeffs, bool& bigger) {
  fmpz_poly_t poly;
  fmpz_poly_init(poly);
  fmpz_t val;
  fmpz_init(val);
  for(size_t i = 0; i < coeffs.size(); i++) {
    fmpz_set_mpz(val, coeffs[i].get_mpz_t());
    fmpz_poly_set_coeff_fmpz(poly, i, val);
  }
  fmpz_clear(val);

  fmpz_poly_factor_t facs;
  fmpz_poly_factor_init(facs);
  fmpz_poly_factor(facs, poly);
  bool irred = facs->num == 1 && facs->exp[0] == 1;
  fmpz_poly_factor_clear(facs);
  if(!irred) {
    fmpz_poly_clear(poly);
    throw invalid_argument("polynomial is reducible");
  }

  slong deg = fmpz_poly_degree(poly);
  qqbar_ptr roots = _qqbar_vec_init(deg);
  qqbar_roots_fmpz_poly(roots, poly, QQBAR_ROOTS_IRREDUCIBLE);

  int pos = 0;
  qqbar_ptr alpha = nullptr;
  for(slong i = 0; i < deg; ++i) {
    if(qqbar_is_real(roots + i) && qqbar_sgn_re(roots + i) > 0) {
      alpha = roots + i;
      ++pos;
    }
  }
  bigger = false;
  if(pos == 1) {
    for(slong i = 0; i < deg; ++i) {
      bigger |= qqbar_cmpabs(roots + i, alpha) > 0;
    }
  }

  _qqbar_vec_clear(roots, deg);
  fmpz_poly_clear(poly);
  return pos;
}

// cor 3.13
bool cyclic_check(const vector<string>& coeffs, vector<mpz_class>& monic, string& reason) {
  vector<mpq_class> poly = read_poly(coeffs);
  if(poly.size() < 2) {
    throw invalid_argument("need degree at least one");
  }

  mpz_class lcm;
  bool bigger = false;
  int pos = inspect(clear_denoms(poly, lcm), bigger);

  for(const auto& c : poly) {
    mpq_class scaled = c / poly.back();
    if(scaled.get_den() != 1) {
      reason = "alpha is not an algebraic integer";
      return false;
    }
    monic.push_back(scaled.get_num());
  }

  if(pos != 1) {
    reason = "alpha is not the only positive root";
    return false;
  }
  if(bigger) {
    reason = "a conjugate has modulus larger than alpha";
    return false;
  }
  reason = "alpha is a weak Perron algebraic integer and the only positive root";
  return true;
}

// section 5.2, in the basis 1, alpha, ..., alpha^(d-1): v_(n+d) = -sum a_i v_(n+i), and k is the first power lying in the span of the earlier ones
// monic really has to be monic and minimal, cyclic_check guarantees both
Points cyclic_gens(const vector<mpz_class>& monic) {
  size_t d = monic.size() - 1;
  Points powers;
  for(size_t i = 0; i < d; ++i) {
    Point unit(d);
    unit[i] = 1;
    powers.push_back(std::move(unit));
  }


  for(size_t s = d;; ++s) {
    Point next(d);
    for(size_t i = 0; i < d; ++i) {
      for(size_t j = 0; j < d; ++j) {
        next[j] -= monic[i]*powers[s-d+i][j];
      }
    }

    if(contains(powers, next)) {
      return powers;
    }
    powers.push_back(std::move(next));
  }
}


Result decide_quad(const vector<Quad>& alphas, const mpz_class& D, bool krull) {
  Field field{mpq_class(D)};

  Result res;
  res.property = krull ? "Krull" : "Mori";
  if(!quad_check(field, alphas, res.reason)) {
    res.reason += "; not finitely generated";
    return res;
  }

  vector<Quad> gens = gen_set(field, alphas);
  for(const auto& x : gens) {
    res.gens.push_back(field.text(x));
  }
  if(!krull) {
    res.holds = true;
    res.reason += "; Mori is the same as finite generation here";
    return res;
  }

  mpz_class common;
  Points lattice = to_lattice(gens, common);
  krull_step(res, lattice, [&](const Point& point) {
    return field.text({mpq_class(point[0]) / common, mpq_class(point[1]) / common});
  });
  return res;
}

Result decide_cyclic(const vector<string>& poly, bool krull) {
  vector<mpz_class> monic;

  Result res;
  res.property = krull ? "Krull" : "Mori";
  if(!cyclic_check(poly, monic, res.reason)) {
    res.reason += "; not finitely generated";
    return res;
  }

  Points gens = cyclic_gens(monic);
  for(size_t i = 0; i < gens.size(); ++i) {
    string out = "alpha^" + std::to_string(i) + " = [";
    for(size_t j = 0; j < gens[i].size(); j++) {
      if(j) {
        out += ", ";
      }
      out += gens[i][j].get_str();
    }
    res.gens.push_back(out + "]");
  }
  if(!krull) {
    res.holds = true;
    res.reason += "; Mori is the same as finite generation";
    return res;
  }

  krull_step(res, gens, [](const Point& point) {
    string out = "[";
    for(size_t i = 0; i < point.size(); i++) {
      if(i) {
        out += ", ";
      }
      out += point[i].get_str();
    }
    return out + "]";
  });
  return res;
}

// one gen is the cyclic case and works in any deg
// several only work when they all live in one quadratic field, which is as far as thm 4.2 reaches
// TODO deg 3 and up is open, so there is nothing to call here yet
Result decide(const vector<Generator>& given, bool krull) {
  if(given.empty()) {
    throw invalid_argument("need at least one generator");
  }
  if(given.size() == 1) {
    if(!given.front().root.empty()) {
      throw invalid_argument("@ does nothing with one generator");
    }
    return decide_cyclic(given.front().poly, krull);
  }

  vector<Quad> alphas;
  mpz_class D = 0;
  for(const auto& gen : given) {
    mpz_class d;
    alphas.push_back(parse_root(gen, d));
    if(d != 0) {
      if(D != 0 && D != d) {
        throw invalid_argument("generators lie in different quadratic fields");
      }
      D = d;
    }
  }
  if(D == 0) {
    throw invalid_argument("generators only span Q");
  }

  return decide_quad(alphas, D, krull);
}

} // namespace

Result mori(const vector<Generator>& gens) {
  return decide(gens, false);
}
Result krull(const vector<Generator>& gens) {
  return decide(gens, true);
}

} // namespace morikrull
