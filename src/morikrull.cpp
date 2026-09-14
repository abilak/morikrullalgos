#include "morikrull.hpp"

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace morikrull {
namespace {

using Int = boost::multiprecision::cpp_int;
using Vec = std::vector<Int>;

Int abs(Int n) { return n < 0 ? -n : n; }

Int gcd(Int a, Int b)
{
  a = abs(a);
  b = abs(b);
  while (b != 0) {
    Int r = a % b;
    a = b;
    b = r;
  }
  return a;
}

Int lcm(const Int& a, const Int& b)
{
  if (a == 0 || b == 0) return 0;
  return abs(a / gcd(a, b) * b);
}

struct Rat {
  Int n = 0;
  Int d = 1;

  Rat() = default;
  Rat(long long value) : n(value) {}
  Rat(Int numerator, Int denominator = 1) : n(std::move(numerator)), d(std::move(denominator))
  {
    normalize();
  }

  explicit Rat(const std::string& text)
  {
    auto slash = text.find('/');
    if (slash == std::string::npos) n = Int(text);
    else {
      n = Int(text.substr(0, slash));
      d = Int(text.substr(slash + 1));
    }
    normalize();
  }

  void normalize()
  {
    if (d == 0) throw std::invalid_argument("a rational denominator is zero");
    if (d < 0) {
      n = -n;
      d = -d;
    }
    Int g = gcd(n, d);
    if (g != 0) {
      n /= g;
      d /= g;
    }
  }

  std::string str() const
  {
    if (d == 1) return n.convert_to<std::string>();
    return n.convert_to<std::string>() + "/" + d.convert_to<std::string>();
  }
};

Rat operator+(const Rat& a, const Rat& b) { return {a.n * b.d + b.n * a.d, a.d * b.d}; }
Rat operator-(const Rat& a, const Rat& b) { return {a.n * b.d - b.n * a.d, a.d * b.d}; }
Rat operator-(const Rat& a) { return {-a.n, a.d}; }
Rat operator*(const Rat& a, const Rat& b) { return {a.n * b.n, a.d * b.d}; }
Rat operator/(const Rat& a, const Rat& b)
{
  if (b.n == 0) throw std::domain_error("division by zero");
  return {a.n * b.d, a.d * b.n};
}
bool operator==(const Rat& a, const Rat& b) { return a.n == b.n && a.d == b.d; }
bool operator!=(const Rat& a, const Rat& b) { return !(a == b); }
bool operator<(const Rat& a, const Rat& b) { return a.n * b.d < b.n * a.d; }
bool operator>(const Rat& a, const Rat& b) { return b < a; }

int sign(const Rat& x) { return (x.n > 0) - (x.n < 0); }

struct RatLess {
  bool operator()(const Rat& a, const Rat& b) const { return a < b; }
};

struct VecLess {
  bool operator()(const Vec& a, const Vec& b) const
  {
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
  }
};

Int dot(const Vec& a, const Vec& b)
{
  Int total = 0;
  for (size_t i = 0; i < a.size(); ++i) total += a[i] * b[i];
  return total;
}

bool is_zero(const Vec& a)
{
  return std::all_of(a.begin(), a.end(), [](const Int& x) { return x == 0; });
}

Vec subtract_multiple(const Vec& a, const Vec& b, const Int& n)
{
  Vec out(a.size());
  for (size_t i = 0; i < a.size(); ++i) out[i] = a[i] - n * b[i];
  return out;
}

Int determinant(std::vector<Vec> rows)
{
  const size_t n = rows.size();
  if (n == 0) return 1;
  for (const auto& row : rows)
    if (row.size() != n) throw std::invalid_argument("determinant of a nonsquare matrix");

  Int previous = 1;
  int parity = 1;
  for (size_t k = 0; k + 1 < n; ++k) {
    size_t pivot = k;
    while (pivot < n && rows[pivot][k] == 0) ++pivot;
    if (pivot == n) return 0;
    if (pivot != k) {
      std::swap(rows[pivot], rows[k]);
      parity = -parity;
    }
    Int value = rows[k][k];
    for (size_t i = k + 1; i < n; ++i) {
      for (size_t j = k + 1; j < n; ++j)
        rows[i][j] = (rows[i][j] * value - rows[i][k] * rows[k][j]) / previous;
    }
    previous = value;
  }
  return parity * rows[n - 1][n - 1];
}

Int det_columns(const std::vector<Vec>& columns)
{
  if (columns.empty()) return 1;
  std::vector<Vec> rows(columns.size(), Vec(columns.size()));
  for (size_t j = 0; j < columns.size(); ++j)
    for (size_t i = 0; i < columns.size(); ++i) rows[i][j] = columns[j][i];
  return determinant(std::move(rows));
}

int matrix_rank(const std::vector<Vec>& columns)
{
  if (columns.empty()) return 0;
  size_t m = columns.front().size();
  size_t n = columns.size();
  std::vector<std::vector<Rat>> a(m, std::vector<Rat>(n));
  for (size_t j = 0; j < n; ++j)
    for (size_t i = 0; i < m; ++i) a[i][j] = Rat(columns[j][i]);

  size_t row = 0;
  for (size_t col = 0; col < n && row < m; ++col) {
    size_t pivot = row;
    while (pivot < m && a[pivot][col].n == 0) ++pivot;
    if (pivot == m) continue;
    std::swap(a[pivot], a[row]);
    Rat p = a[row][col];
    for (size_t j = col; j < n; ++j) a[row][j] = a[row][j] / p;
    for (size_t i = 0; i < m; ++i) {
      if (i == row || a[i][col].n == 0) continue;
      Rat q = a[i][col];
      for (size_t j = col; j < n; ++j) a[i][j] = a[i][j] - q * a[row][j];
    }
    ++row;
  }
  return static_cast<int>(row);
}

std::vector<size_t> independent_rows(const std::vector<Vec>& columns)
{
  const size_t r = columns.size();
  const size_t d = columns.front().size();
  std::vector<size_t> pick;
  std::function<bool(size_t)> search = [&](size_t start) {
    if (pick.size() == r) {
      std::vector<Vec> square(r, Vec(r));
      for (size_t i = 0; i < r; ++i)
        for (size_t j = 0; j < r; ++j) square[i][j] = columns[j][pick[i]];
      return determinant(std::move(square)) != 0;
    }
    for (size_t i = start; i < d; ++i) {
      pick.push_back(i);
      if (search(i + 1)) return true;
      pick.pop_back();
    }
    return false;
  };
  if (!search(0)) throw std::logic_error("could not find independent matrix rows");
  return pick;
}

bool solve_nonnegative(const std::vector<Vec>& columns, const Vec& target)
{
  if (columns.empty()) return is_zero(target);
  auto rows = independent_rows(columns);
  size_t r = columns.size();
  std::vector<Vec> square(r, Vec(r));
  for (size_t i = 0; i < r; ++i)
    for (size_t j = 0; j < r; ++j) square[i][j] = columns[j][rows[i]];

  Int den = determinant(square);
  std::vector<Int> coefficient(r);
  for (size_t j = 0; j < r; ++j) {
    auto changed = square;
    for (size_t i = 0; i < r; ++i) changed[i][j] = target[rows[i]];
    Int num = determinant(std::move(changed));
    if (num % den != 0) return false;
    coefficient[j] = num / den;
    if (coefficient[j] < 0) return false;
  }

  for (size_t i = 0; i < target.size(); ++i) {
    Int sum = 0;
    for (size_t j = 0; j < r; ++j) sum += coefficient[j] * columns[j][i];
    if (sum != target[i]) return false;
  }
  return true;
}

bool affine_contains(std::vector<Vec> generators, const Vec& target, const Vec& weight)
{
  generators.erase(std::remove_if(generators.begin(), generators.end(), is_zero), generators.end());
  if (generators.empty()) return is_zero(target);
  if (target.size() != weight.size()) throw std::invalid_argument("weight dimension mismatch");

  Int target_weight = dot(weight, target);
  if (target_weight < 0) return false;
  for (const auto& g : generators)
    if (dot(weight, g) <= 0) throw std::logic_error("the cone weight is not strictly positive");

  std::vector<Vec> basis;
  std::vector<Vec> rest;
  for (const auto& g : generators) {
    auto trial = basis;
    trial.push_back(g);
    if (matrix_rank(trial) > matrix_rank(basis)) basis.push_back(g);
    else rest.push_back(g);
  }

  std::sort(rest.begin(), rest.end(), [&](const Vec& a, const Vec& b) {
    return dot(weight, a) > dot(weight, b);
  });

  std::function<bool(size_t, const Vec&)> visit = [&](size_t at, const Vec& residual) {
    Int remaining_weight = dot(weight, residual);
    if (remaining_weight < 0) return false;
    if (at == rest.size()) return solve_nonnegative(basis, residual);

    Int step = dot(weight, rest[at]);
    Int limit = remaining_weight / step;
    for (Int c = 0; c <= limit; ++c)
      if (visit(at + 1, subtract_multiple(residual, rest[at], c))) return true;
    return false;
  };
  return visit(0, target);
}

std::string vector_string(const Vec& v)
{
  std::ostringstream out;
  out << '[';
  for (size_t i = 0; i < v.size(); ++i) {
    if (i) out << ", ";
    out << v[i];
  }
  out << ']';
  return out.str();
}

Int integer_sqrt(const Int& n)
{
  if (n < 0) throw std::invalid_argument("square root of a negative integer");
  if (n < 2) return n;
  unsigned bits = boost::multiprecision::msb(n) + 1;
  Int x = Int(1) << ((bits + 1) / 2);
  for (;;) {
    Int y = (x + n / x) / 2;
    if (y >= x) return x;
    x = y;
  }
}

struct Quadratic {
  Rat a;
  Rat b;
  long long D = 0;
};

Quadratic operator*(const Quadratic& x, const Quadratic& y)
{
  Rat d(x.D);
  return {x.a * y.a + d * x.b * y.b, x.a * y.b + x.b * y.a, x.D};
}

struct QuadraticLess {
  bool operator()(const Quadratic& x, const Quadratic& y) const
  {
    if (x.a != y.a) return RatLess{}(x.a, y.a);
    return RatLess{}(x.b, y.b);
  }
};

int sign(const Quadratic& x)
{
  int sa = sign(x.a);
  int sb = sign(x.b);
  if (sa == 0) return sb;
  if (sb == 0 || sa == sb) return sa;
  Rat left = x.a * x.a;
  Rat right = Rat(x.D) * x.b * x.b;
  if (left == right) return 0;
  if (sa > 0) return left > right ? 1 : -1;
  return right > left ? 1 : -1;
}

Quadratic conjugate(const Quadratic& x) { return {x.a, -x.b, x.D}; }

bool algebraic_integer(const Quadratic& x)
{
  if (x.b.n == 0) return x.a.d == 1;
  Rat trace = Rat(2) * x.a;
  Rat norm = x.a * x.a - Rat(x.D) * x.b * x.b;
  return trace.d == 1 && norm.d == 1;
}

bool weak_perron(const Quadratic& x)
{
  if (!algebraic_integer(x) || sign(x) <= 0) return false;
  Quadratic c = conjugate(x);
  if (sign(c) < 0) c = {-c.a, -c.b, c.D};
  Quadratic difference{x.a - c.a, x.b - c.b, x.D};
  return sign(difference) >= 0;
}

std::string quadratic_string(const Quadratic& x)
{
  std::ostringstream out;
  out << x.a.str();
  if (x.b.n >= 0) out << " + ";
  else out << " - ";
  Rat magnitude = x.b.n >= 0 ? x.b : -x.b;
  out << magnitude.str() << "*sqrt(" << x.D << ')';
  return out.str();
}

Vec quadratic_weight(long long D, const std::vector<Vec>& values)
{
  Int q = 1;
  for (int attempt = 0; attempt < 500; ++attempt) {
    Int p = integer_sqrt(Int(D) * q * q);
    Vec weight{q, p};
    bool good = true;
    for (const auto& value : values) good = good && dot(weight, value) > 0;
    if (good) return weight;
    q *= 10;
  }
  throw std::runtime_error("could not construct a positive quadratic cone functional");
}

std::pair<std::vector<Vec>, Vec> scale_quadratics(
  const std::vector<Quadratic>& generators, const Quadratic& target)
{
  Int common = 1;
  for (const auto& g : generators) {
    common = lcm(common, g.a.d);
    common = lcm(common, g.b.d);
  }
  common = lcm(common, target.a.d);
  common = lcm(common, target.b.d);

  std::vector<Vec> scaled;
  for (const auto& g : generators)
    scaled.push_back({g.a.n * (common / g.a.d), g.b.n * (common / g.b.d)});
  Vec t{target.a.n * (common / target.a.d), target.b.n * (common / target.b.d)};
  return {scaled, t};
}

bool quadratic_contains(const std::vector<Quadratic>& generators, const Quadratic& target)
{
  auto [scaled, t] = scale_quadratics(generators, target);
  auto values = scaled;
  values.push_back(t);
  return affine_contains(std::move(scaled), t, quadratic_weight(target.D, values));
}

void monomials_of_degree(const std::vector<Quadratic>& alpha, size_t at, int degree,
  const Quadratic& value, std::set<Quadratic, QuadraticLess>& out)
{
  if (at == alpha.size()) {
    if (degree == 0) out.insert(value);
    return;
  }
  Quadratic power{Rat(1), Rat(0), value.D};
  for (int e = 0; e <= degree; ++e) {
    monomials_of_degree(alpha, at + 1, degree - e, value * power, out);
    power = power * alpha[at];
  }
}

std::vector<Quadratic> quadratic_generating_set(const std::vector<Quadratic>& alpha)
{
  long long D = alpha.front().D;
  std::vector<Quadratic> basis{{Rat(1), Rat(0), D}};

  for (int degree = 1;; ++degree) {
    std::set<Quadratic, QuadraticLess> level;
    monomials_of_degree(alpha, 0, degree, {Rat(1), Rat(0), D}, level);
    for (const auto& x : level)
      if (!quadratic_contains(basis, x)) basis.push_back(x);

    bool closed = true;
    for (size_t i = 0; i < basis.size() && closed; ++i)
      for (const auto& a : alpha)
        if (!quadratic_contains(basis, basis[i] * a)) {
          closed = false;
          break;
        }
    if (closed) return basis;
  }
}

using Poly = std::vector<Rat>;

void trim(Poly& p)
{
  while (p.size() > 1 && p.back().n == 0) p.pop_back();
}

Poly derivative(const Poly& p)
{
  if (p.size() <= 1) return {Rat(0)};
  Poly out(p.size() - 1);
  for (size_t i = 1; i < p.size(); ++i) out[i - 1] = Rat(Int(i)) * p[i];
  trim(out);
  return out;
}

std::pair<Poly, Poly> divide(const Poly& numerator, const Poly& denominator)
{
  if (denominator.size() == 1 && denominator[0].n == 0)
    throw std::domain_error("polynomial division by zero");
  Poly remainder = numerator;
  trim(remainder);
  Poly quotient(remainder.size() >= denominator.size()
    ? remainder.size() - denominator.size() + 1 : 1, Rat(0));

  while (!(remainder.size() == 1 && remainder[0].n == 0)
      && remainder.size() >= denominator.size()) {
    size_t shift = remainder.size() - denominator.size();
    Rat q = remainder.back() / denominator.back();
    quotient[shift] = quotient[shift] + q;
    for (size_t i = 0; i < denominator.size(); ++i)
      remainder[i + shift] = remainder[i + shift] - q * denominator[i];
    trim(remainder);
  }
  trim(quotient);
  return {quotient, remainder};
}

Rat evaluate(const Poly& p, const Rat& x)
{
  Rat value(0);
  for (auto i = p.rbegin(); i != p.rend(); ++i) value = value * x + *i;
  return value;
}

std::vector<Poly> sturm_sequence(Poly p)
{
  trim(p);
  std::vector<Poly> sequence{p, derivative(p)};
  while (!(sequence.back().size() == 1 && sequence.back()[0].n == 0)) {
    auto remainder = divide(sequence[sequence.size() - 2], sequence.back()).second;
    for (Rat& x : remainder) x = -x;
    trim(remainder);
    if (remainder.size() == 1 && remainder[0].n == 0) break;
    sequence.push_back(std::move(remainder));
  }
  return sequence;
}

int variations_at(const std::vector<Poly>& sturm, const Rat& x)
{
  int changes = 0;
  int previous = 0;
  for (const auto& p : sturm) {
    int current = sign(evaluate(p, x));
    if (current == 0) continue;
    if (previous != 0 && current != previous) ++changes;
    previous = current;
  }
  return changes;
}

int roots_between(const std::vector<Poly>& sturm, const Rat& left, const Rat& right)
{
  return variations_at(sturm, left) - variations_at(sturm, right);
}

Int root_bound(const Poly& monic)
{
  Rat maximum(0);
  for (size_t i = 0; i + 1 < monic.size(); ++i) {
    Rat value = monic[i].n < 0 ? -monic[i] : monic[i];
    if (value > maximum) maximum = value;
  }
  Int ceiling = maximum.n / maximum.d;
  if (maximum.n % maximum.d != 0) ++ceiling;
  return ceiling + 1;
}

struct Interval {
  Rat lo;
  Rat hi;
};

Interval isolate_positive_root(const Poly& polynomial, const std::vector<Poly>& sturm)
{
  Interval interval{Rat(0), Rat(root_bound(polynomial))};
  if (roots_between(sturm, interval.lo, interval.hi) != 1)
    throw std::logic_error("expected one positive root");

  for (int step = 0; step < 400; ++step) {
    Rat middle = (interval.lo + interval.hi) / Rat(2);
    Rat value = evaluate(polynomial, middle);
    if (value.n == 0) return {middle, middle};
    int left = roots_between(sturm, interval.lo, middle);
    if (left == 1) interval.hi = middle;
    else interval.lo = middle;
  }
  return interval;
}

Poly monic(Poly p)
{
  trim(p);
  if (p.back().n == 0) throw std::invalid_argument("zero polynomial");
  Rat leading = p.back();
  for (Rat& x : p) x = x / leading;
  return p;
}

Poly polynomial_gcd(Poly a, Poly b)
{
  trim(a);
  trim(b);
  while (!(b.size() == 1 && b[0].n == 0)) {
    Poly remainder = divide(a, b).second;
    a = std::move(b);
    b = std::move(remainder);
  }
  return monic(std::move(a));
}

Poly squarefree_part(const Poly& p)
{
  Poly common = polynomial_gcd(p, derivative(p));
  return monic(divide(p, common).first);
}

using Matrix = std::vector<Vec>;

Matrix multiply(const Matrix& a, const Matrix& b)
{
  size_t n = a.size();
  Matrix product(n, Vec(n));
  for (size_t i = 0; i < n; ++i)
    for (size_t k = 0; k < n; ++k)
      if (a[i][k] != 0)
        for (size_t j = 0; j < n; ++j) product[i][j] += a[i][k] * b[k][j];
  return product;
}

Int trace(const Matrix& a)
{
  Int value = 0;
  for (size_t i = 0; i < a.size(); ++i) value += a[i][i];
  return value;
}

Poly pair_product_polynomial(const std::vector<Int>& polynomial)
{
  size_t d = polynomial.size() - 1;
  Matrix companion(d, Vec(d));
  for (size_t i = 1; i < d; ++i) companion[i][i - 1] = 1;
  for (size_t i = 0; i < d; ++i) companion[i][d - 1] = -polynomial[i];

  size_t n = d * d;
  std::vector<Int> power_traces(n + 1);
  Matrix power = companion;
  for (size_t k = 1; k <= n; ++k) {
    Int t = trace(power);
    power_traces[k] = t * t;
    if (k != n) power = multiply(power, companion);
  }

  // Newton identities for the characteristic polynomial of C tensor C. Its
  // eigenvalues are exactly beta_i*beta_j for the conjugates beta_i of alpha.
  std::vector<Int> characteristic(n + 1);
  characteristic[0] = 1;
  for (size_t k = 1; k <= n; ++k) {
    Int sum = 0;
    for (size_t i = 1; i <= k; ++i) sum += characteristic[k - i] * power_traces[i];
    if (sum % Int(k) != 0) throw std::logic_error("Newton identity was not integral");
    characteristic[k] = -sum / Int(k);
  }

  Poly result(n + 1);
  for (size_t i = 0; i <= n; ++i) result[i] = Rat(characteristic[n - i]);
  return squarefree_part(result);
}

Poly compose_square(const Poly& p)
{
  Poly result(2 * (p.size() - 1) + 1, Rat(0));
  for (size_t i = 0; i < p.size(); ++i) result[2 * i] = p[i];
  trim(result);
  return result;
}

void refine_root_interval(const Poly& p, const std::vector<Poly>& sturm, Interval& interval)
{
  if (interval.lo == interval.hi) return;
  Rat middle = (interval.lo + interval.hi) / Rat(2);
  if (evaluate(p, middle).n == 0) {
    interval = {middle, middle};
    return;
  }
  if (roots_between(sturm, interval.lo, middle) == 1) interval.hi = middle;
  else interval.lo = middle;
}

int sign_at_alpha_squared(const Poly& q, const Poly& minimal,
  const std::vector<Poly>& minimal_sturm, Interval& alpha)
{
  Poly composed = compose_square(q);
  Poly remainder = divide(composed, minimal).second;
  trim(remainder);
  if (remainder.size() == 1 && remainder[0].n == 0) return 0;
  if (alpha.lo == alpha.hi) return sign(evaluate(composed, alpha.lo));

  Poly roots = squarefree_part(composed);
  auto root_sturm = sturm_sequence(roots);
  for (int attempt = 0; attempt < 10000; ++attempt) {
    Rat left_value = evaluate(composed, alpha.lo);
    Rat right_value = evaluate(composed, alpha.hi);
    if (left_value.n != 0 && right_value.n != 0
        && roots_between(root_sturm, alpha.lo, alpha.hi) == 0) {
      return sign(evaluate(composed, (alpha.lo + alpha.hi) / Rat(2)));
    }
    refine_root_interval(minimal, minimal_sturm, alpha);
  }
  throw std::runtime_error("could not separate an algebraic sign");
}

int variations_at_alpha_squared(const std::vector<Poly>& sequence,
  const Poly& minimal, const std::vector<Poly>& minimal_sturm, Interval& alpha)
{
  int changes = 0;
  int previous = 0;
  for (const auto& p : sequence) {
    int current = sign_at_alpha_squared(p, minimal, minimal_sturm, alpha);
    if (current == 0) continue;
    if (previous != 0 && current != previous) ++changes;
    previous = current;
  }
  return changes;
}

int variations_at_positive_infinity(const std::vector<Poly>& sequence)
{
  int changes = 0;
  int previous = 0;
  for (const auto& p : sequence) {
    int current = sign(p.back());
    if (current == 0) continue;
    if (previous != 0 && current != previous) ++changes;
    previous = current;
  }
  return changes;
}

bool conjugate_outside_alpha(const Poly& minimal, const std::vector<Int>& integer_polynomial,
  const std::vector<Poly>& minimal_sturm, Interval alpha)
{
  Poly products = pair_product_polynomial(integer_polynomial);
  auto product_sturm = sturm_sequence(products);
  int at_alpha = variations_at_alpha_squared(product_sturm, minimal, minimal_sturm, alpha);
  int at_infinity = variations_at_positive_infinity(product_sturm);
  return at_alpha - at_infinity > 0;
}

struct CyclicData {
  bool finite = false;
  std::string reason;
  Poly polynomial;
  std::vector<Int> integer_polynomial;
  Interval alpha;
};

CyclicData cyclic_finite_generation(const std::vector<std::string>& input)
{
  if (input.size() < 2) throw std::invalid_argument("the polynomial must have positive degree");
  Poly p;
  for (const auto& coefficient : input) p.emplace_back(coefficient);
  trim(p);
  if (p.size() < 2 || p.back().n == 0) throw std::invalid_argument("invalid polynomial");

  Rat leading = p.back();
  for (Rat& coefficient : p) coefficient = coefficient / leading;

  CyclicData data;
  data.polynomial = p;
  for (const Rat& coefficient : p) {
    if (coefficient.d != 1) {
      data.reason = "alpha is not an algebraic integer (its monic minimal polynomial is not integral)";
      return data;
    }
    data.integer_polynomial.push_back(coefficient.n);
  }

  if (p[0].n == 0 && p.size() > 2)
    throw std::invalid_argument("the supplied polynomial is not a minimal polynomial");
  auto sturm = sturm_sequence(p);
  Int bound = root_bound(p);
  int positive_roots = roots_between(sturm, Rat(0), Rat(bound));
  if (positive_roots != 1) {
    data.reason = "the minimal polynomial does not have exactly one positive root";
    return data;
  }
  data.alpha = isolate_positive_root(p, sturm);
  if (conjugate_outside_alpha(p, data.integer_polynomial, sturm, data.alpha)) {
    data.reason = "a conjugate has modulus larger than alpha";
    return data;
  }

  data.finite = true;
  data.reason = "alpha is a weak Perron algebraic integer and is the only positive conjugate";
  return data;
}

Vec cyclic_weight(const Poly& polynomial, const Interval& original, const std::vector<Vec>& values)
{
  if (values.empty()) return {};
  Interval interval = original;
  size_t d = values.front().size();
  auto sturm = sturm_sequence(polynomial);

  for (int attempt = 0; attempt < 2000; ++attempt) {
    Rat r = (interval.lo + interval.hi) / Rat(2);
    std::vector<Rat> powers(d, Rat(1));
    for (size_t i = 1; i < d; ++i) powers[i] = powers[i - 1] * r;

    bool good = true;
    for (const auto& v : values) {
      Rat sum(0);
      for (size_t i = 0; i < d; ++i) sum = sum + Rat(v[i]) * powers[i];
      good = good && sum.n > 0;
    }
    if (good) {
      Int denominator = 1;
      for (const auto& power : powers) denominator = lcm(denominator, power.d);
      Vec weight(d);
      for (size_t i = 0; i < d; ++i)
        weight[i] = powers[i].n * (denominator / powers[i].d);
      return weight;
    }

    if (interval.lo == interval.hi) break;
    if (evaluate(polynomial, r).n == 0) break;
    if (roots_between(sturm, interval.lo, r) == 1) interval.hi = r;
    else interval.lo = r;
  }
  throw std::runtime_error("could not construct a positive cyclic cone functional");
}

std::vector<Vec> cyclic_generating_vectors(const CyclicData& data)
{
  size_t d = data.integer_polynomial.size() - 1;
  std::vector<Vec> powers;
  for (size_t i = 0; i < d; ++i) {
    Vec unit(d);
    unit[i] = 1;
    powers.push_back(std::move(unit));
  }

  for (size_t s = d;; ++s) {
    Vec next(d);
    for (size_t i = 0; i < d; ++i)
      for (size_t j = 0; j < d; ++j)
        next[j] -= data.integer_polynomial[i] * powers[s - d + i][j];

    auto values = powers;
    values.push_back(next);
    Vec weight = cyclic_weight(data.polynomial, data.alpha, values);
    if (affine_contains(powers, next, weight)) return powers;
    powers.push_back(std::move(next));
  }
}

Vec primitive(Vec v)
{
  Int g = 0;
  for (const Int& x : v) g = gcd(g, x);
  if (g == 0) return v;
  for (Int& x : v) x /= g;
  return v;
}

bool in_fundamental_parallelepiped(const std::vector<Vec>& columns, const Vec& point)
{
  Int determinant_value = det_columns(columns);
  Int den = abs(determinant_value);
  int s = determinant_value > 0 ? 1 : -1;
  for (size_t j = 0; j < columns.size(); ++j) {
    auto changed = columns;
    changed[j] = point;
    Int numerator = det_columns(changed) * s;
    if (numerator < 0 || numerator >= den) return false;
  }
  return true;
}

void enumerate_box(const Vec& lo, const Vec& hi, size_t at, Vec& point,
  const std::function<void(const Vec&)>& visit)
{
  if (at == point.size()) {
    visit(point);
    return;
  }
  for (Int x = lo[at]; x <= hi[at]; ++x) {
    point[at] = x;
    enumerate_box(lo, hi, at + 1, point, visit);
  }
}

std::vector<Vec> cone_lattice_candidates(const std::vector<Vec>& input)
{
  if (input.empty()) return {};
  size_t d = input.front().size();
  std::set<Vec, VecLess> rays;
  for (auto v : input) if (!is_zero(v)) rays.insert(primitive(std::move(v)));
  std::vector<Vec> generators(rays.begin(), rays.end());
  if (matrix_rank(generators) != static_cast<int>(d))
    throw std::invalid_argument("the cone is not full-dimensional in its quotient lattice");

  std::set<Vec, VecLess> candidates(rays.begin(), rays.end());
  std::vector<size_t> chosen;
  std::function<void(size_t)> combinations = [&](size_t start) {
    if (chosen.size() == d) {
      std::vector<Vec> columns;
      for (size_t i : chosen) columns.push_back(generators[i]);
      if (det_columns(columns) == 0) return;

      Vec lo(d), hi(d), point(d);
      for (size_t row = 0; row < d; ++row) {
        for (const auto& column : columns) {
          if (column[row] < 0) lo[row] += column[row];
          if (column[row] > 0) hi[row] += column[row];
        }
      }
      enumerate_box(lo, hi, 0, point, [&](const Vec& x) {
        if (!is_zero(x) && in_fundamental_parallelepiped(columns, x)) candidates.insert(x);
      });
      return;
    }
    for (size_t i = start; i < generators.size(); ++i) {
      chosen.push_back(i);
      combinations(i + 1);
      chosen.pop_back();
    }
  };
  combinations(0);
  return {candidates.begin(), candidates.end()};
}

std::vector<Vec> hilbert_basis(std::vector<Vec> candidates, const Vec& weight)
{
  std::sort(candidates.begin(), candidates.end(), [&](const Vec& a, const Vec& b) {
    Int wa = dot(weight, a);
    Int wb = dot(weight, b);
    if (wa != wb) return wa < wb;
    return VecLess{}(a, b);
  });

  std::vector<Vec> hilbert;
  for (const auto& point : candidates)
    if (!affine_contains(hilbert, point, weight)) hilbert.push_back(point);
  return hilbert;
}

struct ExtendedGcd {
  Int g;
  Int x;
  Int y;
};

ExtendedGcd extended_gcd(Int a, Int b)
{
  Int sa = a < 0 ? -1 : 1;
  Int sb = b < 0 ? -1 : 1;
  a = abs(a);
  b = abs(b);
  Int old_r = a, r = b;
  Int old_s = 1, s = 0;
  Int old_t = 0, t = 1;
  while (r != 0) {
    Int q = old_r / r;
    Int next = old_r - q * r;
    old_r = r;
    r = next;
    next = old_s - q * s;
    old_s = s;
    s = next;
    next = old_t - q * t;
    old_t = t;
    t = next;
  }
  return {old_r, old_s * sa, old_t * sb};
}

struct QuadraticLattice {
  Int denominator;
  std::array<Vec, 2> basis;
  std::vector<Vec> generators;
};

QuadraticLattice quadratic_lattice(const std::vector<Quadratic>& generators)
{
  Int denominator = 1;
  for (const auto& x : generators) {
    denominator = lcm(denominator, x.a.d);
    denominator = lcm(denominator, x.b.d);
  }

  std::vector<Vec> scaled;
  for (const auto& x : generators)
    scaled.push_back({x.a.n * (denominator / x.a.d), x.b.n * (denominator / x.b.d)});
  if (matrix_rank(scaled) != 2) throw std::invalid_argument("the generators do not span a quadratic field");

  Int gy = 0;
  Int x_for_gy = 0;
  for (const auto& v : scaled) {
    ExtendedGcd e = extended_gcd(gy, v[1]);
    x_for_gy = e.x * x_for_gy + e.y * v[0];
    gy = e.g;
  }

  Int determinant_gcd = 0;
  for (size_t i = 0; i < scaled.size(); ++i)
    for (size_t j = i + 1; j < scaled.size(); ++j)
      determinant_gcd = gcd(determinant_gcd,
        scaled[i][0] * scaled[j][1] - scaled[i][1] * scaled[j][0]);
  determinant_gcd = abs(determinant_gcd);
  if (gy == 0 || determinant_gcd == 0 || determinant_gcd % gy != 0)
    throw std::logic_error("could not construct the quadratic quotient lattice");

  Int hx = determinant_gcd / gy;
  x_for_gy %= hx;
  if (x_for_gy < 0) x_for_gy += hx;
  std::array<Vec, 2> basis{Vec{hx, 0}, Vec{x_for_gy, gy}};

  std::vector<Vec> coordinates;
  for (const auto& v : scaled) {
    Int second = v[1] / gy;
    Int numerator = v[0] - second * x_for_gy;
    if (v[1] % gy != 0 || numerator % hx != 0)
      throw std::logic_error("quadratic lattice basis does not contain a generator");
    coordinates.push_back({numerator / hx, second});
  }
  return {denominator, basis, coordinates};
}

Vec from_lattice_coordinates(const std::array<Vec, 2>& basis, const Vec& x)
{
  return {
    basis[0][0] * x[0] + basis[1][0] * x[1],
    basis[0][1] * x[0] + basis[1][1] * x[1]
  };
}

Vec quadratic_lattice_weight(long long D, const std::array<Vec, 2>& basis,
  const std::vector<Vec>& coordinates)
{
  std::vector<Vec> physical;
  for (const auto& x : coordinates) physical.push_back(from_lattice_coordinates(basis, x));
  Vec w = quadratic_weight(D, physical);
  return {dot(w, basis[0]), dot(w, basis[1])};
}

std::string quadratic_lattice_string(const QuadraticLattice& lattice, const Vec& x, long long D)
{
  Vec physical = from_lattice_coordinates(lattice.basis, x);
  Quadratic value{
    Rat(physical[0], lattice.denominator),
    Rat(physical[1], lattice.denominator),
    D
  };
  return quadratic_string(value);
}

bool squarefree(long long n)
{
  if (n <= 1) return false;
  for (long long p = 2; p <= n / p; ++p) {
    long long square = p * p;
    if (n % square == 0) return false;
  }
  return true;
}

struct QuadraticData {
  bool finite = false;
  std::string reason;
  std::vector<Quadratic> alpha;
};

QuadraticData quadratic_finite_generation(long long D,
  const std::vector<QuadraticGenerator>& input)
{
  if (!squarefree(D)) throw std::invalid_argument("D must be a positive squarefree integer greater than one");
  if (input.empty()) throw std::invalid_argument("at least one generator is required");

  QuadraticData data;
  bool spans_quadratic_field = false;
  bool negative_conjugate = false;
  for (const auto& generator : input) {
    Quadratic x{Rat(generator.a), Rat(generator.b), D};
    if (sign(x) <= 0) throw std::invalid_argument("every quadratic generator must be positive");
    spans_quadratic_field = spans_quadratic_field || x.b.n != 0;
    if (!algebraic_integer(x)) {
      data.reason = quadratic_string(x) + " is not an algebraic integer";
      return data;
    }
    if (!weak_perron(x)) {
      data.reason = quadratic_string(x) + " is not weak Perron";
      return data;
    }
    negative_conjugate = negative_conjugate || sign(conjugate(x)) < 0;
    data.alpha.push_back(std::move(x));
  }
  if (!spans_quadratic_field)
    throw std::invalid_argument("the supplied generators span Q, not a quadratic field");
  if (!negative_conjugate) {
    data.reason = "none of the generators has a negative conjugate";
    return data;
  }

  data.finite = true;
  data.reason = "all generators are weak Perron algebraic integers and one has a negative conjugate";
  return data;
}

std::vector<std::string> quadratic_generator_strings(const std::vector<Quadratic>& generators)
{
  std::vector<std::string> out;
  for (const auto& x : generators) out.push_back(quadratic_string(x));
  return out;
}

std::vector<std::string> cyclic_generator_strings(const std::vector<Vec>& generators)
{
  std::vector<std::string> out;
  for (size_t i = 0; i < generators.size(); ++i)
    out.push_back("alpha^" + std::to_string(i) + " = " + vector_string(generators[i]));
  return out;
}

Answer quadratic_answer(long long D, const std::vector<QuadraticGenerator>& input, bool krull)
{
  QuadraticData data = quadratic_finite_generation(D, input);
  Answer answer;
  answer.property = krull ? "Krull" : "Mori";
  if (!data.finite) {
    answer.reason = data.reason + "; the monoid is not finitely generated";
    return answer;
  }

  std::vector<Quadratic> generators = quadratic_generating_set(data.alpha);
  answer.generators = quadratic_generator_strings(generators);
  if (!krull) {
    answer.value = true;
    answer.reason = data.reason + "; Mori is equivalent to finite generation here";
    return answer;
  }

  QuadraticLattice lattice = quadratic_lattice(generators);
  auto candidates = cone_lattice_candidates(lattice.generators);
  Vec weight = quadratic_lattice_weight(D, lattice.basis, candidates);
  auto hilbert = hilbert_basis(std::move(candidates), weight);
  for (const auto& h : hilbert)
    answer.hilbert_generators.push_back(quadratic_lattice_string(lattice, h, D));

  for (const auto& h : hilbert) {
    if (!affine_contains(lattice.generators, h, weight)) {
      answer.reason = "the finite generating set is not saturated in its quotient lattice";
      answer.hole = quadratic_lattice_string(lattice, h, D);
      return answer;
    }
  }
  answer.value = true;
  answer.reason = "the monoid is finitely generated and every Hilbert-basis element is in it";
  return answer;
}

Answer cyclic_answer(const std::vector<std::string>& polynomial, bool krull)
{
  CyclicData data = cyclic_finite_generation(polynomial);
  Answer answer;
  answer.property = krull ? "Krull" : "Mori";
  if (!data.finite) {
    answer.reason = data.reason + "; the monoid is not finitely generated";
    return answer;
  }

  auto generators = cyclic_generating_vectors(data);
  answer.generators = cyclic_generator_strings(generators);
  if (!krull) {
    answer.value = true;
    answer.reason = data.reason + "; cyclic Mori is equivalent to finite generation";
    return answer;
  }

  auto candidates = cone_lattice_candidates(generators);
  Vec weight = cyclic_weight(data.polynomial, data.alpha, candidates);
  auto hilbert = hilbert_basis(std::move(candidates), weight);
  for (const auto& h : hilbert) answer.hilbert_generators.push_back(vector_string(h));

  for (const auto& h : hilbert) {
    if (!affine_contains(generators, h, weight)) {
      answer.reason = "the affine monoid is not saturated in Z[alpha]";
      answer.hole = vector_string(h);
      return answer;
    }
  }
  answer.value = true;
  answer.reason = "the monoid is finitely generated and every Hilbert-basis element belongs to it";
  return answer;
}

} // namespace

Answer quadratic_mori(long long D, const std::vector<QuadraticGenerator>& generators)
{
  return quadratic_answer(D, generators, false);
}

Answer quadratic_krull(long long D, const std::vector<QuadraticGenerator>& generators)
{
  return quadratic_answer(D, generators, true);
}

Answer cyclic_mori(const std::vector<std::string>& polynomial)
{
  return cyclic_answer(polynomial, false);
}

Answer cyclic_krull(const std::vector<std::string>& polynomial)
{
  return cyclic_answer(polynomial, true);
}

} // namespace morikrull
