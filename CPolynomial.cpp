#ifndef __PROGTEST__
#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cfloat>
#include <climits>
#include <cmath>
#include <compare>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <span>
#include <sstream>
#include <string>
#include <vector>
#endif /* __PROGTEST__ */

// keep this dummy version if you do not implement a real manipulator
std::ios_base& (*poly_var(const std::string& name))(std::ios_base& x) {
  return [](std::ios_base& ios) -> std::ios_base& { return ios; };
}

class CPolynomial {
 public:
  CPolynomial() : m_Coefficients(1, 0.0) {}

  CPolynomial(const std::vector<double>& Coefficients) {
    this->m_Coefficients = Coefficients;
    normalize();
  }

  CPolynomial(const CPolynomial& oth) = default;

  ~CPolynomial() = default;

  CPolynomial& operator=(const CPolynomial& oth) = default;

  CPolynomial operator*(double scalar) const;

  CPolynomial operator*(const CPolynomial& oth) const;

  CPolynomial& operator*=(double scalar);

  CPolynomial& operator*=(const CPolynomial& oth);

  bool operator==(const CPolynomial& oth) const;

  bool operator!=(const CPolynomial& oth) const;

  double operator[](size_t index) const;

  double& operator[](size_t index);

  double operator()(double x) const;

  bool operator!() const;

  explicit operator bool() const;

  unsigned degree() const;

  friend std::ostream& operator<<(std::ostream& os,
                                  const CPolynomial& polynomial);

  // getCoefficients() is declared in public as used in dumpMatch()
  // I don't really know how is dumpMatch() implemented in ProgTest,
  // so will keep getCoefficients() here.

  std::vector<double> getCoefficients() const;

 private:
  std::vector<double> m_Coefficients;
  using complex_t = std::complex<double>;  // Used in FFT for sorting complex
                                           // values of n^th roots.
  static const int DEG_UP_BOUND = 200;     // Upper bound of degree() for which
                                        // linear multiplication is available.

  void normalize();

  // FFT Algorithm implementation based on Math Paper at
  // https://www.cs.cmu.edu/~15451-f22/lectures/lec25-fft.pdf
  // - - - - - - - - - - - - - - - - - - - - - -
  // Code based on sample at
  // https://www.geeksforgeeks.org/fast-fourier-transformation-poynomial-multiplication/

  static void fft(std::vector<complex_t>& a, bool invert);

  static std::vector<double> multiply_by_fft(const std::vector<double>& a,
                                             const std::vector<double>& b);
};

CPolynomial CPolynomial::operator*(double scalar) const {
  std::vector<double> result = m_Coefficients;
  for (double& coef : result) {
    coef *= scalar;
  }
  return CPolynomial(result);
}

CPolynomial CPolynomial::operator*(const CPolynomial& oth) const {
  if (degree() > DEG_UP_BOUND && oth.degree() > DEG_UP_BOUND) {
    return CPolynomial(
        multiply_by_fft(this->m_Coefficients, oth.m_Coefficients));
  } else {
    std::vector<double> result(
        m_Coefficients.size() + oth.m_Coefficients.size() - 1, 0.0);
    for (size_t i = 0; i < m_Coefficients.size(); ++i) {
      for (size_t j = 0; j < oth.m_Coefficients.size(); ++j) {
        result[i + j] += m_Coefficients[i] * oth.m_Coefficients[j];
      }
    }
    return CPolynomial(result);
  }
}

bool CPolynomial::operator==(const CPolynomial& oth) const {
  if (degree() != oth.degree()) return false;
  // Use rel. and abs. tolerance to compare floating point nums,
  // allowing a little numerical unprecision after FFT or  *=.

  for (size_t i = 0; i < m_Coefficients.size(); ++i) {
    if (!(std::fabs(m_Coefficients[i] - oth.m_Coefficients[i]) <=
              std::max(std::fabs(m_Coefficients[i]),
                       std::fabs(oth.m_Coefficients[i])) *
                  1e-4 ||
          std::fabs(m_Coefficients[i] - oth.m_Coefficients[i]) < 1e-10))
      return false;
  }
  return true;
}

CPolynomial& CPolynomial::operator*=(double scalar) {
  for (double& c : m_Coefficients) {
    c *= scalar;
  }
  normalize();
  return *this;
}

CPolynomial& CPolynomial::operator*=(const CPolynomial& oth) {
  *this = *this * oth;
  return *this;
}

bool CPolynomial::operator!=(const CPolynomial& oth) const {
  return !(*this == oth);
}

double CPolynomial::operator[](size_t index) const {
  if (index < m_Coefficients.size()) return m_Coefficients[index];
  return 0.0;
}

double& CPolynomial::operator[](size_t index) {
  // Resize if out of bound index.
  // Then can write on index > size(), the rest will be zeroes

  if (index >= m_Coefficients.size()) {
    m_Coefficients.resize(index + 1, 0.0);
  }
  return m_Coefficients[index];
}

double CPolynomial::operator()(double x) const {
  double result = 0.0;
  double power = 1.0;
  for (double c : m_Coefficients) {
    result += c * power;
    power *= x;
  }
  return result;
}

bool CPolynomial::operator!() const { return !static_cast<bool>(*this); }

CPolynomial::operator bool() const {
  return std::any_of(m_Coefficients.begin(), m_Coefficients.end(),
                     [](double coef) { return coef != 0.0; });
}

unsigned CPolynomial::degree() const {
  if (m_Coefficients.empty()) return 0;
  for (size_t i = m_Coefficients.size(); i > 0; --i) {
    if (std::fabs(m_Coefficients[i - 1]) > 1e-10) {
      return i - 1;
    }
  }
  return 0;
}

std::ostream& operator<<(std::ostream& os, const CPolynomial& polynomial) {
  bool first = true;
  // We set a flag which indicates if we have printed any non-zero term yet.
  // If we didn't -> polynomial is zero -> print zero;

  for (int i = polynomial.m_Coefficients.size() - 1; i >= 0; --i) {
    double coef = polynomial.m_Coefficients[i];
    if (coef == 0) continue;

    if (!first) {
      os << (coef > 0 ? " + " : " - ");
    } else {
      if (coef < 0) os << "- ";
    }

    coef = std::abs(coef);
    if (coef != 1 && i == 0) os << coef;
    if (coef == 1 && i != 0) os << "x^" << i;
    if (coef != 1 && i != 0) os << coef << "*x^" << i;

    first = false;
  }
  if (first) os << "0";
  return os;
}

std::vector<double> CPolynomial::getCoefficients() const {
  return m_Coefficients;
}

void CPolynomial::normalize() {
  // Removing until no zeroes remain, if empty -> just set to zero.

  while (!m_Coefficients.empty() && std::fabs(m_Coefficients.back()) < 1e-10) {
    m_Coefficients.pop_back();
  }
  if (m_Coefficients.empty()) {
    m_Coefficients.push_back(0.0);
  }
}

void CPolynomial::fft(std::vector<complex_t>& a, bool invert) {
  size_t n = a.size();
  if (n == 1) return;

  std::vector<complex_t> a0(n / 2), a1(n / 2);
  for (size_t i = 0; 2 * i < n; i++) {
    a0[i] = a[i * 2];
    a1[i] = a[i * 2 + 1];
  }

  fft(a0, invert);
  fft(a1, invert);

  double angle = 2 * acos(-1) / n * (invert ? -1 : 1);
  complex_t w(1), wn(cos(angle), sin(angle));

  for (size_t i = 0; 2 * i < n; i++) {
    a[i] = a0[i] + w * a1[i];
    a[i + n / 2] = a0[i] - w * a1[i];

    if (invert) {
      a[i] /= 2;
      a[i + n / 2] /= 2;
    }
    w *= wn;
  }
}

std::vector<double> CPolynomial::multiply_by_fft(const std::vector<double>& a,
                                                 const std::vector<double>& b) {
  std::vector<complex_t> fa(a.begin(), a.end()), fb(b.begin(), b.end());
  size_t finalSize = 1;
  while (finalSize < a.size() + b.size()) {
    finalSize *= 2;
  }

  fa.resize(finalSize);
  fb.resize(finalSize);

  fft(fa, false);
  fft(fb, false);

  for (size_t i = 0; i < finalSize; i++) {
    fa[i] *= fb[i];
  }

  fft(fa, true);

  std::vector<double> result(finalSize);
  for (size_t i = 0; i < finalSize; i++) {
    result[i] = round(fa[i].real() * 1e9) / 1e9;
  }

  while (!result.empty() && abs(result.back()) < 1e-9) {
    result.pop_back();
  }

  return result;
}

#ifndef __PROGTEST__
bool smallDiff(double a, double b) {
  return std::fabs(a - b) <= std::max(std::fabs(a), std::fabs(b)) * 1e-4;
}

bool dumpMatch(const CPolynomial& x, const std::vector<double>& ref) {
  return x.getCoefficients() == ref;
}

void test0() {
  CPolynomial a, b, c;
  std::ostringstream out, tmp;
  a[0] = -10;
  a[1] = 3.5;
  a[3] = 1;
  assert(smallDiff(a(2), 5));
  out.str("");
  out << a;
  assert(out.str() == "x^3 + 3.5*x^1 - 10");
  c = a * -2;
  assert(c.degree() == 3 &&
         dumpMatch(c, std::vector<double>{20.0, -7.0, -0.0, -2.0}));

  out.str("");
  out << c;
  assert(out.str() == "- 2*x^3 - 7*x^1 + 20");
  out.str("");
  out << b;
  assert(out.str() == "0");
  b[5] = -1;
  b[2] = 3;
  out.str("");
  out << b;
  assert(out.str() == "- x^5 + 3*x^2");
  c = a * b;
  assert(c.degree() == 8 &&
         dumpMatch(c, std::vector<double>{-0.0, -0.0, -30.0, 10.5, -0.0, 13.0,
                                          -3.5, 0.0, -1.0}));

  out.str("");
  out << c;
  assert(out.str() == "- x^8 - 3.5*x^6 + 13*x^5 + 10.5*x^3 - 30*x^2");
  a *= 5;
  assert(a.degree() == 3 &&
         dumpMatch(a, std::vector<double>{-50.0, 17.5, 0.0, 5.0}));

  a *= b;
  assert(a.degree() == 8 &&
         dumpMatch(a, std::vector<double>{0.0, 0.0, -150.0, 52.5, -0.0, 65.0,
                                          -17.5, -0.0, -5.0}));

  assert(a != b);
  b[5] = 0;
  assert(static_cast<bool>(b));
  assert(!!b);
  b[2] = 0;
  assert(!(a == b));
  a *= 0;
  assert(a.degree() == 0 && dumpMatch(a, std::vector<double>{0.0}));

  assert(a == b);
  assert(!static_cast<bool>(b));
  assert(!b);
}

void test1_base() {
  CPolynomial p0;
  assert(p0.degree() == 0);
  assert(p0(0) == 0.0);
  assert(!p0);
  assert(!!p0 == false);

  CPolynomial p1;
  p1[0] = 5.0;
  assert(p1.degree() == 0);
  assert(p1(10) == 5.0);
  assert(p1(0) == 5.0);
  assert(p1 == CPolynomial(std::vector<double>{5.0}));

  CPolynomial p2;
  p2[1] = 2.0;
  p2[0] = 1.0;
  assert(p2(3.0) == 7.0);
}

void test2_borderline() {
  CPolynomial p3({0.0, 1.0, 0.0, 2.0, 0.0, 0.0});
  assert(p3.degree() == 3);
  assert(p3[3] == 2.0);
  assert(p3[5] == 0.0);
  assert(p3[10] == 0.0);

  CPolynomial p4;
  p4[100] = 1.0;
  assert(p4.degree() == 100);
  assert(p4[99] == 0.0);

  CPolynomial p5({1, 2, 3});
  CPolynomial zero;
  CPolynomial prod = p5 * zero;
  assert(prod == CPolynomial({0.0}));
  assert(prod.degree() == 0);

  CPolynomial p6({1.0000000001, 2.0});
  CPolynomial p7({1.0000000002, 2.0});
  assert(p6 == p7);
}

void test3_large() {
  const size_t DEG = 100000000;
  CPolynomial p;
  p[DEG] = 1.0;
  assert(p.degree() == DEG);

  const size_t DEG2 = 500000000;
  CPolynomial c;
  for (size_t i = 0; i <= DEG2; ++i) c[i] = 1.0;

  assert(c.degree() == DEG2);
  assert(std::abs(c(1.0) - (DEG2 + 1)) < 1e-6);

  CPolynomial s;
  for (size_t i = 0; i < 100000; i += 1000) s[i] = 2.0;

  for (size_t i = 0; i < 100000; i++) {
    if (i % 1000 == 0)
      assert(s[i] == 2.0);
    else
      assert(s[i] == 0.0);
  }

  CPolynomial bigZero;
  for (size_t i = 0; i < 100000; ++i) bigZero[i] = 0.0;

  assert(bigZero.degree() == 0);
  assert(bigZero[9999] < 1e-10);
  assert(!bigZero);

  CPolynomial poly;
  for (size_t i = 0; i < 10000; ++i) poly[i] = 1.0;

  CPolynomial scaled = poly * 1000.0;
  for (size_t i = 0; i < 10000; ++i)
    assert(std::abs(scaled[i] - 1000.0) < 1e-9);

  CPolynomial a, b;
  for (size_t i = 0; i < 5000; ++i) {
    a[i] = std::sin(i);
    b[i] = std::sin(i);
  }
  assert(a == b);

  b[4999] += 1e-12;
  assert(a == b);

  b[4999] += 1e-2;
  assert(a != b);

  CPolynomial z;
  for (size_t i = 0; i < 1000; ++i) z[i] = 1.0;

  double sum = 0.0;
  for (size_t i = 0; i < 1000; ++i) sum += std::pow(2.0, i);
  assert(std::abs(z(2.0) - sum) < 1e-6);
}

int main() {
  test0();
  test1_base();
  test2_borderline();
  test3_large();
  std::cout << "All tests passed" << std::endl;
  // // bonus - manipulators

  // out.str("");
  // out << poly_var("y") << c;
  // assert(out.str() == "- y^8 - 3.5*y^6 + 13*y^5 + 10.5*y^3 - 30*y^2");
  // out.str("");
  // tmp << poly_var("abc");
  // out.copyfmt(tmp);
  // out << c;
  // assert(out.str() == "- abc^8 - 3.5*abc^6 + 13*abc^5 + 10.5*abc^3 -
  // 30*abc^2");
  return EXIT_SUCCESS;
}
#endif /* __PROGTEST__ */
