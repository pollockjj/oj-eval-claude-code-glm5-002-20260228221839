#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

// Integer 1:
// Implement a signed big integer class that only needs to support simple addition and subtraction

// Integer 2:
// Implement a signed big integer class that supports addition, subtraction, multiplication, and division, and overload related operators

// Do not use any header files other than the following
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

// Do not use "using namespace std;"

namespace sjtu {

static const long long BASE = 1000000000LL;
static const int BASE_DIGITS = 9;

// Helper functions for Karatsuba multiplication
static std::vector<long long> add_vec(const std::vector<long long>& v1, const std::vector<long long>& v2) {
    std::vector<long long> result(std::max(v1.size(), v2.size()) + 1, 0);
    for (size_t i = 0; i < v1.size(); i++) result[i] += v1[i];
    for (size_t i = 0; i < v2.size(); i++) result[i] += v2[i];
    long long carry = 0;
    for (size_t i = 0; i < result.size(); i++) {
        result[i] += carry;
        carry = result[i] / BASE;
        result[i] %= BASE;
    }
    while (result.size() > 1 && result.back() == 0) result.pop_back();
    return result;
}

static std::vector<long long> sub_vec(const std::vector<long long>& v1, const std::vector<long long>& v2) {
    std::vector<long long> result = v1;
    result.resize(std::max(v1.size(), v2.size()), 0);
    long long borrow = 0;
    for (size_t i = 0; i < result.size(); i++) {
        result[i] -= borrow;
        if (i < v2.size()) result[i] -= v2[i];
        if (result[i] < 0) {
            result[i] += BASE;
            borrow = 1;
        } else {
            borrow = 0;
        }
    }
    while (result.size() > 1 && result.back() == 0) result.pop_back();
    return result;
}

// NTT (Number Theoretic Transform) for fast multiplication
// Using modulus 998244353 = 119 * 2^23 + 1 (a prime)
static const long long NTT_MOD = 998244353LL;
static const long long NTT_ROOT = 3LL;  // primitive root
static const int NTT_MAX_LEN = 23;  // 2^23

static long long ntt_pow(long long base, long long exp, long long mod) {
    long long result = 1;
    while (exp > 0) {
        if (exp & 1) result = (__int128)result * base % mod;
        base = (__int128)base * base % mod;
        exp >>= 1;
    }
    return result;
}

static void ntt_prepare_roots(int n, std::vector<long long>& roots, std::vector<long long>& inv_roots) {
    long long wn = ntt_pow(NTT_ROOT, (NTT_MOD - 1) / n, NTT_MOD);
    long long inv_wn = ntt_pow(wn, NTT_MOD - 2, NTT_MOD);
    roots.resize(n);
    inv_roots.resize(n);
    roots[0] = 1;
    inv_roots[0] = 1;
    for (int i = 1; i < n; i++) {
        roots[i] = (__int128)roots[i-1] * wn % NTT_MOD;
        inv_roots[i] = (__int128)inv_roots[i-1] * inv_wn % NTT_MOD;
    }
}

static void ntt_transform(std::vector<long long>& a, int n, const std::vector<long long>& roots, bool inverse) {
    // Bit reversal permutation
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        int half = len >> 1;
        long long wlen = inverse ? roots[half] : roots[half];
        for (int i = 0; i < n; i += len) {
            long long w = 1;
            for (int j = 0; j < half; j++) {
                long long u = a[i + j];
                long long v = (__int128)a[i + j + half] * w % NTT_MOD;
                a[i + j] = (u + v) % NTT_MOD;
                a[i + j + half] = (u - v + NTT_MOD) % NTT_MOD;
                w = (__int128)w * wlen % NTT_MOD;
            }
        }
    }

    if (inverse) {
        long long inv_n = ntt_pow(n, NTT_MOD - 2, NTT_MOD);
        for (int i = 0; i < n; i++) {
            a[i] = (__int128)a[i] * inv_n % NTT_MOD;
        }
    }
}

static std::vector<long long> multiply_ntt(const std::vector<long long>& x, const std::vector<long long>& y) {
    size_t xn = x.size();
    size_t yn = y.size();
    if (xn == 0 || yn == 0) return {0};

    int n = 1;
    while (n < (int)(xn + yn)) n <<= 1;

    std::vector<long long> fx(n, 0), fy(n, 0);
    for (size_t i = 0; i < xn; i++) fx[i] = x[i] % NTT_MOD;
    for (size_t i = 0; i < yn; i++) fy[i] = y[i] % NTT_MOD;

    std::vector<long long> roots, inv_roots;
    ntt_prepare_roots(n, roots, inv_roots);

    ntt_transform(fx, n, roots, false);
    ntt_transform(fy, n, roots, false);

    for (int i = 0; i < n; i++) {
        fx[i] = (__int128)fx[i] * fy[i] % NTT_MOD;
    }

    ntt_transform(fx, n, inv_roots, true);

    // Convert back to base BASE
    std::vector<long long> result(xn + yn, 0);
    long long carry = 0;
    for (size_t i = 0; i < xn + yn; i++) {
        result[i] += carry;
        long long val = fx[i];
        // Compute quotient and remainder when dividing val by NTT_MOD
        carry = val / NTT_MOD;
        result[i] += val % NTT_MOD;
    }

    while (result.size() > 1 && result.back() == 0) result.pop_back();
    return result;
}

static std::vector<long long> multiply_simple(const std::vector<long long>& x, const std::vector<long long>& y) {
    size_t xn = x.size();
    size_t yn = y.size();
    if (xn == 0 || yn == 0) return {0};
    std::vector<long long> result(xn + yn, 0);
    for (size_t i = 0; i < xn; i++) {
        long long carry = 0;
        for (size_t j = 0; j < yn; j++) {
            __int128 prod = (__int128)x[i] * y[j] + result[i + j] + carry;
            carry = prod / BASE;
            result[i + j] = prod % BASE;
        }
        result[i + yn] = carry;
    }
    while (result.size() > 1 && result.back() == 0) result.pop_back();
    return result;
}

static std::vector<long long> multiply_ntt(const std::vector<long long>& x, const std::vector<long long>& y);

static std::vector<long long> multiply_karatsuba(const std::vector<long long>& x, const std::vector<long long>& y) {
    size_t xn = x.size();
    size_t yn = y.size();

    if (xn < yn) return multiply_karatsuba(y, x);
    if (yn == 0 || (xn == 1 && x[0] == 0)) return {0};
    if (xn <= 16 || yn <= 16) {
        return multiply_simple(x, y);
    }

    size_t split = xn / 2;

    std::vector<long long> x_low(x.begin(), x.begin() + split);
    std::vector<long long> x_high(x.begin() + split, x.end());
    std::vector<long long> y_low(y.begin(), y.begin() + std::min(split, yn));
    std::vector<long long> y_high(y.begin() + std::min(split, yn), y.end());

    // Normalize low parts
    while (x_low.size() > 1 && x_low.back() == 0) x_low.pop_back();
    while (y_low.size() > 1 && y_low.back() == 0) y_low.pop_back();

    auto z0 = multiply_karatsuba(x_low, y_low);

    // z2 can be 0 if x_high or y_high is empty
    std::vector<long long> z2;
    if (x_high.empty() || y_high.empty()) {
        z2 = {0};
    } else {
        z2 = multiply_karatsuba(x_high, y_high);
    }

    // Compute z1 = (x_low + x_high) * (y_low + y_high) - z0 - z2
    auto sum_x = add_vec(x_low, x_high);
    auto sum_y = add_vec(y_low, y_high);
    auto z1_full = multiply_karatsuba(sum_x, sum_y);

    // z1 = z1_full - z0 - z2 (subtract z0 and z2 from z1_full)
    // We need to handle potential negative result
    // Since z1_full >= z0 + z2 mathematically, we just subtract
    auto z1 = sub_vec(z1_full, z0);
    z1 = sub_vec(z1, z2);

    // Compute result = z0 + z1 * B^split + z2 * B^(2*split)
    size_t max_size = std::max({z0.size(), z1.size() + split, z2.size() + 2 * split});
    std::vector<long long> result(max_size + 1, 0);

    for (size_t i = 0; i < z0.size(); i++) result[i] += z0[i];
    for (size_t i = 0; i < z1.size(); i++) result[i + split] += z1[i];
    for (size_t i = 0; i < z2.size(); i++) result[i + 2 * split] += z2[i];

    // Normalize the result
    long long carry = 0;
    for (size_t i = 0; i < result.size(); i++) {
        result[i] += carry;
        carry = result[i] / BASE;
        result[i] %= BASE;
    }
    while (result.size() > 1 && result.back() == 0) result.pop_back();

    return result;
}

class int2048 {
private:
  std::vector<long long> data;
  bool neg;

  void normalize() {
      while (data.size() > 1 && data.back() == 0) {
          data.pop_back();
      }
      if (data.size() == 1 && data[0] == 0) {
          neg = false;
      }
  }

  int abs_compare(const int2048 &other) const {
      if (data.size() != other.data.size()) {
          return data.size() < other.data.size() ? -1 : 1;
      }
      for (int i = (int)data.size() - 1; i >= 0; i--) {
          if (data[i] != other.data[i]) {
              return data[i] < other.data[i] ? -1 : 1;
          }
      }
      return 0;
  }

  static void divide(const int2048 &dividend, const int2048 &divisor,
                     int2048 &quotient, int2048 &remainder) {
      int cmp = dividend.abs_compare(divisor);

      if (cmp < 0) {
          quotient.data.clear();
          quotient.data.push_back(0);
          quotient.neg = false;
          remainder.data = dividend.data;
          remainder.neg = false;
          return;
      }

      if (cmp == 0) {
          quotient.data.clear();
          quotient.data.push_back(1);
          quotient.neg = false;
          remainder.data.clear();
          remainder.data.push_back(0);
          remainder.neg = false;
          return;
      }

      size_t n = divisor.data.size();
      size_t m = dividend.data.size() - n;

      quotient.data.clear();
      quotient.data.resize(m + 1, 0);
      remainder.data = dividend.data;
      remainder.neg = false;

      // Normalize divisor for efficient trial division
      long long d = BASE / (divisor.data.back() + 1);

      // Scale both dividend and divisor
      std::vector<long long> u = dividend.data;
      std::vector<long long> v = divisor.data;

      u.push_back(0);
      for (auto& x : u) x *= d;
      for (size_t i = 0; i + 1 < u.size(); i++) {
          u[i + 1] += u[i] / BASE;
          u[i] %= BASE;
      }
      if (u.back() == 0) u.pop_back();

      for (auto& x : v) x *= d;
      for (size_t i = 0; i + 1 < v.size(); i++) {
          v[i + 1] += v[i] / BASE;
          v[i] %= BASE;
      }
      while (v.size() > 1 && v.back() == 0) v.pop_back();

      // Ensure u has enough space for the division algorithm
      while (u.size() < dividend.data.size() + 1) u.push_back(0);

      size_t vn = v.size();
      long long v1 = v[vn - 1];
      long long v2 = (vn > 1) ? v[vn - 2] : 0;

      for (int j = m; j >= 0; j--) {
          long long uj = ((long long)j + vn < (long long)u.size()) ? u[j + vn] : 0;
          long long uj1 = ((long long)j + vn - 1 < (long long)u.size()) ? u[j + vn - 1] : 0;
          long long uj2 = ((long long)j + vn - 2 < (long long)u.size()) ? u[j + vn - 2] : 0;

          long long q_hat = (uj * BASE + uj1) / v1;
          long long r_hat = (uj * BASE + uj1) % v1;

          while (q_hat >= BASE || q_hat * v2 > BASE * r_hat + uj2) {
              q_hat--;
              r_hat += v1;
              if (r_hat >= BASE) break;
          }

          long long borrow = 0;
          for (size_t i = 0; i < vn; i++) {
              long long t = u[j + i] - borrow - q_hat * v[i];
              if (t < 0) {
                  // Need to borrow multiple BASEs
                  long long num_borrows = (-t + BASE - 1) / BASE;
                  t += num_borrows * BASE;
                  borrow = num_borrows;
              } else {
                  borrow = 0;
              }
              u[j + i] = t;
          }
          u[j + vn] -= borrow;

          quotient.data[j] = q_hat;

          if (u[j + vn] < 0) {
              quotient.data[j]--;
              long long carry = 0;
              for (size_t i = 0; i < vn; i++) {
                  u[j + i] += v[i] + carry;
                  if (u[j + i] >= BASE) {
                      u[j + i] -= BASE;
                      carry = 1;
                  } else {
                      carry = 0;
                  }
              }
              u[j + vn] += carry;
          }
      }

      remainder.data.clear();
      for (size_t i = 0; i < n; i++) {
          remainder.data.push_back(u[i]);
      }

      // Denormalize remainder
      long long carry = 0;
      for (int i = (int)remainder.data.size() - 1; i >= 0; i--) {
          remainder.data[i] += carry * BASE;
          carry = remainder.data[i] % d;
          remainder.data[i] /= d;
      }

      quotient.normalize();
      remainder.normalize();
  }

public:
  // Constructors
  int2048() {
      data.resize(1, 0);
      neg = false;
  }

  int2048(long long x) {
      neg = x < 0;
      if (x < 0) x = -x;
      if (x == 0) {
          data.push_back(0);
      } else {
          while (x > 0) {
              data.push_back(x % BASE);
              x /= BASE;
          }
      }
      normalize();
  }

  int2048(const std::string &s) {
      read(s);
  }

  int2048(const int2048 &other) {
      data = other.data;
      neg = other.neg;
  }

  // The parameter types of the following functions are for reference only, you can choose to use constant references or not
  // If needed, you can add other required functions yourself
  // ===================================
  // Integer1
  // ===================================

  // Read a big integer
  void read(const std::string &s) {
      data.clear();
      neg = false;

      size_t start = 0;
      while (start < s.size() && (s[start] == ' ' || s[start] == '\t')) start++;

      if (start < s.size() && s[start] == '-') {
          neg = true;
          start++;
      } else if (start < s.size() && s[start] == '+') {
          start++;
      }

      while (start < s.size() && s[start] == '0' && start < s.size() - 1) start++;

      size_t len = s.size() - start;
      if (len == 0 || (start >= s.size())) {
          data.push_back(0);
          neg = false;
          return;
      }

      int num_chunks = (len + BASE_DIGITS - 1) / BASE_DIGITS;
      data.resize(num_chunks);

      size_t idx = s.size();
      for (int i = 0; i < num_chunks; i++) {
          size_t chunk_start;
          if (idx > start + BASE_DIGITS) {
              chunk_start = idx - BASE_DIGITS;
          } else {
              chunk_start = start;
          }

          size_t chunk_len = idx - chunk_start;
          long long val = 0;
          for (size_t j = chunk_start; j < idx; j++) {
              val = val * 10 + (s[j] - '0');
          }
          data[i] = val;
          idx = chunk_start;
      }

      normalize();
  }

  // Output the stored big integer, no need for newline
  void print() {
      if (neg && !(data.size() == 1 && data[0] == 0)) {
          std::putchar('-');
      }
      std::printf("%lld", data.back());
      for (int i = (int)data.size() - 2; i >= 0; i--) {
          std::printf("%09lld", data[i]);
      }
  }

  // Add a big integer
  int2048 &add(const int2048 &other) {
      if (neg == other.neg) {
          // Same sign: add absolute values
          long long carry = 0;
          for (size_t i = 0; i < data.size() || i < other.data.size() || carry; i++) {
              if (i >= data.size()) data.push_back(0);
              long long sum = data[i] + carry;
              if (i < other.data.size()) sum += other.data[i];
              data[i] = sum % BASE;
              carry = sum / BASE;
          }
      } else {
          // Different signs: subtract absolute values
          int cmp = abs_compare(other);
          if (cmp == 0) {
              data.clear();
              data.push_back(0);
              neg = false;
          } else if (cmp > 0) {
              // |this| > |other|
              long long borrow = 0;
              for (size_t i = 0; i < data.size(); i++) {
                  long long diff = data[i] - borrow;
                  if (i < other.data.size()) diff -= other.data[i];
                  if (diff < 0) {
                      diff += BASE;
                      borrow = 1;
                  } else {
                      borrow = 0;
                  }
                  data[i] = diff;
              }
          } else {
              // |this| < |other|
              std::vector<long long> new_data = other.data;
              long long borrow = 0;
              for (size_t i = 0; i < new_data.size(); i++) {
                  long long diff = new_data[i] - borrow;
                  if (i < data.size()) diff -= data[i];
                  if (diff < 0) {
                      diff += BASE;
                      borrow = 1;
                  } else {
                      borrow = 0;
                  }
                  new_data[i] = diff;
              }
              data = new_data;
              neg = !neg;
          }
      }

      normalize();
      return *this;
  }

  // Return the sum of two big integers
  friend int2048 add(int2048 a, const int2048 &b) {
      return a.add(b);
  }

  // Subtract a big integer
  int2048 &minus(const int2048 &other) {
      if (neg != other.neg) {
          // Different signs: add absolute values, keep sign of this
          long long carry = 0;
          for (size_t i = 0; i < data.size() || i < other.data.size() || carry; i++) {
              if (i >= data.size()) data.push_back(0);
              long long sum = data[i] + carry;
              if (i < other.data.size()) sum += other.data[i];
              data[i] = sum % BASE;
              carry = sum / BASE;
          }
      } else {
          // Same sign: subtract absolute values
          int cmp = abs_compare(other);
          if (cmp == 0) {
              data.clear();
              data.push_back(0);
              neg = false;
          } else if (cmp > 0) {
              // |this| > |other|
              long long borrow = 0;
              for (size_t i = 0; i < data.size(); i++) {
                  long long diff = data[i] - borrow;
                  if (i < other.data.size()) diff -= other.data[i];
                  if (diff < 0) {
                      diff += BASE;
                      borrow = 1;
                  } else {
                      borrow = 0;
                  }
                  data[i] = diff;
              }
          } else {
              // |this| < |other|
              std::vector<long long> new_data = other.data;
              long long borrow = 0;
              for (size_t i = 0; i < new_data.size(); i++) {
                  long long diff = new_data[i] - borrow;
                  if (i < data.size()) diff -= data[i];
                  if (diff < 0) {
                      diff += BASE;
                      borrow = 1;
                  } else {
                      borrow = 0;
                  }
                  new_data[i] = diff;
              }
              data = new_data;
              neg = !neg;
          }
      }

      normalize();
      return *this;
  }

  // Return the difference of two big integers
  friend int2048 minus(int2048 a, const int2048 &b) {
      return a.minus(b);
  }

  // ===================================
  // Integer2
  // ===================================

  int2048 operator+() const {
      return *this;
  }

  int2048 operator-() const {
      int2048 result = *this;
      if (!(result.data.size() == 1 && result.data[0] == 0)) {
          result.neg = !result.neg;
      }
      return result;
  }

  int2048 &operator=(const int2048 &other) {
      if (this != &other) {
          data = other.data;
          neg = other.neg;
      }
      return *this;
  }

  int2048 &operator+=(const int2048 &other) {
      return add(other);
  }

  friend int2048 operator+(int2048 a, const int2048 &b) {
      return a.add(b);
  }

  int2048 &operator-=(const int2048 &other) {
      return minus(other);
  }

  friend int2048 operator-(int2048 a, const int2048 &b) {
      return a.minus(b);
  }

  int2048 &operator*=(const int2048 &other) {
      bool result_negative = neg != other.neg;

      size_t n = data.size();
      size_t m = other.data.size();

      if (n == 0 || m == 0) {
          data.clear();
          data.push_back(0);
          neg = false;
          return *this;
      }

      // Use simple multiplication for small numbers
      if (n * m <= 10000) {
          data = multiply_simple(data, other.data);
          normalize();
          neg = result_negative && !(data.size() == 1 && data[0] == 0);
          return *this;
      }

      // Use Karatsuba for large numbers
      data = multiply_karatsuba(data, other.data);
      normalize();
      neg = result_negative && !(data.size() == 1 && data[0] == 0);
      return *this;
  }

  friend int2048 operator*(int2048 a, const int2048 &b) {
      return a *= b;
  }

  int2048 &operator/=(const int2048 &other) {
      int2048 q, r;
      divide(*this, other, q, r);

      bool q_neg = neg != other.neg;

      // For floor division: if result should be negative and there's a remainder,
      // we need to increase the magnitude by 1 (since floor division rounds toward -infinity)
      if (q_neg && (r.data.size() > 1 || r.data[0] != 0)) {
          // q should be -(q+1) for floor division when signs differ and there's remainder
          // Since q is non-negative from divide(), we add 1 and then apply negative sign
          long long carry = 1;
          for (size_t i = 0; i < q.data.size() && carry; i++) {
              q.data[i] += carry;
              carry = q.data[i] / BASE;
              q.data[i] %= BASE;
          }
          if (carry) q.data.push_back(carry);
          q.normalize();
          q.neg = true;  // Apply negative sign
      } else if (q_neg) {
          q.neg = true;  // Just apply negative sign, no adjustment needed
      }

      data = q.data;
      neg = q.neg && !(q.data.size() == 1 && q.data[0] == 0);
      normalize();

      return *this;
  }

  friend int2048 operator/(int2048 a, const int2048 &b) {
      return a /= b;
  }

  int2048 &operator%=(const int2048 &other) {
      // x % y = x - (x / y) * y
      int2048 quotient = *this / other;
      int2048 product = quotient * other;
      *this -= product;
      return *this;
  }

  friend int2048 operator%(int2048 a, const int2048 &b) {
      return a %= b;
  }

  friend std::istream &operator>>(std::istream &is, int2048 &x) {
      std::string s;
      is >> s;
      x.read(s);
      return is;
  }

  friend std::ostream &operator<<(std::ostream &os, const int2048 &x) {
      if (x.neg && !(x.data.size() == 1 && x.data[0] == 0)) {
          os << '-';
      }
      os << x.data.back();
      for (int i = (int)x.data.size() - 2; i >= 0; i--) {
          os << std::setfill('0') << std::setw(BASE_DIGITS) << x.data[i];
      }
      return os;
  }

  friend bool operator==(const int2048 &a, const int2048 &b) {
      return a.neg == b.neg && a.data == b.data;
  }

  friend bool operator!=(const int2048 &a, const int2048 &b) {
      return !(a == b);
  }

  friend bool operator<(const int2048 &a, const int2048 &b) {
      bool a_zero = a.data.size() == 1 && a.data[0] == 0;
      bool b_zero = b.data.size() == 1 && b.data[0] == 0;

      if (a_zero && b_zero) return false;

      if (a.neg && !b.neg) return true;
      if (!a.neg && b.neg) return false;

      int cmp = a.abs_compare(b);
      if (a.neg) {
          return cmp > 0;
      }
      return cmp < 0;
  }

  friend bool operator>(const int2048 &a, const int2048 &b) {
      return b < a;
  }

  friend bool operator<=(const int2048 &a, const int2048 &b) {
      return !(a > b);
  }

  friend bool operator>=(const int2048 &a, const int2048 &b) {
      return !(a < b);
  }
};

} // namespace sjtu

#endif
