// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_PrintDouble.hpp"
#include "Teuchos_BigUInt.hpp"

#include <cstring>

namespace Teuchos {

namespace {

int ndigits_for(std::uint32_t x) {
  int n = 0;
  while (x) {
    ++n;
    x /= 10;
  }
  return n;
}

}

/**

\note This function is based on the journal article:

Burger, Robert G., and R. Kent Dybvig.
"Printing floating-point numbers quickly and accurately."
ACM SIGPLAN Notices. Vol. 31. No. 5. ACM, 1996.

Burger and Dybvig describe an algorithm for printing
a decimal representation of a double-precision floating-point
number using the minimum number of digits such that
the value is bitwise exactly preserved when scanned back in.

Their work is meant to improve on the work of Steele and White
in terms of runtime performance.
It requires high-precision integer arithmetic, hence the
Teuchos::BigUInt class was developed to support this implementation.

Variables in this code are meant to closely reflect the variables
in Burger and Dybvig's work.

Note that there is an error in the original paper which is corrected
in our code:
On page 4, step 4 of the Integer Arithmetic procedure, stopping
condition (2) has the wrong sign, and should be:
$r_n + m_n^+ > s_n$

We extend their work slightly to include an intelligent choice of
leading/trailing zeros versus scientific notation, choosing whichever
minimizes the number of characters.
We also wrote this implementation such that the resulting string
cannot be confused with an integer, by adding decimal points
(e.g. "5." instead of "5").
*/

void print_double(std::ostream& os, double v) {
  char buffer[64];
  constexpr std::uint64_t one = 1;
  constexpr std::uint64_t p = 53;
  std::uint64_t pun;
  std::memcpy(&pun, &v, sizeof(v));
  auto sign = pun >> 63;
  pun -= sign << 63;
  auto be = pun >> (p - 1);
  pun -= be << (p - 1);
  auto m = pun;
  int bp = 0;
  if (be == 2047) {
    if (m == 0) {
      if (sign) buffer[bp++] = '-';
      buffer[bp++] = 'i';
      buffer[bp++] = 'n';
      buffer[bp++] = 'f';
    } else {
      buffer[bp++] = 'n';
      buffer[bp++] = 'a';
      buffer[bp++] = 'n';
    }
  } else {
    if (sign) buffer[bp++] = '-';
    std::uint64_t f;
    if (be == 0) {
      f = m;
    } else {
      f = m + (one << (p - 1));
    }
    auto e = std::int64_t(be) - 1075;
    BigUInt<34> r, s, mp, mm;
    if (e >= 0) {
      if (f != (one << (p - 1))) {
        r = BigUInt<34>(f);
        r <<= (e + 1);
        s = BigUInt<34>(2);
        mp = BigUInt<34>(1);
        mp <<= e;
        mm = BigUInt<34>(1);
        mm <<= e;
      } else {
        r = BigUInt<34>(f);
        r <<= (e + 2);
        s = BigUInt<34>(2 * 2);
        mp = BigUInt<34>(1);
        mp <<= (e + 1);
        mm = BigUInt<34>(1);
        mm <<= e;
      }
    } else {
      if ((be == 0) || (f != (one << (p - 1)))) {
        r = BigUInt<34>(f);
        r <<= 1;
        s = BigUInt<34>(1);
        s <<= (1 - e);
        mp = BigUInt<34>(1);
        mm = BigUInt<34>(1);
      } else {
        r = BigUInt<34>(f);
        r <<= 2;
        s = BigUInt<34>(1);
        s <<= (2 - e);
        mp = BigUInt<34>(2);
        mm = BigUInt<34>(1);
      }
    }
    std::int32_t k = 0;
    BigUInt<34> B_k{1};
    auto r_p_mp = r + mp;
    auto r_p_mp_comp = comp(r_p_mp, s);
    if (r_p_mp_comp == 0) {
    } else if (r_p_mp_comp == 1) {
      while (r_p_mp > (s * B_k)) {
        ++k, B_k *= 10;
      }
    } else {
      while ((r_p_mp * B_k) < s) {
        --k, B_k *= 10;
      }
      ++k;
      B_k = B_k / 10;
    }
    if (k >= 0) {
      s = s * B_k;
    } else {
      r = r * B_k;
      mp = mp * B_k;
      mm = mm * B_k;
    }
    char last_d = '0';
    int n;
    for (n = 0; true; ++n) {
      auto r_x_10 = r;
      r_x_10 *= 10;
      auto d_np1 = r_x_10 / s;
      auto cond1 = r < mm;
      auto cond2 = (r + mp) > s;
      if (cond1 && cond2) {
        r <<= 1;
        if (r < s) {
          buffer[bp++] = last_d;
        } else {
          buffer[bp++] = last_d + 1;
        }
        break;
      } else if (cond1) {
        buffer[bp++] = last_d;
        break;
      } else if (cond2) {
        buffer[bp++] = last_d + 1;
        break;
      } else {
        if (n) buffer[bp++] = last_d;
        r = r_x_10;
        r -= (s * d_np1);
        mp *= 10;
        mm *= 10;
        last_d = char(d_np1[0]) + '0';
      }
    }
    if (v == 0.0) {
      k = 1;
      ++n;
    }
    int dot_pos = -1;
    bool do_scientific = false;
    if (0 <= k && k <= n) {
      // dot is touching significant digits
      dot_pos = k;
    } else if (k < 0) {
      auto nchars_sci = ndigits_for(-k + 1) + 2;
      if (n > 1) nchars_sci += 1; // add a dot to scientific notation if more than one digit
      if (nchars_sci < (-k + 1)) {
        // scientific notation requires fewer chars than trailing zeros
        if (n > 1) dot_pos = 1;
        do_scientific = true;
      } else {
        // trailing zeros are no more chars than scientific
        for (int i = 0; i < n; ++i) {
          buffer[bp + (-k) - i - 1] = buffer[bp - i - 1];
        }
        for (int i = 0; i < -k; ++i) {
          buffer[bp - n + i] = '0';
        }
        dot_pos = bp - n;
        bp += -k;
        n += -k;
      }
    } else if (k > n) {
      auto nchars_sci = ndigits_for(k - 1) + 1;
      if (n > 1) nchars_sci += 1; // add a dot to scientific notation if more than one digit
      if (nchars_sci < ((k-n)+1)) {
        // scientific notation requires fewer chars than trailing zeros
        if (n > 1) dot_pos = 1;
        do_scientific = true;
      } else {
        // trailing zeros are no more chars than scientific
        for (; n < k; ++n) buffer[bp++] = '0';
        dot_pos = n;
      }
    }
    if (dot_pos != -1) {
      for (int i = 0; i < (n - dot_pos) && i < bp; ++i) buffer[bp - i] = buffer[bp - i - 1];
      buffer[bp - n + dot_pos] = '.';
      ++bp;
    }
    if (do_scientific) {
      buffer[bp++] = 'e';
      auto decimal_exponent = (k - 1);
      if (decimal_exponent < 0) {
        buffer[bp++] = '-';
        decimal_exponent = -decimal_exponent;
      }
      int ne;
      for (ne = 0; decimal_exponent; ++ne) {
        buffer[bp++] = char(decimal_exponent % 10) + '0';
        decimal_exponent /= 10;
      }
      for (int i = 0; i < ne / 2; ++i) {
        auto tmp = buffer[bp - ne + i];
        buffer[bp - ne + i] = buffer[bp - i - 1];
        buffer[bp - i - 1] = tmp;
      }
    }
  }
  buffer[bp] = '\0';
  os << buffer;
}

}
