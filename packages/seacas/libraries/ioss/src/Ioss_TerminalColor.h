#ifndef IOSS__TRMCLR_H__
#define IOSS__TRMCLR_H__

// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <cstdint>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>

// ## Example

// ```c++
// #include "trmclr.hpp"
// #include <iostream>

// int main()
// {
//     trmclr::Style fancyStyle(trmclr::Background::LIGHT_BLUE    |
//                              trmclr::Foreground::WHITE         |
//                              trmclr::Attribute::UNDERLINED     |
//                              trmclr::Attribute::BOLD);

//     trmclr::Style basicStyle(trmclr::Attribute::DEFAULT);

//     std::cout << fancyStyle << "Hello "
//               << basicStyle << "World!\n";

//     return 0;
// }

// /*
// Note you can also do things like:
// auto bold = [](trmclr::Style style) { return trmclr::Style(style | trmclr::Attribute::BOLD); };
// */
// `

namespace Ioss {
  namespace trmclr {

    struct Style
    {
      explicit Style(uint32_t value) : _value(value) {}

      operator uint32_t() const { return _value; }

      uint32_t _value;
    };

    enum StyleTypes { FOREGROUND, ATTRIBUTE, BACKGROUND, N_STYLE_TYPES };

    static const uint32_t STYLE_SHIFT = std::numeric_limits<uint32_t>::digits / N_STYLE_TYPES;

    struct Attribute
    {
      static const uint32_t SHIFT = STYLE_SHIFT * ATTRIBUTE;

      enum {
        DEFAULT    = 0x001 << SHIFT,
        BOLD       = 0x002 << SHIFT,
        DIM        = 0x004 << SHIFT,
        UNDERLINED = 0x010 << SHIFT,
        BLINK      = 0x020 << SHIFT,
        REVERSE    = 0x080 << SHIFT,
        HIDDEN     = 0x100 << SHIFT
      };
    };

    struct Foreground
    {
      static const uint32_t SHIFT = STYLE_SHIFT * FOREGROUND;

      enum {
        BLACK         = 30 << SHIFT,
        RED           = 31 << SHIFT,
        GREEN         = 32 << SHIFT,
        YELLOW        = 33 << SHIFT,
        BLUE          = 34 << SHIFT,
        MAGENTA       = 35 << SHIFT,
        CYAN          = 36 << SHIFT,
        LIGHT_GRAY    = 37 << SHIFT,
        DEFAULT       = 39 << SHIFT,
        DARK_GRAY     = 90 << SHIFT,
        LIGHT_RED     = 91 << SHIFT,
        LIGHT_GREEN   = 92 << SHIFT,
        LIGHT_YELLOW  = 93 << SHIFT,
        LIGHT_BLUE    = 94 << SHIFT,
        LIGHT_MAGENTA = 95 << SHIFT,
        LIGHT_CYAN    = 96 << SHIFT,
        WHITE         = 97 << SHIFT
      };
    };

    struct Background
    {
      static const uint32_t SHIFT = STYLE_SHIFT * BACKGROUND;

      enum {
        BLACK         = 40 << SHIFT,
        RED           = 41 << SHIFT,
        GREEN         = 42 << SHIFT,
        YELLOW        = 43 << SHIFT,
        BLUE          = 44 << SHIFT,
        MAGENTA       = 45 << SHIFT,
        CYAN          = 46 << SHIFT,
        LIGHT_GRAY    = 47 << SHIFT,
        DEFAULT       = 49 << SHIFT,
        DARK_GRAY     = 100 << SHIFT,
        LIGHT_RED     = 101 << SHIFT,
        LIGHT_GREEN   = 102 << SHIFT,
        LIGHT_YELLOW  = 103 << SHIFT,
        LIGHT_BLUE    = 104 << SHIFT,
        LIGHT_MAGENTA = 105 << SHIFT,
        LIGHT_CYAN    = 106 << SHIFT,
        WHITE         = 107 << SHIFT
      };
    };

    static Style black(Ioss::trmclr::Foreground::BLACK);
    static Style red(Ioss::trmclr::Foreground::RED);
    static Style green(Ioss::trmclr::Foreground::GREEN);
    static Style yellow(Ioss::trmclr::Foreground::YELLOW);
    static Style blue(Ioss::trmclr::Foreground::BLUE);
    static Style magenta(Ioss::trmclr::Foreground::MAGENTA);
    static Style cyan(Ioss::trmclr::Foreground::CYAN);
    static Style normal(Ioss::trmclr::Attribute::DEFAULT);

    inline std::ostream &operator<<(std::ostream &os, const Style &style)
    {
      const uint32_t base    = 1 << STYLE_SHIFT;
      uint32_t       encoded = style / base;
      uint32_t       decoded = style % base;

      os << "\x1B[" << (decoded ? decoded : Foreground::DEFAULT >> Foreground::SHIFT);

      decoded = encoded % base;

      for (uint32_t i = 0; decoded != 0; decoded >>= 1, i++) {
        if (decoded & 1) {
          os << ";" << i;
        }
      }

      encoded = encoded / base;
      decoded = encoded % base;

      os << ";" << (decoded ? decoded : Background::DEFAULT >> Background::SHIFT) << "m";

      return os;
    }

  } // namespace trmclr
} // namespace Ioss
#endif // end __TRMCLR_H__
