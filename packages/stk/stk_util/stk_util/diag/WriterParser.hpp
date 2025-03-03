// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
// 

#ifndef STK_UTIL_Diag_WriterParser_h
#define STK_UTIL_Diag_WriterParser_h

#include "stk_util/diag/Option.hpp"  // for OptionMaskParser, OptionMaskParser::Mask
#include <string>                    // for string

namespace stk {
namespace diag {

///
/// @addtogroup DiagWriterDetail
/// @{
///

/**
 * @brief Class <b>WriterParser</b> implements a parser a Writer PrintMask string.
 *
 */
class WriterParser : public OptionMaskParser
{
public:
  /**
   * @brief Typedef <b>Mask</b> bring the OptionMaskParser Mask definition into this
   * namespace.
   *
   */
  typedef OptionMaskParser::Mask  Mask;

public:
  /**
   * @brief Creates a new <b>WriterParser</b> instance containing the lowerest
   * level PrintMask names.
   *
   */
  WriterParser();

  /**
   * @brief Member function <b>parse</b> returns the mask which results from parsing the
   * <b>mask_string</b>.
   *
   * @param mask_string    a <b>std::string</b> const reference to the string to be
   *        parsed.
   *
   * @return      a <b>Mask</b> value of the result from parsing the mask
   *        string.
   */
  Mask parse(const char *mask_string) const override;

  /**
   * @brief Member function <b>parseArg</b> parses the argument and its argument
   * values.
   *
   * @param name    a <b>std::string</b> const reference to the argument
   *        name.
   *
   * @param arg      a <b>std::string</b> const reference to the argument
   *        values.
   */
  virtual void parseArg(const std::string &name, const std::string &arg) const override;
};

///
/// @}
///

} // namespace diag
} // namespace stk

namespace sierra {
namespace Diag {

typedef stk::diag::WriterParser WriterParser;

} // namespace Diag
} // namespace sierra

#endif // STK_UTIL_Diag_WriterParser_h
