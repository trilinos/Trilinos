/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/environment/FormatMemorySize.hpp>
#include <boost/lexical_cast.hpp>       // for lexical_cast
#include <sstream>                      // for basic_stringbuf<>::int_type, etc



namespace stk {

std::string
formatMemorySize(
  double                size)
{
  std::string           result;
  
  static const double kb = 1024.0;
  // static const double mb = kb * kb;
  // static const double gb = kb * kb * kb;

  if (size < 0.0) {
    result = "-";
    size = -size;
  }

  // output size in kilo bytes
  result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size / kb));
  result += " KB";
  // if (size < kb) {
  //   // output size in bytes
  //   result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size));
  //   result += " B";
  // }
  // else if (size < mb) {
  //   // output size in kilo bytes
  //   result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size / kb));
  //   result += " KB";
  // }
  // else if (size < gb) {
  //   // output size in mega bytes
  //   result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size / mb));
  //   result += " MB";
  // }
  // else {
  //   // everything else output in giga bytes
  //   result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size / gb));
  //   result += " GB";
  // }
  
  return result;
}


std::string
formatMemorySize(
  MemorySize            size)
{
  std::string           result;
  
  static const MemorySize kb = 1024;
  // static const MemorySize mb = kb * kb;
  // static const MemorySize gb = kb * kb * kb;

  // output size in kilo bytes
  result = boost::lexical_cast<std::string>(size / kb);
  result += " KB";
  
  // if (size < kb) {
  //   // output size in bytes
  //   result = boost::lexical_cast<std::string>(size);
  //   result += " B";
  // }
  // else if (size < mb) {
  //   // output size in kilo bytes
  //   result = boost::lexical_cast<std::string>(size / kb);
  //   result += " KB";
  // }
  // else if (size < gb) {
  //   // output size in mega bytes
  //   result = boost::lexical_cast<std::string>(size / mb);
  //   result += " MB";
  // }
  // else {
  //   // everything else output in giga bytes
  //   result = boost::lexical_cast<std::string>(size / gb);
  //   result += " GB";
  // }
  
  return result;
}

} // namespace stk
