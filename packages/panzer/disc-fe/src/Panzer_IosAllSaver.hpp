// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

// *******************************************************************
// This file contains a copy of one function from the io state saver
// code from boost.  Modified to use the panzer namespace and remove
// locale support. Boost copyright is below.
// *******************************************************************

//  Copyright 2002, 2005 Daryle Walker.  Use, modification, and distribution
//  are subject to the Boost Software License, Version 1.0.  (See accompanying
//  file LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/io/> for the library's home page.

// ******************************************************************* 
// ******************************************************************* 

#ifndef PANZER_IOS_ALL_SAVER_HPP
#define PANZER_IOS_ALL_SAVER_HPP

#include "PanzerDiscFE_config.hpp"
#include <iomanip>
#include <iosfwd>  // for std::char_traits (declaration)

namespace panzer {

  template < typename Ch, class Tr = ::std::char_traits<Ch> >
  class basic_ios_all_saver;

  typedef basic_ios_all_saver<char>            ios_all_saver;
  typedef basic_ios_all_saver<wchar_t>        wios_all_saver;

  template < typename Ch, class Tr >
  class basic_ios_all_saver
  {
  public:
    typedef ::std::basic_ios<Ch, Tr>  state_type;

    explicit  basic_ios_all_saver( state_type &s )
      : s_save_( s ), a1_save_( s.flags() ), a2_save_( s.precision() )
      , a3_save_( s.width() ), a4_save_( s.rdstate() )
      , a5_save_( s.exceptions() ), a6_save_( s.tie() )
      , a7_save_( s.rdbuf() ), a8_save_( s.fill() )
    {}

    ~basic_ios_all_saver()
    { this->restore(); }

    void  restore()
    {
      s_save_.fill( a8_save_ );
      s_save_.rdbuf( a7_save_ );
      s_save_.tie( a6_save_ );
      s_save_.exceptions( a5_save_ );
      s_save_.clear( a4_save_ );
      s_save_.width( a3_save_ );
      s_save_.precision( a2_save_ );
      s_save_.flags( a1_save_ );
    }

  private:
    state_type &                            s_save_;
    typename state_type::fmtflags const     a1_save_;
    ::std::streamsize const                 a2_save_;
    ::std::streamsize const                 a3_save_;
    typename state_type::iostate const      a4_save_;
    typename state_type::iostate const      a5_save_;
    ::std::basic_ostream<Ch, Tr> * const    a6_save_;
    ::std::basic_streambuf<Ch, Tr> * const  a7_save_;
    typename state_type::char_type const    a8_save_;

    basic_ios_all_saver& operator=(const basic_ios_all_saver&);
  };

}

#endif
