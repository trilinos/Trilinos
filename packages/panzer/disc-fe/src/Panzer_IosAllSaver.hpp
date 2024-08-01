// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
