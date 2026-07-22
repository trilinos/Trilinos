// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_O_BLACK_HOLE_STREAM_H
#define TEUCHOS_O_BLACK_HOLE_STREAM_H

#include "Teuchos_basic_oblackholestream.hpp"

namespace Teuchos {
  /** \brief .
   * \ingroup teuchos_outputting_grp
   */
	typedef basic_oblackholestream<char,std::char_traits<char> >   oblackholestream;
}

#endif // TEUCHOS_O_BLACK_HOLE_STREAM_H
