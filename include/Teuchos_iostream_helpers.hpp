// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_IOSTREAM_HELPERS_HPP
#define TEUCHOS_IOSTREAM_HELPERS_HPP


/*! \file Teuchos_VerbosityLevel.hpp
    \brief .
*/


#include "Teuchos_Assert.hpp"


namespace Teuchos {


template<class EnumType>
std::istream& enumIstreamExtractor(std::istream& std_is, EnumType& enum_value)
{
  int intval;
  std_is >> intval;
  enum_value = static_cast<EnumType>(intval);
  return std_is;
}


} // namespace Teuchos


#define TEUCHOS_ENUM_INPUT_STREAM_OPERATOR(ENUMTYPE) \
inline \
std::istream& operator>>(std::istream& std_is, ENUMTYPE& enum_value) \
{ return Teuchos::enumIstreamExtractor(std_is, enum_value); }


#endif // TEUCHOS_IOSTREAM_HELPERS_HPP
