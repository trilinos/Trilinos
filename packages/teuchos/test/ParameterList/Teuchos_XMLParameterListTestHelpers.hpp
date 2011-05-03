// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_XML_PARAMETER_LIST_TEST_HELPERS_HPP
#define TEUCHOS_XML_PARAMETER_LIST_TEST_HELPERS_HPP


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DependencySheet.hpp"


namespace Teuchos {


/** \brief Write a parameter list to xml and then read that xml back in via
 * a string. The intent of this function is to be used for testing purposes.
 *
 * \param paramList [in] Contains the parameters and sublists that will be
 * written out and then read back in.
 *
 * \return The read in parameter list.
 * \ingroup XML
 */
TEUCHOS_LIB_DLL_EXPORT
RCP<ParameterList> writeThenReadPL(ParameterList& myList);


/** \brief Write a parameter list to xml and then read that xml back in via
 * a string. The intent of this function is to be used for testing purposes.
 *
 * \param paramList [in] Contains the parameters and sublists that will be
 * written out and then read back in.
 *
 * \param depSheetIn [in] The Dependency Sheet from which Dependencies should be
 * should be written.
 * \param depSheetOut [out] The Dependency Sheet into which Dependencies should 
 * be placed once read.
 * 
 * \return The read in parameter list.
 * \ingroup XML
 */
TEUCHOS_LIB_DLL_EXPORT
RCP<ParameterList> writeThenReadPL(ParameterList& myList, RCP<DependencySheet> depSheetIn,
  RCP<DependencySheet> depSheetOut);


} // namespace Teuchos


#endif // TEUCHOS_XML_PARAMETER_LIST_TEST_HELPERS_HPP
