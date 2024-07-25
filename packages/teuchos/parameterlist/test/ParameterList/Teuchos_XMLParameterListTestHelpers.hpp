// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
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
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
RCP<ParameterList> writeThenReadPL(ParameterList& myList, RCP<DependencySheet> depSheetIn,
  RCP<DependencySheet> depSheetOut);


} // namespace Teuchos


#endif // TEUCHOS_XML_PARAMETER_LIST_TEST_HELPERS_HPP
