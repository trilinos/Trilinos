// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_UNIT_TEST_BASE_HPP
#define TEUCHOS_UNIT_TEST_BASE_HPP


/*! \file Teuchos_UnitTestBase.hpp
    \brief Unit testing support.
*/


#include "Teuchos_Describable.hpp"
#include "Teuchos_FancyOStream.hpp"


namespace Teuchos {


/** \brief Unit test base class. */
class TEUCHOSCORE_LIB_DLL_EXPORT UnitTestBase : public Describable {
public:

  /** \brief . */
  UnitTestBase(const std::string groupName, std::string testName);

  /** \brief . */
  bool runUnitTest(FancyOStream &out) const;

  /** \brief . */
  virtual std::string unitTestFile() const = 0;

  /** \brief . */
  virtual long int unitTestFileLineNumber() const = 0;

protected:

  /** \brief . */
  virtual void runUnitTestImpl(FancyOStream &out, bool &success) const = 0;

};


} // namespace Teuchos


#endif  // TEUCHOS_UNIT_TEST_BASE_HPP
