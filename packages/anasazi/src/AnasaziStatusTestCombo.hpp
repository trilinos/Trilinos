
// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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
//

#ifndef ANASAZI_STATUS_TEST_COMBO_HPP
#define ANASAZI_STATUS_TEST_COMBO_HPP

/*!
  \file AnasaziStatusTestCombo.hpp
  \brief Status test for forming logical combinations of other status tests.
*/


#include "AnasaziTypes.hpp"
#include "AnasaziStatusTest.hpp"
#include "Teuchos_Array.hpp"

  /*! 
    \class Anasazi::StatusTestCombo
    \brief Status test for forming logical combinations of other status tests.
    
    Test types include OR, AND, SEQOR and SEQAND.  The OR and AND tests
    evaluate all of the tests, in the order they were passed to the
    StatusTestCombo.  The SEQOR and SEQAND run only the tests necessary to
    determine the final outcome, short-circuiting on the first test that
    conclusively decides the outcome. More formally, SEQAND runs the tests in
    the order they were given to the StatusTestCombo class and stops after the
    first test that evaluates Undefined. SEQOR run the tests in the order they
    were given to the StatusTestCombo class and stops after the first test that
    evaluates Passed.
  */

namespace Anasazi {


template <class ScalarType, class MV, class OP>
class StatusTestCombo : public StatusTest<ScalarType,MV,OP> {

 private:
   typedef Teuchos::Array< Teuchos::RefCountPtr< StatusTest<ScalarType,MV,OP> > > STPArray;

 public:

 //!  \brief Enumerated type to list the types of StatusTestCombo combo types.
 enum ComboType
   {
     OR,           /*!< Logical OR which evaluates all tests */
     AND,          /*!< Logical AND which evaluates all tests */
     SEQOR,        /*!< Short-circuited logical OR */
     SEQAND        /*!< Short-circuited logical AND */
   };


#ifndef DOXYGEN_SHOULD_SKIP_THIS

  typedef Teuchos::Array< Teuchos::RefCountPtr< StatusTest<ScalarType,MV,OP> > > t_arr;
  typedef std::vector< Teuchos::RefCountPtr< StatusTest<ScalarType,MV,OP> > > st_vector;
  typedef typename st_vector::iterator                 iterator;
  typedef typename st_vector::const_iterator     const_iterator;

#endif // DOXYGEN_SHOULD_SKIP_THIS

  //! @name Constructors/destructors
  //@{ 

  //! Constructor
  //! \brief Default constructor has no tests and initializes to ComboType OR.
  StatusTestCombo() : _state(Undefined) {}

  //! Constructor
  //! \brief Constructor specifying the ComboType and the tests.
  StatusTestCombo(ComboType type, Teuchos::Array< Teuchos::RefCountPtr< StatusTest<ScalarType,MV,OP> > > tests) :
    _state(Undefined), 
    _type(type)
  {
    setTests(tests);
  };

  //! Destructor
  virtual ~StatusTestCombo() {};
  //@}

  //! @name Status methods
  //@{ 
  /*! Check status as defined by test.
    
    \return TestStatus indicating whether the test passed or failed.
  */
  TestStatus checkStatus( Eigensolver<ScalarType,MV,OP>* solver );

  //! Return the result of the most recent checkStatus call.
  TestStatus getStatus() const {
    return _state;
  }
  //@}

  //! @name Accessor methods
  //@{ 

  /*! \brief Set the maximum number of iterations.
   *  This also resets the test status to Undefined.
   */
  void setComboType(ComboType type) {
    _type = type;
    _state = Undefined;
  }

  //! Get the maximum number of iterations.
  ComboType getComboType() const {return _type;}

  /*! \brief Set the tests
   *  This also resets the test status to Undefined.
   */
  void setTests(Teuchos::Array<Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > > tests) {
    _tests = tests;
    _state = Undefined;
  }

  //! Get the tests
  Teuchos::Array<Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > > getTests() const {return _tests;}

  /*! \brief Add a test to the combination.
   *
   *  This also resets the test status to Undefined.
   */
  void addTest(Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > test) {
    _tests.push_back(test);
    _state = Undefined;
  }

  /*! \brief Removes a test from the combination, if it exists in the tester.
   *
   * This also resets the test status to Undefined, if a test was removed.
   */
  void removeTest(const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > &test);

  //@}

  //! @name Reset methods
  //@{ 
  //! \brief Informs the status test that it should reset its internal configuration to the uninitialized state.
  /*! The StatusTestCombo class has no internal state, but children classes might, so this method will call
     reset() on all child status tests. It also resets the test status to Undefined.
  */
  void reset();

  //! \brief Clears the results of the last status test.
  /*! This should be distinguished from the reset() method, as it only clears the cached result from the last 
   * status test, so that a call to getStatus() will return Undefined. This is necessary for the SEQOR and SEQAND
   * tests in the StatusTestCombo class, which may short circuit and not evaluate all of the StatusTests contained
   * in them. It also resets the test status to Undefined.
  */
  void clearStatus();

  //@}

  //! @name Print methods
  //@{ 
  
  //! Output formatted description of stopping test to output stream.
  ostream& print(ostream& os, int indent = 0) const;
 
  //@}
  private:

  TestStatus evalOR(Eigensolver<ScalarType,MV,OP>* solver);
  TestStatus evalAND(Eigensolver<ScalarType,MV,OP>* solver);
  TestStatus evalSEQOR(Eigensolver<ScalarType,MV,OP>* solver);
  TestStatus evalSEQAND(Eigensolver<ScalarType,MV,OP>* solver);

  TestStatus _state;
  ComboType _type;
  STPArray _tests;

};


template <class ScalarType, class MV, class OP>
void StatusTestCombo<ScalarType,MV,OP>::removeTest(const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > &test) 
{
  typename STPArray::iterator iter1;
  iter1 = find(_tests.begin(),_tests.end(),test);
  if (iter1 != _tests.end()) {
    _tests.erase(iter1);
    _state = Undefined;
  }
}


template <class ScalarType, class MV, class OP>
TestStatus StatusTestCombo<ScalarType,MV,OP>::checkStatus( Eigensolver<ScalarType,MV,OP>* solver ) {
  clearStatus();
  switch (_type) {
    case OR:
      _state = evalOR(solver);
      break;
    case AND:
      _state = evalAND(solver);
      break;
    case SEQOR:
      _state = evalSEQOR(solver);
      break;
    case SEQAND:
      _state = evalSEQAND(solver);
      break;
  }
  return _state;
}


template <class ScalarType, class MV, class OP>
void StatusTestCombo<ScalarType,MV,OP>::reset() {
  _state = Undefined;
  for (iterator i=_tests.begin(); i != _tests.end(); i++) {
    (*i)->reset();
  }
}

template <class ScalarType, class MV, class OP>
void StatusTestCombo<ScalarType,MV,OP>::clearStatus() {
  _state = Undefined;
  for (iterator i=_tests.begin(); i != _tests.end(); i++) {
    (*i)->clearStatus();
  }
}

template <class ScalarType, class MV, class OP>
ostream& StatusTestCombo<ScalarType,MV,OP>::print(ostream& os, int indent) const {
  string ind(indent,' ');
  os << ind << "ComboTest: ";
  switch (_state) {
  case Passed:
    os << "Passed" << endl;
    break;
  case Failed:
    os << "Failed" << endl;
    break;
  case Undefined:
    os << "Undefined" << endl;
    break;
  }
  // print children, with extra indention
  for (const_iterator i=_tests.begin(); i != _tests.end(); i++) {
    (*i)->print(os,indent+2);
  }
  return os;
}

template <class ScalarType, class MV, class OP>
TestStatus StatusTestCombo<ScalarType,MV,OP>::evalOR( Eigensolver<ScalarType,MV,OP>* solver ) {
  _state = Failed;
  for (iterator i=_tests.begin(); i != _tests.end(); i++) {
    TestStatus r = (*i)->checkStatus(solver);
    if (r == Passed) {
      _state = Passed;
    }
    else if (r != Failed) {
      TEST_FOR_EXCEPTION(true,StatusTestError,
                         "Anasazi::StatusTestCombo::evalOR(): child test gave invalid return");
    }
  }
  return _state;
}

template <class ScalarType, class MV, class OP>
TestStatus StatusTestCombo<ScalarType,MV,OP>::evalSEQOR( Eigensolver<ScalarType,MV,OP>* solver ) {
  _state = Failed;
  for (iterator i=_tests.begin(); i != _tests.end(); i++) {
    TestStatus r = (*i)->checkStatus(solver);
    if (r == Passed) {
      _state = Passed;
      break;
    }
    else if (r != Failed) {
      TEST_FOR_EXCEPTION(true,StatusTestError,
                         "Anasazi::StatusTestCombo::evalSEQOR(): child test gave invalid return");
    }
  }
  return _state;
}

template <class ScalarType, class MV, class OP>
TestStatus StatusTestCombo<ScalarType,MV,OP>::evalAND( Eigensolver<ScalarType,MV,OP>* solver ) {
  _state = Passed;
  for (iterator i=_tests.begin(); i != _tests.end(); i++) {
    TestStatus r = (*i)->checkStatus(solver);
    if (r == Failed) {
      _state = Failed;
    }
    else if (r != Passed) {
      TEST_FOR_EXCEPTION(true,StatusTestError,
                         "Anasazi::StatusTestCombo::evalAND(): child test gave invalid return");
    }
  }
  return _state;
}

template <class ScalarType, class MV, class OP>
TestStatus StatusTestCombo<ScalarType,MV,OP>::evalSEQAND( Eigensolver<ScalarType,MV,OP>* solver ) {
  _state = Passed;
  for (iterator i=_tests.begin(); i != _tests.end(); i++) {
    TestStatus r = (*i)->checkStatus(solver);
    if (r == Failed) {
      _state = Failed;
      break;
    }
    else if (r != Passed) {
      TEST_FOR_EXCEPTION(true,StatusTestError,
                         "Anasazi::StatusTestCombo::evalAND(): child test gave invalid return");
    }
  }
  return _state;
}



} // end of Anasazi namespace

#endif /* ANASAZI_STATUS_TEST_COMBO_HPP */
