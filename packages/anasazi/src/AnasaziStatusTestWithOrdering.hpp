// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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

#ifndef ANASAZI_STATUS_TEST_ORDEREDRESNORM_HPP
#define ANASAZI_STATUS_TEST_ORDEREDRESNORM_HPP

/*!
  \file AnasaziStatusTestWithOrdering.hpp
  \brief A status test for testing the norm of the eigenvectors residuals along with a 
         set of auxiliary eigenvalues.
*/


#include "AnasaziStatusTest.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_LAPACK.hpp"

  /*! 
    \class Anasazi::StatusTestWithOrdering 
    
    \brief A status test for testing the norm of the eigenvectors residuals
    along with a set of auxiliary eigenvalues. 
    
    The test evaluates to ::Passed when then the most significant of the
    eigenvalues all have a residual below a certain threshhold. The purpose of
    the test is to not only test convergence for some number of eigenvalues,
    but to test convergence for the correct ones.
    
    In addition to specifying the tolerance, the user may specify:
    <ul>
      <li> the norm to be used: 2-norm or OrthoManager::norm() or getRitzRes2Norms()
      <li> the scale: absolute or relative to magnitude of Ritz value 
      <li> the quorum: the number of vectors required for the test to 
           evaluate as ::Passed.
    </ul>

    Finally, the user must specify the Anasazi::SortManager used for deciding
    significance. 
  */

namespace Anasazi {


template <class ScalarType, class MV, class OP>
class StatusTestWithOrdering : public StatusTest<ScalarType,MV,OP> {

 private:
  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
  typedef Teuchos::ScalarTraits<MagnitudeType> MT;

 public:

  //! @name Constructors/destructors
  //@{ 

  //! Constructor
  StatusTestWithOrdering(Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test, Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > sorter, int quorum = -1);

  //! Destructor
  virtual ~StatusTestWithOrdering() {};
  //@}

  //! @name Status methods
  //@{ 
  /*! Check status as defined by test.

    \return TestStatus indicating whether the test passed or failed.
  */
  TestStatus checkStatus( Eigensolver<ScalarType,MV,OP>* solver );

  //! Return the result of the most recent checkStatus call, or undefined if it has not been run.
  TestStatus getStatus() const { return state_; }

  //! Get the indices for the vectors that passed the test.
  /*!
   * Non-negative indices correspond to passing vectors from the constituent status test. 
   * Negative entries correspond to auxilliary values, where the first auxilliary value
   * is indexed by -NumAuxVals, the second by -NumAuxVals+1, and so forth.
   */
  std::vector<int> whichVecs() const {
    return ind_;
  }

  //! Get the number of vectors that passed the test.
  int howMany() const {
    return ind_.size();
  }

  //@}

  //! @name Accessor methods
  //@{ 

  /*! \brief Set quorum.
   *
   *  Setting quorum to -1 signifies that all residuals from the solver must meet the tolerance.
   *  This also resets the test status to ::Undefined.
   */
  void setQuorum(int quorum) {
    state_ = Undefined;
    quorum_ = quorum;
  }

  /*! \brief Get quorum.
   */
  int getQuorum() const {
    return quorum_;
  }

  //@}

  //! @name Reset methods
  //@{ 
  //! Informs the status test that it should reset its internal configuration to the uninitialized state.
  /*! This is necessary for the case when the status test is being reused by another solver or for another
    eigenvalue problem. The status test may have information that pertains to a particular problem or solver 
    state. The internal information will be reset back to the uninitialized state. The user specified information 
    that the convergence test uses will remain.
  */
  void reset() { 
    ind_.resize(0);
    state_ = Undefined;
    test_->reset();
  }

  //! Clears the results of the last status test.
  /*! This should be distinguished from the reset() method, as it only clears the cached result from the last 
   * status test, so that a call to getStatus() will return ::Undefined. This is necessary for the SEQOR and SEQAND
   * tests in the StatusTestCombo class, which may short circuit and not evaluate all of the StatusTests contained
   * in them.
  */
  void clearStatus() {
    ind_.resize(0);
    state_ = Undefined;
    test_->clearStatus();
  }

  /*! \brief Set the auxiliary eigenvalues.
   *
   *  This routine sets only the real part of the auxiliary eigenvalues; the imaginary part is set to zero. This routine also resets the state to ::Undefined.
   */
  void setAuxVals(const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &vals) {
    rvals_ = vals;
    ivals_.resize(rvals_.size(),MT::zero());
    state_ = Undefined;
  }

  /*! \brief Set the auxiliary eigenvalues.
   *
   *  This routine sets both the real and imaginary parts of the auxiliary eigenvalues. This routine also resets the state to ::Undefined.
   */
  void setAuxVals(const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &rvals, const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &ivals) {
    rvals_ = rvals;
    ivals_ = ivals;
    state_ = Undefined;
  }

  /*! \brief Get the auxiliary eigenvalues.
   *
   *  This routine gets the real and imaginary parts of the auxiliary eigenvalues.
   */
  void getAuxVals(std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &rvals, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &ivals) const {
    rvals = rvals_;
    ivals = ivals_;
  }

  //@}

  //! @name Print methods
  //@{ 
  
  //! Output formatted description of stopping test to output stream.
  std::ostream& print(std::ostream& os, int indent = 0) const;
 
  //@}
  private:
    TestStatus state_;
    std::vector<int> ind_;
    int quorum_;
    std::vector<MagnitudeType> rvals_, ivals_;
    Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > sorter_;
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test_;
};


template <class ScalarType, class MV, class OP>
StatusTestWithOrdering<ScalarType,MV,OP>::StatusTestWithOrdering(Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test, Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > sorter, int quorum)
  : state_(Undefined), ind_(0), quorum_(quorum), rvals_(0), ivals_(0), sorter_(sorter), test_(test)
{
  TEUCHOS_TEST_FOR_EXCEPTION(sorter_ == Teuchos::null, StatusTestError, "StatusTestWithOrdering::constructor() was passed null pointer for constituent SortManager.");
  TEUCHOS_TEST_FOR_EXCEPTION(test_ == Teuchos::null, StatusTestError, "StatusTestWithOrdering::constructor() was passed null pointer for constituent StatusTest.");
}

template <class ScalarType, class MV, class OP>
TestStatus StatusTestWithOrdering<ScalarType,MV,OP>::checkStatus( Eigensolver<ScalarType,MV,OP>* solver ) {

  // Algorithm
  // we PASS iff the "most signficant" values of "all values" PASS
  // "most significant" is measured by sorter
  // auxilliary values are automatically PASSED
  // constituent status test determines who else passes
  // "all values" mean {auxilliary values} UNION {solver's ritz values}
  //
  // test_->checkStatus()                   calls the constituent status test
  // cwhch = test_->whichVecs()             gets the passing indices from the constituent test
  // solval = solver->getRitzValues()       gets the solver's ritz values
  // allval = {solval auxval}               contains all values
  // sorter_->sort(allval,perm)             sort all values (we just need the perm vector)
  //
  // allpass = {cwhch numsolval+1 ... numsolval+numAux}
  // mostsig = {perm[1] ... perm[quorum]}
  // whichVecs = mostsig INTERSECT allpass
  // if size(whichVecs) >= quorum,
  //    PASS
  // else
  //    FAIL
  // 
  // finish: this needs to be better tested and revisited for efficiency.

  // call the constituent solver manager
  test_->checkStatus(solver);
  std::vector<int> cwhch( test_->whichVecs() );

  // get the ritzvalues from solver
  std::vector<Value<ScalarType> > solval = solver->getRitzValues();
  int numsolval = solval.size();
  std::cout << "Number of approximate Ritz values from solver: " << numsolval << std::endl;
  int numauxval = rvals_.size();
  std::cout << "Number of auxiliary values from solver: " << numauxval << std::endl;
  int numallval = numsolval + numauxval;
  std::cout << "Number of all Ritz values (solver + aux): " << numallval << std::endl;

  if (numallval == 0) {
    ind_.resize(0);
    return Failed;
  }

  // make space for all values, real and imaginary parts
  std::vector<MagnitudeType> allvalr(numallval), allvali(numallval);
  // separate the real and imaginary parts of solver ritz values
  for (int i=0; i<numsolval; ++i) {
    allvalr[i] = solval[i].realpart;
    allvali[i] = solval[i].imagpart;
  }
  // put the auxiliary values at the end of the solver ritz values
  std::copy(rvals_.begin(),rvals_.end(),allvalr.begin()+numsolval);
  std::copy(ivals_.begin(),ivals_.end(),allvali.begin()+numsolval);

  // sort all values
  std::vector<int> perm(numallval);
  sorter_->sort(allvalr,allvali,Teuchos::rcpFromRef(perm),numallval);

  std::cout << "After sorting allval: " << std::endl;
  for (int i=0; i<numallval; i++)
  {
    std::cout << allvalr[i] << " + i " << allvali[i] << std::endl;
  } 

  // make the set of passing values: allpass = {cwhch -1 ... -numauxval}
  std::vector<int> allpass(cwhch.size() + numauxval);
  std::copy(cwhch.begin(),cwhch.end(),allpass.begin());
  for (int i=0; i<numauxval; i++) {
    allpass[cwhch.size()+i] = -(i+1);
  }

  // make list, with length quorum, of most significant values, if there are that many
  int numsig = quorum_ < numallval ? quorum_ : numallval;
  // int numsig = cwhch.size() + numauxval;
  std::cout << "Number of most significant values (quorum_ < numallvall ? quorum_ : numallvall : " << numsig << std::endl;
  std::vector<int> mostsig(numsig);
  for (int i=0; i<numsig; ++i) {
    mostsig[i] = perm[i];
    // if perm[i] >= numsolval, then it corresponds to the perm[i]-numsolval aux val
    // the first aux val gets index -numauxval, the second -numauxval+1, and so forth
    if (mostsig[i] >= numsolval) {
      mostsig[i] = mostsig[i]-numsolval-numauxval;
    }
  }

  // who passed?
  // to use set_intersection, ind_ must have room for the result, which will have size() <= min( allpass.size(), mostsig.size() )
  // also, allpass and mostsig must be in ascending order; neither are
  // lastly, the return from set_intersection points to the last element in the intersection, which tells us how big the intersection is
  ind_.resize(numsig);
  std::vector<int>::iterator end;
  std::sort(mostsig.begin(),mostsig.end());
  std::sort(allpass.begin(),allpass.end());
  //for (int i=0; i<(int)cwhch.size(); ++i) 
  //{
  //  std::cout << "allpass[" << i << "] = " << allpass[i] << std::endl;
  //  std::cout << "mostsig[" << i << "] = " << mostsig[i] << std::endl;
  //}
  end = std::set_intersection(mostsig.begin(),mostsig.end(),allpass.begin(),allpass.end(),ind_.begin());
  std::cout << "Number of most significant values that have passed: " << (int)(end-ind_.begin()) << std::endl;
  ind_.resize(end - ind_.begin());
  for (int i=0; i<(int)ind_.size(); ++i) 
    std::cout << "ind_[" << i << "] = " << ind_[i] << std::endl;

  // did we pass, overall
  if (ind_.size() >= (unsigned int)quorum_) {
    state_ = Passed;
  }
  else {
    state_ = Failed;
  }
  return state_;
}


template <class ScalarType, class MV, class OP>
std::ostream& StatusTestWithOrdering<ScalarType,MV,OP>::print(std::ostream& os, int indent) const {
  // build indent string
  std::string ind(indent,' ');
  // print header
  os << ind << "- StatusTestWithOrdering: ";
  switch (state_) {
  case Passed:
    os << "Passed" << std::endl;
    break;
  case Failed:
    os << "Failed" << std::endl;
    break;
  case Undefined:
    os << "Undefined" << std::endl;
    break;
  }
  // print parameters, namely, quorum
  os << ind << "  Quorum: " << quorum_ << std::endl;
  // print any auxilliary values
  os << ind << "  Auxiliary values: ";
  if (rvals_.size() > 0) {
    for (unsigned int i=0; i<rvals_.size(); i++) {
      os << "(" << rvals_[i] << ", " << ivals_[i] << ")  ";
    }
    os << std::endl;
  }
  else {
    os << "[empty]" << std::endl;
  }
  // print which vectors have passed
  if (state_ != Undefined) {
    os << ind << "  Which vectors: ";
    if (ind_.size() > 0) {
      for (unsigned int i=0; i<ind_.size(); i++) os << ind_[i] << " ";
      os << std::endl;
    }
    else {
      os << "[empty]" << std::endl;
    }
  }
  // print constituent test
  test_->print(os,indent+2);
  return os;
}


} // end of Anasazi namespace

#endif /* ANASAZI_STATUS_TEST_ORDEREDRESNORM_HPP */
