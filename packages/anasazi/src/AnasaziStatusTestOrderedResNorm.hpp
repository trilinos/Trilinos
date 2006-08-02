
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

#ifndef ANASAZI_STATUS_TEST_ORDEREDRESNORM_HPP
#define ANASAZI_STATUS_TEST_ORDEREDRESNORM_HPP

/*!
  \file AnasaziStatusTestOrderedResNorm.hpp
  \brief A status test for testing the norm of the eigenvectors residuals along with a 
         set of auxiliary eigenvalues.
*/


#include "AnasaziStatusTest.hpp"
#include "Teuchos_ScalarTraits.hpp"

  /*! 
    \class Anasazi::StatusTestOrderedResNorm 
    
    \brief A status test for testing the norm of the eigenvectors residuals
    along with a set of auxiliary eigenvalues. 
    
    The test evaluates to ::Passed when then the most significant of the
    eigenvalues all have a residual below a certain threshhold. The purpose of
    the test is to not only test convergence for some number of eigenvalues,
    but to test convergence for the correct ones.
    
    In addition to specifying the tolerance, the user may specify:
    <ul>
      <li> the norm to be used: 2-norm or OrthoManager::norm()
      <li> the scale: absolute or relative to magnitude of Ritz value 
      <li> the quorum: the number of vectors required for the test to 
           evaluate as ::Passed.
    </ul>

    Finally, the user must specify the Anasazi::SortManager used for deciding
    significance. 
  */

namespace Anasazi {


template <class ScalarType, class MV, class OP>
class StatusTestOrderedResNorm : public StatusTest<ScalarType,MV,OP> {

 private:
  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
  typedef Teuchos::ScalarTraits<MagnitudeType> MT;
  typedef Teuchos::ScalarTraits<ScalarType>    SCT;

 public:
  //! @name Constructors/destructors
  //@{ 

  //! Constructor
  StatusTestOrderedResNorm(Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > sorter, MagnitudeType tol, int quorum = -1, bool use2Norm = false, bool scaled = true);

  //! Destructor
  virtual ~StatusTestOrderedResNorm() {};
  //@}

  //! @name Status methods
  //@{ 
  /*! Check status as defined by test.
    \return TestStatus indicating whether the test passed or failed.
  */
  TestStatus checkStatus( Eigensolver<ScalarType,MV,OP>* solver );

  //! Return the result of the most recent checkStatus call.
  TestStatus getStatus() const { return _state; }
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
    _state = Undefined;
  }

  //! Clears the results of the last status test.
  /*! This should be distinguished from the reset() method, as it only clears the cached result from the last 
   * status test, so that a call to getStatus() will return ::Undefined. This is necessary for the SEQOR and SEQAND
   * tests in the StatusTestCombo class, which may short circuit and not evaluate all of the StatusTests contained
   * in them.
  */
  void clearStatus() {
    _state = Undefined;
  }

  //! Set the auxiliary eigenvalues.
  /*! This routine also resets the state to ::Undefined.
   */
  void setAuxVals(const std::vector<MagnitudeType> &vals) {
    _vals = vals;
    _state = Undefined;
  }

  //! Get the auxiliary eigenvalues and their residuals.
  void getAuxVals(std::vector<MagnitudeType> &vals) const {
    vals = _vals;
  }

  //@}

  //! @name Accessor methods
  //@{ 

  /*! \brief Set tolerance.
   *  This also resets the test status to ::Undefined.
   */
  void setTolerance(MagnitudeType tol) {
    TEST_FOR_EXCEPTION(tol <= MT::zero(), StatusTestError, "StatusTestOrderedResNorm: test tolerance must be strictly positive.");
    _state = Undefined;
    _tol = tol;
  }

  //! Get tolerance.
  MagnitudeType getTolerance() {return _tol;}

  /*! \brief Set the norm. If true, the test uses the 2-norm. Otherwise, it uses the norm defined by the OrthoManager.
   *  This also resets the test status to ::Undefined.
   */
  void setNorm(bool use2Norm) {
    _state = Undefined;
    _use2Norm = use2Norm;
  }

  //! Returns true if the test is using the 2-norm. Returns false if using the norm defined by the OrthoManager.
  bool uses2Norm() {return _use2Norm;}

  /*! \brief Instruct test to scale norms by eigenvalue estimates (relative scale).
   *  This also resets the test status to ::Undefined.
   */
  void setScale(bool relscale) {
    _state = Undefined;
    _scaled = relscale;
  }

  //! Returns true if the test scales the norms by the eigenvalue estimates (relative scale).
  bool getScale() {return _scaled;}

  //! Get the indices for the vectors that passed the test.
  std::vector<int> whichVecs() {
    return _ind;
  }

  //! Get the number of vectors that passed the test.
  int howMany() {
    return _ind.size();
  }

  //@}

  //! @name Print methods
  //@{ 
  
  //! Output formatted description of stopping test to output stream.
  ostream& print(ostream& os, int indent = 0) const;
 
  //@}
  private:
    TestStatus _state;
    MagnitudeType _tol;
    std::vector<int> _ind;
    int _quorum;
    bool _use2Norm, _scaled;
    std::vector<MagnitudeType> _vals;
    Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > _sorter;
};


template <class ScalarType, class MV, class OP>
StatusTestOrderedResNorm<ScalarType,MV,OP>::StatusTestOrderedResNorm(Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > sorter, MagnitudeType tol, int quorum, bool use2Norm, bool scaled)
  : _state(Undefined), _quorum(quorum), _use2Norm(use2Norm), _scaled(scaled), _sorter(sorter) 
{
  TEST_FOR_EXCEPTION(_sorter == Teuchos::null, StatusTestError, "StatusTestOrderedResNorm::constructor() was passed null pointer for SortManager.");
  setTolerance(tol); 
}

template <class ScalarType, class MV, class OP>
TestStatus StatusTestOrderedResNorm<ScalarType,MV,OP>::checkStatus( Eigensolver<ScalarType,MV,OP>* solver ) {

  int bs = solver->getBlockSize();
  int numaux = _vals.size();
  int num = bs + numaux;
  std::vector<MagnitudeType> res; 
  std::vector<MagnitudeType> evals;

  // get the eigenvalues
  evals = solver->getRitzValues();
  evals.resize(solver->getBlockSize());
  // put the auxiliary values in the vectors as well
  evals.insert(evals.end(),_vals.begin(),_vals.end());

  // get the eigenvector residuals norms (using the appropriate norm)
  if (_use2Norm) {
    res = solver->getRes2Norms();
  }
  else {
    res = solver->getResNorms();
  }
  // if appropriate, scale the norms by the magnitude of the eigenvalue estimate
  if (_scaled) {
    for (int i=0; i<bs; i++) {
      if ( SCT::magnitude(evals[i]) != SCT::zero() ) {
        res[i] /= SCT::magnitude(evals[i]);
      }
    }
  }
  // add -1 residuals for the auxiliary values (-1 < _tol)
  res.insert(res.end(),numaux,-MT::one());

  // sort the eigenvalues; SortManager takes a vector of ScalarType, so promote our eigenvalues
  // we don't actually need the sorted eigenvalues; just the permutation vector
  std::vector<ScalarType> evals_st(evals.size());
  copy(evals.begin(),evals.end(),evals_st.begin());
  std::vector<int> perm(num,-1);
  _sorter->sort(solver,num,&evals_st[0],&perm);
  // we don't need to the eigenvalues anymore, so don't bother moving the sorted values back to the MagnitudeType vector

  // sort the residuals
  std::vector<MagnitudeType> oldres = res;
  // apply the sorting to the residuals and original indices
  for (int i=0; i<num; i++) {
    res[i] = oldres[perm[i]];
  }

  // indices: [0,bs) are from solver, [bs,bs+numaux) are from _vals
  _ind.resize(num);

  // test the norms: we want res [0,quorum) to be <= tol
  int have = 0;
  int need = (_quorum == -1) ? num : _quorum;
  for (int i=0; i<need; i++) {
    TEST_FOR_EXCEPTION( SCT::isnaninf(res[i]), StatusTestError, "StatusTestOrderedResNorm::checkStatus(): residual norm is nan or inf" );
    if (res[i] < _tol) {
      _ind[have] = perm[i];
      have++;
    }
  }
  _ind.resize(have);
  _state = (have >= need) ? Passed : Failed;
  return _state;
}


template <class ScalarType, class MV, class OP>
ostream& StatusTestOrderedResNorm<ScalarType,MV,OP>::print(ostream& os, int indent) const {
  string ind(indent,' ');
  os << ind << "- StatusTestOrderedResNorm: ";
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
  os << ind << "(Tolerance,Use2Norm,Scaled,Quorum): " 
            << "(" << _tol 
            << "," << (_use2Norm ? "true" : "false") 
            << "," << (_scaled   ? "true" : "false")
            << "," << _quorum 
            << ")" << endl;
  os << ind << "Auxiliary values: ";
  if (_vals.size() > 0) {
    for (unsigned int i=0; i<_vals.size(); i++) {
      os << _vals[i] << ", ";
    }
    os << endl;
  }
  else {
    os << "[empty]" << endl;
  }

  if (_state != Undefined) {
    os << ind << "Which vectors: ";
    if (_ind.size() > 0) {
      for (unsigned int i=0; i<_ind.size(); i++) os << _ind[i] << " ";
      os << endl;
    }
    else {
      os << "[empty]" << endl;
    }
  }
  return os;
}


} // end of Anasazi namespace

#endif /* ANASAZI_STATUS_TEST_ORDEREDRESNORM_HPP */
