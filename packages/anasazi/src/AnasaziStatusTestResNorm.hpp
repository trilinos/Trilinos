
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

#ifndef ANASAZI_STATUS_TEST_RESNORM_HPP
#define ANASAZI_STATUS_TEST_RESNORM_HPP

/*!
  \file AnasaziStatusTestResNorm.hpp
  \brief A status test for testing the norm of the eigenvectors residuals.
*/


#include "AnasaziStatusTest.hpp"
#include "Teuchos_ScalarTraits.hpp"

  /*! 
    \class Anasazi::StatusTestResNorm
    \brief A status test for testing the norm of the eigenvectors residuals.
    
    Anasazi::StatusTestResNorm was designed to be used as a test for
    convergence. The tester compares the norms of the residual vectors against
    a user specified tolerance. 
    
    In addition to specifying the tolerance, the user may specify:
    <ul>
      <li> the norm to be used: 2-norm or OrthoManager::norm()
      <li> the scale: absolute or relative to magnitude of Ritz value 
      <li> the quorum: the number of vectors required to trigger the test
    </ul>
  */

namespace Anasazi {


template <class ScalarType, class MV, class OP>
class StatusTestResNorm : public StatusTest<ScalarType,MV,OP> {

  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

 public:
  //@{ \name Constructors/destructors.

  //! Constructor
  StatusTestResNorm(MagnitudeType tol, int quorum = -1, bool use2Norm = false, bool scaled = true);

  //! Destructor
  virtual ~StatusTestResNorm() {};
  //@}

  //@{ \name Status methods
  /*! Check status as defined by test.
    \return TestStatus indicating whether the test passed or failed.
  */
  TestStatus checkStatus( Eigensolver<ScalarType,MV,OP>* solver );

  //! Return the result of the most recent checkStatus call.
  TestStatus getStatus() const { return _state; }
  //@}

  //@{ \name Reset methods
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
   * status test, so that a call to getStatus() will return Undefined. This is necessary for the SEQOR and SEQAND
   * tests in the StatusTestCombo class, which may short circuit and not evaluate all of the StatusTests contained
   * in them.
  */
  void clearStatus() {
    _state = Undefined;
  }

  //@}

  //@{ \name Accessor methods

  /*! \brief Set quorum.
   *
   *  Setting quorum to -1 signifies that all residuals from the solver must meet the tolerance.
   *  This also resets the test status to Undefined.
   */
  void setQuorum(int quorum) {
    _state = Undefined;
    _quorum = quorum;
  }

  /*! \brief Get quorum.
   */
  int getQuorum() {
    return _quorum;
  }

  /*! \brief Set tolerance.
   *  This also resets the test status to Undefined.
   */
  void setTolerance(MagnitudeType tol) {
    _state = Undefined;
    _tol = tol;
  }

  //! Get tolerance.
  MagnitudeType getTolerance() {return _tol;}

  /*! \brief Set the norm. If true, the test uses the 2-norm. Otherwise, it uses the norm defined by the OrthoManager.
   *  This also resets the test status to Undefined.
   */
  void setNorm(bool use2Norm) {
    _state = Undefined;
    _use2Norm = use2Norm;
  }

  //! Returns true if the test is using the 2-norm. Returns false if using the norm defined by the OrthoManager.
  bool uses2Norm() {return _use2Norm;}

  /*! \brief Instruct test to scale norms by eigenvalue estimates (relative scale).
   *  This also resets the test status to Undefined.
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

  //@{ \name Print methods
  
  //! Output formatted description of stopping test to output stream.
  ostream& print(ostream& os, int indent = 0) const;
 
  //@}
  private:
    TestStatus _state;
    MagnitudeType _tol;
    std::vector<int> _ind;
    int _quorum;
    bool _use2Norm, _scaled;
};


template <class ScalarType, class MV, class OP>
StatusTestResNorm<ScalarType,MV,OP>::StatusTestResNorm(MagnitudeType tol, int quorum, bool use2Norm, bool scaled)
  : _state(Undefined), _tol(tol), _quorum(quorum), _use2Norm(use2Norm), _scaled(scaled) {}

template <class ScalarType, class MV, class OP>
TestStatus StatusTestResNorm<ScalarType,MV,OP>::checkStatus( Eigensolver<ScalarType,MV,OP>* solver ) {
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  
  std::vector<MagnitudeType> res;

  // get the eigenvector residuals norms (using the appropriate norm)
  if (_use2Norm) {
    res = solver->getRes2Norms();
  }
  else {
    res = solver->getResNorms();
  }

  // if appropriate, scale the norms by the magnitude of the eigenvalue estimate
  if (_scaled) {
    std::vector<MagnitudeType> vals = solver->getEigenvalues();
    for (unsigned int i=0; i<res.size(); i++) {
      if ( vals[i] != SCT::zero() ) {
        res[i] /= vals[i];
      }
    }
  }

  // test the norms
  int have = 0;
  _ind.resize(res.size());
  for (unsigned int i=0; i<res.size(); i++) {
    TEST_FOR_EXCEPTION( SCT::isnaninf(res[i]), StatusTestError, "StatusTestResNorm::checkStatus(): residual norm is nan or inf" );
    if (res[i] < _tol) {
      _ind[have] = i;
      have++;
    }
  }
  _ind.resize(have);
  int need = (_quorum == -1) ? res.size() : _quorum;
  _state = (have >= need) ? Passed : Failed;
  return _state;
}


template <class ScalarType, class MV, class OP>
ostream& StatusTestResNorm<ScalarType,MV,OP>::print(ostream& os, int indent) const {
  string ind(indent,' ');
  if (_state == Undefined) {
    os << ind << "ResNorm <= " << _tol << " undefined." << endl;
  }
  else {
    os << ind << "ResNorm <= " << _tol << " for " << _ind.size() 
       << " vectors. Status " << (_state == Passed ? "Passed" : "Failed") << endl;
  }
  return os;
}


} // end of Anasazi namespace

#endif /* ANASAZI_STATUS_TEST_RESNORM_HPP */
