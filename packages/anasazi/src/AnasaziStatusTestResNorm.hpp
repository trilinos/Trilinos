
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
      <li> the norm to be used: 2-norm or OrthoManager::norm() or getRitzRes2Norms()
      <li> the scale: absolute or relative to magnitude of Ritz value 
      <li> the quorum: the number of vectors required for the test to 
           evaluate as ::Passed.
    </ul>
  */

namespace Anasazi {


template <class ScalarType, class MV, class OP>
class StatusTestResNorm : public StatusTest<ScalarType,MV,OP> {

  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

 public:

  //! @name Enums
  //@{

  /*! \enum ResType 
      \brief Enumerated type used to specify which residual norm used by this status test.
   */
  enum ResType {
    RES_ORTH,
    RES_2NORM,
    RITZRES_2NORM
  };

  //@}

  //! @name Constructors/destructors
  //@{ 

  //! Constructor
  StatusTestResNorm(MagnitudeType tol, int quorum = -1, ResType whichNorm = RES_ORTH, bool scaled = true);

  //! Destructor
  virtual ~StatusTestResNorm() {};
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

  //@}

  //! @name Accessor methods
  //@{ 

  /*! \brief Set quorum.
   *
   *  Setting quorum to -1 signifies that all residuals from the solver must meet the tolerance.
   *  This also resets the test status to ::Undefined.
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
   *  This also resets the test status to ::Undefined.
   */
  void setTolerance(MagnitudeType tol) {
    _state = Undefined;
    _tol = tol;
  }

  //! Get tolerance.
  MagnitudeType getTolerance() {return _tol;}

  /*! \brief Set the residual norm to be used by the status test.
   *
   *  This also resets the test status to ::Undefined.
   */
  void setWhichNorm(ResType whichNorm) {
    _state = Undefined;
    _whichNorm = whichNorm;
  }

  //! Return the residual norm used by the status test.
  ResType getWhichNorm() {return _whichNorm;}

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
    bool _scaled;
    ResType _whichNorm;
};


template <class ScalarType, class MV, class OP>
StatusTestResNorm<ScalarType,MV,OP>::StatusTestResNorm(MagnitudeType tol, int quorum, ResType whichNorm, bool scaled)
  : _state(Undefined), _tol(tol), _quorum(quorum), _scaled(scaled), _whichNorm(whichNorm) {}

template <class ScalarType, class MV, class OP>
TestStatus StatusTestResNorm<ScalarType,MV,OP>::checkStatus( Eigensolver<ScalarType,MV,OP>* solver ) {
  typedef Teuchos::ScalarTraits<MagnitudeType> MT;
  
  std::vector<MagnitudeType> res;

  // get the eigenvector/ritz residuals norms (using the appropriate norm)
  // get the eigenvalues/ritzvalues and ritz index as well
  std::vector<MagnitudeType> vals = solver->getRitzValues();
  std::vector<int> ind = solver->getRitzIndex();
  switch (_whichNorm) {
    case RES_2NORM:
      res = solver->getRes2Norms();
      // we want only the ritz values corresponding to our eigenvector residuals
      vals.resize(res.size());
      break;
    case RES_ORTH:
      res = solver->getResNorms();
      // we want only the ritz values corresponding to our eigenvector residuals
      vals.resize(res.size());
      break;
    case RITZRES_2NORM:
      res = solver->getRitzRes2Norms();
      break;
  }

  // if appropriate, scale the norms by the magnitude of the eigenvalue estimate
  if (_scaled) {
    Teuchos::LAPACK<int,MagnitudeType> lapack;

    for (unsigned int i=0; i<res.size(); i++) {
      MagnitudeType tmp;
      if ( ind[i] == 0 ) {
        // real ritz value
        tmp = MT::magnitude(vals[i]);
      }
      else if ( ind[i] == +1 ) {
        // positive part of complex ritz value
        tmp = lapack.LAPY2( vals[i], vals[i+1] );
      }
      else if ( ind[i] == -1 ) {
        // negative part of complex ritz value
        tmp = lapack.LAPY2( vals[i-1], vals[i] );
      }
      else {
        TEST_FOR_EXCEPTION( true, std::logic_error, "Anasazi::StatusTestOrderedResNorm::checkStatus(): invalid Ritz index returned from solver." );
      }
      // scale by the newly computed magnitude of the ritz values
      if ( tmp != MT::zero() ) {
        res[i] /= tmp;
      }
    }
  }

  // test the norms
  int have = 0;
  _ind.resize(res.size());
  for (unsigned int i=0; i<res.size(); i++) {
    TEST_FOR_EXCEPTION( MT::isnaninf(res[i]), StatusTestError, "StatusTestResNorm::checkStatus(): residual norm is nan or inf" );
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
  os << ind << "- StatusTestResNorm: ";
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
  os << ind << "(Tolerance,WhichNorm,Scaled,Quorum): " 
            << "(" << _tol;
  switch (_whichNorm) {
  case RES_ORTH:
    os << ",RES_ORTH";
    break;
  case RES_2NORM:
    os << ",RES_2NORM";
    break;
  case RITZRES_2NORM:
    os << ",RITZRES_2NORM";
    break;
  }
  os        << "," << (_scaled   ? "true" : "false")
            << "," << _quorum 
            << ")" << endl;

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

#endif /* ANASAZI_STATUS_TEST_RESNORM_HPP */
