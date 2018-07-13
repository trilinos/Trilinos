// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ANASAZI_EIGENSOLVER_HPP
#define ANASAZI_EIGENSOLVER_HPP

/*! \file AnasaziEigensolver.hpp
    \brief Pure virtual base class which describes the basic interface to the iterative eigensolver.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEigensolverDecl.hpp"
#include "AnasaziStatusTestDecl.hpp"

#include "AnasaziEigenproblem.hpp"
#include "AnasaziSortManager.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziOrthoManager.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"


namespace Anasazi {

template<class ScalarType, class MV, class OP>
class Eigensolver {

  public:

  //! @name Constructors/Destructor
  //@{ 

  //! Default Constructor.
  Eigensolver() {};

  //! Basic Constructor.
  /*! This constructor, implemented by all Anasazi eigensolvers, takes an Anasazi::Eigenproblem,
    Anasazi::SortManager, Anasazi::OutputManager, and Teuchos::ParameterList as input.  These
    four arguments are sufficient enough for constructing any Anasazi::Eigensolver object.
  */
  Eigensolver( const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem, 
               const Teuchos::RCP<SortManager<ScalarType> >        &sorter,
               const Teuchos::RCP<OutputManager<ScalarType> >      &printer,
               const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >   &tester,
               const Teuchos::RCP<OrthoManager<ScalarType,MV> >    &ortho,
               Teuchos::ParameterList &params );

  //! Destructor.
  virtual ~Eigensolver() {};
  //@}


  //! @name Solver methods
  //@{ 
  
  /*! \brief This method performs eigensolvers iterations until the status test
    indicates the need to stop or an error occurs (in which case, an exception is thrown).
  */
  virtual void iterate() = 0;

  /*! \brief Initialize the solver with the initial vectors from the eigenproblem
   *  or random data.
   */
  virtual void initialize() = 0;

  //@}

    
    //! @name Status methods
  //@{ 

  //! \brief Get the current iteration count.
  virtual int getNumIters() const = 0;

  //! \brief Reset the iteration count.
  virtual void resetNumIters() = 0;

  /*! \brief Get the Ritz vectors from the previous iteration. These are indexed using getRitzIndex().
   *
   * For a description of the indexing scheme, see getRitzIndex().
   */
  virtual Teuchos::RCP<const MV> getRitzVectors() = 0;

  //! \brief Get the Ritz values from the previous iteration.
  virtual std::vector<Value<ScalarType> > getRitzValues() = 0;

  /*! \brief Get the index used for indexing the compressed storage used for Ritz vectors for real, non-Hermitian problems. 
   *
   *  index has length numVecs, where each entry is 0, +1, or -1. These have the following interpretation:
   *     - index[i] == 0: signifies that the corresponding eigenvector is stored as the i column of Evecs. This will usually be the 
   *       case when ScalarType is complex, an eigenproblem is Hermitian, or a real, non-Hermitian eigenproblem has a real eigenvector.
   *     - index[i] == +1: signifies that the corresponding eigenvector is stored in two vectors: the real part in the i column of Evecs and the <i><b>positive</b></i> imaginary part in the i+1 column of Evecs.
   *     - index[i] == -1: signifies that the corresponding eigenvector is stored in two vectors: the real part in the i-1 column of Evecs and the <i><b>negative</b></i> imaginary part in the i column of Evecs
   */
  virtual std::vector<int> getRitzIndex() = 0;

  //! \brief Get the current residual norms
  /*! \return A vector of length blockSize containing the norms of the residuals, 
      according to the orthogonalization manager norm() method.
   */
  virtual std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getResNorms() = 0;

  //! Get the current residual 2-norms
  //! \return A vector of length blockSize containing the 2-norms of the residuals. 
  virtual std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRes2Norms() = 0;

  //! Get the 2-norms of the Ritz residuals.
  //! \return A vector of length blockSize containing the 2-norms of the Ritz residuals.
  virtual std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRitzRes2Norms() = 0;

  //! Get the dimension of the search subspace used to generate the current eigenvectors and eigenvalues.
  virtual int getCurSubspaceDim() const = 0;

  //! Get the maximum dimension allocated for the search subspace.
  virtual int getMaxSubspaceDim() const = 0;

  //@}


  
    //! @name Accessor methods
  //@{ 

  //! Set a new StatusTest for the solver.
  virtual void setStatusTest(Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test) = 0;

  //! Get the current StatusTest used by the solver.
  virtual Teuchos::RCP<StatusTest<ScalarType,MV,OP> > getStatusTest() const = 0;

  //! Get a constant reference to the eigenvalue problem.
  virtual const Eigenproblem<ScalarType,MV,OP>& getProblem() const = 0;

  //! Get the blocksize to be used by the iterative solver in solving this eigenproblem.
  virtual int getBlockSize() const = 0;
  
  //! \brief Set the blocksize to be used by the iterative solver in solving this eigenproblem.
  virtual void setBlockSize(int blockSize) = 0;

  //! Set the auxiliary vectors for the solver.
  virtual void setAuxVecs(const Teuchos::Array<Teuchos::RCP<const MV> > &auxvecs) = 0;

  //! Get the auxiliary vectors for the solver.
  virtual Teuchos::Array<Teuchos::RCP<const MV> > getAuxVecs() const = 0;

  //! States whether the solver has been initialized or not.
  virtual bool isInitialized() const = 0;

  //@}

    //! @name Output methods
  //@{ 

  //! This method requests that the solver print out its current status to screen.
  virtual void currentStatus(std::ostream &os) = 0;

  //@}
  
};

} // end Anasazi namespace

#endif /* ANASAZI_EIGENSOLVER_HPP */
