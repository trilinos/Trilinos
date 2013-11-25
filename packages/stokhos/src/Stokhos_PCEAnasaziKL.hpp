// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_PCE_ANASAZI_KL_HPP
#define STOKHOS_PCE_ANASAZI_KL_HPP

#include "Stokhos_ConfigDefs.h"
#ifdef HAVE_STOKHOS_ANASAZI

#include "Teuchos_ParameterList.hpp"
#include "Stokhos_PCECovarianceOp.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziEpetraAdapter.hpp"

namespace Stokhos {

  //! Compute a KL expansion from a PCE expansion using Anasazi
  /*!
   * This class computes a KL expansion from a given PCE expansion by solving
   * the eigenvalue problem C*z = lambda*z using Anasazi, where C is the 
   * covariance operator of the expansion.
   */
  class PCEAnasaziKL {
  public:

    //! Constructor
    PCEAnasaziKL(const Stokhos::VectorOrthogPoly<Epetra_Vector>& X_poly,
		 int num_KL_);

    //! Constructor with block-vector X
    PCEAnasaziKL(const Teuchos::RCP<const EpetraExt::BlockVector>& X,
		 const Stokhos::OrthogPolyBasis<int,double>& basis,
		 int num_KL_);

    //! Constructor with multi-vector X
    PCEAnasaziKL(const Teuchos::RCP<const Epetra_MultiVector>& X,
		 const Stokhos::OrthogPolyBasis<int,double>& basis,
		 int num_KL_);
    
    //! Destructor
    virtual ~PCEAnasaziKL() {}

    //! Get default parameters
    Teuchos::ParameterList getDefaultParams() const;

    //! Compute KL expansion
    bool computeKL(Teuchos::ParameterList& anasazi_params);

    //! Get KL eigenvalues
    Teuchos::Array<double> getEigenvalues() const;

    //! Get KL eigenvectors
    Teuchos::RCP<Epetra_MultiVector> getEigenvectors() const;

  private:

    //! Prohibit copying
    PCEAnasaziKL(const PCEAnasaziKL&);

    //! Prohibit copying
    PCEAnasaziKL& operator=(const PCEAnasaziKL&);

  protected:

    typedef double ScalarType;
    typedef Teuchos::ScalarTraits<ScalarType>          SCT;
    typedef SCT::magnitudeType               MagnitudeType;
    typedef Epetra_MultiVector                          MV;
    typedef Epetra_Operator                             OP;
    typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
    typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;

    //! Covariance operator
    Teuchos::RCP<Stokhos::PCECovarianceOp> covOp;

    //! Number of KL terms
    int num_KL;

    //! Eigen problem
    Teuchos::RCP<Anasazi::BasicEigenproblem<ScalarType,MV,OP> > anasazi_problem;

    // Eigen solution
    Anasazi::Eigensolution<ScalarType,MV> sol;

  }; // PCEAnasaziKL

} // namespace Stokhos

#endif // HAVE_STOKHOS_ANASAZI

#endif // STOKHOS_PCE_ANASAZI_KL_HPP
