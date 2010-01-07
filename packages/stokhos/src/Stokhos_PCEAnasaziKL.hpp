// $Id: Stokhos_MatrixFreeEpetraOp.hpp,v 1.7 2009/09/14 18:35:48 etphipp Exp $ 
// $Source: /space/CVS/Trilinos/packages/stokhos/src/Stokhos_MatrixFreeEpetraOp.hpp,v $ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
