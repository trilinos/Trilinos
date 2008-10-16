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

#ifndef BELOS_TPETRA_ADAPTER_HPP
#define BELOS_TPETRA_ADAPTER_HPP

/*! \file BelosTpetraAdapter.hpp
    \brief Provides several interfaces between Belos virtual classes and Tpetra concrete classes.
*/

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"

namespace Belos {
 
  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::MultiVecTraits for Tpetra::MultiVector.
  //
  ////////////////////////////////////////////////////////////////////

  template<class Scalar>
  class MultiVecTraits<Scalar, Tpetra::MultiVector<int,Scalar> >
  {
  public:

    static Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > Clone( const Tpetra::MultiVector<int,Scalar>& mv, const int numvecs )
    { TEST_FOR_EXCEPT(true); }

    static Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > CloneCopy( const Tpetra::MultiVector<int,Scalar>& mv )
    { TEST_FOR_EXCEPT(true); }

    static Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > CloneCopy( const Tpetra::MultiVector<int,Scalar>& mv, const std::vector<int>& index )
    { TEST_FOR_EXCEPT(true); }

    static Teuchos::RCP<Tpetra::MultiVector<int,Scalar> > CloneView( Tpetra::MultiVector<int,Scalar>& mv, const std::vector<int>& index )
    {  TEST_FOR_EXCEPT(true); }

    static Teuchos::RCP<const Tpetra::MultiVector<int,Scalar> > CloneView( const Tpetra::MultiVector<int,Scalar>& mv, const std::vector<int>& index )
    { TEST_FOR_EXCEPT(true); }

    static int GetVecLength( const Tpetra::MultiVector<int,Scalar>& mv )
    { TEST_FOR_EXCEPT(true); }

    static int GetNumberVecs( const Tpetra::MultiVector<int,Scalar>& mv )
    { TEST_FOR_EXCEPT(true); }

    static void MvTimesMatAddMv( const Scalar alpha, const Tpetra::MultiVector<int,Scalar>& A, 
                                 const Teuchos::SerialDenseMatrix<int,Scalar>& B, 
                                 const Scalar beta, Tpetra::MultiVector<int,Scalar>& mv )
    { TEST_FOR_EXCEPT(true); }

    static void MvAddMv( const Scalar alpha, const Tpetra::MultiVector<int,Scalar>& A, const Scalar beta, const Tpetra::MultiVector<int,Scalar>& B, Tpetra::MultiVector<int,Scalar>& mv )
    { TEST_FOR_EXCEPT(true); }

    static void MvScale ( Tpetra::MultiVector<int,Scalar>& mv, Scalar alpha )
    { TEST_FOR_EXCEPT(true); }

    static void MvScale ( Tpetra::MultiVector<int,Scalar>& mv, const std::vector<Scalar>& alpha )
    { TEST_FOR_EXCEPT(true); }

    static void MvTransMv( const Scalar alpha, const Tpetra::MultiVector<int,Scalar>& A, const Tpetra::MultiVector<int,Scalar>& mv, Teuchos::SerialDenseMatrix<int,Scalar>& B )
    { TEST_FOR_EXCEPT(true); }

    static void MvDot( const Tpetra::MultiVector<int,Scalar>& mv, const Tpetra::MultiVector<int,Scalar>& A, std::vector<Scalar>& b )
    { TEST_FOR_EXCEPT(true); }

    static void MvNorm( const Tpetra::MultiVector<int,Scalar>& mv, std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec, NormType type=TwoNorm )
    { TEST_FOR_EXCEPT(true); }

    static void SetBlock( const Tpetra::MultiVector<int,Scalar>& A, const std::vector<int>& index, Tpetra::MultiVector<int,Scalar>& mv )
    { TEST_FOR_EXCEPT(true); }

    static void MvRandom( Tpetra::MultiVector<int,Scalar>& mv )
    { TEST_FOR_EXCEPT(true); }

    static void MvInit( Tpetra::MultiVector<int,Scalar>& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    { TEST_FOR_EXCEPT(true); }

    static void MvPrint( const Tpetra::MultiVector<int,Scalar>& mv, std::ostream& os )
    { TEST_FOR_EXCEPT(true); }
  };        

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::OperatorTraits for Tpetra::Operator.
  //
  ////////////////////////////////////////////////////////////////////

  template <class Scalar> 
  class OperatorTraits < Scalar, Tpetra::MultiVector<int,Scalar>, Tpetra::Operator<int,Scalar> >
  {
  public:
    static void Apply ( const Tpetra::Operator<int,Scalar>& Op, 
                        const Tpetra::MultiVector<int,Scalar>& X,
                              Tpetra::MultiVector<int,Scalar>& Y,
                        ETrans trans=NOTRANS )
    { TEST_FOR_EXCEPT(true); }
  };

} // end of Belos namespace 

#endif 
// end of file BELOS_TPETRA_ADAPTER_HPP
