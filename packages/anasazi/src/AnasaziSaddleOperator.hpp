// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file AnasaziSaddleOperator.hpp
 *  \brief An operator of the form [A Y; Y' 0] where A is a sparse matrix and Y a multivector.
 *
 *  Used by TraceMin to solve the saddle point problem.
*/

#ifndef ANASAZI_SADDLE_OPERATOR_HPP
#define ANASAZI_SADDLE_OPERATOR_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziSaddleContainer.hpp"
#include "AnasaziTraceMinRitzOp.hpp"

#include "Teuchos_SerialDenseSolver.hpp"

using Teuchos::RCP;

enum PrecType {NO_PREC, NONSYM, BD_PREC, HSS_PREC};

namespace Anasazi {
namespace Experimental {

template <class ScalarType, class MV, class OP>
class SaddleOperator : public TraceMinOp<ScalarType,SaddleContainer<ScalarType,MV>,OP>
{
  typedef Anasazi::MultiVecTraits<ScalarType,MV>       MVT;
  typedef Teuchos::SerialDenseMatrix<int,ScalarType>   SerialDenseMatrix;

public:
  // Default constructor
  SaddleOperator( ) { };
  SaddleOperator( const Teuchos::RCP<OP> A, const Teuchos::RCP<const MV> B, PrecType pt=NO_PREC, const ScalarType alpha=1. );

  // Applies the saddle point operator to a "multivector"
  void Apply(const SaddleContainer<ScalarType,MV>& X, SaddleContainer<ScalarType,MV>& Y) const;

  void removeIndices(const std::vector<int>& indicesToRemove) { A_->removeIndices(indicesToRemove); };

private:
  // A is the 1-1 block, and B the 1-2 block
  Teuchos::RCP<OP> A_;
  Teuchos::RCP<const MV> B_;
  Teuchos::RCP<SerialDenseMatrix> Schur_;
  PrecType pt_;
  ScalarType alpha_;
};



// Default constructor
template <class ScalarType, class MV, class OP>
SaddleOperator<ScalarType, MV, OP>::SaddleOperator( const Teuchos::RCP<OP> A, const Teuchos::RCP<const MV> B, PrecType pt, const ScalarType alpha )
{
  // Get a pointer to A and B
  A_ = A;
  B_ = B;
  pt_ = pt;
  alpha_ = alpha;

  if(pt == BD_PREC)
  {
    // Form the Schur complement
    int nvecs = MVT::GetNumberVecs(*B);
    Teuchos::RCP<MV> AinvB = MVT::Clone(*B,nvecs);
    Schur_ = rcp(new SerialDenseMatrix(nvecs,nvecs));

    A_->Apply(*B_,*AinvB);

    MVT::MvTransMv(1., *B_, *AinvB, *Schur_);
  }
}

// Applies the saddle point operator to a "multivector"
template <class ScalarType, class MV, class OP>
void SaddleOperator<ScalarType, MV, OP>::Apply(const SaddleContainer<ScalarType,MV>& X, SaddleContainer<ScalarType,MV>& Y) const
{
  RCP<SerialDenseMatrix> Xlower = X.getLower();
  RCP<SerialDenseMatrix> Ylower = Y.getLower();

  if(pt_ == NO_PREC)
  {
    // trans does literally nothing, because the operator is symmetric
    // Y.bottom = B'X.top
    MVT::MvTransMv(1., *B_, *(X.upper_), *Ylower);
  
    // Y.top = A*X.top+B*X.bottom
    A_->Apply(*(X.upper_), *(Y.upper_));
    MVT::MvTimesMatAddMv(1., *B_, *Xlower, 1., *(Y.upper_));
  }
  else if(pt_ == NONSYM)
  {
    // Y.bottom = -B'X.top
	MVT::MvTransMv(-1., *B_, *(X.upper_), *Ylower);
  
    // Y.top = A*X.top+B*X.bottom
    A_->Apply(*(X.upper_), *(Y.upper_));
    MVT::MvTimesMatAddMv(1., *B_, *Xlower, 1., *(Y.upper_));
  }
  else if(pt_ == BD_PREC)
  {
    Teuchos::SerialDenseSolver<int,ScalarType> MySolver;

    // Solve A Y.X = X.X
    A_->Apply(*(X.upper_),*(Y.upper_));

    // So, let me tell you a funny story about how the SerialDenseSolver destroys the original matrix...
    Teuchos::RCP<SerialDenseMatrix> localSchur = Teuchos::rcp(new SerialDenseMatrix(*Schur_));

    // Solve the small system
    MySolver.setMatrix(localSchur);
    MySolver.setVectors(Ylower, Xlower);
    MySolver.solve();
  }
  // Hermitian-Skew Hermitian splitting has some extra requirements
  // We need B'B = I, which is true for standard eigenvalue problems, but not generalized
  // We also need to use gmres, because our operator is no longer symmetric
  else if(pt_ == HSS_PREC)
  {
//    std::cout << "applying preconditioner to";
//    X.MvPrint(std::cout);

    // Solve (H + alpha I) Y1 = X
    // 1.  Apply preconditioner
    A_->Apply(*(X.upper_),*(Y.upper_));
    // 2. Scale by 1/alpha
    *Ylower = *Xlower; 
    Ylower->scale(1./alpha_);
	
//    std::cout << "H preconditioning produced";
//	Y.setLower(Ylower);
//    Y.MvPrint(std::cout);

    // Solve (S + alpha I) Y = Y1
    // 1.  Y_lower = (B' Y1_upper + alpha Y1_lower) / (1 + alpha^2)
    Teuchos::RCP<SerialDenseMatrix> Y1_lower = Teuchos::rcp(new SerialDenseMatrix(*Ylower));
    MVT::MvTransMv(1,*B_,*(Y.upper_),*Ylower);
//	std::cout << "Y'b1 " << *Ylower;
    Y1_lower->scale(alpha_);
//	std::cout << "alpha b2 " << *Y1_lower;
    *Ylower += *Y1_lower;
//	std::cout << "alpha b2 + Y'b1 " << *Ylower;
    Ylower->scale(1/(1+alpha_*alpha_));
    // 2.  Y_upper = (Y1_upper - B Y_lower) / alpha
    MVT::MvTimesMatAddMv(-1/alpha_,*B_,*Ylower,1/alpha_,*(Y.upper_));

//    std::cout << "preconditioning produced";
//	Y.setLower(Ylower);
//    Y.MvPrint(std::cout);
  }
  else
  {
    std::cout << "Not a valid preconditioner type\n";
  }

  Y.setLower(Ylower);
  
//  std::cout << "result of applying operator";
//  Y.MvPrint(std::cout);
}

} // End namespace Experimental

template<class ScalarType, class MV, class OP>
class OperatorTraits<ScalarType, Experimental::SaddleContainer<ScalarType,MV>, Experimental::SaddleOperator<ScalarType,MV,OP> >
{
public:
  static void Apply( const Experimental::SaddleOperator<ScalarType,MV,OP>& Op, 
                     const Experimental::SaddleContainer<ScalarType,MV>& x, 
                     Experimental::SaddleContainer<ScalarType,MV>& y)
    { Op.Apply( x, y); };
};

} // end namespace Anasazi

#ifdef HAVE_ANASAZI_BELOS
namespace Belos {

template<class ScalarType, class MV, class OP>
class OperatorTraits<ScalarType, Anasazi::Experimental::SaddleContainer<ScalarType,MV>, Anasazi::Experimental::SaddleOperator<ScalarType,MV,OP> >
{
public:
  static void Apply( const Anasazi::Experimental::SaddleOperator<ScalarType,MV,OP>& Op, 
                     const Anasazi::Experimental::SaddleContainer<ScalarType,MV>& x, 
                     Anasazi::Experimental::SaddleContainer<ScalarType,MV>& y)
    { Op.Apply( x, y); };
};

} // end namespace Belos
#endif

#endif
