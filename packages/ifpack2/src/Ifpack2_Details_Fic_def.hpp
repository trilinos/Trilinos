// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// @file Ifpack2_fic_def.hpp

#ifndef __IFPACK2_FIC_DEF_HPP__ 
#define __IFPACK2_FIC_DEF_HPP__ 

#include "Ifpack2_Details_Fic_decl.hpp"
#include "Ifpack2_Details_CrsArrays.hpp"
#include <Kokkos_Timer.hpp>
#include <shylu_fastic.hpp>

namespace Ifpack2
{
namespace Details
{

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Fic<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Fic(Teuchos::RCP<const TRowMatrix> A) :
  FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>(A) {}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
int Fic<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getSweeps() const
{
  return localPrec_->getNFact();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
std::string Fic<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getSpTrsvType() const
{
  return localPrec_->getSpTrsvType();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
int Fic<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getNTrisol() const
{
  return localPrec_->getNTrisol();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Fic<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
checkLocalIC() const
{
  localPrec_->checkIC();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Fic<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
initLocalPrec()
{
  auto nRows = this->mat_->getLocalNumRows();
  auto& p = this->params_;
  localPrec_ = Teuchos::rcp(new LocalFIC(this->localRowPtrs_, this->localColInds_, this->localValues_, nRows, (p.sptrsv_algo != FastILU::SpTRSV::Fast),
                                         p.nFact, p.nTrisol, p.level, p.omega, p.shift, p.guessFlag ? 1 : 0, p.blockSizeILU, p.blockSize));
  localPrec_->initialize();
  this->initTime_ = localPrec_->getInitializeTime();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Fic<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
computeLocalPrec()
{
  //update values in local prec (until compute(), values aren't needed)
  localPrec_->setValues(this->localValues_);
  localPrec_->compute();
  this->computeTime_ = localPrec_->getComputeTime();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Fic<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
applyLocalPrec(ImplScalarArray x, ImplScalarArray y) const
{
  localPrec_->apply(x, y);
  //since this may be applied to multiple vectors, add to applyTime_ instead of setting it
  this->applyTime_ += localPrec_->getApplyTime();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
std::string Fic<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getName() const
{
  return "Fic";
}

#define IFPACK2_DETAILS_FIC_INSTANT(S, L, G, N) \
template class Ifpack2::Details::Fic<S, L, G, N>;

} //namespace Details
} //namespace Ifpack2

#endif

