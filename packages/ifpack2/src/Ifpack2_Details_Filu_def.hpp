// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// @file Ifpack2_filu_def.hpp

#ifndef __IFPACK2_FILU_DEF_HPP__ 
#define __IFPACK2_FILU_DEF_HPP__ 

#include "Ifpack2_Details_Filu_decl.hpp"
#include "Ifpack2_Details_CrsArrays.hpp"
#include "Ifpack2_Details_getCrsMatrix.hpp"
#include <Kokkos_Timer.hpp>
#include <shylu_fastilu.hpp>

namespace Ifpack2
{
namespace Details
{

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, bool BlockCrsEnabled>
Filu<Scalar, LocalOrdinal, GlobalOrdinal, Node, BlockCrsEnabled>::
Filu(Teuchos::RCP<const TRowMatrix> A) :
  FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>(A) {}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, bool BlockCrsEnabled>
int Filu<Scalar, LocalOrdinal, GlobalOrdinal, Node, BlockCrsEnabled>::
getSweeps() const
{
  return localPrec_->getNFact();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, bool BlockCrsEnabled>
std::string Filu<Scalar, LocalOrdinal, GlobalOrdinal, Node, BlockCrsEnabled>::
getSpTrsvType() const
{
  return localPrec_->getSpTrsvType();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, bool BlockCrsEnabled>
int Filu<Scalar, LocalOrdinal, GlobalOrdinal, Node, BlockCrsEnabled>::
getNTrisol() const
{
  return localPrec_->getNTrisol();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, bool BlockCrsEnabled>
void Filu<Scalar, LocalOrdinal, GlobalOrdinal, Node, BlockCrsEnabled>::
checkLocalILU() const
{
  localPrec_->checkILU();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, bool BlockCrsEnabled>
void Filu<Scalar, LocalOrdinal, GlobalOrdinal, Node, BlockCrsEnabled>::
checkLocalIC() const
{
  localPrec_->checkIC();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, bool BlockCrsEnabled>
void Filu<Scalar, LocalOrdinal, GlobalOrdinal, Node, BlockCrsEnabled>::
initLocalPrec()
{
  auto nRows = this->mat_->getLocalNumRows();
  auto& p = this->params_;
  auto matCrs = Ifpack2::Details::getCrsMatrix(this->mat_);

  if (p.blockCrsSize > 1 && !BlockCrsEnabled) {
    throw std::runtime_error("Must use prec type FAST_ILU_B if you want blockCrs support");
  }

  bool skipSortMatrix = !matCrs.is_null() && matCrs->getCrsGraph()->isSorted() &&
                        !p.use_metis;
  localPrec_ =
    Teuchos::rcp(new LocalFILU(skipSortMatrix, this->localRowPtrs_, this->localColInds_, this->localValues_, nRows, p.sptrsv_algo,
                               p.nFact, p.nTrisol, p.level, p.omega, p.shift, p.guessFlag ? 1 : 0, p.blockSizeILU, p.blockSize,
                               p.blockCrsSize));

  #ifdef HAVE_IFPACK2_METIS
  if (p.use_metis) {
    localPrec_->setMetisPerm(this->metis_perm_, this->metis_iperm_);
  }
  #endif

  localPrec_->initialize();
  this->initTime_ = localPrec_->getInitializeTime();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, bool BlockCrsEnabled>
void Filu<Scalar, LocalOrdinal, GlobalOrdinal, Node, BlockCrsEnabled>::
computeLocalPrec()
{
  //update values in local prec (until compute(), values aren't needed)
  localPrec_->setValues(this->localValues_);
  localPrec_->compute();
  this->computeTime_ = localPrec_->getComputeTime();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, bool BlockCrsEnabled>
void Filu<Scalar, LocalOrdinal, GlobalOrdinal, Node, BlockCrsEnabled>::
applyLocalPrec(ImplScalarArray x, ImplScalarArray y) const
{
  localPrec_->apply(x, y);

  //since this may be applied to multiple vectors, add to applyTime_ instead of setting it
  this->applyTime_ += localPrec_->getApplyTime();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, bool BlockCrsEnabled>
std::string Filu<Scalar, LocalOrdinal, GlobalOrdinal, Node, BlockCrsEnabled>::
getName() const
{
  return "Filu";
}

#define IFPACK2_DETAILS_FILU_INSTANT(S, L, G, N) \
template class Ifpack2::Details::Filu<S, L, G, N,false>; \
template class Ifpack2::Details::Filu<S, L, G, N,true>;

} //namespace Details
} //namespace Ifpack2

#endif

