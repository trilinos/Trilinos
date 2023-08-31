/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

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
  const int blockCrsSize = p.blockCrs ? p.blockCrsSize : 1;
  localPrec_ = Teuchos::rcp(new LocalFILUB(skipSortMatrix, this->localRowPtrs_, this->localColInds_, this->localValues_, nRows, p.sptrsv_algo,
                                          p.nFact, p.nTrisol, p.level, p.omega, p.shift, p.guessFlag ? 1 : 0, p.blockSizeILU, p.blockSize,
                                          blockCrsSize));
  localPrec2_ = Teuchos::rcp(new LocalFILU(skipSortMatrix, this->localRowPtrs2_, this->localColInds2_, this->localValues2_, nRows*blockCrsSize, p.sptrsv_algo,
                                           p.nFact, p.nTrisol, p.level, p.omega, p.shift, p.guessFlag ? 1 : 0, p.blockSizeILU, p.blockSize));
  localPrec3_ = Teuchos::rcp(new LocalFILUOrig(skipSortMatrix, this->localRowPtrs2_, this->localColInds2_, this->localValues2_, nRows*blockCrsSize, p.sptrsv_algo,
                                               p.nFact, p.nTrisol, p.level, p.omega, p.shift, p.guessFlag ? 1 : 0, p.blockSizeILU, p.blockSize));

  #ifdef HAVE_IFPACK2_METIS
  if (p.use_metis) {
    localPrec_->setMetisPerm(this->metis_perm_, this->metis_iperm_);
  }
  #endif

  localPrec_->initialize();
  localPrec2_->initialize();
  localPrec3_->initialize();
  localPrec2_->verify(*localPrec3_, "initialization_unblocked_vs_orig", true);
  localPrec_->verify(*localPrec2_, "initialization_blocked_vs_unblocked", true);
  this->initTime_ = localPrec_->getInitializeTime();

  std::cout << "JGF initialization of blocked took: " << localPrec_->getInitializeTime() << std::endl;
  std::cout << "JGF initialization of unblocked took: " << localPrec2_->getInitializeTime() << std::endl;
  std::cout << "JGF initialization of orig took: " << localPrec3_->getInitializeTime() << std::endl;
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, bool BlockCrsEnabled>
void Filu<Scalar, LocalOrdinal, GlobalOrdinal, Node, BlockCrsEnabled>::
computeLocalPrec()
{
  //update values in local prec (until compute(), values aren't needed)
  localPrec_->setValues(this->localValues_);
  localPrec2_->setValues(this->localValues2_);
  localPrec3_->setValues(this->localValues2_);
  localPrec_->compute();
  localPrec2_->compute();
  localPrec3_->compute();
  this->computeTime_ = localPrec_->getComputeTime();

  localPrec2_->verify(*localPrec3_, "compute_unblocked_vs_orig");
  localPrec_->verify(*localPrec2_, "compute_blocked_vs_unblocked");

  std::cout << "JGF compute of blocked took: " << localPrec_->getComputeTime() << std::endl;
  std::cout << "JGF compute of unblocked took: " << localPrec2_->getComputeTime() << std::endl;
  std::cout << "JGF compute of orig took: " << localPrec3_->getComputeTime() << std::endl;
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, bool BlockCrsEnabled>
void Filu<Scalar, LocalOrdinal, GlobalOrdinal, Node, BlockCrsEnabled>::
applyLocalPrec(ImplScalarArray x, ImplScalarArray y) const
{
  ImplScalarArray x1("x1", x.extent(0)), x2("x2", x.extent(0)),
    y1("y1", y.extent(0)), y2("y2", y.extent(0));
  Kokkos::deep_copy(x1, x);
  Kokkos::deep_copy(x2, x);
  Kokkos::deep_copy(y1, y);
  Kokkos::deep_copy(y2, y);

  localPrec_->apply(x, y);
  localPrec2_->apply(x1, y1);
  localPrec3_->apply(x2, y2);

  localPrec2_->verify(*localPrec3_, "apply_unblocked_vs_orig");
  localPrec_->verify(*localPrec2_, "apply_blocked_vs_unblocked");
  compare_views(x, x1, "x1");
  compare_views(x, x2, "x2");
  compare_views(y, y1, "y1");
  compare_views(y, y2, "y2");

  //since this may be applied to multiple vectors, add to applyTime_ instead of setting it
  this->applyTime_ += localPrec_->getApplyTime();

  std::cout << "JGF apply of blocked took: " << localPrec_->getApplyTime() << std::endl;
  std::cout << "JGF apply of unblocked took: " << localPrec2_->getApplyTime() << std::endl;
  std::cout << "JGF apply of orig took: " << localPrec3_->getApplyTime() << std::endl;
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

