/*
//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ************************************************************************
//@HEADER
*/

#include "BelosTpetraAdapter.hpp"

#ifdef HAVE_BELOS_TPETRA_TIMERS

#include <Kokkos_SerialNode.hpp>
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::SerialNode   > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::SerialNode   > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::SerialNode   > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::SerialNode   > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::SerialNode   > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::SerialNode   > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::SerialNode   > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::SerialNode   > >::mvTimesMatAddMvTimer_ = Teuchos::null;

#ifdef HAVE_KOKKOSCLASSIC_TBB
#include <Kokkos_TBBNode.hpp>
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::TBBNode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::TBBNode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::TBBNode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::TBBNode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::TBBNode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::TBBNode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::TBBNode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::TBBNode      > >::mvTransMvTimer_ = Teuchos::null;
#endif

#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#include <Kokkos_TPINode.hpp>
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::TPINode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::TPINode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::TPINode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::TPINode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::TPINode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::TPINode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::TPINode      > >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::TPINode      > >::mvTimesMatAddMvTimer_ = Teuchos::null;
#endif

#ifdef HAVE_KOKKOSCLASSIC_THRUST
#include <Kokkos_ThrustGPUNode.hpp>
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::ThrustGPUNode> >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<float, Tpetra::MultiVector<float, int, int, Kokkos::ThrustGPUNode> >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::ThrustGPUNode> >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<double, Tpetra::MultiVector<double, int, int, Kokkos::ThrustGPUNode> >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::ThrustGPUNode> >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<float>, Tpetra::MultiVector<std::complex<float>, int, int, Kokkos::ThrustGPUNode> >::mvTimesMatAddMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::ThrustGPUNode> >::mvTransMvTimer_ = Teuchos::null;
template <> Teuchos::RCP<Teuchos::Time> Belos::MultiVecTraits<std::complex<double>, Tpetra::MultiVector<std::complex<double>, int, int, Kokkos::ThrustGPUNode> >::mvTimesMatAddMvTimer_ = Teuchos::null;
#endif

#endif
