// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_Array.hpp>
#include "Xpetra_EpetraCrsMatrix.hpp"

namespace Xpetra {

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
    (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
  template class EpetraCrsMatrixT<int, Xpetra::EpetraNode >;
#endif

#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraCrsMatrixT<int, Kokkos::Compat::KokkosSerialWrapperNode >;
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraCrsMatrixT<int, Kokkos::Compat::KokkosThreadsWrapperNode>;
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraCrsMatrixT<int, Kokkos::Compat::KokkosOpenMPWrapperNode >;
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Kokkos::Compat::KokkosCudaWrapperNode default_node_type;
template class EpetraCrsMatrixT<int, default_node_type >;
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Kokkos::Compat::KokkosHIPWrapperNode default_node_type;
template class EpetraCrsMatrixT<int, default_node_type >;
#endif
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
    (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
  template class EpetraCrsMatrixT<long long, Xpetra::EpetraNode >;
#endif

#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraCrsMatrixT<long long, Kokkos::Compat::KokkosSerialWrapperNode >;
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraCrsMatrixT<long long, Kokkos::Compat::KokkosThreadsWrapperNode>;
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraCrsMatrixT<long long, Kokkos::Compat::KokkosOpenMPWrapperNode >;
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Kokkos::Compat::KokkosCudaWrapperNode default_node_type;
template class EpetraCrsMatrixT<long long, default_node_type >;
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Kokkos::Compat::KokkosHIPWrapperNode default_node_type;
template class EpetraCrsMatrixT<long long, default_node_type >;
#endif
#endif

}
