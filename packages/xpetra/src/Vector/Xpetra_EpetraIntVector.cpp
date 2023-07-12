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
#include "Xpetra_EpetraIntVector.hpp"
#include "Xpetra_EpetraImport.hpp"
#include "Xpetra_EpetraExport.hpp"

namespace Xpetra {


  // TODO: move that elsewhere
  template<class GlobalOrdinal, class Node>
  Epetra_IntVector & toEpetra(Vector<int, int, GlobalOrdinal,Node> &x) {
    XPETRA_DYNAMIC_CAST(      EpetraIntVectorT<GlobalOrdinal XPETRA_COMMA Node>, x, tX, "toEpetra");
    return *tX.getEpetra_IntVector();
  }

  template<class GlobalOrdinal, class Node>
  const Epetra_IntVector & toEpetra(const Vector<int, int, GlobalOrdinal, Node> &x) {
    XPETRA_DYNAMIC_CAST(const EpetraIntVectorT<GlobalOrdinal XPETRA_COMMA Node>, x, tX, "toEpetra");
    return *tX.getEpetra_IntVector();
  }
  //

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
    (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraIntVectorT<int, Xpetra::EpetraNode >;
template Epetra_IntVector & toEpetra<int,Xpetra::EpetraNode >(Vector<int, int, int, Xpetra::EpetraNode> &);
template const Epetra_IntVector & toEpetra<int, Xpetra::EpetraNode >(const Vector<int, int, int, Xpetra::EpetraNode> &);
#endif
#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraIntVectorT<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode >;
template Epetra_IntVector & toEpetra<int,Tpetra::KokkosCompat::KokkosSerialWrapperNode >(Vector<int, int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
template const Epetra_IntVector & toEpetra<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode >(const Vector<int, int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraIntVectorT<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template Epetra_IntVector & toEpetra<int,Tpetra::KokkosCompat::KokkosThreadsWrapperNode >(Vector<int, int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
template const Epetra_IntVector & toEpetra<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode >(const Vector<int, int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraIntVectorT<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode >;
template Epetra_IntVector & toEpetra<int,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode >(Vector<int, int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
template const Epetra_IntVector & toEpetra<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode >(const Vector<int, int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
template class EpetraIntVectorT<int, default_node_type >;
template Epetra_IntVector & toEpetra<int,default_node_type >(Vector<int, int, int, default_node_type> &);
template const Epetra_IntVector & toEpetra<int,default_node_type >(const Vector<int, int, int, default_node_type> &);
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
template class EpetraIntVectorT<int, default_node_type >;
template Epetra_IntVector & toEpetra<int,default_node_type >(Vector<int, int, int, default_node_type> &);
template const Epetra_IntVector & toEpetra<int,default_node_type >(const Vector<int, int, int, default_node_type> &);
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraIntVectorT<int, default_node_type >;
template Epetra_IntVector & toEpetra<int,default_node_type >(Vector<int, int, int, default_node_type> &);
template const Epetra_IntVector & toEpetra<int,default_node_type >(const Vector<int, int, int, default_node_type> &);
#endif // HAVE_XPETRA_TPETRA
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
    (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraIntVectorT<long long, Xpetra::EpetraNode >;
template Epetra_IntVector & toEpetra<long long,Xpetra::EpetraNode >(Vector<int, int, long long, Xpetra::EpetraNode> &);
template const Epetra_IntVector & toEpetra<long long, Xpetra::EpetraNode >(const Vector<int, int, long long, Xpetra::EpetraNode> &);
#endif
#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraIntVectorT<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode >;
template Epetra_IntVector & toEpetra<long long,Tpetra::KokkosCompat::KokkosSerialWrapperNode >(Vector<int, int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
template const Epetra_IntVector & toEpetra<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode >(const Vector<int, int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraIntVectorT<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template Epetra_IntVector & toEpetra<long long,Tpetra::KokkosCompat::KokkosThreadsWrapperNode >(Vector<int, int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
template const Epetra_IntVector & toEpetra<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode >(const Vector<int, int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraIntVectorT<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode >;
template Epetra_IntVector & toEpetra<long long,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode >(Vector<int, int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
template const Epetra_IntVector & toEpetra<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode >(const Vector<int, int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
template class EpetraIntVectorT<long long, default_node_type >;
template Epetra_IntVector & toEpetra<long long,default_node_type >(Vector<int, int, long long, default_node_type> &);
template const Epetra_IntVector & toEpetra<long long,default_node_type >(const Vector<int, int, long long, default_node_type> &);
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
template class EpetraIntVectorT<long long, default_node_type >;
template Epetra_IntVector & toEpetra<long long,default_node_type >(Vector<int, int, long long, default_node_type> &);
template const Epetra_IntVector & toEpetra<long long,default_node_type >(const Vector<int, int, long long, default_node_type> &);
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraIntVectorT<long long, default_node_type >;
template Epetra_IntVector & toEpetra<long long,default_node_type >(Vector<int, int, long long, default_node_type> &);
template const Epetra_IntVector & toEpetra<long long,default_node_type >(const Vector<int, int, long long, default_node_type> &);
#endif // HAVE_XPETRA_TPETRA
#endif

} // namespace Xpetra
