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
#include "Xpetra_EpetraExport.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

template <class GlobalOrdinal, class Node>
RCP<const Export<int, GlobalOrdinal, Node> > toXpetra(const Epetra_Export *exp) {
  if (exp != NULL) {
    RCP<const Epetra_Export> eexp = rcp(new Epetra_Export(*exp));  // NOTE: non consitent: return pointer, take ref
    return rcp(new Xpetra::EpetraExportT<GlobalOrdinal, Node>(eexp));
  }

  return Teuchos::null;
}

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraExportT<int, Xpetra::EpetraNode>;
template RCP<const Export<int, int, Xpetra::EpetraNode> > toXpetra<int, Xpetra::EpetraNode>(const Epetra_Export *);
#endif

#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraExportT<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>;
template RCP<const Export<int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraExportT<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template RCP<const Export<int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraExportT<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
template RCP<const Export<int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
template class EpetraExportT<int, default_node_type>;
template RCP<const Export<int, int, default_node_type> > toXpetra<int, default_node_type>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
template class EpetraExportT<int, default_node_type>;
template RCP<const Export<int, int, default_node_type> > toXpetra<int, default_node_type>(const Epetra_Export *);
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraExportT<int, default_node_type>;
template RCP<const Export<int, int, default_node_type> > toXpetra<int, default_node_type>(const Epetra_Export *);
#endif  // HAVE_XPETRA_TPETRA
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraExportT<long long, Xpetra::EpetraNode>;
template RCP<const Export<int, long long, Xpetra::EpetraNode> > toXpetra<long long, Xpetra::EpetraNode>(const Epetra_Export *);
#endif

#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraExportT<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>;
template RCP<const Export<int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraExportT<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template RCP<const Export<int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraExportT<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
template RCP<const Export<int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
template class EpetraExportT<long long, default_node_type>;
template RCP<const Export<int, long long, default_node_type> > toXpetra<long long, default_node_type>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
template class EpetraExportT<long long, default_node_type>;
template RCP<const Export<int, long long, default_node_type> > toXpetra<long long, default_node_type>(const Epetra_Export *);
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraExportT<long long, default_node_type>;
template RCP<const Export<int, long long, default_node_type> > toXpetra<long long, default_node_type>(const Epetra_Export *);
#endif  // HAVE_XPETRA_TPETRA
#endif

}  // namespace Xpetra
