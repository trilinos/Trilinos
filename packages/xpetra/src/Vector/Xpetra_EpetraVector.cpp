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
#include "Xpetra_EpetraVector.hpp"

//TODO: replace double -> Scalar etc.

namespace Xpetra {

  template<class EpetraGlobalOrdinal, class Node>
  EpetraVectorT<EpetraGlobalOrdinal, Node>::EpetraVectorT(const Teuchos::RCP<const Map<int,GlobalOrdinal,Node> > &map, bool zeroOut) : EpetraMultiVectorT<GlobalOrdinal, Node>(map,1,zeroOut) { }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraVectorT<EpetraGlobalOrdinal,Node>::replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { XPETRA_MONITOR("EpetraVectorT::replaceGlobalValue"); this->EpetraMultiVectorT<GlobalOrdinal,Node>::getEpetra_MultiVector()->ReplaceGlobalValue(globalRow, 0, value); }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraVectorT<EpetraGlobalOrdinal,Node>::sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { XPETRA_MONITOR("EpetraVectorT::sumIntoGlobalValue");this->EpetraMultiVectorT<GlobalOrdinal,Node>::getEpetra_MultiVector()->SumIntoGlobalValue(globalRow, 0, value); }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraVectorT<EpetraGlobalOrdinal,Node>::replaceLocalValue(LocalOrdinal myRow, const Scalar &value) { XPETRA_MONITOR("EpetraVectorT::replaceLocalValue");this->EpetraMultiVectorT<GlobalOrdinal,Node>::getEpetra_MultiVector()->ReplaceMyValue(myRow, 0, value); }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraVectorT<EpetraGlobalOrdinal,Node>::sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value) { XPETRA_MONITOR("EpetraVectorT::sumIntoLocalValue");this->EpetraMultiVectorT<GlobalOrdinal,Node>::getEpetra_MultiVector()->SumIntoMyValue(myRow, 0, value); }

  template<class EpetraGlobalOrdinal, class Node>
  double EpetraVectorT<EpetraGlobalOrdinal,Node>::dot(const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &a) const {
    XPETRA_MONITOR("EpetraVectorT::dot");

    XPETRA_DYNAMIC_CAST(const EpetraVectorT, a, tA, "This Xpetra::EpetraVectorT method only accept Xpetra::EpetraVectorT as input arguments.");
    //      return getEpetra_Vector()->Dot(*tA.getEpetra_Vector());

    // other way: use the MultiVector Dot instead of VectorDot:
    double r;
    this->EpetraMultiVectorT<GlobalOrdinal,Node>::getEpetra_MultiVector()->Epetra_MultiVector::Dot(*tA.getEpetra_MultiVector(), &r);
    return r;
  }

  template<class EpetraGlobalOrdinal, class Node>
  Teuchos::ScalarTraits<double>::magnitudeType EpetraVectorT<EpetraGlobalOrdinal,Node>::norm1() const { XPETRA_MONITOR("EpetraVectorT::norm1"); double r; this->EpetraMultiVectorT<GlobalOrdinal,Node>::getEpetra_MultiVector()->Norm1(&r); return r; }

  template<class EpetraGlobalOrdinal, class Node>
  Teuchos::ScalarTraits<double>::magnitudeType EpetraVectorT<EpetraGlobalOrdinal,Node>::norm2() const { XPETRA_MONITOR("EpetraVectorT::norm2"); double r; this->EpetraMultiVectorT<GlobalOrdinal,Node>::getEpetra_MultiVector()->Norm2(&r); return r; }

  template<class EpetraGlobalOrdinal, class Node>
  Teuchos::ScalarTraits<double>::magnitudeType EpetraVectorT<EpetraGlobalOrdinal,Node>::normInf() const { XPETRA_MONITOR("EpetraVectorT::normInf"); double r; this->EpetraMultiVectorT<GlobalOrdinal,Node>::getEpetra_MultiVector()->NormInf(&r); return r; }

  template<class EpetraGlobalOrdinal, class Node>
  double EpetraVectorT<EpetraGlobalOrdinal,Node>::meanValue() const {
    XPETRA_MONITOR("EpetraVectorT::meanValue");
    double r;
    this->EpetraMultiVectorT<GlobalOrdinal,Node>::getEpetra_MultiVector()->MeanValue(&r);
    return r;
  }

    template<class EpetraGlobalOrdinal, class Node>
    std::string EpetraVectorT<EpetraGlobalOrdinal,Node>::description() const {
    XPETRA_MONITOR("EpetraVectorT::description");
    // This implementation come from Epetra_Vector_def.hpp (without modification)
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{length="<<this->getGlobalLength()
        << "}";
    return oss.str();
  }

    template<class EpetraGlobalOrdinal, class Node>
    void EpetraVectorT<EpetraGlobalOrdinal,Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
      XPETRA_MONITOR("EpetraVectorT::describe");

      if (verbLevel > Teuchos::VERB_NONE) {
        getEpetra_Vector()->Print (out);
      }
    }

  template<class EpetraGlobalOrdinal, class Node>
  EpetraVectorT<EpetraGlobalOrdinal,Node>::EpetraVectorT(const RCP<Epetra_MultiVector> &mv, size_t j)
    : EpetraMultiVectorT<GlobalOrdinal,Node>(rcp((*mv)(j), false)), // view of the vector number j. false == I do not own the data.
      internalRefToBaseMV_(mv)                 // keep an internal reference to the initial MultiVector to avoid desallocation of the view.
  {
    // The view of the internal data of 'mv' is only valid until the destruction of 'mv'.
    // The new vector hold an internal reference to 'mv' in order to keep the view valid after disappearance of 'mv' references in user code.
    // This implements the logic of subArray rcp (as required by the Tpetra interface).
  }

  // TODO: move that elsewhere
  template<class GlobalOrdinal, class Node>
  Epetra_Vector & toEpetra(Vector<double, int, GlobalOrdinal,Node> &x) {
    XPETRA_DYNAMIC_CAST(      EpetraVectorT<GlobalOrdinal COMMA Node>, x, tX, "toEpetra");
    return *tX.getEpetra_Vector();
  }

  template<class GlobalOrdinal, class Node>
  const Epetra_Vector & toEpetra(const Vector<double, int, GlobalOrdinal, Node> &x) {
    XPETRA_DYNAMIC_CAST(const EpetraVectorT<GlobalOrdinal COMMA Node>, x, tX, "toEpetra");
    return *tX.getEpetra_Vector();
  }
  //

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES

#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_SERIAL
template class EpetraVectorT<int, Kokkos::Compat::KokkosSerialWrapperNode >;
//template RCP<Vector<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode > > toXpetra<int, Kokkos::Compat::KokkosSerialWrapperNode>(RCP<Epetra_Vector>);
template Epetra_Vector & toEpetra<int,Kokkos::Compat::KokkosSerialWrapperNode>(Vector<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode> &);
template const Epetra_Vector & toEpetra<int, Kokkos::Compat::KokkosSerialWrapperNode >(const Vector<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode > &);
#endif
#ifdef HAVE_XPETRA_PTHREAD
template class EpetraVectorT<int, Kokkos::Compat::KokkosThreadsWrapperNode>;
//template RCP<Vector<double, int, int, Kokkos::Compat::KokkosThreadsWrapperNode > > toXpetra<int, Kokkos::Compat::KokkosThreadsWrapperNode>(RCP<Epetra_Vector>);
template Epetra_Vector & toEpetra<int,Kokkos::Compat::KokkosThreadsWrapperNode>(Vector<double, int, int,Kokkos::Compat::KokkosThreadsWrapperNode> &);
template const Epetra_Vector & toEpetra<int, Kokkos::Compat::KokkosThreadsWrapperNode >(const Vector<double, int, int, Kokkos::Compat::KokkosThreadsWrapperNode > &);
#endif
#ifdef HAVE_XPETRA_OPENMP
template class EpetraVectorT<int, Kokkos::Compat::KokkosOpenMPWrapperNode >;
//template RCP<Vector<double, int, int, Kokkos::Compat::KokkosOpenMPWrapperNode > > toXpetra<int, Kokkos::Compat::KokkosOpenMPWrapperNode>(RCP<Epetra_Vector>);
template Epetra_Vector & toEpetra<int,Kokkos::Compat::KokkosOpenMPWrapperNode>(Vector<double, int, int,Kokkos::Compat::KokkosOpenMPWrapperNode> &);
template const Epetra_Vector & toEpetra<int, Kokkos::Compat::KokkosOpenMPWrapperNode >(const Vector<double, int, int, Kokkos::Compat::KokkosOpenMPWrapperNode > &);
#endif
#ifdef HAVE_XPETRA_CUDA
typedef Kokkos::Compat::KokkosCudaWrapperNode default_node_type;
template class EpetraVectorT<int, default_node_type >;
//template RCP<Vector<double, int, int, default_node_type > toXpetra<int, default_node_type>(RCP<Epetra_Vector>);
template Epetra_Vector & toEpetra<int,default_node_type >(Vector<double, int, int,default_node_type> &);
template const Epetra_Vector & toEpetra<int, default_node_type >(const Vector<double, int, int, default_node_type > &);
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Kokkos::Compat::KokkosSerialWrapperNode default_node_type;
template class EpetraVectorT<int, default_node_type >;
//template RCP<Vector<double, int, int, default_node_type > toXpetra<int, default_node_type>(RCP<Epetra_Vector>);
template Epetra_Vector & toEpetra<int,default_node_type >(Vector<double, int, int,default_node_type> &);
template const Epetra_Vector & toEpetra<int, default_node_type >(const Vector<double, int, int, default_node_type > &);
#endif // HAVE_XPETRA_TPETRA
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_SERIAL
template class EpetraVectorT<long long, Kokkos::Compat::KokkosSerialWrapperNode >;
//template RCP<Vector<double, int, long long, Kokkos::Compat::KokkosSerialWrapperNode > > toXpetra<long long, Kokkos::Compat::KokkosSerialWrapperNode>(RCP<Epetra_Vector>);
template Epetra_Vector & toEpetra<long long,Kokkos::Compat::KokkosSerialWrapperNode>(Vector<double, int, long long, Kokkos::Compat::KokkosSerialWrapperNode> &);
template const Epetra_Vector & toEpetra<long long, Kokkos::Compat::KokkosSerialWrapperNode >(const Vector<double, int, long long, Kokkos::Compat::KokkosSerialWrapperNode > &);
#endif
#ifdef HAVE_XPETRA_PTHREAD
template class EpetraVectorT<long long, Kokkos::Compat::KokkosThreadsWrapperNode>;
//template RCP<Vector<double, int, long long, Kokkos::Compat::KokkosThreadsWrapperNode > > toXpetra<long long, Kokkos::Compat::KokkosThreadsWrapperNode>(RCP<Epetra_Vector>);
template Epetra_Vector & toEpetra<long long,Kokkos::Compat::KokkosThreadsWrapperNode>(Vector<double, int, long long, Kokkos::Compat::KokkosThreadsWrapperNode> &);
template const Epetra_Vector & toEpetra<long long, Kokkos::Compat::KokkosThreadsWrapperNode >(const Vector<double, int, long long, Kokkos::Compat::KokkosThreadsWrapperNode > &);
#endif
#ifdef HAVE_XPETRA_OPENMP
template class EpetraVectorT<long long, Kokkos::Compat::KokkosOpenMPWrapperNode >;
//template RCP<Vector<double, int, long long, Kokkos::Compat::KokkosOpenMPWrapperNode > > toXpetra<long long, Kokkos::Compat::KokkosOpenMPWrapperNode>(RCP<Epetra_Vector>);
template Epetra_Vector & toEpetra<long long,Kokkos::Compat::KokkosOpenMPWrapperNode>(Vector<double, int, long long, Kokkos::Compat::KokkosOpenMPWrapperNode> &);
template const Epetra_Vector & toEpetra<long long, Kokkos::Compat::KokkosOpenMPWrapperNode >(const Vector<double, int, long long, Kokkos::Compat::KokkosOpenMPWrapperNode > &);
#endif
#ifdef HAVE_XPETRA_CUDA
typedef Kokkos::Compat::KokkosCudaWrapperNode default_node_type;
template class EpetraVectorT<long long, default_node_type >;
//template RCP<Vector<double, int, long long, default_node_type > toXpetra<long long, default_node_type>(RCP<Epetra_Vector>);
template Epetra_Vector & toEpetra<long long,default_node_type >(Vector<double, int, long long, default_node_type> &);
template const Epetra_Vector & toEpetra<long long, default_node_type >(const Vector<double, int, long long, default_node_type > &);
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Kokkos::Compat::KokkosSerialWrapperNode default_node_type;
template class EpetraVectorT<long long, default_node_type >;
//template RCP<Vector<double, int, long long, default_node_type > toXpetra<int, default_node_type>(RCP<Epetra_Vector>);
template Epetra_Vector & toEpetra<long long,default_node_type >(Vector<double, int, long long,default_node_type> &);
template const Epetra_Vector & toEpetra<long long, default_node_type >(const Vector<double, int, long long, default_node_type > &);
#endif // HAVE_XPETRA_TPETRA
#endif

}
