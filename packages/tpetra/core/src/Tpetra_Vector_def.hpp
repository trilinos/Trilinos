// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_VECTOR_DEF_HPP
#define TPETRA_VECTOR_DEF_HPP

/// \file Tpetra_Vector_def.hpp
/// \brief Definition of the Tpetra::Vector class
///
/// If you want to use Tpetra::Vector, include "Tpetra_Vector.hpp" (a
/// file which CMake generates and installs for you).  If you only
/// want the declaration of Tpetra::Vector, include
/// "Tpetra_Vector_decl.hpp".

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "KokkosCompat_View.hpp"
#include "KokkosBlas1_nrm2w_squared.hpp"

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Vector ()
    : base_type ()
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const bool zeroOut)
    : base_type (map, 1, zeroOut)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Vector (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source)
    : base_type (source)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Vector (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source,
          const Teuchos::DataAccess copyOrView)
    : base_type (source, copyOrView)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const Teuchos::ArrayView<const Scalar>& values)
    : base_type (map, values, values.size (), 1)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const dual_view_type& view)
    : base_type (map, view)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const dual_view_type& view,
          const dual_view_type& origView)
    : base_type (map, view, origView)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Vector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
          const size_t j)
    : base_type (X, j)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  ~Vector ()
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceGlobalValue (const GlobalOrdinal globalRow, const Scalar& value) const {
    this->base_type::replaceGlobalValue (globalRow, 0, value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sumIntoGlobalValue (const GlobalOrdinal globalRow,
                      const Scalar& value,
                      const bool atomic) const
  {
    this->base_type::sumIntoGlobalValue (globalRow, 0, value, atomic);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceLocalValue (const LocalOrdinal myRow, const Scalar& value) const {
    this->base_type::replaceLocalValue (myRow, 0, value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sumIntoLocalValue (const LocalOrdinal globalRow,
                     const Scalar& value,
                     const bool atomic) const
  {
    this->base_type::sumIntoLocalValue (globalRow, 0, value, atomic);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  get1dCopy (const Teuchos::ArrayView<Scalar>& A) const {
    const size_t lda = this->getLocalLength ();
    this->get1dCopy (A, lda);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dot_type
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  dot (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y) const
  {
    dot_type result;
    this->dot (y, Teuchos::arrayView (&result, 1));
    return result;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Scalar
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  meanValue () const
  {
    Scalar mean;
    this->meanValue (Teuchos::arrayView (&mean, 1));
    return mean;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  norm1 () const
  {
    mag_type norm;
    this->norm1 (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  norm2 () const
  {
    mag_type norm;
    this->norm2 (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  normInf () const
  {
    mag_type norm;
    this->normInf (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type
  TPETRA_DEPRECATED
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  normWeighted (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& weights) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef Kokkos::Details::ArithTraits<impl_scalar_type> ATS;
    typedef Kokkos::Details::ArithTraits<mag_type> ATM;
    typedef Kokkos::View<mag_type, device_type> norm_view_type; // just one
    const char tfecfFuncName[] = "normWeighted: ";

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! this->getMap ()->isCompatible (*weights.getMap ()), std::runtime_error,
      "Vectors do not have compatible Maps:" << std::endl
      << "this->getMap(): " << std::endl << *this->getMap()
      << "weights.getMap(): " << std::endl << *weights.getMap() << std::endl);
#else
    const size_t lclNumRows = this->getLocalLength ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      lclNumRows != weights.getLocalLength (), std::runtime_error,
      "Vectors do not have the same local length.");
#endif // HAVE_TPETRA_DEBUG

    norm_view_type lclNrm ("lclNrm"); // local norm
    mag_type gblNrm = ATM::zero (); // return value

    auto X_lcl = this->template getLocalView<device_type> ();
    auto W_lcl = this->template getLocalView<device_type> ();
    KokkosBlas::nrm2w_squared (lclNrm,
                               subview (X_lcl, ALL (), 0),
                               subview (W_lcl, ALL (), 0));
    const mag_type OneOverN =
      ATM::one () / static_cast<mag_type> (this->getGlobalLength ());
    RCP<const Comm<int> > comm = this->getMap ().is_null () ?
      Teuchos::null : this->getMap ()->getComm ();

    if (! comm.is_null () && this->isDistributed ()) {
      // Assume that MPI can access device memory.
      reduceAll<int, mag_type> (*comm, REDUCE_SUM, 1, lclNrm.data (),
                                &gblNrm);
      gblNrm = ATM::sqrt (gblNrm * OneOverN);
    }
    else {
      auto lclNrm_h = Kokkos::create_mirror_view (lclNrm);
      Kokkos::deep_copy (lclNrm_h, lclNrm);
      gblNrm = ATM::sqrt (ATS::magnitude (lclNrm_h()) * OneOverN);
    }

    return gblNrm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  description () const
  {
    return this->descriptionImpl ("Tpetra::Vector");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    this->describeImpl (out, "Tpetra::Vector", verbLevel);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  createCopy (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& src)
  {
    // The 2-argument copy constructor with second argument =
    // Teuchos::Copy does a deep copy of its input.
    Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> dst (src, Teuchos::Copy);

    // The Kokkos refactor version of Vector has view semantics, so
    // returning the Vector directly, rather than through RCP, only
    // does a shallow copy.
    return dst;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  offsetView (const Teuchos::RCP<const map_type>& subMap,
              const size_t offset) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::rcp;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> V;

    const size_t newNumRows = subMap->getNodeNumElements ();
    const bool tooManyElts = newNumRows + offset > this->getOrigNumLocalRows ();
    if (tooManyElts) {
      const int myRank = this->getMap ()->getComm ()->getRank ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows + offset > this->getLocalLength (), std::runtime_error,
        "Tpetra::Vector::offsetView(NonConst): Invalid input Map.  The input "
        "Map owns " << newNumRows << " entries on process " << myRank << ".  "
        "offset = " << offset << ".  Yet, the Vector contains only "
        << this->getOrigNumLocalRows () << " rows on this process.");
    }

    const std::pair<size_t, size_t> offsetPair (offset, offset + newNumRows);
    // Need 'this->' to get view_ and origView_ from parent class.
    return rcp (new V (subMap,
                       subview (this->view_, offsetPair, ALL ()),
                       this->origView_));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  offsetViewNonConst (const Teuchos::RCP<const map_type>& subMap,
                      const size_t offset)
  {
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> V;
    return Teuchos::rcp_const_cast<V> (this->offsetView (subMap, offset));
  }

} // namespace Tpetra

/// \macro TPETRA_VECTOR_INSTANT
/// \brief Explicit instantiation macro for Tpetra::Vector.
///
/// \warning This is ONLY for Trilinos developers and expert users.
///   If you don't know what explicit template instantiation is, you
///   shouldn't even THINK about using this macro.
///
/// \warning This macro must be invoked within the Tpetra namespace!
#define TPETRA_VECTOR_INSTANT(SCALAR,LO,GO,NODE) \
  template class Vector< SCALAR , LO , GO , NODE >; \
  template Vector< SCALAR , LO , GO , NODE > createCopy (const Vector< SCALAR , LO , GO , NODE >& src);

#endif // TPETRA_VECTOR_DEF_HPP
