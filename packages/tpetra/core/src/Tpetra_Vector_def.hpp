// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_VECTOR_DEF_HPP
#define TPETRA_VECTOR_DEF_HPP

/// \file Tpetra_Vector_def.hpp
/// \brief Definition of the Tpetra::Vector class

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "KokkosCompat_View.hpp"
#include "KokkosBlas1_nrm2w_squared.hpp"
#include "Teuchos_CommHelpers.hpp"

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
  Vector (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source,
          const Teuchos::DataAccess copyOrView)
    : base_type (source, copyOrView)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Vector (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source,
          const Teuchos::RCP<const map_type>& map,
          const local_ordinal_type rowOffset) :
    base_type (source, map, rowOffset)
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
  Vector (const Teuchos::RCP<const map_type>& map,
          const wrapped_dual_view_type& view)
    : base_type (map, view)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Vector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
          const size_t j)
    : base_type (X, j)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceGlobalValue (const GlobalOrdinal globalRow, const Scalar& value) {
    this->base_type::replaceGlobalValue (globalRow, 0, value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sumIntoGlobalValue (const GlobalOrdinal globalRow,
                      const Scalar& value,
                      const bool atomic)
  {
    this->base_type::sumIntoGlobalValue (globalRow, 0, value, atomic);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceLocalValue (const LocalOrdinal myRow, const Scalar& value) {
    this->base_type::replaceLocalValue (myRow, 0, value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sumIntoLocalValue (const LocalOrdinal globalRow,
                     const Scalar& value,
                     const bool atomic)
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
  Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  offsetView (const Teuchos::RCP<const map_type>& subMap,
              const size_t offset) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::rcp;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> V;

    const size_t newNumRows = subMap->getLocalNumElements ();
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

    // Need 'this->' to get view_ from parent class.
    return rcp (new V (subMap, wrapped_dual_view_type(this->view_,Kokkos::pair<int,int>(offset,offset+newNumRows),ALL())));
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  createCopy (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& src)
  {
    using vec_type = Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    // The 2-argument copy constructor with second argument =
    // Teuchos::Copy does a deep copy of its input.
    vec_type dst (src, Teuchos::Copy);

    // The Kokkos refactor version of Vector has view semantics, so
    // returning the Vector directly, rather than through RCP, only
    // does a shallow copy.
    return dst;
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
