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

#ifndef TPETRANEW_VECTOR_DEF_HPP
#define TPETRANEW_VECTOR_DEF_HPP

/// \file TpetraNew_Vector_def.hpp
/// \brief Definition of the TpetraNew::Vector class
///
/// If you want to use TpetraNew::Vector, include
/// "TpetraNew_Vector.hpp" (a file which CMake generates and installs
/// for you).  If you only want the declaration of TpetraNew::Vector,
/// include "TpetraNew_Vector_decl.hpp".

#include "TpetraNew_Map.hpp"
#include "TpetraNew_MultiVector.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "KokkosCompat_View.hpp"
#include "Teuchos_CommHelpers.hpp"

namespace TpetraNew {

  template <class Scalar>
  Vector<Scalar>::
  Vector ()
    : base_type ()
  {}

  template <class Scalar>
  Vector<Scalar>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const bool zeroOut)
    : base_type (map, 1, zeroOut)
  {}

  template <class Scalar>
  Vector<Scalar>::
  Vector (const Vector<Scalar>& source)
    : base_type (source)
  {}

  template <class Scalar>
  Vector<Scalar>::
  Vector (const Vector<Scalar>& source,
          const Teuchos::DataAccess copyOrView)
    : base_type (source, copyOrView)
  {}

  template <class Scalar>
  Vector<Scalar>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const Teuchos::ArrayView<const Scalar>& values)
    : base_type (map, values, values.size (), 1)
  {}

  template <class Scalar>
  Vector<Scalar>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const dual_view_type& view)
    : base_type (map, view)
  {}

  template <class Scalar>
  Vector<Scalar>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const dual_view_type& view,
          const dual_view_type& origView)
    : base_type (map, view, origView)
  {}

  template <class Scalar>
  Vector<Scalar>::
  Vector (const MultiVector<Scalar>& X,
          const size_t j)
    : base_type (X, j)
  {}

  template <class Scalar>
  Vector<Scalar>::
  ~Vector ()
  {}

  template <class Scalar>
  void
  Vector<Scalar>::
  replaceGlobalValue (const global_ordinal_type globalRow,
		      const Scalar& value) const
  {
    this->base_type::replaceGlobalValue (globalRow, 0, value);
  }

  template <class Scalar>
  void
  Vector<Scalar>::
  sumIntoGlobalValue (const global_ordinal_type globalRow,
                      const Scalar& value,
                      const bool atomic) const
  {
    this->base_type::sumIntoGlobalValue (globalRow, 0, value, atomic);
  }

  template <class Scalar>
  void
  Vector<Scalar>::
  replaceLocalValue (const local_ordinal_type myRow,
		     const Scalar& value) const
  {
    this->base_type::replaceLocalValue (myRow, 0, value);
  }

  template <class Scalar>
  void
  Vector<Scalar>::
  sumIntoLocalValue (const local_ordinal_type globalRow,
                     const Scalar& value,
                     const bool atomic) const
  {
    this->base_type::sumIntoLocalValue (globalRow, 0, value, atomic);
  }

  template <class Scalar>
  void
  Vector<Scalar>::
  get1dCopy (const Teuchos::ArrayView<Scalar>& A) const
  {
    const auto lda = this->getLocalLength ();
    this->get1dCopy (A, lda);
  }

  template <class Scalar>
  typename Vector<Scalar>::dot_type
  Vector<Scalar>::
  dot (const Vector<Scalar>& y) const
  {
    dot_type result;
    this->dot (y, Teuchos::arrayView (&result, 1));
    return result;
  }

  template <class Scalar>
  Scalar
  Vector<Scalar>::
  meanValue () const
  {
    Scalar mean;
    this->meanValue (Teuchos::arrayView (&mean, 1));
    return mean;
  }

  template <class Scalar>
  typename Vector<Scalar>::mag_type
  Vector<Scalar>::
  norm1 () const
  {
    mag_type norm;
    this->norm1 (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar>
  typename Vector<Scalar>::mag_type
  Vector<Scalar>::
  norm2 () const
  {
    mag_type norm;
    this->norm2 (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar>
  typename Vector<Scalar>::mag_type
  Vector<Scalar>::
  normInf () const
  {
    mag_type norm;
    this->normInf (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar>
  std::string Vector<Scalar>::
  description () const
  {
    return this->descriptionImpl ("TpetraNew::Vector");
  }

  template <class Scalar>
  void Vector<Scalar>::
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    this->describeImpl (out, "TpetraNew::Vector", verbLevel);
  }

  template <class Scalar>
  Teuchos::RCP<const Vector<Scalar> >
  Vector<Scalar>::
  offsetView (const Teuchos::RCP<const map_type>& subMap,
              const size_t offset) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::rcp;
    typedef Vector<Scalar> V;
    using LO = local_ordinal_type;

    const LO newNumRows = subMap->getMyNumIndices ();
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

    const std::pair<size_t, size_t> offsetPair (offset, offset + size_t (newNumRows));
    // Need 'this->' to get view_ and origView_ from parent class.
    return rcp (new V (subMap,
                       subview (this->view_, offsetPair, ALL ()),
                       this->origView_));
  }

  template <class Scalar>
  Teuchos::RCP<Vector<Scalar> >
  Vector<Scalar>::
  offsetViewNonConst (const Teuchos::RCP<const map_type>& subMap,
                      const size_t offset)
  {
    using Teuchos::rcp_const_cast;
    return rcp_const_cast<Vector<Scalar> > (this->offsetView (subMap, offset));
  }

  template <class Scalar>
  Vector<Scalar>
  createCopy (const Vector<Scalar>& src)
  {
    // The 2-argument copy constructor with second argument =
    // Teuchos::Copy does a deep copy of its input.  Returning by
    // value is fine, because Vector has view semantics, just like
    // Kokkos::View.
    return Vector<Scalar> (src, Teuchos::Copy);
  }

} // namespace TpetraNew

/// \macro TPETRANEW_VECTOR_INSTANT
/// \brief Explicit instantiation macro for TpetraNew::Vector.
///
/// \warning This is ONLY for Trilinos developers and expert users.
///   If you don't know what explicit template instantiation is, you
///   shouldn't even THINK about using this macro.
///
/// \warning This macro must be invoked within the TpetraNew namespace!
#define TPETRANEW_VECTOR_INSTANT( SCALAR ) \
  template class Vector< SCALAR >; \
  template Vector< SCALAR > createCopy (const Vector< SCALAR >& src);

#endif // TPETRANEW_VECTOR_DEF_HPP
