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
#include "KokkosCompat_View.hpp"
#include "Kokkos_Blas1_MV.hpp"

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const bool zeroOut)
    : base_type (map, 1, zeroOut)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  Vector (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>& source)
    : base_type (source)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  Vector (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>& source,
          const Teuchos::DataAccess copyOrView)
    : base_type (source, copyOrView)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const Teuchos::ArrayView<const Scalar>& values)
    : base_type (map, values, values.size (), 1)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const dual_view_type& view)
    : base_type (map, view)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const dual_view_type& view,
          const dual_view_type& origView)
    : base_type (map, view, origView)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  ~Vector ()
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  replaceGlobalValue (const GlobalOrdinal globalRow, const Scalar& value) const {
    this->base_type::replaceGlobalValue (globalRow, 0, value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  sumIntoGlobalValue (const GlobalOrdinal globalRow,
                      const Scalar& value,
                      const bool atomic) const
  {
    this->base_type::sumIntoGlobalValue (globalRow, 0, value, atomic);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  replaceLocalValue (const LocalOrdinal myRow, const Scalar& value) const {
    this->base_type::replaceLocalValue (myRow, 0, value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  sumIntoLocalValue (const LocalOrdinal globalRow,
                     const Scalar& value,
                     const bool atomic) const
  {
    this->base_type::sumIntoLocalValue (globalRow, 0, value, atomic);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  get1dCopy (const Teuchos::ArrayView<Scalar>& A) const {
    const size_t lda = this->getLocalLength ();
    this->get1dCopy (A, lda);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  typename Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::dot_type
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  dot (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>& y) const
  {
    dot_type result;
    this->dot (y, Teuchos::arrayView (&result, 1));
    return result;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Scalar
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  meanValue () const
  {
    Scalar mean;
    this->meanValue (Teuchos::arrayView (&mean, 1));
    return mean;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  typename Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::mag_type
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  norm1 () const
  {
    mag_type norm;
    this->norm1 (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  typename Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::mag_type
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  norm2 () const
  {
    mag_type norm;
    this->norm2 (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  typename Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::mag_type
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  normInf () const
  {
    mag_type norm;
    this->normInf (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  typename Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::mag_type
  TPETRA_DEPRECATED
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  normWeighted (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>& weights) const
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
      reduceAll<int, mag_type> (*comm, REDUCE_SUM, 1, lclNrm.ptr_on_device (),
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  std::string Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  description () const
  {
    using Teuchos::TypeNameTraits;

    std::ostringstream out;
    out << "\"Tpetra::Vector\": {";
    out << "Template parameters: {Scalar: " << TypeNameTraits<Scalar>::name ()
        << ", LocalOrdinal: " << TypeNameTraits<LocalOrdinal>::name ()
        << ", GlobalOrdinal: " << TypeNameTraits<GlobalOrdinal>::name ()
        << ", Node" << Node::name ()
        << "}, ";
    if (this->getObjectLabel () != "") {
      out << "Label: \"" << this->getObjectLabel () << "\", ";
    }
    out << "Global length: " << this->getGlobalLength ();
    out << "}";

    return out.str ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> V;

    const Teuchos::EVerbosityLevel vl =
      (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;
    const Teuchos::Comm<int>& comm = * (this->getMap ()->getComm ());
    const int myImageID = comm.getRank ();
    const int numImages = comm.getSize ();

    size_t width = 1;
    for (size_t dec=10; dec<this->getGlobalLength(); dec *= 10) {
      ++width;
    }
    Teuchos::OSTab tab(out);
    if (vl != VERB_NONE) {
      // VERB_LOW and higher prints description()
      if (myImageID == 0) out << this->description() << std::endl;
      for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
        if (myImageID == imageCtr) {
          if (vl != VERB_LOW) {
            // VERB_MEDIUM and higher prints getLocalLength()
            out << "Process " << setw(width) << myImageID << ":" << endl;
            Teuchos::OSTab tab1 (out);
            const size_t lclNumRows = this->getLocalLength ();

            out << "Local length: " << lclNumRows << endl;
            if (vl != VERB_MEDIUM) {
              // VERB_HIGH and higher prints isConstantStride() and stride()
              if (vl == VERB_EXTREME && lclNumRows > 0) {
                // VERB_EXTREME prints values

                // We have to be able to access the data on host in
                // order to print it.
                //
                // FIXME (mfh 06 Mar 2015) For now, just sync to host.
                // At some point, we might like to check whether the
                // host execution space can access device memory, so
                // that we can avoid the sync
                const_cast<V*> (this)->template sync<Kokkos::HostSpace> ();
                auto X_host = this->template getLocalView<Kokkos::HostSpace> ();
                for (size_t i = 0; i < lclNumRows; ++i) {
                  out << setw(width) << this->getMap ()->getGlobalElement (i)
                      << ": " << X_host(i,0) << endl;
                }
              }
            }
            else {
              out << endl;
            }
          }
        }
        comm.barrier ();
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>
  createCopy (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>& src)
  {
    // The 2-argument copy constructor with second argument =
    // Teuchos::Copy does a deep copy of its input.
    Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> dst (src, Teuchos::Copy);

    // The Kokkos refactor version of Vector has view semantics, so
    // returning the Vector directly, rather than through RCP, only
    // does a shallow copy.
    return dst;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> >
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  offsetView (const Teuchos::RCP<const map_type>& subMap,
              const size_t offset) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::rcp;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> V;

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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> >
  Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>::
  offsetViewNonConst (const Teuchos::RCP<const map_type>& subMap,
                      const size_t offset)
  {
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic> V;
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
