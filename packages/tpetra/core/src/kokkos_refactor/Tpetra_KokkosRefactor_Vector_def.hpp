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

#ifndef TPETRA_KOKKOS_REFACTOR_VECTOR_DEF_HPP
#define TPETRA_KOKKOS_REFACTOR_VECTOR_DEF_HPP

#include <Kokkos_DefaultArithmetic.hpp>
#include <Kokkos_NodeTrace.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_Vector_decl.hpp"
#endif
#include <KokkosCompat_View.hpp>
#include <Kokkos_MV.hpp>
#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const bool zeroOut)
    : base_type (map, 1, zeroOut)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  Vector (const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>& source)
    : base_type (source)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  Vector (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>& source,
          const Teuchos::DataAccess copyOrView)
    : base_type (source, copyOrView)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const Teuchos::ArrayView<const Scalar>& values)
    : base_type (map, values, values.size (), 1)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const dual_view_type& view)
    : base_type (map, view)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const dual_view_type& view,
          const dual_view_type& origView)
    : base_type (map, view, origView)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  ~Vector ()
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  replaceGlobalValue (GlobalOrdinal globalRow, const Scalar& value) {
    this->base_type::replaceGlobalValue (globalRow, 0, value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  sumIntoGlobalValue (GlobalOrdinal globalRow, const Scalar& value) {
    this->base_type::sumIntoGlobalValue (globalRow, 0, value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  replaceLocalValue (LocalOrdinal myRow, const Scalar& value) {
    this->base_type::replaceLocalValue (myRow, 0, value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  sumIntoLocalValue (LocalOrdinal myRow, const Scalar& value) {
    this->base_type::sumIntoLocalValue (myRow, 0, value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  get1dCopy (const Teuchos::ArrayView<Scalar>& A) const {
    const size_t lda = this->getLocalLength ();
    this->get1dCopy (A, lda);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  typename Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::dot_type
  Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  dot (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>& y) const
  {
    using Kokkos::subview;
    using Teuchos::outArg;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;

    const size_t lclNumRows = this->getLocalLength ();
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! this->getMap ()->isCompatible (* (y.getMap ())), std::runtime_error,
      "Tpetra::Vector::dots: Vectors do not have compatible Maps:" << std::endl
      << "this->getMap(): " << std::endl << * (this->getMap ())
      << "y.getMap(): " << std::endl << * (y.getMap ()) << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION(
      lclNumRows != y.getLocalLength (), std::runtime_error,
      "Tpetra::Vector::dots: Vectors do not have the same local length.");
#endif // HAVE_TPETRA_DEBUG

    typedef Kokkos::View<impl_scalar_type*,
      typename dual_view_type::array_layout,
      DeviceType> col_view_type;

    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    // FIXME (mfh 04 Jan 2015) We should actually prefer executing
    // this kernel where the data were most recently updated.  The
    // problem is that *this and y may have been most recently updated
    // in different memory spaces.  What should we do then?  Perhaps
    // the best option would be to allocate a temporary in the
    // "capacity" memory space and do the work there.  Or, we could
    // ignore the "const" appellation and just sync *this.
    col_view_type X_0 = subview<col_view_type> (this->view_.d_view, rowRng, 0);
    col_view_type Y_0 = subview<col_view_type> (y.view_.d_view, rowRng, 0);

    dot_type gbldot = Kokkos::V_Dot (X_0, Y_0, lclNumRows);
    if (this->isDistributed ()) {
      Teuchos::RCP<const map_type> map = this->getMap ();
      Teuchos::RCP<const Teuchos::Comm<int> > comm =
        map.is_null () ? Teuchos::null : map->getComm ();
      if (! comm.is_null ()) {
        dot_type lcldot = gbldot;
        reduceAll<int, dot_type> (*comm, REDUCE_SUM, lcldot, outArg (gbldot));
      }
    }
    return gbldot;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Scalar
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  meanValue () const
  {
    Scalar mean;
    this->meanValue (Teuchos::arrayView (&mean, 1));
    return mean;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  typename Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::mag_type
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  norm1 () const
  {
    mag_type norm;
    this->norm1 (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  typename Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::mag_type
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  norm2 () const
  {
    mag_type norm;
    this->norm2 (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  typename Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::mag_type
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  normInf () const
  {
    mag_type norm;
    this->normInf (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  typename Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::mag_type
  Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  normWeighted (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>& weights) const
  {
    mag_type norm;
    this->normWeighted (weights, Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  std::string Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  description () const
  {
    using Teuchos::TypeNameTraits;

    std::ostringstream out;
    out << "\"Tpetra::Vector\": {";
    out << "Template parameters: {Scalar: " << TypeNameTraits<Scalar>::name ()
        << ", LocalOrdinal: " << TypeNameTraits<LocalOrdinal>::name ()
        << ", GlobalOrdinal: " << TypeNameTraits<GlobalOrdinal>::name ()
        << ", Node" << Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>::name ()
        << "}, ";
    if (this->getObjectLabel () != "") {
      out << "Label: \"" << this->getObjectLabel () << "\", ";
    }
    out << "Global length: " << this->getGlobalLength ();
    out << "}";

    return out.str ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  describe (Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
  {
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;

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
            out << "node " << setw(width) << myImageID << ": local length=" << this->getLocalLength() << endl;
            if (vl != VERB_MEDIUM) {
              // VERB_HIGH and higher prints isConstantStride() and stride()
              if (vl == VERB_EXTREME && this->getLocalLength() > 0) {
                /*RCP<Node> node = this->lclMV_.getNode();
                KOKKOS_NODE_TRACE("Vector::describe()")
                ArrayRCP<const Scalar> myview = node->template viewBuffer<Scalar>(
                                                               this->getLocalLength(),
                                                               MVT::getValues(this->lclMV_) );
                // VERB_EXTREME prints values
                for (size_t i=0; i<this->getLocalLength(); ++i) {
                  out << setw(width) << this->getMap()->getGlobalElement(i)
                      << ": "
                      << myview[i] << endl;
                }
                myview = Teuchos::null;*/
              }
            }
            else {
              out << endl;
            }
          }
        }
        comm.barrier();
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Vector<Scalar, LocalOrdinal, GlobalOrdinal,
         Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
  createCopy (const Vector<Scalar, LocalOrdinal, GlobalOrdinal,
                           Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& src)
  {
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> NT;

    // The 2-argument copy constructor with second argument =
    // Teuchos::Copy does a deep copy of its input.
    Vector<Scalar, LocalOrdinal, GlobalOrdinal, NT> dst (src, Teuchos::Copy);

    // The Kokkos refactor version of Vector has view semantics, so
    // returning the Vector directly, rather than through RCP, only
    // does a shallow copy.
    return dst;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Vector<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false> >
  Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  offsetView (const Teuchos::RCP<const map_type>& subMap,
              const size_t offset) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::rcp;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > V;

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
                       subview<dual_view_type> (this->view_, offsetPair, ALL ()),
                       this->origView_));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<Vector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false> >
  Vector<Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  offsetViewNonConst (const Teuchos::RCP<const map_type>& subMap,
                      const size_t offset)
  {
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > V;
    return Teuchos::rcp_const_cast<V> (this->offsetView (subMap, offset));
  }

} // namespace Tpetra


#endif // TPETRA_VECTOR_DEF_HPP
