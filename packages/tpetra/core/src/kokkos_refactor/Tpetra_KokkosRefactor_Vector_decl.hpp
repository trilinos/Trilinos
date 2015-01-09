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

#ifndef TPETRA_KOKKOS_REFACTOR_VECTOR_DECL_HPP
#define TPETRA_KOKKOS_REFACTOR_VECTOR_DECL_HPP

#include "Tpetra_ConfigDefs.hpp"
#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

namespace Tpetra {

/// \brief Partial specialization of Vector for Kokkos refactor Node types.
///
/// \tparam Scalar The type of each entry of the vector.  (You can use
///  real-valued or complex-valued types here, unlike in Epetra, where
///  the scalar type is always \c double.)
/// \tparam LocalOrdinal The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GlobalOrdinal The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam DeviceType The Kokkos execution space type.
///
/// This class inherits from MultiVector, and has the same template
/// parameters.  A Vector is a special case of a MultiVector that has
/// only one vector (column).  It may be used wherever a MultiVector
/// may be used.  Please see the documentation of MultiVector for more
/// details.
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class DeviceType>
class Vector<Scalar,
             LocalOrdinal,
             GlobalOrdinal,
             Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>,
             false> :
   public MultiVector<Scalar,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
{
private:
  // need this so that MultiVector::operator() can call Vector's
  // private view constructor
  friend class MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >;
  typedef MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > base_type;

public:
  //! \name Typedefs to facilitate template metaprogramming
  //@{

  /// \brief This class' first template parameter; the type of each
  ///   entry in the Vector.
  typedef Scalar scalar_type;
  /// \brief The type used internally in place of \c Scalar.
  ///
  /// Some \c Scalar types might not work with Kokkos on all execution
  /// spaces, due to missing CUDA device macros or volatile overloads.
  /// The C++ standard type std::complex<T> has this problem.  To fix
  /// this, we replace std::complex<T> values internally with the
  /// (usually) bitwise identical type Kokkos::complex<T>.  The latter
  /// is the \c impl_scalar_type corresponding to \c Scalar =
  /// std::complex.
  typedef typename base_type::impl_scalar_type impl_scalar_type;
  //! This class' second template parameter; the type of local indices.
  typedef LocalOrdinal local_ordinal_type;
  //! This class' third template parameter; the type of global indices.
  typedef GlobalOrdinal global_ordinal_type;
  //! The Kokkos Node type.
  typedef typename base_type::node_type node_type;

  /// \brief Type of an inner ("dot") product result.
  ///
  /// This is usually the same as <tt>impl_scalar_type</tt>, but may
  /// differ if <tt>impl_scalar_type</tt> is e.g., an uncertainty
  /// quantification type from the Stokhos package.
  typedef typename base_type::dot_type dot_type;

  /// \brief Type of a norm result.
  ///
  /// This is usually the same as the type of the magnitude (absolute
  /// value) of <tt>impl_scalar_type</tt>, but may differ if
  /// <tt>impl_scalar_type</tt> is e.g., an uncertainty quantification
  /// type from the Stokhos package.
  typedef typename base_type::mag_type mag_type;

  //! Kokkos::DualView specialization used by this class.
  typedef typename base_type::dual_view_type dual_view_type;

  //! The type of the Map specialization used by this class.
  typedef typename base_type::map_type map_type;

  //@}
  //! \name Constructors and destructor
  //@{

  /// \brief Basic constructor.
  ///
  /// \param map [in] The Vector's Map.  The Map describes the
  ///   distribution of rows over process(es) in the Map's
  ///   communicator.
  //
  /// \param zeroOut [in] If true (the default), require that all the
  ///   Vector's entries be zero on return.  If false, the Vector's
  ///   entries have undefined values on return, and must be set
  ///   explicitly.
  explicit Vector (const Teuchos::RCP<const map_type>& map,
                   const bool zeroOut = true);

  /// \brief Copy constructor (always a shallow copy).
  ///
  /// In this, the Kokkos refactor version of Tpetra, the "copy
  /// constructor" does a shallow copy.  Use the nonmember function
  /// deep_copy() to do a deep copy from one existing Vector to
  /// another, and use the two-argument copy constructor below (with
  /// copyOrView=Teuchos::Copy) to create a Vector which is a deep
  /// copy of an existing Vector.
  Vector (const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>& source);

  /// \brief Copy constructor (shallow or deep copy).
  ///
  /// \param source [in] The Vector to copy.
  ///
  /// \param copyOrView [in] If Teuchos::View, return a shallow copy
  ///   (a view) of \c source.  If Teuchos::Copy, return a deep copy
  ///   of \c source.  Regardless, the result has "view semantics."
  ///   This means that copy construction or assignment (operator=)
  ///   with the resulting object will always do a shallow copy, and
  ///   will transmit view semantics to the result of the shallow
  ///   copy.
  Vector (const Vector<Scalar, LocalOrdinal, GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>& source,
          const Teuchos::DataAccess copyOrView);

  //! \brief Set vector values from an existing array (copy)
  Vector (const Teuchos::RCP<const map_type>& map,
          const Teuchos::ArrayView<const Scalar>& A);

  /// \brief Expert mode constructor, that takes a Kokkos::DualView of
  ///   the Vector's data, and returns a Vector that views those data.
  ///
  /// \warning This constructor is only for expert users.  We make
  ///   no promises about backwards compatibility for this interface.
  ///   It may change or go away at any time.
  ///
  /// See the documentation of the MultiVector (parent class)
  /// constructor that takes the same arguments.
  ///
  /// \param map [in] Map describing the distribution of rows.
  /// \param view [in] View of the data (shallow copy).
  Vector (const Teuchos::RCP<const map_type>& map,
          const dual_view_type& view);

  /// \brief Expert mode constructor, that takes a Kokkos::DualView of
  ///   the Vector's data and the "original" Kokkos::DualView of the
  ///   data, and returns a Vector that views those data.
  ///
  /// \warning This constructor is only for expert users.  We make
  ///   no promises about backwards compatibility for this interface.
  ///   It may change or go away at any time.
  ///
  /// See the documentation of the MultiVector (parent class)
  /// constructor that takes the same arguments.
  ///
  /// \param map [in] Map describing the distribution of rows.
  /// \param view [in] View of the data (shallow copy).
  /// \param origView [in] "Original" view of the data (shallow copy).
  Vector (const Teuchos::RCP<const map_type>& map,
          const dual_view_type& view,
          const dual_view_type& origView);

  //! Destructor.
  virtual ~Vector ();

  //@}
  //! \name Clone method
  //@{

  /// \brief Return a deep copy of <tt>*this</tt> with a different
  ///   Node type (and therefore a different Device type).
  /// \tparam Node2 The returned Vector's Node type.
  ///
  /// \param node2 [in] The returned Vector's Kokkos Node instance.
  template <class Node2>
  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node2, Node2::classic> >
  clone (const Teuchos::RCP<Node2>& node2);

  //@}
  //! @name Post-construction modification routines
  //@{

  //! Replace current value at the specified location with specified value.
  /** \pre \c globalRow must be a valid global element on this node, according to the row map.
   */
  void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value);

  //! Adds specified value to existing value at the specified location.
  /** \pre \c globalRow must be a valid global element on this node, according to the row map.
   */
  void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value);

  //! Replace current value at the specified location with specified values.
  /** \pre \c localRow must be a valid local element on this node, according to the row map.
   */
  void replaceLocalValue(LocalOrdinal myRow, const Scalar &value);

  //! Adds specified value to existing value at the specified location.
  /** \pre \c localRow must be a valid local element on this node, according to the row map.
   */
  void sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value);

  //@}

  //! @name Extraction methods
  //@{

  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>::get1dCopy; // overloading, not hiding
  //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
  void get1dCopy (const Teuchos::ArrayView<Scalar>& A) const;

  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>::getDataNonConst; // overloading, not hiding
  //! View of the local values of this vector.
  Teuchos::ArrayRCP<Scalar> getDataNonConst () { return getDataNonConst (0); }

  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>::getData; // overloading, not hiding
  //! Const view of the local values of this vector.
  Teuchos::ArrayRCP<const Scalar> getData () const { return getData (0); }

  Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false> >
  offsetView (const Teuchos::RCP<const map_type>& subMap,
              const size_t offset) const;

  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false> >
  offsetViewNonConst (const Teuchos::RCP<const map_type>& subMap,
                      const size_t offset);
  //@}
  //! @name Mathematical methods
  //@{

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type>::dot; // overloading, not hiding
  //! Computes dot product of this Vector against input Vector x.
  dot_type dot (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>& y) const;

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type,false>::norm1; // overloading, not hiding
  //! Return 1-norm of this Vector.
  mag_type norm1() const;

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type,false>::norm2; // overloading, not hiding
  //! Compute 2-norm of this Vector.
  mag_type norm2() const;

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type,false>::normInf; // overloading, not hiding
  //! Compute Inf-norm of this Vector.
  mag_type normInf() const;

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type,false>::normWeighted; // overloading, not hiding
  //! Compute Weighted 2-norm (RMS Norm) of this Vector.
  mag_type
  normWeighted (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>& weights) const;

  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type, false>::meanValue; // overloading, not hiding
  //! Compute mean (average) value of this Vector.
  Scalar meanValue () const;

  //@}
  //! @name Implementation of the Teuchos::Describable interface
  //@{

  //! A simple one-line description of this object.
  virtual std::string description() const;

  //! Print the object with some verbosity level to a FancyOStream.
  virtual void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;
  //@}

protected:
  template <class S,class LO,class GO,class N>
  friend Teuchos::RCP< Vector<S,LO,GO,N> >
  createVectorFromView(const RCP<const Map<LO,GO,N> > &,const ArrayRCP<S> &);

  typedef KokkosClassic::MultiVector<Scalar,node_type> KMV;
  typedef KokkosClassic::DefaultArithmetic<KMV>   MVT;
}; // class Vector


/// \brief Nonmember function to create a Vector with view semantics
///   using user-allocated data.
///
/// \warning This function is DEPRECATED and exists only for backwards
///   compatibility with the "classic" version of Tpetra.
/// \relatesalso Vector
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
createVectorFromView (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &map,
                      const ArrayRCP<Scalar> &view)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error, "Tpetra::createVectorFromView: "
    "Not implemented for Node = KokkosDeviceWrapperNode");
}

/// \brief Return a deep copy of the Vector \c src.
///
/// This is the preferred way to make a deep copy of a Vector.  The
/// new Kokkos refactor implementations of Tpetra objects have view
/// semantics, which means that the copy constructor and assignment
/// operator (operator=) make shallow copies.  To make a deep copy,
/// call the nonmember function createCopy().
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
createCopy (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& src);

} // namespace Tpetra

#endif // TPETRA_KOKKOS_REFACTOR_VECTOR_DECL_HPP
