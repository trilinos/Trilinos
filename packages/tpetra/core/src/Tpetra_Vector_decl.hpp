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

#ifndef TPETRA_VECTOR_DECL_HPP
#define TPETRA_VECTOR_DECL_HPP

/// \file Tpetra_Vector_decl.hpp
/// \brief Declaration of the Tpetra::Vector class
///
/// If you want to use Tpetra::Vector, include "Tpetra_Vector.hpp" (a
/// file which CMake generates and installs for you).  If you only
/// want the declaration of Tpetra::Vector, include this file
/// (Tpetra_Vector_decl.hpp).

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_MultiVector_decl.hpp"

namespace Tpetra {

/// \class Vector
/// \brief A distributed dense vector.
///
/// \tparam Scalar The type of each entry of the vector.  (You can use
///  real-valued or complex-valued types here, unlike in Epetra, where
///  the scalar type is always \c double.)
/// \tparam LocalOrdinal The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GlobalOrdinal The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam Node The Kokkos Node type.  See the documentation of Map
///   for requirements.
///
/// This class inherits from MultiVector, and has the same template
/// parameters.  A Vector is a special case of a MultiVector that has
/// only one vector (column).  It may be used wherever a MultiVector
/// may be used.  Please see the documentation of MultiVector for more
/// details.
template <class Scalar = ::Tpetra::Details::DefaultTypes::scalar_type,
          class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
          class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
          class Node = ::Tpetra::Details::DefaultTypes::node_type>
class Vector : public MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{
private:
  friend class MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> base_type;

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
  //! The Kokkos device type.
  typedef typename Node::device_type device_type;

  //! The Kokkos Node type.
  typedef Node node_type;

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

  //! Default constructor: makes a Vector with no rows or columns.
  Vector ();

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
  Vector (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source);

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
  Vector (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source,
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

  /// \brief Create a Vector that views a single column of the input
  ///   MultiVector.
  ///
  /// \param X [in] Input MultiVector to view (in possibly nonconst fashion).
  /// \param j [in] The column of X to view.
  Vector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
          const size_t j);

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
  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node2> >
  clone (const Teuchos::RCP<Node2>& node2);

  //@}
  //! @name Post-construction modification routines
  //@{

  //! Replace current value at the specified location with specified value.
  /** \pre \c globalRow must be a valid global element on this node, according to the row map.
   */
  void replaceGlobalValue (const GlobalOrdinal globalRow, const Scalar& value) const;

  /// \brief Add value to existing value, using global (row) index.
  ///
  /// Add the given value to the existing value at row \c globalRow
  /// (a global index).
  ///
  /// This method affects the host memory version of the data.  If the
  /// \c DeviceType template parameter is a device that has two memory
  /// spaces, and you want to modify the non-host version of the data,
  /// you must access the device View directly by calling
  /// getLocalView().  Please see modify(), sync(), and the discussion
  /// of DualView semantics elsewhere in the documentation.
  ///
  /// \param globalRow [in] Global row index of the entry to modify.
  ///   This <i>must</i> be a valid global row index on the calling
  ///   process with respect to the Vector's Map.
  /// \param value [in] Incoming value to add to the entry.
  /// \param atomic [in] Whether to use an atomic update.  If this
  ///   class' execution space is not Kokkos::Serial, then this is
  ///   true by default, else it is false by default.
  void
  sumIntoGlobalValue (const GlobalOrdinal globalRow,
                      const Scalar& value,
                      const bool atomic = base_type::useAtomicUpdatesByDefault) const;

  //! Replace current value at the specified location with specified values.
  /** \pre \c localRow must be a valid local element on this node, according to the row map.
   */
  void replaceLocalValue (const LocalOrdinal myRow, const Scalar& value) const;

  /// \brief Add \c value to existing value, using local (row) index.
  ///
  /// Add the given value to the existing value at row \c localRow (a
  /// local index).
  ///
  /// This method affects the host memory version of the data.  If the
  /// \c DeviceType template parameter is a device that has two memory
  /// spaces, and you want to modify the non-host version of the data,
  /// you must access the device View directly by calling
  /// getLocalView().  Please see modify(), sync(), and the discussion
  /// of DualView semantics elsewhere in the documentation.
  ///
  /// \param localRow [in] Local row index of the entry to modify.
  /// \param value [in] Incoming value to add to the entry.
  /// \param atomic [in] Whether to use an atomic update.  If this
  ///   class' execution space is not Kokkos::Serial, then this is
  ///   true by default, else it is false by default.
  void
  sumIntoLocalValue (const LocalOrdinal myRow,
                     const Scalar& value,
                     const bool atomic = base_type::useAtomicUpdatesByDefault) const;

  //@}

  //! @name Extraction methods
  //@{

  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get1dCopy; // overloading, not hiding
  //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
  void get1dCopy (const Teuchos::ArrayView<Scalar>& A) const;

  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDataNonConst; // overloading, not hiding
  //! View of the local values of this vector.
  Teuchos::ArrayRCP<Scalar> getDataNonConst () {
    return getDataNonConst (0);
  }

  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getData; // overloading, not hiding
  //! Const view of the local values of this vector.
  Teuchos::ArrayRCP<const Scalar> getData () const {
    return getData (0);
  }

  Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  offsetView (const Teuchos::RCP<const map_type>& subMap,
              const size_t offset) const;

  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  offsetViewNonConst (const Teuchos::RCP<const map_type>& subMap,
                      const size_t offset);
  //@}
  //! @name Mathematical methods
  //@{

  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dot; // overloading, not hiding
  //! Return the dot product of this Vector and the input Vector x.
  dot_type dot (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y) const;

  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::norm1; // overloading, not hiding
  //! Return the one-norm of this Vector.
  mag_type norm1() const;

  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::norm2; // overloading, not hiding
  //! Return the two-norm of this Vector.
  mag_type norm2() const;

  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::normInf; // overloading, not hiding
  //! Return the infinity-norm of this Vector.
  mag_type normInf() const;

  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::normWeighted; // overloading, not hiding
  /// \brief Compute Weighted 2-norm (RMS Norm) of this Vector.
  ///
  /// \warning This method is DEPRECATED.
  mag_type TPETRA_DEPRECATED
  normWeighted (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& weights) const;

  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::meanValue; // overloading, not hiding
  //! Compute mean (average) value of this Vector.
  Scalar meanValue () const;

  //@}
  //! @name Implementation of the Teuchos::Describable interface
  //@{

  //! Return a one-line description of this object.
  virtual std::string description () const;

  /// \brief Describe this object in a human-readable way to the
  ///   given output stream.
  ///
  /// You must call this method as a collective over all processes
  /// in this object's communicator.
  ///
  /// \param out [out] Output stream to which to write.  Only
  ///   Process 0 in this object's communicator may write to the
  ///   output stream.
  ///
  /// \param verbLevel [in] Verbosity level.  This also controls
  ///   whether this method does any communication.  At verbosity
  ///   levels higher (greater) than Teuchos::VERB_LOW, this method
  ///   may behave as a collective over the object's communicator.
  ///
  /// Teuchos::FancyOStream wraps std::ostream.  It adds features
  /// like tab levels.  If you just want to wrap std::cout, try
  /// this:
  /// \code
  /// auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::out));
  /// \endcode
  virtual void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const;
  //@}
}; // class Vector


/// \brief Return a deep copy of the given Vector.
/// \relatesalso Vector
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>
createCopy (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& src);

/// \brief Nonmember Vector "constructor": Create a Vector from a given Map.
/// \relatesalso Vector
///
/// \param map [in] Map describing the distribution of rows of the
///   resulting Vector.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
createVector (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& map)
{
  return rcp (new Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map));
}

} // namespace Tpetra

#endif // TPETRA_VECTOR_DECL_HPP
