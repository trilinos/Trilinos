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

#ifndef TPETRANEW_VECTOR_DECL_HPP
#define TPETRANEW_VECTOR_DECL_HPP

/// \file TpetraNew_Vector_decl.hpp
/// \brief Declaration of the TpetraNew::Vector class
///
/// If you want to use TpetraNew::Vector, include
/// "TpetraNew_Vector.hpp" (a file which CMake generates and installs
/// for you).  If you only want the declaration of TpetraNew::Vector,
/// include this file (TpetraNew_Vector_decl.hpp).

#include "TpetraNew_Vector_fwd.hpp"
#include "TpetraNew_MultiVector_decl.hpp"

namespace TpetraNew {

/// \class Vector
/// \brief A distributed dense vector.
///
/// \tparam Scalar The type of each entry of the vector.  (You can use
///  real-valued or complex-valued types here, unlike in Epetra, where
///  the scalar type is always \c double.)
///
/// This class inherits from MultiVector, and has the same template
/// parameters.  A Vector is a special case of a MultiVector that has
/// only one vector (column).  It may be used wherever a MultiVector
/// may be used.  Please see the documentation of MultiVector for more
/// details.
template <class Scalar>
class Vector : public MultiVector<Scalar> {
private:
  using base_type = MultiVector<Scalar>;

public:
  //! \name Typedefs to facilitate template metaprogramming
  //@{

  /// \brief This class' template parameter; the type of each entry in
  ///   the Vector.
  using scalar_type = Scalar;
  /// \brief The type used internally in place of \c Scalar.
  ///
  /// Some \c Scalar types might not work with Kokkos on all execution
  /// spaces, due to missing CUDA device macros or volatile overloads.
  /// The C++ standard type <tt>std::complex<T></tt> has this problem.
  /// To fix this, we replace <tt>std::complex<T></tt> values
  /// internally with the bitwise identical type
  /// <tt>Kokkos::complex<T></tt>.  The latter is the
  /// <tt>impl_scalar_type</tt> corresponding to <tt>Scalar =
  /// std::complex<T></tt>.
  using impl_scalar_type = typename base_type::impl_scalar_type;
  //! The type of local indices.
  using local_ordinal_type =
    ::Tpetra::Details::DefaultTypes::local_ordinal_type;
  //! The type of global indices.
  using global_ordinal_type =
    ::Tpetra::Details::DefaultTypes::global_ordinal_type;
  //! The Kokkos execution space.
  using execution_space =
    ::Tpetra::Details::DefaultTypes::execution_space;
  //! The Kokkos memory space.
  using memory_space = ::Tpetra::Details::DefaultTypes::memory_space;
  //! The Kokkos device type.
  using device_type = ::Tpetra::Details::DefaultTypes::device_type;

  /// \brief Type of an inner ("dot") product result.
  ///
  /// This is usually the same as <tt>impl_scalar_type</tt>, but may
  /// differ if <tt>impl_scalar_type</tt> is e.g., an uncertainty
  /// quantification type from the Stokhos package.
  using dot_type = typename base_type::dot_type;

  /// \brief Type of a norm result.
  ///
  /// This is usually the same as the type of the magnitude (absolute
  /// value) of <tt>impl_scalar_type</tt>, but may differ if
  /// <tt>impl_scalar_type</tt> is e.g., an uncertainty quantification
  /// type from the Stokhos package.
  using mag_type = typename base_type::mag_type;

  //! Kokkos::DualView specialization used by this class.
  using dual_view_type = typename base_type::dual_view_type;

  //! The type of the Map specialization used by this class.
  using map_type = typename base_type::map_type;

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
  Vector (const Vector<Scalar>& source);

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
  Vector (const Vector<Scalar>& source,
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
  Vector (const MultiVector<Scalar>& X,
          const size_t j);

  //! Destructor.
  virtual ~Vector ();

  //@}
  //! @name Post-construction modification routines
  //@{

  /// \brief Replace current value at the specified location with
  ///   specified value.
  ///
  /// \pre <tt>this->getMap()->isMyGlobalIndex(globalRow)</tt> is
  ///   <tt>true</tt>.
  void
  replaceGlobalValue (const global_ordinal_type globalRow,
		      const Scalar& value) const;

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
  sumIntoGlobalValue (const global_ordinal_type globalRow,
                      const Scalar& value,
                      const bool atomic = base_type::useAtomicUpdatesByDefault) const;

  /// \brief Replace current value at the specified location with
  ///   specified value.
  ///
  /// \pre <tt>this->getMap()->isMyLocalIndex(myRow)</tt> is
  ///   <tt>true</tt>.
  void
  replaceLocalValue (const local_ordinal_type myRow,
		     const Scalar& value) const;

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
  sumIntoLocalValue (const local_ordinal_type myRow,
                     const Scalar& value,
                     const bool atomic = base_type::useAtomicUpdatesByDefault) const;

  //@}

  //! @name Extraction methods
  //@{

  using MultiVector<Scalar>::get1dCopy; // overloading, not hiding
  //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
  void get1dCopy (const Teuchos::ArrayView<Scalar>& A) const;

  using MultiVector<Scalar>::getDataNonConst; // overloading, not hiding
  //! View of the local values of this vector.
  Teuchos::ArrayRCP<Scalar> getDataNonConst () {
    return getDataNonConst (0);
  }

  using MultiVector<Scalar>::getData; // overloading, not hiding
  //! Const view of the local values of this vector.
  Teuchos::ArrayRCP<const Scalar> getData () const {
    return getData (0);
  }

  Teuchos::RCP<const Vector<Scalar> >
  offsetView (const Teuchos::RCP<const map_type>& subMap,
              const size_t offset) const;

  Teuchos::RCP<Vector<Scalar> >
  offsetViewNonConst (const Teuchos::RCP<const map_type>& subMap,
                      const size_t offset);
  //@}
  //! @name Mathematical methods
  //@{

  using MultiVector<Scalar>::dot; // overloading, not hiding
  //! Return the dot product of this Vector and the input Vector x.
  dot_type dot (const Vector<Scalar>& y) const;

  using MultiVector<Scalar>::norm1; // overloading, not hiding
  //! Return the one-norm of this Vector.
  mag_type norm1() const;

  using MultiVector<Scalar>::norm2; // overloading, not hiding
  //! Return the two-norm of this Vector.
  mag_type norm2() const;

  using MultiVector<Scalar>::normInf; // overloading, not hiding
  //! Return the infinity-norm of this Vector.
  mag_type normInf() const;

  using MultiVector<Scalar>::meanValue; // overloading, not hiding
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
template<class Scalar>
Vector<Scalar>
createCopy (const Vector<Scalar>& src);

/// \brief Nonmember Vector "constructor": Create a Vector from a given Map.
/// \relatesalso Vector
///
/// \param map [in] Map describing the distribution of rows of the
///   resulting Vector.
template <class Scalar>
Teuchos::RCP<Vector<Scalar> >
createVector (const Teuchos::RCP<const Map>& map)
{
  return Teuchos::rcp (new Vector<Scalar> (map));
}

} // namespace TpetraNEW

#endif // TPETRANEW_VECTOR_DECL_HPP
