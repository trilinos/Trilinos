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

#ifndef TPETRA_FEMULTIVECTOR_DECL_HPP
#define TPETRA_FEMULTIVECTOR_DECL_HPP

/// \file Tpetra_FEMultiVector_decl.hpp
/// \brief Declaration of the Tpetra::MultiVector class
///
/// If you want to use Tpetra::MultiVector, include "Tpetra_MultiVector.hpp"
/// (a file which CMake generates and installs for you).  If you only want
/// the declaration of Tpetra::MultiVector, include this file
/// (Tpetra_MultiVector_decl.hpp).
///

//#include "Tpetra_DistObject.hpp"
//#include "Tpetra_Map_decl.hpp"
#include "Tpetra_MultiVector_decl.hpp"
//#include "Teuchos_Import_decl.hpp"
//#include "Kokkos_DualView.hpp"
//#include "Teuchos_BLAS_types.hpp"
//#include "Teuchos_DataAccess.hpp"
//#include "Teuchos_Range1D.hpp"

//#include "Kokkos_ArithTraits.hpp"
//#include "Kokkos_InnerProductSpaceTraits.hpp"
//#include "Tpetra_KokkosRefactor_Details_MultiVectorLocalDeepCopy.hpp"
//#include <type_traits>

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of MultiVector (declared later in this file)
  template<class S, class LO, class GO, class N> class FEMultiVector;
#endif // DOXYGEN_SHOULD_SKIP_THIS

  // CMS: Removed all of the non-member stuff... should this be kept?


  template <class Scalar = ::Tpetra::Details::DefaultTypes::scalar_type,
            class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
            class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
            class Node = ::Tpetra::Details::DefaultTypes::node_type>
  class FEMultiVector :
    public MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  {

  public:
    //! @name Typedefs to facilitate template metaprogramming.
    //@{

    //! This class' first template parameter; the Scalar type.
    typedef Scalar scalar_type;
    //! This class' second template parameter; the type of local indices.
    typedef LocalOrdinal local_ordinal_type;
    //! This class' third template parameter; the type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! This class' fourth template parameter; the Kokkos Node type.
    typedef Node node_type;

    //! The dual_view_type picked up from MultiVector
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type dual_view_type;

    //! The type of the Map specialization used by this class.
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::map_type map_type;

    //! Grab impl_scalar_type from superclass
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::impl_scalar_type impl_scalar_type;

    //! Grab dot_type from superclass
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dot_type dot_type;

    //! Grab mag_type from superclass
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type mag_type;

    //! Grab device_type from superclass
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::device_type device_type;

    //! Grab execution_space from superclass
    typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::execution_space execution_space;

    //@}
    //! @name Constructors and destructor
    //@{
    /// \brief Basic constuctor.
    // CMS - A map AND an imported need to be arguments because in serial, the importer will be null
    FEMultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & map,
                  const Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> >& importer,
                  const size_t numVecs,
                  const bool zeroOut = true);

     // CMS - I do not feel that the ArrayView constructors need to be implemented.
     // I am a more open to implementing variants of the dual_view constructors, but are those used?


     // CMS - These are the constructors to MultiVector that I do not want to implement here
   private:

     // WCMCLEN: @CMS: shouldn't this be MultiVector -- not sure how that affects the stuff in _def.hpp
     //! The type of the base class of this class.
     typedef DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> base_type;

     //! Default constructor: makes a MultiVector with no rows or columns.
     FEMultiVector ();

     /// \brief Basic constuctor.
     FEMultiVector (const Teuchos::RCP<const map_type>& map,
                    const size_t numVecs,
                    const bool zeroOut = true);

     /// \brief Create multivector by copying two-dimensional array of local data.
     FEMultiVector (const Teuchos::RCP<const map_type>& map,
                    const Teuchos::ArrayView<const Scalar>& A,
                    const size_t LDA,
                    const size_t NumVectors);

     /// \brief Create multivector by copying array of views of local data.
     FEMultiVector (const Teuchos::RCP<const map_type>& map,
                    const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> >&ArrayOfPtrs,
                    const size_t NumVectors);

     /// \brief Constructor, that takes a Kokkos::DualView of the
     FEMultiVector (const Teuchos::RCP<const map_type>& map,
                    const dual_view_type& view);

     /// \brief Constructor, that takes a Kokkos::View of the
     ///   MultiVector's data (living in the Device's memory space),
     ///   and returns a MultiVector that views those data.
     FEMultiVector (const Teuchos::RCP<const map_type>& map,
                    const typename dual_view_type::t_dev& d_view);

     /// \brief Expert mode constructor, that takes a Kokkos::DualView
     ///   of the MultiVector's data and the "original"
     ///   Kokkos::DualView of the data, and returns a MultiVector that
     ///   views those data.
     FEMultiVector (const Teuchos::RCP<const map_type>& map,
                    const dual_view_type& view,
                    const dual_view_type& origView);

    /// \brief Expert mode constructor for noncontiguous views.
    FEMultiVector (const Teuchos::RCP<const map_type>& map,
                   const dual_view_type& view,
                   const Teuchos::ArrayView<const size_t>& whichVectors);

    /// \brief Expert mode constructor for noncontiguous views, with
    ///   original view.
    FEMultiVector (const Teuchos::RCP<const map_type>& map,
                   const dual_view_type& view,
                   const dual_view_type& origView,
                   const Teuchos::ArrayView<const size_t>& whichVectors);

    /// \brief "Offset view" constructor; make a view of a contiguous
    ///   subset of rows on each process.
    /// WCMCLEN - no equivalent C'tor in MultiVector -- was this deprecated?
    /*
    FEMultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                   const map_type& subMap,
                   const size_t offset = 0);
    */


     /// \brief Return a deep copy of this MultiVector, with a
    ///   different Node type.
    template <class Node2>
    Teuchos::RCP<FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node2> >
    clone (const Teuchos::RCP<Node2>& node2) const;

    //! Shallow copy: assign \c source to \c *this.
    MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
    operator= (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source);

    /// \brief Replace the underlying Map in place.
    // CMS -  We will need some analogue of this guy --- it normally gets used for communicator restriction, which is a use case that matters.
    void replaceMap (const Teuchos::RCP<const map_type>& map);


  public:

     //! Copy constructor (shallow copy!).
     FEMultiVector (const FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source);

     /// \brief Copy constructor, with option to do deep or shallow copy.
     ///
     /// The Kokkos refactor version of Tpetra, unlike the "classic"
     /// version, always has view semantics.  Thus, copyOrView =
     /// Teuchos::View has no effect, and copyOrView = Teuchos::Copy
     /// does not mark this MultiVector as having copy semantics.
     /// However, copyOrView = Teuchos::Copy will make the resulting
     /// MultiVector a deep copy of the input MultiVector.
     FEMultiVector (const FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source,
                  const Teuchos::DataAccess copyOrView);


   protected:

     /// \brief Single-column subview constructor, for derived classes ONLY.
     ///
     /// \param X [in] Input MultiVector to view (in possibly nonconst fashion).
     /// \param j [in] The column of X to view.
     FEMultiVector (const FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                  const size_t j);

   public:
    /// \brief Return a deep copy of this MultiVector, with a
    ///   different Node type.
    // CMS - cloning to a *FEMultiVector* is OK, cloning to a MultiVector is not
    ///
    /// \param node2 [in/out] The new Node type.
    ///
    /// \warning We prefer that you use Tpetra::deep_copy (see below)
    ///   rather than this method.  This method will go away at some
    ///   point.
    template <class Node2>
    Teuchos::RCP<FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node2> >
    clone (const Teuchos::RCP<Node2>& node2) const;

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~FEMultiVector ();

    //@}
    //! @name Post-construction modification routines
    //@{

  public:

    // CMS - FEMultiVector Specific Routines

    enum FEWhichActive
    {
      FE_ACTIVE_TARGET,
      FE_ACTIVE_SOURCE
    };


    //! All routines now return results the Vector corresponding to the Importer's 'source' map
    // WARNING: This routine is not thread-safe.
    void activateSourceMultiVector();
    //! All routines now return results the Vector corresponding to the Importer's 'target' map
    // WARNING: This routine is not thread-safe.
    void activateTargetMultiVector();

    //! Returns true if the source multivector is 'active'
    void isSourceMultiVectorActive() const;
    //! Returns true if the target multivector is 'active'
    void isTargetMultiVectorActive() const;

    //! Returns the source MultiVector
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getSourceMultiVector();
    //! Returns the target MultiVector
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getTargetMultiVector();
    //! Returns the "active" target multivector
    MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & getActiveMultiVector();

    //! Migrate data from the source to the target map
    void doSourceToTarget();

    //! Migrate data from the target to the source map
    // Since this is non-unique -> unique, we need a combine mode.
    // CMS - Should there be a default?
    void doTargetToSource(const CombineMode CM);

    // ! Calls endFill()
    void globalAssemble() {endFill();}

    //! Calls doTargetToSource() and then activateSourceMultiVector()
    //CMS - Should add a sanity check to make sure I start in Target mode (if we're not in serial)

    void endFill() {
      doTargetToSource(Tpetra::ADD);
      activateSourceMultiVector();
    }
    //! Activates the target map mode
    // CMS
    void beginFill() {
      activateTargetMultiVector();
    }


  protected:
    // CMS - The actual data that we'd store
    // We'd rely on "this" to keep the *target* view (aka the bigger one)
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > sourceMultiVector_;
    FEWhichActive activeMultiVector_;

    /// \brief Whether sumIntoLocalValue and sumIntoGlobalValue should
    ///   use atomic updates by default.
    ///
    /// \warning This is an implementation detail.
    static const bool useAtomicUpdatesByDefault =
#ifdef KOKKOS_HAVE_SERIAL
      ! std::is_same<execution_space, Kokkos::Serial>::value;
#else
      true;
#endif // KOKKOS_HAVE_SERIAL


  public:
    /// \brief Replace value in host memory, using global row index.
    ///
    /// Replace the current value at row \c gblRow (a global index)
    /// and column \c col with the given value.  The column index is
    /// zero based.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.  You are responsible for calling modify() and
    /// sync(), if needed; this method doesn't do that.
    ///
    /// This method does not have an "atomic" option like
    /// sumIntoGlobalValue.  This is deliberate.  Replacement is not
    /// commutative, unlike += (modulo rounding error).  Concurrent
    /// calls to replaceGlobalValue on different threads that modify
    /// the same entry/ies have undefined results.  (It's not just
    /// that one thread might win; it's that the value might get
    /// messed up.)
    ///
    /// \param gblRow [in] Global row index of the entry to modify.
    ///   This <i>must</i> be a valid global row index on the calling
    ///   process with respect to the MultiVector's Map.
    /// \param col [in] Column index of the entry to modify.
    /// \param value [in] Incoming value to add to the entry.
    void
    replaceGlobalValue (const GlobalOrdinal gblRow,
                        const size_t col,
                        const impl_scalar_type& value) const;

    /// \brief Like the above replaceGlobalValue, but only enabled if
    ///   T differs from impl_scalar_type.
    ///
    /// This method only exists if its template parameter \c T and
    /// impl_scalar_type differ, and if it is syntactically possible
    /// to convert \c T to impl_scalar_type.  This method is mainly
    /// useful for backwards compatibility, when the Scalar template
    /// parameter differs from impl_scalar_type.  That is commonly
    /// only the case when Scalar is std::complex<U> for some type U.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling getLocalView().
    /// Please see modify(), sync(), and the discussion of DualView
    /// semantics elsewhere in the documentation.  You are responsible
    /// for calling modify() and sync(), if needed; this method
    /// doesn't do that.
    ///
    /// This method does not have an "atomic" option like
    /// sumIntoGlobalValue.  This is deliberate.  Replacement is not
    /// commutative, unlike += (modulo rounding error).  Concurrent
    /// calls to replaceGlobalValue on different threads that modify
    /// the same entry/ies have undefined results.  (It's not just
    /// that one thread might win; it's that the value might get
    /// messed up.)
    ///
    /// \param gblRow [in] Global row index of the entry to modify.
    ///   This <i>must</i> be a valid global row index on the calling
    ///   process with respect to the MultiVector's Map.
    /// \param col [in] Column index of the entry to modify.
    /// \param value [in] Incoming value to add to the entry.
    template<typename T>
    typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
    replaceGlobalValue (GlobalOrdinal globalRow,
                        size_t col,
                        const T& value) const;

    /// \brief Update (+=) a value in host memory, using global row index.
    ///
    /// Add the given value to the existing value at row \c gblRow (a
    /// global index) and column \c col.  The column index is zero
    /// based.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.  You are responsible for calling modify() and
    /// sync(), if needed; this method doesn't do that.
    ///
    /// \param gblRow [in] Global row index of the entry to modify.
    ///   This <i>must</i> be a valid global row index on the calling
    ///   process with respect to the MultiVector's Map.
    /// \param col [in] Column index of the entry to modify.
    /// \param value [in] Incoming value to add to the entry.
    /// \param atomic [in] Whether to use an atomic update.  If this
    ///   class' execution space is not Kokkos::Serial, then this is
    ///   true by default, else it is false by default.
    void
    sumIntoGlobalValue (const GlobalOrdinal gblRow,
                        const size_t col,
                        const impl_scalar_type& value,
                        const bool atomic = useAtomicUpdatesByDefault) const;

    /// \brief Like the above sumIntoGlobalValue, but only enabled if
    ///   T differs from impl_scalar_type.
    ///
    /// This method only exists if its template parameter \c T and
    /// impl_scalar_type differ, and if it is syntactically possible
    /// to convert \c T to impl_scalar_type.  This method is mainly
    /// useful for backwards compatibility, when the Scalar template
    /// parameter differs from impl_scalar_type.  That is commonly
    /// only the case when Scalar is std::complex<U> for some type U.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.  You are responsible for calling modify() and
    /// sync(), if needed; this method doesn't do that.
    ///
    /// \param gblRow [in] Global row index of the entry to modify.
    ///   This <i>must</i> be a valid global row index on the calling
    ///   process with respect to the MultiVector's Map.
    /// \param col [in] Column index of the entry to modify.
    /// \param val [in] Incoming value to add to the entry.
    /// \param atomic [in] Whether to use an atomic update.  If this
    ///   class' execution space is not Kokkos::Serial, then this is
    ///   true by default, else it is false by default.
    template<typename T>
    typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
    sumIntoGlobalValue (const GlobalOrdinal gblRow,
                        const size_t col,
                        const T& val,
                        const bool atomic = useAtomicUpdatesByDefault) const;

    /// \brief Replace value in host memory, using local (row) index.
    ///
    /// Replace the current value at row \c lclRow (a local index) and
    /// column \c col with the given value.  The column index is zero
    /// based.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.  You are responsible for calling modify() and
    /// sync(), if needed; this method doesn't do that.
    ///
    /// This method does not have an "atomic" option like
    /// sumIntoLocalValue.  This is deliberate.  Replacement is not
    /// commutative, unlike += (modulo rounding error).  Concurrent
    /// calls to replaceLocalValue on different threads that modify
    /// the same entry/ies have undefined results.  (It's not just
    /// that one thread might win; it's that the value might get
    /// messed up.)
    ///
    /// \param lclRow [in] Local row index of the entry to modify.
    ///   Must be a valid local index in this MultiVector's Map on the
    ///   calling process.
    /// \param col [in] Column index of the entry to modify.
    /// \param value [in] Incoming value to add to the entry.
    void
    replaceLocalValue (const LocalOrdinal lclRow,
                       const size_t col,
                       const impl_scalar_type& value) const;

    /// \brief Like the above replaceLocalValue, but only enabled if
    ///   T differs from impl_scalar_type.
    ///
    /// This method only exists if its template parameter \c T and
    /// impl_scalar_type differ, and if it is syntactically possible
    /// to convert \c T to impl_scalar_type.  This method is mainly
    /// useful for backwards compatibility, when the Scalar template
    /// parameter differs from impl_scalar_type.  That is commonly
    /// only the case when Scalar is std::complex<U> for some type U.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.  You are responsible for calling modify() and
    /// sync(), if needed; this method doesn't do that.
    ///
    /// This method does not have an "atomic" option like
    /// sumIntoLocalValue.  This is deliberate.  Replacement is not
    /// commutative, unlike += (modulo rounding error).  Concurrent
    /// calls to replaceLocalValue on different threads that modify
    /// the same entry/ies have undefined results.  (It's not just
    /// that one thread might win; it's that the value might get
    /// messed up.)
    ///
    /// \param lclRow [in] Local row index of the entry to modify.
    ///   Must be a valid local index in this MultiVector's Map on the
    ///   calling process.
    /// \param col [in] Column index of the entry to modify.
    /// \param val [in] Incoming value to add to the entry.
    template<typename T>
    typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
    replaceLocalValue (const LocalOrdinal lclRow,
                       const size_t col,
                       const T& val) const;

    /// \brief Update (+=) a value in host memory, using local row index.
    ///
    /// Add the given value to the existing value at row \c localRow
    /// (a local index) and column \c col.  The column index is zero
    /// based.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.  You are responsible for calling modify() and
    /// sync(), if needed; this method doesn't do that.
    ///
    /// \param lclRow [in] Local row index of the entry to modify.
    ///   Must be a valid local index in this MultiVector's Map on the
    ///   calling process.
    /// \param col [in] Column index of the entry to modify.
    /// \param val [in] Incoming value to add to the entry.
    /// \param atomic [in] Whether to use an atomic update.  If this
    ///   class' execution space is not Kokkos::Serial, then this is
    ///   true by default, else it is false by default.
    void
    sumIntoLocalValue (const LocalOrdinal lclRow,
                       const size_t col,
                       const impl_scalar_type& val,
                       const bool atomic = useAtomicUpdatesByDefault) const;

    /// \brief Like the above sumIntoLocalValue, but only enabled if
    ///   T differs from impl_scalar_type.
    ///
    /// This method only exists if its template parameter \c T and
    /// impl_scalar_type differ, and if it is syntactically possible
    /// to convert \c T to impl_scalar_type.  This method is mainly
    /// useful for backwards compatibility, when the Scalar template
    /// parameter differs from impl_scalar_type.  That is commonly
    /// only the case when Scalar is std::complex<U> for some type U.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.  You are responsible for calling modify() and
    /// sync(), if needed; this method doesn't do that.
    ///
    /// \param lclRow [in] Local row index of the entry to modify.
    /// \param col [in] Column index of the entry to modify.
    /// \param val [in] Incoming value to add to the entry.
    /// \param atomic [in] Whether to use an atomic update.  If this
    ///   class' execution space is not Kokkos::Serial, then this is
    ///   true by default, else it is false by default.
    template<typename T>
    typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
    sumIntoLocalValue (const LocalOrdinal lclRow,
                       const size_t col,
                       const T& val,
                       const bool atomic = useAtomicUpdatesByDefault) const;

    //! Set all values in the multivector with the given value.
    void putScalar (const Scalar& value);

    /// \brief Set all values in the multivector with the given value.
    ///
    /// This method only exists if its template parameter \c T and
    /// impl_scalar_type differ, and if it is syntactically possible
    /// to convert \c T to impl_scalar_type.  This method is mainly
    /// useful for backwards compatibility, when the Scalar template
    /// parameter differs from impl_scalar_type.  That is commonly
    /// only the case when Scalar is std::complex<U> for some type U.
    template<typename T>
    typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
    putScalar (const T& value);

    /// \brief Set all values in the multivector to pseudorandom numbers.
    ///
    /// \note Do not expect repeatable results.
    /// \note Behavior of this method may or may not depend on
    ///   external use of the C library routines srand() and rand().
    ///   In particular, setting the seed there may not affect it
    ///   here.
    /// \warning This method does <i>not</i> promise to use a
    ///   distributed-memory parallel pseudorandom number generator.
    ///   Corresponding values on different processes might be
    ///   correlated.  It also does not promise to use a high-quality
    ///   pseudorandom number generator within each process.
    void randomize();

    /// \brief Set all values in the multivector to pseudorandom
    ///   numbers in the given range.
    ///
    /// \note Do not expect repeatable results.
    /// \note Behavior of this method may or may not depend on
    ///   external use of the C library routines srand() and rand().
    ///   In particular, setting the seed there may not affect it
    ///   here.
    /// \warning This method does <i>not</i> promise to use a
    ///   distributed-memory parallel pseudorandom number generator.
    ///   Corresponding values on different processes might be
    ///   correlated.  It also does not promise to use a high-quality
    ///   pseudorandom number generator within each process.
    void randomize (const Scalar& minVal, const Scalar& maxVal);

    /// \brief Replace the underlying Map in place.
    void replaceMap (const Teuchos::RCP<const map_type>& map);

    /// \brief Sum values of a locally replicated multivector across all processes.
    ///
    /// \warning This method may only be called for locally replicated
    ///   MultiVectors.
    ///
    /// \pre isDistributed() == false
    void reduce();

    //! Shallow copy: assign \c source to \c *this.
    FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
    operator= (const FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source);

    //@}

    //! @name Get a copy or view of a subset of rows and/or columns
    ///
    /// The following methods get either a (deep) copy or a view
    /// (shallow copy) of a subset of rows and/or columns of the
    /// MultiVector.  They return one of the following:
    ///
    /// <ul>
    /// <li> Another MultiVector </li>
    /// <li> A Kokkos::View or Kokkos::DualView </li>
    /// <li> A Teuchos::ArrayRCP (see the Teuchos Memory Management Classes) </li>
    /// </ul>
    ///
    /// We prefer use of Kokkos classes to Teuchos Memory Management
    /// Classes.  In particular, Teuchos::ArrayRCP reference counts
    /// are not thread safe, while Kokkos::View (and Kokkos::DualView)
    /// reference counts are thread safe.
    ///
    /// Not all of these methods are valid for a particular
    /// MultiVector. For instance, calling a method that accesses a
    /// view of the data in a 1-D format (i.e., get1dView) requires
    /// that the target MultiVector have constant stride.
    ///
    /// This category of methods also includes sync(), modify(), and
    /// getLocalView(), which help MultiVector implement DualView
    /// semantics.
    //@{

    //! Return a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subCopy (const Teuchos::Range1D& colRng) const;

    //! Return a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subCopy (const Teuchos::ArrayView<const size_t>& cols) const;

    //! Return a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subView (const Teuchos::Range1D& colRng) const;

    //! Return a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subView (const Teuchos::ArrayView<const size_t>& cols) const;

    //! Return a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subViewNonConst (const Teuchos::Range1D& colRng);

    //! Return a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subViewNonConst (const Teuchos::ArrayView<const size_t>& cols);

    /// \brief Return a const view of a subset of rows.
    ///
    /// Return a const (nonmodifiable) view of this MultiVector
    /// consisting of a subset of the rows, as specified by an offset
    /// and a subset Map of this MultiVector's current row Map.  If
    /// you want X1 or X2 to be nonconst (modifiable) views, use
    /// offsetViewNonConst() with the same arguments.  "View" means
    /// "alias": if the original (this) MultiVector's data change, the
    /// view will see the changed data.
    ///
    /// \param subMap [in] The row Map for the new MultiVector.  This
    ///   must be a subset Map of this MultiVector's row Map.
    /// \param offset [in] The local row offset at which to start the view.
    ///
    /// Suppose that you have a MultiVector X, and you want to view X,
    /// on all processes in X's (MPI) communicator, as split into two
    /// row blocks X1 and X2.  One could express this in Matlab
    /// notation as X = [X1; X2], except that here, X1 and X2 are
    /// views into X, rather than copies of X's data.  This method
    /// assumes that the <i>local</i> indices of X1 and X2 are each
    /// contiguous, and that the local indices of X2 follow those of
    /// X1.  If that is not the case, you cannot use views to divide X
    /// into blocks like this; you must instead use the Import or
    /// Export functionality, which copies the relevant rows of X.
    ///
    /// Here is how you would construct the views X1 and X2.
    /// \code
    /// using Teuchos::RCP;
    /// typedef Tpetra::Map<LO, GO, Node> map_type;
    /// typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
    ///
    /// MV X (...); // the input MultiVector
    ///  //... fill X with data ...
    ///
    /// // Map that on each process in X's communicator,
    /// // contains the global indices of the rows of X1.
    /// RCP<const map_type> map1 (new map_type (...));
    /// // Map that on each process in X's communicator,
    /// // contains the global indices of the rows of X2.
    /// RCP<const map_type> map2 (new map_type (...));
    ///
    /// // Create the first view X1.  The second argument, the offset,
    /// // is the index of the local row at which to start the view.
    /// // X1 is the topmost block of X, so the offset is zero.
    /// RCP<const MV> X1 = X.offsetView (map1, 0);
    ///
    /// // Create the second view X2.  X2 is directly below X1 in X,
    /// // so the offset is the local number of rows in X1.  This is
    /// // the same as the local number of entries in map1.
    /// RCP<const MV> X2 = X.offsetView (map2, X1->getLocalLength ());
    /// \endcode
    ///
    /// It is legal, in the above example, for X1 or X2 to have zero
    /// local rows on any or all process(es).  In that case, the
    /// corresponding Map must have zero local entries on that / those
    /// process(es).  In particular, if X2 has zero local rows on a
    /// process, then the corresponding offset on that process would
    /// be the number of local rows in X (and therefore in X1) on that
    /// process.  This is the only case in which the sum of the local
    /// number of entries in \c subMap (in this case, zero) and the \c
    /// offset may equal the number of local entries in
    /// <tt>*this</tt>.
    Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    offsetView (const Teuchos::RCP<const map_type>& subMap,
                const size_t offset) const;

    /// \brief Return a nonconst view of a subset of rows.
    ///
    /// Return a nonconst (modifiable) view of this MultiVector
    /// consisting of a subset of the rows, as specified by an offset
    /// and a subset Map of this MultiVector's current row Map.  If
    /// you want X1 or X2 to be const (nonmodifiable) views, use
    /// offsetView() with the same arguments.  "View" means "alias":
    /// if the original (this) MultiVector's data change, the view
    /// will see the changed data, and if the view's data change, the
    /// original MultiVector will see the changed data.
    ///
    /// \param subMap [in] The row Map for the new MultiVector.  This
    ///   must be a subset Map of this MultiVector's row Map.
    /// \param offset [in] The local row offset at which to start the view.
    ///
    /// See the documentation of offsetView() for a code example and
    /// an explanation of edge cases.
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    offsetViewNonConst (const Teuchos::RCP<const map_type>& subMap,
                        const size_t offset);

    //! Return a Vector which is a const view of column j.
    Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    getVector (const size_t j) const;

    //! Return a Vector which is a nonconst view of column j.
    Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    getVectorNonConst (const size_t j);

    //! Const view of the local values in a particular vector of this multivector.
    Teuchos::ArrayRCP<const Scalar> getData (size_t j) const;

    //! View of the local values in a particular vector of this multivector.
    Teuchos::ArrayRCP<Scalar> getDataNonConst (size_t j);

    /// \brief Fill the given array with a copy of this multivector's
    ///   local values.
    ///
    /// \param A [out] View of the array to fill.  We consider A as a
    ///   matrix with column-major storage.
    ///
    /// \param LDA [in] Leading dimension of the matrix A.
    void
    get1dCopy (const Teuchos::ArrayView<Scalar>& A,
               const size_t LDA) const;

    /// \brief Fill the given array with a copy of this multivector's
    ///   local values.
    ///
    /// \param ArrayOfPtrs [out] Array of arrays, one for each column
    ///   of the multivector.  On output, we fill ArrayOfPtrs[j] with
    ///   the data for column j of this multivector.
    void
    get2dCopy (const Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> >& ArrayOfPtrs) const;

    /// \brief Const persisting (1-D) view of this multivector's local values.
    ///
    /// This method assumes that the columns of the multivector are
    /// stored contiguously.  If not, this method throws
    /// std::runtime_error.
    Teuchos::ArrayRCP<const Scalar> get1dView () const;

    //! Return const persisting pointers to values.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > get2dView () const;

    /// \brief Nonconst persisting (1-D) view of this multivector's local values.
    ///
    /// This method assumes that the columns of the multivector are
    /// stored contiguously.  If not, this method throws
    /// std::runtime_error.
    Teuchos::ArrayRCP<Scalar> get1dViewNonConst ();

    //! Return non-const persisting pointers to values.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > get2dViewNonConst ();

    /// \brief Get the Kokkos::DualView which implements local storage.
    ///
    /// \warning This method is scheduled for DEPRECATION.
    ///
    /// \warning This method is ONLY for expert developers.  Its
    ///   interface may change or it may disappear at any time.
    ///
    /// Instead of getting the Kokkos::DualView, we highly recommend
    /// calling the templated getLocalView() method, that returns a
    /// Kokkos::View of the MultiVector's data in a given memory
    /// space.  Since that MultiVector itself implements DualView
    /// semantics, it's much better to use MultiVector's interface to
    /// do "DualView things," like calling modify(), need_sync(), and
    /// sync().
    dual_view_type getDualView () const;

    /// \brief Update data on device or host only if data in the other
    ///   space has been marked as modified.
    ///
    /// If \c TargetDeviceType is the same as this MultiVector's
    /// device type, then copy data from host to device.  Otherwise,
    /// copy data from device to host.  In either case, only copy if
    /// the source of the copy has been modified.
    ///
    /// This is a one-way synchronization only.  If the target of the
    /// copy has been modified, this operation will discard those
    /// modifications.  It will also reset both device and host modified
    /// flags.
    ///
    /// \note This method doesn't know on its own whether you modified
    ///   the data in either memory space.  You must manually mark the
    ///   MultiVector as modified in the space in which you modified
    ///   it, by calling the modify() method with the appropriate
    ///   template parameter.
    template<class TargetDeviceType>
    void sync ();

    //! Whether this MultiVector needs synchronization to the given space.
    template<class TargetDeviceType>
    bool need_sync () const;

    /// \brief Mark data as modified on the given device \c TargetDeviceType.
    ///
    /// If \c TargetDeviceType is the same as this MultiVector's
    /// device type, then mark the device's data as modified.
    /// Otherwise, mark the host's data as modified.
    template<class TargetDeviceType>
    void modify ();

    /// \brief Return a view of the local data on a specific device.
    /// \tparam TargetDeviceType The Kokkos Device type whose data to return.
    ///
    /// Please don't be afraid of the if_c expression in the return
    /// value's type.  That just tells the method what the return type
    /// should be: dual_view_type::t_dev if the \c TargetDeviceType
    /// template parameter matches this Tpetra object's device type,
    /// else dual_view_type::t_host.
    ///
    /// For example, suppose you create a Tpetra::MultiVector for the
    /// Kokkos::Cuda device, like this:
    /// \code
    /// typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Cuda> > node_type;
    /// typedef Tpetra::Map<int, int, node_type> map_type;
    /// typedef Tpetra::MultiVector<float, int, int, node_type> mv_type;
    ///
    /// RCP<const map_type> map = ...;
    /// mv_type DV (map, 3);
    /// \endcode
    /// If you want to get the CUDA device Kokkos::View, do this:
    /// \code
    /// typedef typename mv_type::dual_view_type dual_view_type;
    /// typedef typename dual_view_type::t_dev device_view_type;
    /// device_view_type cudaView = DV.getLocalView<Kokkos::Cuda> ();
    /// \endcode
    /// and if you want to get the host mirror of that View, do this:
    /// \code
    /// typedef typename dual_view_type::host_mirror_space host_execution_space;
    /// typedef typename dual_view_type::t_host host_view_type;
    /// host_view_type hostView = DV.getLocalView<host_execution_space> ();
    /// \endcode
    template<class TargetDeviceType>
    typename Kokkos::Impl::if_c<
      std::is_same<
        typename device_type::memory_space,
        typename TargetDeviceType::memory_space>::value,
      typename dual_view_type::t_dev,
      typename dual_view_type::t_host>::type
    getLocalView () const;

    //@}
    //! @name Mathematical methods
    //@{

    /// \brief Compute the dot product of each corresponding pair of
    ///   vectors (columns) in A and B.
    ///
    /// The "dot product" is the standard Euclidean inner product.  If
    /// the type of entries of the vectors (impl_scalar_type) is complex,
    /// then A is transposed, not <tt>*this</tt>.  For example, if x
    /// and y each have one column, then <tt>x.dot (y, dots)</tt>
    /// computes \f$y^* x = \bar{y}^T x = \sum_i \bar{y}_i \cdot x_i\f$.
    ///
    /// \pre <tt>*this</tt> and A have the same number of columns (vectors).
    /// \pre \c dots has at least as many entries as the number of columns in A.
    ///
    /// \post <tt>dots[j] == (this->getVector[j])->dot (* (A.getVector[j]))</tt>
    void
    dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Teuchos::ArrayView<dot_type>& dots) const;

    /// \brief Compute the dot product of each corresponding pair of
    ///   vectors (columns) in A and B.
    /// \tparam T The output type of the dot products.
    ///
    /// This method only exists if dot_type and T are different types.
    /// For example, if impl_scalar_type and dot_type differ, then this
    /// method ensures backwards compatibility with the previous
    /// interface (that returned dot products as impl_scalar_type rather
    /// than as dot_type).  The complicated \c enable_if expression
    /// just ensures that the method only exists if dot_type and T are
    /// different types; the method still returns \c void, as above.
    template <typename T>
    typename std::enable_if< ! (std::is_same<dot_type, T>::value), void >::type
    dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Teuchos::ArrayView<T> &dots) const;

    //! Like the above dot() overload, but for std::vector output.
    template <typename T>
    typename std::enable_if< ! (std::is_same<dot_type, T>::value), void >::type
    dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, std::vector<T>& dots) const;

    /// \brief Compute the dot product of each corresponding pair of
    ///   vectors (columns) in A and B, storing the result in a device
    ///   View.
    ///
    /// The "dot product" is the standard Euclidean inner product.  If
    /// the type of entries of the vectors (impl_scalar_type) is complex,
    /// then A is transposed, not <tt>*this</tt>.  For example, if x
    /// and y each have one column, then <tt>x.dot (y, dots)</tt>
    /// computes \f$y^* x = \bar{y}^T x = \sum_i \bar{y}_i \cdot x_i\f$.
    ///
    /// \param A [in] MultiVector with which to dot \c *this.
    /// \param dots [out] Device View with getNumVectors() entries.
    ///
    /// \pre <tt>this->getNumVectors () == A.getNumVectors ()</tt>
    /// \pre <tt>dots.dimension_0 () == A.getNumVectors ()</tt>
    ///
    /// \post <tt>dots(j) == (this->getVector[j])->dot (* (A.getVector[j]))</tt>
    void
    dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Kokkos::View<dot_type*, Kokkos::HostSpace>& norms) const;

#if 1  // WCMCLEN
    template<class ViewType>
    void
    dot (typename std::enable_if<std::is_same<typename ViewType::value_type,dot_type>::value &&
                                 std::is_same<typename ViewType::memory_space,typename device_type::memory_space>::value,
                                 const FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>::type& A,
         const ViewType& dots) const;
#endif

#if 0  // WCMCLEN
// TODO: WCMCLEN - This wasn't compiling...
    /// \brief Compute the dot product of each corresponding pair of
    ///   vectors (columns) in A and B, storing the result in a device
    ///   view.
    /// \tparam T The output type of the dot products.
    ///
    /// This method only exists if dot_type and T are different types.
    /// For example, if Scalar and dot_type differ, then this method
    /// ensures backwards compatibility with the previous interface
    /// (that returned dot products as Scalar rather than as
    /// dot_type).  The complicated \c enable_if expression just
    /// ensures that the method only exists if dot_type and T are
    /// different types; the method still returns \c void, as above.
    template <typename T>
    typename std::enable_if< ! (std::is_same<dot_type, T>::value), void >::type
    dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Kokkos::View<T*, device_type>& dots) const;
#endif

    //! Put element-wise absolute values of input Multi-vector in target: A = abs(this)
    void abs (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A);

    //! Put element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    void reciprocal (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A);

    /// \brief Scale in place: <tt>this = alpha*this</tt>.
    ///
    /// Replace this MultiVector with alpha times this MultiVector.
    /// This method will always multiply, even if alpha is zero.  That
    /// means, for example, that if \c *this contains NaN entries
    /// before calling this method, the NaN entries will remain after
    /// this method finishes.
    void scale (const Scalar& alpha);

    /// \brief Scale each column in place: <tt>this[j] = alpha[j]*this[j]</tt>.
    ///
    /// Replace each column j of this MultiVector with
    /// <tt>alpha[j]</tt> times the current column j of this
    /// MultiVector.  This method will always multiply, even if all
    /// the entries of alpha are zero.  That means, for example, that
    /// if \c *this contains NaN entries before calling this method,
    /// the NaN entries will remain after this method finishes.
    void scale (const Teuchos::ArrayView<const Scalar>& alpha);

    /// \brief Scale each column in place: <tt>this[j] = alpha[j]*this[j]</tt>.
    ///
    /// Replace each column j of this MultiVector with
    /// <tt>alpha[j]</tt> times the current column j of this
    /// MultiVector.  This method will always multiply, even if all
    /// the entries of alpha are zero.  That means, for example, that
    /// if \c *this contains NaN entries before calling this method,
    /// the NaN entries will remain after this method finishes.
    void scale (const Kokkos::View<const impl_scalar_type*, device_type>& alpha);

    /// \brief Scale in place: <tt>this = alpha * A</tt>.
    ///
    /// Replace this MultiVector with scaled values of A.  This method
    /// will always multiply, even if alpha is zero.  That means, for
    /// example, that if \c *this contains NaN entries before calling
    /// this method, the NaN entries will remain after this method
    /// finishes.  It is legal for the input A to alias this
    /// MultiVector.
    // TODO: Should A be a MultiVector or FEMultiVector?  (WCMCLEN)
    void
    scale (const Scalar& alpha,
           const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A);

    /// \brief Update: <tt>this = beta*this + alpha*A</tt>.
    ///
    /// Update this MultiVector with scaled values of A.  If beta is
    /// zero, overwrite \c *this unconditionally, even if it contains
    /// NaN entries.  It is legal for the input A to alias this
    /// MultiVector.
    // TODO: Should A be a MultiVector or FEMultiVector?  (WCMCLEN)
    void
    update (const Scalar& alpha,
            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
            const Scalar& beta);

    /// \brief Update: <tt>this = gamma*this + alpha*A + beta*B</tt>.
    ///
    /// Update this MultiVector with scaled values of A and B.  If
    /// gamma is zero, overwrite \c *this unconditionally, even if it
    /// contains NaN entries.  It is legal for the inputs A or B to
    /// alias this MultiVector.
    // TODO: Should A and B be a MultiVector or FEMultiVector?  (WCMCLEN)
    void
    update (const Scalar& alpha,
            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
            const Scalar& beta,
            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
            const Scalar& gamma);

    /// \brief Compute the one-norm of each vector (column), storing
    ///   the result in a device view.
    ///
    /// \param norms [out] Device View with getNumVectors() entries.
    ///
    /// \pre <tt>norms.dimension_0 () == this->getNumVectors ()</tt>
    /// \post <tt>norms(j) == (this->getVector[j])->norm1 (* (A.getVector[j]))</tt>
    ///
    /// The one-norm of a vector is the sum of the magnitudes of the
    /// vector's entries.  On exit, norms(j) is the one-norm of column
    /// j of this MultiVector.
    ///
    /// We use Kokkos::Details::InnerProductSpaceTraits to define
    /// "magnitude."  See the KokkosKernels package, and that class'
    /// documentation in particular, for details.  This matters only
    /// for "interesting" Scalar types, such as one might find in the
    /// Stokhos package.
    template<class ViewType>
    typename std::enable_if<std::is_same<typename ViewType::value_type,mag_type>::value &&
                            std::is_same<typename ViewType::memory_space,typename device_type::memory_space>::value>::type
    norm1 (const ViewType& norms) const;

    void norm1 (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;

    /// \brief Compute the one-norm of each vector (column), storing
    ///   the result in a device view.
    /// \tparam T The output type of the dot products.
    ///
    /// See the above norm1() method for documentation.
    ///
    /// This method only exists if mag_type and T are different types.
    /// For example, if Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// and mag_type differ, then this method ensures backwards
    /// compatibility with the previous interface (that returned norms
    /// products as Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// rather than as mag_type).  The complicated \c enable_if
    /// expression just ensures that the method only exists if
    /// mag_type and T are different types; the method still returns
    /// \c void, as above.
    template <typename T>
    typename std::enable_if< ! (std::is_same<mag_type, T>::value), void >::type
    norm1 (const Kokkos::View<T*, device_type>& norms) const;

    /// \brief Compute the one-norm of each vector (column).
    ///
    /// See the uppermost norm1() method above for documentation.
    void norm1 (const Teuchos::ArrayView<mag_type>& norms) const;

    /// \brief Compute the one-norm of each vector (column).
    /// \tparam T The output type of the norms.
    ///
    /// See the uppermost norm1() method above for documentation.
    ///
    /// This method only exists if mag_type and T are different types.
    /// For example, if Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// and mag_type differ, then this method ensures backwards
    /// compatibility with the previous interface (that returned norms
    /// as Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// rather than as mag_type).  The complicated \c enable_if
    /// expression just ensures that the method only exists if
    /// mag_type and T are different types; the method still returns
    /// \c void, as above.
    template <typename T>
    typename std::enable_if< ! (std::is_same<mag_type,T>::value), void >::type
    norm1 (const Teuchos::ArrayView<T>& norms) const;

    /// \brief Compute the two-norm of each vector (column), storing
    ///   the result in a device view.
    ///
    /// \param norms [out] Device View with getNumVectors() entries.
    ///
    /// \pre <tt>norms.dimension_0 () == this->getNumVectors ()</tt>
    /// \post <tt>norms(j) == (this->getVector[j])->dot (* (A.getVector[j]))</tt>
    ///
    /// The two-norm of a vector is the standard Euclidean norm, the
    /// square root of the sum of squares of the magnitudes of the
    /// vector's entries.  On exit, norms(k) is the two-norm of column
    /// k of this MultiVector.
    ///
    /// We use Kokkos::Details::InnerProductSpaceTraits to define
    /// "magnitude."  See the KokkosKernels package, and that class'
    /// documentation in particular, for details.  This matters only
    /// for "interesting" Scalar types, such as one might find in the
    /// Stokhos package.
    template<class ViewType>
    typename std::enable_if<std::is_same<typename ViewType::value_type,mag_type>::value &&
                            std::is_same<typename ViewType::memory_space,typename device_type::memory_space>::value>::type
    norm2 (const ViewType& norms) const;

    void norm2 (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;

    /// \brief Compute the two-norm of each vector (column), storing
    ///   the result in a device view.
    ///
    /// See the above norm2() method for documentation.
    ///
    /// This method only exists if mag_type and T are different types.
    /// For example, if Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// and mag_type differ, then this method ensures backwards
    /// compatibility with the previous interface (that returned norms
    /// as Teuchos::ScalarTraits<Scalar>::magnitudeType rather than as
    /// mag_type).  The complicated \c enable_if expression just
    /// ensures that the method only exists if mag_type and T are
    /// different types; the method still returns \c void, as above.
    template<typename T>
    typename std::enable_if< ! (std::is_same<mag_type, T>::value), void >::type
    norm2 (const Kokkos::View<T*, device_type>& norms) const;

    /// \brief Compute the two-norm of each vector (column).
    ///
    /// See the uppermost norm2() method above for documentation.
    void norm2 (const Teuchos::ArrayView<mag_type>& norms) const;

    /// \brief Compute the two-norm of each vector (column).
    /// \tparam T The output type of the norms.
    ///
    /// See the uppermost norm2() method above for documentation.
    ///
    /// This method only exists if mag_type and T are different types.
    /// For example, if Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// and mag_type differ, then this method ensures backwards
    /// compatibility with the previous interface (that returned norms
    /// products as Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// rather than as mag_type).  The complicated \c enable_if
    /// expression just ensures that the method only exists if
    /// mag_type and T are different types; the method still returns
    /// \c void, as above.
    template <typename T>
    typename std::enable_if< ! (std::is_same<mag_type,T>::value), void >::type
    norm2 (const Teuchos::ArrayView<T>& norms) const;

    /// \brief Compute the infinity-norm of each vector (column),
    ///   storing the result in a device view.
    ///
    /// The infinity-norm of a vector is the maximum of the magnitudes
    /// of the vector's entries.  On exit, norms(j) is the
    /// infinity-norm of column j of this MultiVector.
    ///
    /// We use Kokkos::Details::InnerProductSpaceTraits to define
    /// "magnitude."  See the KokkosKernels package, and that class'
    /// documentation in particular, for details.  This matters only
    /// for "interesting" Scalar types, such as one might find in the
    /// Stokhos package.
    template<class ViewType>
    typename std::enable_if<std::is_same<typename ViewType::value_type,mag_type>::value &&
                            std::is_same<typename ViewType::memory_space,typename device_type::memory_space>::value>::type
    normInf (const ViewType& norms) const;

    void normInf (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;

    /// \brief Compute the two-norm of each vector (column), storing
    ///   the result in a device view.
    ///
    /// See the above normInf() method for documentation.
    ///
    /// This method only exists if mag_type and T are different types.
    /// For example, if Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// and mag_type differ, then this method ensures backwards
    /// compatibility with the previous interface (that returned norms
    /// as Teuchos::ScalarTraits<Scalar>::magnitudeType rather than as
    /// mag_type).  The complicated \c enable_if expression just
    /// ensures that the method only exists if mag_type and T are
    /// different types; the method still returns \c void, as above.
    template<typename T>
    typename std::enable_if< ! (std::is_same<mag_type, T>::value), void >::type
    normInf (const Kokkos::View<T*, device_type>& norms) const;

    /// \brief Compute the infinity-norm of each vector (column),
    ///   storing the result in a Teuchos::ArrayView.
    ///
    /// See the uppermost normInf() method above for documentation.
    void normInf (const Teuchos::ArrayView<mag_type>& norms) const;

    /// \brief Compute the infinity-norm of each vector (column),
    ///   storing the result in a Teuchos::ArrayView.
    /// \tparam T The output type of the norms.
    ///
    /// See the uppermost normInf() method above for documentation.
    ///
    /// This method only exists if mag_type and T are different types.
    /// For example, if Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// and mag_type differ, then this method ensures backwards
    /// compatibility with the previous interface (that returned norms
    /// products as Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// rather than as mag_type).  The complicated \c enable_if
    /// expression just ensures that the method only exists if
    /// mag_type and T are different types; the method still returns
    /// \c void, as above.
    template <typename T>
    typename std::enable_if< ! (std::is_same<mag_type,T>::value), void >::type
    normInf (const Teuchos::ArrayView<T>& norms) const;

    /// \brief Compute Weighted 2-norm (RMS Norm) of each column.
    ///
    /// \warning This method has been DEPRECATED.
    ///
    /// The results of this method are undefined for scalar types that
    /// are not floating-point types (e.g., int).
    void TPETRA_DEPRECATED
    normWeighted (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& weights,
                  const Teuchos::ArrayView<mag_type>& norms) const;

    /// \brief Compute the weighted 2-norm (RMS Norm) of each column.
    ///
    /// \warning This method is DEPRECATED.
    ///
    /// The outcome of this routine is undefined for non-floating
    /// point scalar types (e.g., int).
    ///
    /// This method only exists if mag_type and T are different types.
    /// For example, if Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// and mag_type differ, then this method ensures backwards
    /// compatibility with the previous interface (that returned norms
    /// as Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// rather than as mag_type).  The complicated \c enable_if
    /// expression just ensures that the method only exists if
    /// mag_type and T are different types; the method still returns
    /// \c void, as above.
    template <typename T>
    typename std::enable_if< ! (std::is_same<mag_type,T>::value), void >::type
    TPETRA_DEPRECATED
    normWeighted (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& weights,
                  const Teuchos::ArrayView<T>& norms) const;

    /// \brief Compute mean (average) value of each column.
    ///
    /// The outcome of this routine is undefined for non-floating
    /// point scalar types (e.g., int).
    void meanValue (const Teuchos::ArrayView<impl_scalar_type>& means) const;

    template <typename T>
    typename std::enable_if<! std::is_same<impl_scalar_type, T>::value, void>::type
    meanValue (const Teuchos::ArrayView<T>& means) const;

    /// \brief Matrix-matrix multiplication: <tt>this = beta*this + alpha*op(A)*op(B)</tt>.
    ///
    /// If beta is zero, overwrite \c *this unconditionally, even if
    /// it contains NaN entries.  This imitates the semantics of
    /// analogous BLAS routines like DGEMM.
    void
    multiply (Teuchos::ETransp transA,
              Teuchos::ETransp transB,
              const Scalar& alpha,
              const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
              const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
              const Scalar& beta);

    /// \brief Multiply a Vector A elementwise by a MultiVector B.
    ///
    /// Compute <tt>this = scalarThis * this + scalarAB * B @ A</tt>
    /// where <tt>@</tt> denotes element-wise multiplication.  In
    /// pseudocode, if C denotes <tt>*this</tt> MultiVector:
    /// \code
    /// C(i,j) = scalarThis * C(i,j) + scalarAB * B(i,j) * A(i,1);
    /// \endcode
    /// for all rows i and columns j of C.
    ///
    /// B must have the same dimensions as <tt>*this</tt>, while A
    /// must have the same number of rows but a single column.
    ///
    /// We do not require that A, B, and <tt>*this</tt> have
    /// compatible Maps, as long as the number of rows in A, B, and
    /// <tt>*this</tt> on each process is the same.  For example, one
    /// or more of these vectors might have a locally replicated Map,
    /// or a Map with a local communicator (<tt>MPI_COMM_SELF</tt>).
    /// This case may occur in block relaxation algorithms when
    /// applying a diagonal scaling.
    void
    elementWiseMultiply (Scalar scalarAB,
                         const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                         const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                         Scalar scalarThis);
    //@}
    //! @name Attribute access functions
    //@{

    //! Number of columns in the multivector.
    size_t getNumVectors() const;

    //! Local number of rows on the calling process.
    size_t getLocalLength() const;

    //! Global number of rows in the multivector.
    global_size_t getGlobalLength() const;

    /// \brief Stride between columns in the multivector.
    ///
    /// This is only meaningful if \c isConstantStride() returns true.
    ///
    /// \warning This may be different on different processes.
    size_t getStride() const;

    /// \brief Whether this multivector has constant stride between columns.
    ///
    /// \warning This may be different on different processes.
    bool isConstantStride() const;

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! A simple one-line description of this object.
    virtual std::string description() const;

    /// \brief Print the object with the given verbosity level to a FancyOStream.
    ///
    /// \param out [out] Output stream to which to print.  For
    ///   verbosity levels VERB_LOW and lower, only the process with
    ///   rank 0 ("Proc 0") in the MultiVector's communicator prints.
    ///   For verbosity levels strictly higher than VERB_LOW, all
    ///   processes in the communicator need to be able to print to
    ///   the output stream.
    ///
    /// \param verbLevel [in] Verbosity level.  The default verbosity
    ///   (verbLevel=VERB_DEFAULT) is VERB_LOW.
    ///
    /// The amount and content of what this method prints depends on
    /// the verbosity level.  In the list below, each higher level
    /// includes all the content of the previous levels, as well as
    /// its own content.
    ///
    /// - VERB_LOW: Only Proc 0 prints; it prints the same thing as \c
    ///   description().
    /// - VERB_MEDIUM: Each process prints its local length (the
    ///   number of rows that it owns).
    /// - VERB_HIGH: Each process prints whether the multivector has
    ///   constant stride (see \c isConstantStride()), and if so, what
    ///   that stride is.  (Stride may differ on different processes.)
    /// - VERB_EXTREME: Each process prints the values in its local
    ///   part of the multivector.  This will print out as many rows
    ///   of data as the global number of rows in the multivector, so
    ///   beware.
    virtual void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const;
    //@}

    /// \brief Remove processes owning zero rows from the Map and their communicator.
    ///
    /// \warning This method is ONLY for use by experts.  We highly
    ///   recommend using the nonmember function of the same name
    ///   defined in Tpetra_DistObject_decl.hpp.
    ///
    /// \warning We make NO promises of backwards compatibility.
    ///   This method may change or disappear at any time.
    ///
    /// \param newMap [in] This <i>must</i> be the result of calling
    ///   the removeEmptyProcesses() method on the row Map.  If it
    ///   is not, this method's behavior is undefined.  This pointer
    ///   will be null on excluded processes.
    virtual void
    removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap);

    /// \brief Set whether this has copy (copyOrView = Teuchos::Copy)
    ///   or view (copyOrView = Teuchos::View) semantics.
    ///
    /// \warning The Kokkos refactor version of MultiVector
    ///   <i>only</i> implements view semantics.  If you attempt to
    ///   call this method with copyOrView == Teuchos::Copy, it will
    ///   throw std::invalid_argument.
    ///
    /// \warning This method is only for expert use.  It may change or
    ///   disappear at any time.
    void setCopyOrView (const Teuchos::DataAccess copyOrView);

    /// \brief Get whether this has copy (copyOrView = Teuchos::Copy)
    ///   or view (copyOrView = Teuchos::View) semantics.
    ///
    /// \note The Kokkos refactor version of MultiVector <i>only</i>
    ///   implements view semantics.  This is not currently true for
    ///   the "classic" version of MultiVector, though that will
    ///   change in the near future.
    ///
    /// \warning This method is only for expert use.  It may change or
    ///   disappear at any time.
    Teuchos::DataAccess getCopyOrView () const;

    /// \brief Copy the contents of \c src into \c *this (deep copy).
    ///
    /// \param src [in] Source MultiVector (input of the deep copy).
    ///
    /// \pre <tt> ! src.getMap ().is_null () && ! this->getMap ().is_null () </tt>
    /// \pre <tt> src.getMap ()->isCompatible (* (this->getMap ()) </tt>
    ///
    /// \post Any outstanding views of \c src or \c *this remain valid.
    ///
    /// \note To implementers: The postcondition implies that the
    ///   implementation must not reallocate any memory of \c *this,
    ///   or otherwise change its dimensions.  This is <i>not</i> an
    ///   assignment operator; it does not change anything in \c *this
    ///   other than the contents of storage.
    void
    assign (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& src);

  protected:
    template <class DS, class DL, class DG, class DN,
              class SS, class SL, class SG, class SN>
    friend void
    deep_copy (FEMultiVector<DS, DL, DG, DN>& dst,
               const FEMultiVector<SS, SL, SG, SN>& src);


    // WCMCLEN: @CMS did you remove these originally?  They're used in _def, do we need them?
    //          see MultiVector_decl for descriptions of what these do
    mutable dual_view_type view_;
    mutable dual_view_type origView_;
    Teuchos::Array<size_t> whichVectors_;

    //! Input argument for normImpl()
    enum EWhichNorm {
      NORM_ONE,  //<! Use the one-norm
      NORM_TWO,  //<! Use the two-norm
      NORM_INF   //<! Use the infinity-norm
    };

    /// \brief Compute the norm of each vector (column), storing the
    ///   result in a device View.
    ///
    /// This method consolidates all common code between the
    /// infinity-norm, 1-norm, and 2-norm calculations.  On exit,
    /// norms(j) is the norm (of the selected type) of column j of
    /// this MultiVector.
    void
    normImpl (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms,
              const EWhichNorm whichNorm) const;

    //@}
    //! @name Misc. implementation details
    //@{

    /// \brief Implementation of description() for this class, and its
    ///   subclass Vector.
    ///
    /// \param className [in] Name of the class calling this method:
    ///   Either "Tpetra::MultiVector" or "Tpetra::Vector" (no quotes
    ///   in the string, in either case).
    std::string
    descriptionImpl (const std::string& className) const;

    /// \brief Print the calling process' verbose describe()
    ///   information to the returned string.
    ///
    /// This is an implementation detail of describe().
    ///
    /// \param vl [in] Verbosity level with which to print.
    std::string
    localDescribeToString (const Teuchos::EVerbosityLevel vl) const;

    /// \brief Implementation of describe() for this class, and its
    ///   subclass Vector.
    ///
    /// \param out [out] Output stream to which to write.  Only
    ///   Process 0 in this object's communicator may write to the
    ///   output stream.
    ///
    /// \param className [in] Name of the class calling this method.
    ///
    /// \param verbLevel [in] Verbosity level.  This also controls
    ///   whether this method does any communication.  At verbosity
    ///   levels higher (greater) than Teuchos::VERB_LOW, this method
    ///   behaves as a collective over the object's communicator.
    void
    describeImpl (Teuchos::FancyOStream& out,
                  const std::string& className,
                  const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const;

    // Return true if and only if VectorIndex is a valid column index.
    bool vectorIndexOutOfRange (const size_t VectorIndex) const;

    /// \brief Persisting view of j-th column in the given ArrayRCP.
    ///
    /// This method considers isConstantStride().  The ArrayRCP may
    /// correspond either to a compute buffer or a host view.
    template <class T>
    Teuchos::ArrayRCP<T>
    getSubArrayRCP (Teuchos::ArrayRCP<T> arr, size_t j) const;

    //! "Original" number of rows in the (local) data.
    size_t getOrigNumLocalRows () const;

    //! "Original" number of columns in the (local) data.
    size_t getOrigNumLocalCols () const;

    //@}
    //! @name Implementation of Tpetra::DistObject
    //@{

    /// \brief Whether data redistribution between \c sourceObj and this object is legal.
    ///
    /// This method is called in DistObject::doTransfer() to check
    /// whether data redistribution between the two objects is legal.
    virtual bool
    checkSizes (const SrcDistObject& sourceObj);

    //! Number of packets to send per LID
    virtual size_t constantNumberOfPackets () const;

    //! Whether this class implements the old or new interface of DistObject.
    virtual bool useNewInterface () { return true; }

    typedef typename DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>::buffer_device_type buffer_device_type;

    virtual void
    copyAndPermuteNew (const SrcDistObject& sourceObj,
                       const size_t numSameIDs,
                       const Kokkos::DualView<const local_ordinal_type*, device_type>& permuteToLIDs,
                       const Kokkos::DualView<const local_ordinal_type*, device_type>& permuteFromLIDs);

    virtual void
    packAndPrepareNew (const SrcDistObject& sourceObj,
                       const Kokkos::DualView<const local_ordinal_type*, device_type>& exportLIDs,
                       Kokkos::DualView<impl_scalar_type*, buffer_device_type>& exports,
                       const Kokkos::DualView<size_t*, buffer_device_type>& /* numPacketsPerLID */,
                       size_t& constantNumPackets,
                       Distributor& /* distor */);

    virtual void
    unpackAndCombineNew (const Kokkos::DualView<const LocalOrdinal*, device_type>& importLIDs,
                         const Kokkos::DualView<const impl_scalar_type*, buffer_device_type>& imports,
                         const Kokkos::DualView<const size_t*, buffer_device_type>& /* numPacketsPerLID */,
                         const size_t constantNumPackets,
                         Distributor& /* distor */,
                         const CombineMode CM);
    //@}
  }; // class MultiVector


} // namespace Tpetra


namespace Teuchos {

  // Give Teuchos::TypeNameTraits<Tpetra::MultiVector<...> > a
  // human-readable definition.
  template<class SC, class LO, class GO, class NT>
  class TypeNameTraits<Tpetra::FEMultiVector<SC, LO, GO, NT> > {
  public:
    static std::string name () {
      return std::string ("Tpetra::FEMultiVector<") +
        TypeNameTraits<SC>::name () + "," +
        TypeNameTraits<LO>::name () + "," +
        TypeNameTraits<GO>::name () + "," +
        TypeNameTraits<NT>::name () + ">";
    }

    static std::string
    concreteName (const Tpetra::FEMultiVector<SC, LO, GO, NT>&) {
      return name ();
    }
  };
} // namespace Teuchos

#endif // TPETRA_FEMULTIVECTOR_DECL_HPP
