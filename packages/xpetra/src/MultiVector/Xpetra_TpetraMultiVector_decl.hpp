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
#ifndef XPETRA_TPETRAMULTIVECTOR_DECL_HPP
#define XPETRA_TPETRAMULTIVECTOR_DECL_HPP

#include "Xpetra_TpetraConfigDefs.hpp"
#include "Xpetra_MultiVector_decl.hpp"

#include "Xpetra_TpetraMap_decl.hpp"
#include "Xpetra_TpetraImport_decl.hpp"
#include "Xpetra_TpetraExport_decl.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Xpetra_Utils.hpp"

namespace Xpetra {

  // TODO: move that elsewhere
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const Tpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node> & toTpetra(const MultiVector< Scalar,LocalOrdinal, GlobalOrdinal, Node> &);

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Tpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node> & toTpetra(MultiVector< Scalar,LocalOrdinal, GlobalOrdinal, Node> &);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of TpetraVector, needed to prevent circular inclusions
  template<class S, class LO, class GO, class N> class TpetraVector;
#endif


  // Because we aren't including the header...
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node > >       toXpetra(RCP<Tpetra::Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > > vec);

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node > > toXpetra(RCP<const Tpetra::Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > > vec);



  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class TpetraMultiVector
    : public virtual MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node >
  {

    // The following typedef are used by the XPETRA_DYNAMIC_CAST() macro.
    typedef TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraMultiVectorClass;

  public:

    //! @name Constructors and destructor
    //@{

    //! Basic constuctor.
    TpetraMultiVector(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map, size_t NumVectors, bool zeroOut=true);

    //! Copy constructor (performs a deep copy).
    TpetraMultiVector(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &source);

    //! Create multivector by copying two-dimensional array of local data.
    TpetraMultiVector(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map, const Teuchos::ArrayView< const Scalar > &A, size_t LDA, size_t NumVectors);

    //! Create multivector by copying array of views of local data.
    TpetraMultiVector(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map, const Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > &ArrayOfPtrs, size_t NumVectors);

    virtual ~TpetraMultiVector();


    //! @name Post-construction modification routines
    //@{

    //! Replace value, using global (row) index.
    void replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value);


    //! Add value to existing value, using global (row) index.
    void sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value);

    //! Replace value, using local (row) index.
    void replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value);

    //! Add value to existing value, using local (row) index.
    void sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value);

    //! Set all values in the multivector with the given value.
    void putScalar(const Scalar &value);

    //! Sum values of a locally replicated multivector across all processes.
    void reduce();

    //@}

    //! @name Data Copy and View get methods
    //@{

    //! Return a Vector which is a const view of column j.
    Teuchos::RCP< const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > > getVector(size_t j) const;

    //! Return a Vector which is a nonconst view of column j.
    Teuchos::RCP< Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > > getVectorNonConst(size_t j);

    //! Const view of the local values in a particular vector of this multivector.
    Teuchos::ArrayRCP< const Scalar > getData(size_t j) const;

    //! View of the local values in a particular vector of this multivector.
    Teuchos::ArrayRCP< Scalar > getDataNonConst(size_t j);

    //! Fill the given array with a copy of this multivector's local values.
    void get1dCopy(Teuchos::ArrayView< Scalar > A, size_t LDA) const;

    //! Fill the given array with a copy of this multivector's local values.
    void get2dCopy(Teuchos::ArrayView< const Teuchos::ArrayView< Scalar > > ArrayOfPtrs) const;

    //! Const persisting (1-D) view of this multivector's local values.
    Teuchos::ArrayRCP< const Scalar > get1dView() const;

    //! Return const persisting pointers to values.
    Teuchos::ArrayRCP< Teuchos::ArrayRCP< const Scalar > > get2dView() const;

    //! Nonconst persisting (1-D) view of this multivector's local values.
    Teuchos::ArrayRCP< Scalar > get1dViewNonConst();

    //! Return non-const persisting pointers to values.
    Teuchos::ArrayRCP< Teuchos::ArrayRCP< Scalar > > get2dViewNonConst();


    //@}

    //! @name Mathematical methods
    //@{

    //! Compute dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i]).
    void dot(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, const Teuchos::ArrayView< Scalar > &dots) const;

    //! Put element-wise absolute values of input Multi-vector in target: A = abs(this).
    void abs(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A);

    //! Put element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    void reciprocal(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A);

    //! Scale the current values of a multi-vector, this = alpha*this.
    void scale(const Scalar &alpha);

    //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
    void scale(Teuchos::ArrayView< const Scalar > alpha);

    //! Replace multi-vector values with scaled values of A, this = alpha*A.
    void scale(const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A);

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    void update(const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, const Scalar &beta);

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    void update(const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, const Scalar &beta, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &B, const Scalar &gamma);

    //! Compute 1-norm of each vector in multi-vector.
    void norm1(const Teuchos::ArrayView< typename Teuchos::ScalarTraits< Scalar >::magnitudeType > &norms) const;

    //!
    void norm2(const Teuchos::ArrayView< typename Teuchos::ScalarTraits< Scalar >::magnitudeType > &norms) const;

    //! Compute Inf-norm of each vector in multi-vector.
    void normInf(const Teuchos::ArrayView< typename Teuchos::ScalarTraits< Scalar >::magnitudeType > &norms) const;

    //! Compute mean (average) value of each vector in multi-vector. The outcome of this routine is undefined for non-floating point scalar types (e.g., int).
    void meanValue(const Teuchos::ArrayView< Scalar > &means) const;

    //! Matrix-matrix multiplication: this = beta*this + alpha*op(A)*op(B).
    void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &B, const Scalar &beta);

    //@}
    //! @name Attribute access functions
    //@{

    //! Number of columns in the multivector.
    size_t getNumVectors() const;

    //! Local number of rows on the calling process.
    size_t getLocalLength() const;

    //! Global number of rows in the multivector.
    global_size_t getGlobalLength() const;

    // \! Checks to see if the local length, number of vectors and size of Scalar type match
    bool isSameSize(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & vec) const;

    //@}
    //! @name Overridden from Teuchos::Describable
    //@{

    //! A simple one-line description of this object.
    std::string description() const;

    //! Print the object with the given verbosity level to a FancyOStream.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

    //@}

    //! Element-wise multiply of a Vector A with a TpetraMultiVector B.
    void elementWiseMultiply(Scalar scalarAB, const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, Scalar scalarThis); // definition at the end of this file
    //TODO: void elementWiseMultiply(Scalar scalarAB, const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &B, Scalar scalarThis){ vec_->elementWiseMultiply(scalarAB, toTpetra(A), toTpetra(B), scalarThis); }

    //! Set multi-vector values to random numbers.
    void randomize(bool bUseXpetraImplementation = false);

    //{@
    // Implements DistObject interface

    Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> > getMap() const;

    void doImport(const DistObject< Scalar, LocalOrdinal,GlobalOrdinal,Node> &source, const Import<LocalOrdinal,GlobalOrdinal,Node> &importer, CombineMode CM);

    void doExport(const DistObject< Scalar, LocalOrdinal, GlobalOrdinal, Node > &dest, const Import<LocalOrdinal,GlobalOrdinal,Node>& importer, CombineMode CM);

    void doImport(const DistObject< Scalar, LocalOrdinal, GlobalOrdinal, Node > &source, const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter, CombineMode CM);

    void doExport(const DistObject< Scalar, LocalOrdinal, GlobalOrdinal, Node > &dest, const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter, CombineMode CM);

    void replaceMap(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map);

//@}


    //! @name Xpetra specific
    //@{

    //! TpetraMultiVector constructor to wrap a Tpetra::MultiVector object
    TpetraMultiVector(const Teuchos::RCP<Tpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node> > &vec);

    //! Get the underlying Tpetra multivector
    RCP< Tpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node> > getTpetra_MultiVector() const;

    //! Set seed for Random function.
    void setSeed(unsigned int seed);


#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
    typedef typename Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type dual_view_type;
    /// \brief Return an unmanaged non-const view of the local data on a specific device.
    /// \tparam TargetDeviceType The Kokkos Device type whose data to return.
    ///
    /// \warning DO NOT USE THIS FUNCTION! There is no reason why you are working directly
    ///          with the Xpetra::TpetraMultiVector object. To write a code which is independent
    ///          from the underlying linear algebra package you should always use the abstract class,
    ///          i.e. Xpetra::MultiVector!
    ///
    /// \warning Be aware that the view on the multivector data is non-persisting, i.e.
    ///          only valid as long as the multivector does not run of scope!
#if 0
    template<class TargetDeviceType>
    typename Kokkos::Impl::if_c<
      std::is_same<
        typename dual_view_type::t_dev_um::execution_space::memory_space,
        typename TargetDeviceType::memory_space>::value,
      typename dual_view_type::t_dev_um,
      typename dual_view_type::t_host_um>::type
    getLocalView () const;
#endif
    typename dual_view_type::t_host_um getHostLocalView () const;

    typename dual_view_type::t_dev_um
    getDeviceLocalView() const;

#endif


  protected:
    /// \brief Implementation of the assignment operator (operator=);
    ///   does a deep copy.
    virtual void
    assign (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs);

  private:
    //! The Tpetra::MultiVector which this class wraps.
    RCP< Tpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node> > vec_;

  }; // TpetraMultiVector class


  // Things we actually need

  // Things we actually need
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node > > toXpetra(RCP<Tpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > > vec) {
    if (!vec.is_null())
      return rcp(new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node >(vec));

    return Teuchos::null;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node > > toXpetra(RCP<const Tpetra::MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > > vec) {
    if (!vec.is_null())
      return rcp(new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node >(vec));

    return Teuchos::null;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const Tpetra::MultiVector< Scalar,LocalOrdinal, GlobalOrdinal, Node> & toTpetra(const MultiVector< Scalar,LocalOrdinal, GlobalOrdinal, Node> &x) {
    typedef TpetraMultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > TpetraMultiVectorClass;
    XPETRA_DYNAMIC_CAST(const TpetraMultiVectorClass, x, tX, "toTpetra");
    return *tX.getTpetra_MultiVector();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Tpetra::MultiVector< Scalar,LocalOrdinal, GlobalOrdinal, Node> & toTpetra(MultiVector< Scalar,LocalOrdinal, GlobalOrdinal, Node> &x) {
    typedef TpetraMultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > TpetraMultiVectorClass;
    XPETRA_DYNAMIC_CAST(      TpetraMultiVectorClass, x, tX, "toTpetra");
    return *tX.getTpetra_MultiVector();
  }

} // Xpetra namespace

#define XPETRA_TPETRAMULTIVECTOR_SHORT
#endif // XPETRA_TPETRAMULTIVECTOR_HPP

