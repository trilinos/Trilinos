// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_EPETRAINTMULTIVECTOR_HPP
#define XPETRA_EPETRAINTMULTIVECTOR_HPP

#include "Xpetra_EpetraConfigDefs.hpp"

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_EpetraMap.hpp"
#include "Xpetra_EpetraMultiVector.hpp"
#include "Epetra_IntMultiVector.h"

#if defined(XPETRA_ENABLE_DEPRECATED_CODE)
#ifdef __GNUC__
#if defined(Xpetra_SHOW_DEPRECATED_WARNINGS)
#warning "The header file Trilinos/packages/xpetra/src/MultiVector/Xpetra_EpetraIntMultiVector.hpp is deprecated."
#endif
#endif
#else
#error "The header file Trilinos/packages/xpetra/src/MultiVector/Xpetra_EpetraIntMultiVector.hpp is deprecated."
#endif

namespace Xpetra {

// TODO: move that elsewhere
template <class GlobalOrdinal, class Node>
XPETRA_DEPRECATED Epetra_IntMultiVector &toEpetra(MultiVector<int, int, GlobalOrdinal, Node> &);

template <class GlobalOrdinal, class Node>
XPETRA_DEPRECATED const Epetra_IntMultiVector &toEpetra(const MultiVector<int, int, GlobalOrdinal, Node> &);
//

// stub implementation for EpetraIntMultiVectorT
template <class EpetraGlobalOrdinal, class Node>
class XPETRA_DEPRECATED EpetraIntMultiVectorT
  : public MultiVector<int, int, EpetraGlobalOrdinal, Node> {
  typedef int Scalar;
  typedef int LocalOrdinal;
  typedef EpetraGlobalOrdinal GlobalOrdinal;

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Sets all vector entries to zero.
  EpetraIntMultiVectorT(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t NumVectors, bool zeroOut = true) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::EpetraIntMultiVector only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)");
  }

  //! MultiVector copy constructor.
  EpetraIntMultiVectorT(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &source, const Teuchos::DataAccess copyOrView = Teuchos::Copy) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::EpetraIntMultiVector only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)");
  }

  //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
  EpetraIntMultiVectorT(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &ArrayOfPtrs, size_t NumVectors) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::EpetraIntMultiVector only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)");
  }

  //! Destructor.
  ~EpetraIntMultiVectorT(){};

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Initialize all values in a multi-vector with specified value.
  void putScalar(const int &value) {}

  //! Set multi-vector values to random numbers.
  void randomize(bool bUseXpetraImplementation = true) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::randomize");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntMultiVectorT::randomize(): Functionnality not available in Epetra");
  }

  //! Set multi-vector values to random numbers.
  void randomize(const Scalar &minVal, const Scalar &maxVal, bool bUseXpetraImplementation = true) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::randomize");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntMultiVectorT::randomize(): Functionnality not available in Epetra");
  }

  //! Set seed for Random function.
  /** Note: this method does not exist in Tpetra interface. Added for MueLu. */
  void setSeed(unsigned int seed) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::setSeed");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "Xpetra::EpetraIntMultiVectorT::setSeed(): Functionnality not available in Epetra");
  }

  //@}

  //! @name Data Copy and View get methods
  //@{

  //! Return a Vector which is a const view of column j.
  Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getVector(size_t j) const {
    return Teuchos::null;
  }

  //! Return a Vector which is a nonconst view of column j.
  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getVectorNonConst(size_t j) {
    return Teuchos::null;
  }

  //! Const Local vector access function.
  //! View of the local values in a particular vector of this multi-vector.
  Teuchos::ArrayRCP<const int> getData(size_t j) const {
    return Teuchos::ArrayRCP<const int>();
  }

  //! Local vector access function.
  //! View of the local values in a particular vector of this multi-vector.
  Teuchos::ArrayRCP<int> getDataNonConst(size_t j) {
    return Teuchos::ArrayRCP<int>();
  }

  //@}

  //! @name Mathematical methods
  //@{
  //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
  void dot(const MultiVector<int, int, GlobalOrdinal, Node> &A,
           const Teuchos::ArrayView<int> &dots) const {
    TEUCHOS_TEST_FOR_EXCEPTION(-1, Xpetra::Exceptions::NotImplemented,
                               "This function is not implemented in Epetra_IntMultiVector");
  }

  //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
  void abs(const MultiVector<int, int, GlobalOrdinal, Node> &A) {}

  //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
  void reciprocal(const MultiVector<int, int, GlobalOrdinal, Node> &A) {
    TEUCHOS_TEST_FOR_EXCEPTION(-1, Xpetra::Exceptions::NotImplemented,
                               "This function is not implemented in Epetra_IntMultiVector");
  }

  //! Scale the current values of a multi-vector, this = alpha*this.
  void scale(const int &alpha) {}

  //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
  void scale(Teuchos::ArrayView<const int> alpha) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::scale");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "Xpetra::EpetraIntMultiVectorT::scale(): Functionnality not available in Epetra");
  }

  //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
  void update(const int &alpha, const MultiVector<int, int, GlobalOrdinal, Node> &A, const int &beta) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::update");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "Xpetra::EpetraIntMultiVectorT::update(): Functionnality not available in Epetra");
  }

  //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
  void update(const int &alpha, const MultiVector<int, int, GlobalOrdinal, Node> &A, const int &beta, const MultiVector<int, int, GlobalOrdinal, Node> &B, const int &gamma) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::update");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "Xpetra::EpetraIntMultiVectorT::update(): Functionnality not available in Epetra");
  }

  //! Compute 1-norm of each vector in multi-vector.
  void norm1(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::norm1");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "Xpetra::EpetraIntMultiVectorT::norm1(): Functionnality not available in Epetra");
  }

  //! Compute 2-norm of each vector in multi-vector.
  void norm2(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::norm2");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "Xpetra::EpetraIntMultiVectorT::norm2(): Functionnality not available in Epetra");
  }

  //! Compute Inf-norm of each vector in multi-vector.
  void normInf(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::normInf");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "Xpetra::EpetraIntMultiVectorT::normInf(): Functionnality not available in Epetra");
  }

  //! Compute mean (average) value of each vector in multi-vector.
  void meanValue(const Teuchos::ArrayView<int> &means) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::meanValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "Xpetra::EpetraIntMultiVectorT::meanValue(): Functionnality not available in Epetra");
  }

  //! Compute max value of each vector in multi-vector.
  void maxValue(const Teuchos::ArrayView<int> &maxs) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::maxValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "Xpetra::EpetraIntMultiVectorT::maxValue(): Functionnality not available in Epetra");
  }

  //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
  void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const int &alpha, const MultiVector<int, int, GlobalOrdinal, Node> &A, const MultiVector<int, int, GlobalOrdinal, Node> &B, const int &beta) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::multiply");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "Xpetra::EpetraIntMultiVectorT::multiply(): Functionnality not available in Epetra");
  }

  //! Element-wise multiply of a Vector A with a EpetraMultiVector B.
  void elementWiseMultiply(int scalarAB, const Vector<int, int, GlobalOrdinal, Node> &A, const MultiVector<int, int, GlobalOrdinal, Node> &B, int scalarThis) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::elementWiseMultiply");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "Xpetra_EpetraIntMultiVector: elementWiseMultiply not implemented because Epetra_IntMultiVector does not support this operation");
  }

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Replace value, using global (row) index.
  void replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::replaceGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Add value to existing value, using global (row) index.
  void sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::sumIntoGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Replace value, using local (row) index.
  void replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::replaceLocalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Add value to existing value, using local (row) index.
  void sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::sumIntoLocalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //@}

  //! @name Attribute access functions
  //@{

  //! Returns the number of vectors in the multi-vector.
  size_t getNumVectors() const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::getNumVectors");
    return 1;
  }

  //! Returns the local vector length on the calling processor of vectors in the multi-vector.
  size_t getLocalLength() const { return 0; }

  //! Returns the global vector length of vectors in the multi-vector.
  global_size_t getGlobalLength() const { return 0; }

  //! Checks to see if the local length, number of vectors and size of Scalar type match
  bool isSameSize(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &vec) const { return false; }

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const {
    return std::string("");
  }

  //! Print the object with some verbosity level to an FancyOStream object.
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {}

  //@}

  RCP<Epetra_IntMultiVector> getEpetra_IntMultiVector() const { return Teuchos::null; }

  const RCP<const Comm<int> > getComm() const {
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO getComm Epetra MultiVector not implemented");
  }

  // Implementing DistObject
  Teuchos::RCP<const Map<int, GlobalOrdinal, Node> > getMap() const {
    return Teuchos::null;
  }

  void doImport(const DistObject<int, int, GlobalOrdinal, Node> &source,
                const Import<int, GlobalOrdinal, Node> &importer, CombineMode CM) {}

  void doExport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &dest,
                const Import<int, GlobalOrdinal, Node> &importer, CombineMode CM) {}

  void doImport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &source,
                const Export<int, GlobalOrdinal, Node> &exporter, CombineMode CM) {}

  void doExport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &dest,
                const Export<int, GlobalOrdinal, Node> &exporter, CombineMode CM) {}

  void replaceMap(const RCP<const Map<int, GlobalOrdinal, Node> > &map) {
    // do nothing
  }

 protected:
  /// \brief Implementation of the assignment operator (operator=);
  ///   does a deep copy.
  virtual void
  assign(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &rhs) {}

 private:
  //! The Epetra_IntMultiVector which this class wraps.
  // RCP< Epetra_IntMultiVector > vec_;

};  // class EpetraIntMultiVectorT

// specialization on GO=int and Node=Serial
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
template <>
class EpetraIntMultiVectorT<int, EpetraNode>
  : public virtual MultiVector<int, int, int, EpetraNode> {
  typedef int Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef EpetraNode Node;

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Sets all vector entries to zero.
  EpetraIntMultiVectorT(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t NumVectors, bool zeroOut = true) {
    vec_ = rcp(new Epetra_IntMultiVector(toEpetra<GlobalOrdinal, Node>(map), NumVectors, zeroOut));
  }

  //! MultiVector copy constructor.
  EpetraIntMultiVectorT(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &source, const Teuchos::DataAccess copyOrView = Teuchos::Copy) {
    if (copyOrView == Teuchos::Copy)
      vec_ = rcp(new Epetra_IntMultiVector(toEpetra<GlobalOrdinal, Node>(source)));
    else {
      int *indices = new int[getNumVectors()];
      for (size_t i = 0; i < getNumVectors(); i++)
        indices[i] = i;
      vec_ = Teuchos::rcp(new Epetra_IntMultiVector(View, toEpetra<GlobalOrdinal, Node>(source), indices, getNumVectors()));
      delete[] indices;
    }
  }

  //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
  EpetraIntMultiVectorT(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > & /* map */, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > & /* ArrayOfPtrs */, size_t /* NumVectors */) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::EpetraIntMultiVector only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)");
  }

  //! Destructor.
  ~EpetraIntMultiVectorT(){};

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Initialize all values in a multi-vector with specified value.
  void putScalar(const int &value) {
    int ierr = 0;
    ierr     = vec_->PutScalar(value);
    TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, Xpetra::Exceptions::RuntimeError, "Epetra_IntMultiVector::PutScalar returned an error.")
  }

  //! Set multi-vector values to random numbers.
  void randomize(bool /* bUseXpetraImplementation */ = true) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::randomize");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntMultiVectorT::randomize(): Functionnality not available in Epetra");
  }

  //! Set multi-vector values to random numbers.
  void randomize(const Scalar &minVal, const Scalar &maxVal, bool bUseXpetraImplementation = true) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::randomize");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntMultiVectorT::randomize(): Functionnality not available in Epetra");
  }

  //! Set seed for Random function.
  /** Note: this method does not exist in Tpetra interface. Added for MueLu. */
  void setSeed(unsigned int /* seed */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::setSeed");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntMultiVectorT::setSeed(): Functionnality not available in Epetra");
  }

  typedef typename Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type dual_view_type;

  typename dual_view_type::t_host_const_um getHostLocalView(Access::ReadOnlyStruct) const override { return getHostLocalView(Access::ReadWrite); }

  typename dual_view_type::t_dev_const_um getDeviceLocalView(Access::ReadOnlyStruct) const override { return getDeviceLocalView(Access::ReadWrite); }

  typename dual_view_type::t_host_um getHostLocalView(Access::OverwriteAllStruct) const override { return getHostLocalView(Access::ReadWrite); }

  typename dual_view_type::t_dev_um getDeviceLocalView(Access::OverwriteAllStruct) const override { return getDeviceLocalView(Access::ReadWrite); }

  typename dual_view_type::t_host_um getHostLocalView(Access::ReadWriteStruct) const override {
    typedef Kokkos::View<typename dual_view_type::t_host::data_type,
                         Kokkos::LayoutLeft,
                         typename dual_view_type::t_host::device_type,
                         Kokkos::MemoryUnmanaged>
        epetra_view_type;
    // access Epetra multivector data
    Scalar *data = NULL;
    int myLDA;
    vec_->ExtractView(&data, &myLDA);
    int localLength = vec_->MyLength();
    int numVectors  = getNumVectors();

    // create view
    epetra_view_type test                  = epetra_view_type(data, localLength, numVectors);
    typename dual_view_type::t_host_um ret = subview(test, Kokkos::ALL(), Kokkos::ALL());

    return ret;
  }

  typename dual_view_type::t_dev_um getDeviceLocalView(Access::ReadWriteStruct) const override { return getHostLocalView(Access::ReadWrite); }

  //@}

  //! @name Data Copy and View get methods
  //@{

  //! Return a Vector which is a const view of column j.
  Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getVector(size_t /* j */) const {
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Return a Vector which is a nonconst view of column j.
  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getVectorNonConst(size_t /* j */) {
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Const Local vector access function.
  //! View of the local values in a particular vector of this multi-vector.
  Teuchos::ArrayRCP<const int> getData(size_t j) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::getData");

    int **arrayOfPointers;
    vec_->ExtractView(&arrayOfPointers);
    int *data       = arrayOfPointers[j];
    int localLength = vec_->MyLength();

    return ArrayRCP<int>(data, 0, localLength, false);  // not ownership
  }

  //! Local vector access function.
  //! View of the local values in a particular vector of this multi-vector.
  Teuchos::ArrayRCP<int> getDataNonConst(size_t j) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::getDataNonConst");

    int **arrayOfPointers;
    vec_->ExtractView(&arrayOfPointers);
    int *data       = arrayOfPointers[j];
    int localLength = vec_->MyLength();

    return ArrayRCP<int>(data, 0, localLength, false);  // not ownership
  }

  //@}

  //! @name Mathematical methods
  //@{
  //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
  void dot(const MultiVector<int, int, GlobalOrdinal, Node> & /* A */,
           const Teuchos::ArrayView<int> & /* dots */) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::dot");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "This function is not implemented in Epetra_IntMultiVector");
  }

  //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
  void abs(const MultiVector<int, int, GlobalOrdinal, Node> & /* A */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::abs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "This function is not available in Epetra_IntMultiVector");
  }

  //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
  void reciprocal(const MultiVector<int, int, GlobalOrdinal, Node> & /* A */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::reciprocal");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::RuntimeError, "The reciprocal of an IntMultiVector is not defined!");
  }

  //! Scale the current values of a multi-vector, this = alpha*this.
  void scale(const int & /* alpha */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::scale");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
  void scale(Teuchos::ArrayView<const int> /* alpha */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::scale");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
  void update(const int & /* alpha */, const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const int & /* beta */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::update");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
  void update(const int & /* alpha */, const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const int & /* beta */, const MultiVector<int, int, GlobalOrdinal, Node> & /* B */, const int & /* gamma */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::update");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, B, eB, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute 1-norm of each vector in multi-vector.
  void norm1(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> & /* norms */) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::norm1");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute 2-norm of each vector in multi-vector.
  void norm2(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> & /* norms */) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::norm2");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute Inf-norm of each vector in multi-vector.
  void normInf(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> & /* norms */) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::normInf");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute mean (average) value of each vector in multi-vector.
  void meanValue(const Teuchos::ArrayView<int> & /* means */) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::meanValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute max value of each vector in multi-vector.
  void maxValue(const Teuchos::ArrayView<int> & /* maxs */) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::maxValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
  void multiply(Teuchos::ETransp /* transA */, Teuchos::ETransp /* transB */, const int & /* alpha */, const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const MultiVector<int, int, GlobalOrdinal, Node> & /* B */, const int & /* beta */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::multiply");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Not available in Epetra");
  }

  //! Element-wise multiply of a Vector A with a EpetraMultiVector B.
  void elementWiseMultiply(int /* scalarAB */, const Vector<int, int, GlobalOrdinal, Node> & /* A */, const MultiVector<int, int, GlobalOrdinal, Node> & /* B */, int /* scalarThis */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::elementWiseMultiply");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra_EpetraIntMultiVector: elementWiseMultiply not implemented because Epetra_IntMultiVector does not support this operation");
  }

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Replace value, using global (row) index.
  void replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) {
    vec_->ReplaceGlobalValue(globalRow, vectorIndex, value);
  }

  //! Add value to existing value, using global (row) index.
  void sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) {
    vec_->SumIntoGlobalValue(globalRow, vectorIndex, value);
  }

  //! Replace value, using local (row) index.
  void replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) {
    vec_->ReplaceMyValue(myRow, vectorIndex, value);
  }

  //! Add value to existing value, using local (row) index.
  void sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) {
    vec_->SumIntoMyValue(myRow, vectorIndex, value);
  }

  //@}

  //! @name Attribute access functions
  //@{

  //! Returns the number of vectors in the multi-vector.
  size_t getNumVectors() const {
    return vec_->NumVectors();
  }

  //! Returns the local vector length on the calling processor of vectors in the multi-vector.
  size_t getLocalLength() const {
    return vec_->MyLength();
  }

  //! Returns the global vector length of vectors in the multi-vector.
  global_size_t getGlobalLength() const { return vec_->GlobalLength64(); }

  //! Checks to see if the local length, number of vectors and size of Scalar type match
  bool isSameSize(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &vec) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::isSameSize");
    auto vv = toEpetra<GlobalOrdinal, Node>(vec);
    return ((getLocalLength() == Teuchos::as<size_t>(vv.MyLength())) &&
            (getNumVectors() == Teuchos::as<size_t>(vv.NumVectors())));
  }

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::description");

    // This implementation come from Epetra_Vector_def.hpp (without modification)
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{length=" << this->getGlobalLength()
        << "}";
    return oss.str();
  }

  //! Print the object with some verbosity level to an FancyOStream object.
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::describe");

    // This implementation come from Tpetra_Vector_def.hpp (without modification) // JG: true?
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_EXTREME;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_NONE;

    if (verbLevel > Teuchos::VERB_NONE)
      vec_->Print(out);
  }

  //@}

  RCP<Epetra_IntMultiVector> getEpetra_IntMultiVector() const { return vec_; }

  const RCP<const Comm<int> > getComm() const {
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO getComm Epetra MultiVector not implemented");
  }

  // Implementing DistObject
  Teuchos::RCP<const Map<int, GlobalOrdinal, Node> > getMap() const {
    RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(vec_->Map()));
    return rcp(new Xpetra::EpetraMapT<GlobalOrdinal, Node>(map));
  }

  void doImport(const DistObject<int, int, GlobalOrdinal, Node> &source,
                const Import<int, GlobalOrdinal, Node> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraIntMultiVectorT<GlobalOrdinal XPETRA_COMMA Node>, source, tSource, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraIntMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal XPETRA_COMMA Node>, importer, tImporter, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntMultiVector &v = *tSource.getEpetra_IntMultiVector();
    int err                        = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void doExport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &dest,
                const Import<int, GlobalOrdinal, Node> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraIntMultiVectorT<GlobalOrdinal XPETRA_COMMA Node>, dest, tDest, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraIntMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal XPETRA_COMMA Node>, importer, tImporter, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntMultiVector &v = *tDest.getEpetra_IntMultiVector();
    int err                        = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void doImport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &source,
                const Export<int, GlobalOrdinal, Node> &exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraIntMultiVectorT<GlobalOrdinal XPETRA_COMMA Node>, source, tSource, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraIntMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal XPETRA_COMMA Node>, exporter, tExporter, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntMultiVector &v = *tSource.getEpetra_IntMultiVector();
    int err                        = vec_->Import(v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void doExport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &dest,
                const Export<int, GlobalOrdinal, Node> &exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraIntMultiVectorT<GlobalOrdinal XPETRA_COMMA Node>, dest, tDest, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraIntMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal XPETRA_COMMA Node>, exporter, tExporter, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntMultiVector &v = *tDest.getEpetra_IntMultiVector();
    int err                        = vec_->Export(v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void replaceMap(const RCP<const Map<int, GlobalOrdinal, Node> > &map) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::replaceMap");
    int err = 0;
    if (!map.is_null()) {
      err = this->getEpetra_IntMultiVector()->ReplaceMap(toEpetra<GlobalOrdinal, Node>(map));

    } else {
      // Replace map with a dummy map to avoid potential hangs later
      Epetra_SerialComm SComm;
      Epetra_Map NewMap((GlobalOrdinal)vec_->MyLength(), (GlobalOrdinal)vec_->Map().IndexBase64(), SComm);
      err = this->getEpetra_IntMultiVector()->ReplaceMap(NewMap);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

 protected:
  /// \brief Implementation of the assignment operator (operator=);
  ///   does a deep copy.
  virtual void
  assign(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &rhs) {
    typedef EpetraIntMultiVectorT<GlobalOrdinal, Node> this_type;
    const this_type *rhsPtr = dynamic_cast<const this_type *>(&rhs);
    TEUCHOS_TEST_FOR_EXCEPTION(
        rhsPtr == NULL, std::invalid_argument,
        "Xpetra::MultiVector::operator=: "
        "The left-hand side (LHS) of the assignment has a different type than "
        "the right-hand side (RHS).  The LHS has type Xpetra::EpetraIntMultiVectorT "
        "(which means it wraps an Epetra_IntMultiVector), but the RHS has some "
        "other type.  This probably means that the RHS wraps either an "
        "Tpetra::MultiVector, or an Epetra_MultiVector.  Xpetra::MultiVector "
        "does not currently implement assignment from a Tpetra object to an "
        "Epetra object, though this could be added with sufficient interest.");

    RCP<const Epetra_IntMultiVector> rhsImpl = rhsPtr->getEpetra_IntMultiVector();
    RCP<Epetra_IntMultiVector> lhsImpl       = this->getEpetra_IntMultiVector();

    TEUCHOS_TEST_FOR_EXCEPTION(
        rhsImpl.is_null(), std::logic_error,
        "Xpetra::MultiVector::operator= "
        "(in Xpetra::EpetraIntMultiVectorT::assign): *this (the right-hand side of "
        "the assignment) has a null RCP<Epetra_IntMultiVector> inside.  Please "
        "report this bug to the Xpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        lhsImpl.is_null(), std::logic_error,
        "Xpetra::MultiVector::operator= "
        "(in Xpetra::EpetraIntMultiVectorT::assign): The left-hand side of the "
        "assignment has a null RCP<Epetra_IntMultiVector> inside.  Please report "
        "this bug to the Xpetra developers.");

    // Epetra_IntMultiVector's assignment operator does a deep copy.
    *lhsImpl = *rhsImpl;
  }

 private:
  //! The Epetra_IntMultiVector which this class wraps.
  RCP<Epetra_IntMultiVector> vec_;
};
#endif

// specialization on GO=long long and Node=Serial
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
template <>
class EpetraIntMultiVectorT<long long, EpetraNode>
  : public virtual MultiVector<int, int, long long, EpetraNode> {
  typedef int Scalar;
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Sets all vector entries to zero.
  EpetraIntMultiVectorT(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t NumVectors, bool zeroOut = true) {
    vec_ = rcp(new Epetra_IntMultiVector(toEpetra<GlobalOrdinal, Node>(map), NumVectors, zeroOut));
  }

  //! MultiVector copy constructor.
  EpetraIntMultiVectorT(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &source) {
    vec_ = rcp(new Epetra_IntMultiVector(toEpetra<GlobalOrdinal, Node>(source)));
  }

  //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
  EpetraIntMultiVectorT(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > & /* map */, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > & /* ArrayOfPtrs */, size_t /* NumVectors */) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::EpetraIntMultiVector only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)");
  }

  //! Destructor.
  ~EpetraIntMultiVectorT(){};

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Initialize all values in a multi-vector with specified value.
  void putScalar(const int &value) {
    int ierr = 0;
    ierr     = vec_->PutScalar(value);
    TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, Xpetra::Exceptions::RuntimeError, "Epetra_IntMultiVector::PutScalar returns a non zero error.");
  }

  //! Set multi-vector values to random numbers.
  void randomize(bool /* bUseXpetraImplementation */ = true) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::randomize");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntMultiVectorT::randomize(): Functionnality not available in Epetra");
  }

  //! Set multi-vector values to random numbers.
  void randomize(const Scalar &minVal, const Scalar &maxVal, bool bUseXpetraImplementation = true) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::randomize");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntMultiVectorT::randomize(): Functionnality not available in Epetra");
  }

  //! Set seed for Random function.
  /** Note: this method does not exist in Tpetra interface. Added for MueLu. */
  void setSeed(unsigned int /* seed */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::setSeed");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntMultiVectorT::setSeed(): Functionnality not available in Epetra");
  }

  //@}

  //! @name Data Copy and View get methods
  //@{

  //! Return a Vector which is a const view of column j.
  Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getVector(size_t /* j */) const {
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Return a Vector which is a nonconst view of column j.
  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getVectorNonConst(size_t /* j */) {
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Const Local vector access function.
  //! View of the local values in a particular vector of this multi-vector.
  Teuchos::ArrayRCP<const int> getData(size_t j) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::getData");

    int **arrayOfPointers;
    vec_->ExtractView(&arrayOfPointers);
    int *data       = arrayOfPointers[j];
    int localLength = vec_->MyLength();

    return ArrayRCP<int>(data, 0, localLength, false);  // not ownership
  }

  //! Local vector access function.
  //! View of the local values in a particular vector of this multi-vector.
  Teuchos::ArrayRCP<int> getDataNonConst(size_t j) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::getDataNonConst");

    int **arrayOfPointers;
    vec_->ExtractView(&arrayOfPointers);
    int *data       = arrayOfPointers[j];
    int localLength = vec_->MyLength();

    return ArrayRCP<int>(data, 0, localLength, false);  // not ownership
  }

  //@}

  //! @name Mathematical methods
  //@{
  //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
  void dot(const MultiVector<int, int, GlobalOrdinal, Node> & /* A */,
           const Teuchos::ArrayView<int> & /* dots */) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::dot");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
  void abs(const MultiVector<int, int, GlobalOrdinal, Node> & /* A */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::abs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "This function is not available in Epetra_IntMultiVector");
  }

  //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
  void reciprocal(const MultiVector<int, int, GlobalOrdinal, Node> & /* A */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::reciprocal");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented,
                               "This function is not implemented in Epetra_IntMultiVector");
  }

  //! Scale the current values of a multi-vector, this = alpha*this.
  void scale(const int & /* alpha */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::scale");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
  void scale(Teuchos::ArrayView<const int> /* alpha */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::scale");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
  void update(const int & /* alpha */, const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const int & /* beta */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::update");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
  void update(const int & /* alpha */, const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const int & /* beta */, const MultiVector<int, int, GlobalOrdinal, Node> & /* B */, const int & /* gamma */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::update");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, B, eB, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute 1-norm of each vector in multi-vector.
  void norm1(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> & /* norms */) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::norm1");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute 2-norm of each vector in multi-vector.
  void norm2(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> & /* norms */) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::norm2");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute Inf-norm of each vector in multi-vector.
  void normInf(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> & /* norms */) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::normInf");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute mean (average) value of each vector in multi-vector.
  void meanValue(const Teuchos::ArrayView<int> & /* means */) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::meanValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute max value of each vector in multi-vector.
  void maxValue(const Teuchos::ArrayView<int> & /* maxs */) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::maxValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
  void multiply(Teuchos::ETransp /* transA */, Teuchos::ETransp /* transB */, const int & /* alpha */, const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const MultiVector<int, int, GlobalOrdinal, Node> & /* B */, const int & /* beta */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::multiply");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Not available in Epetra");
  }

  //! Element-wise multiply of a Vector A with a EpetraMultiVector B.
  void elementWiseMultiply(int /* scalarAB */, const Vector<int, int, GlobalOrdinal, Node> & /* A */, const MultiVector<int, int, GlobalOrdinal, Node> & /* B */, int /* scalarThis */) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::elementWiseMultiply");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra_EpetraIntMultiVector: elementWiseMultiply not implemented because Epetra_IntMultiVector does not support this operation");
  }

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Replace value, using global (row) index.
  void replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) {
    vec_->ReplaceGlobalValue(globalRow, vectorIndex, value);
  }

  //! Add value to existing value, using global (row) index.
  void sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) {
    vec_->SumIntoGlobalValue(globalRow, vectorIndex, value);
  }

  //! Replace value, using local (row) index.
  void replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) {
    vec_->ReplaceMyValue(myRow, vectorIndex, value);
  }

  //! Add value to existing value, using local (row) index.
  void sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) {
    vec_->SumIntoMyValue(myRow, vectorIndex, value);
  }

  //@}

  //! @name Attribute access functions
  //@{

  //! Returns the number of vectors in the multi-vector.
  size_t getNumVectors() const {
    return vec_->NumVectors();
  }

  //! Returns the local vector length on the calling processor of vectors in the multi-vector.
  size_t getLocalLength() const { return vec_->MyLength(); }

  //! Returns the global vector length of vectors in the multi-vector.
  global_size_t getGlobalLength() const { return vec_->GlobalLength64(); }

  //! Checks to see if the local length, number of vectors and size of Scalar type match
  bool isSameSize(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &vec) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::isSameSize");
    auto vv = toEpetra<GlobalOrdinal, Node>(vec);
    return ((getLocalLength() == Teuchos::as<size_t>(vv.MyLength())) &&
            (getNumVectors() == Teuchos::as<size_t>(vv.NumVectors())));
  }
  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::description");

    // This implementation come from Epetra_Vector_def.hpp (without modification)
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{length=" << this->getGlobalLength()
        << "}";
    return oss.str();
  }

  //! Print the object with some verbosity level to an FancyOStream object.
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {
    XPETRA_MONITOR("EpetraIntMultiVectorT::describe");

    // This implementation come from Tpetra_Vector_def.hpp (without modification) // JG: true?
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_EXTREME;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_NONE;

    if (verbLevel > Teuchos::VERB_NONE)
      vec_->Print(out);
  }

  //@}

  RCP<Epetra_IntMultiVector> getEpetra_IntMultiVector() const { return vec_; }

  const RCP<const Comm<int> > getComm() const {
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO getComm Epetra MultiVector not implemented");
  }

  // Implementing DistObject
  Teuchos::RCP<const Map<int, GlobalOrdinal, Node> > getMap() const {
    RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(vec_->Map()));
    return rcp(new Xpetra::EpetraMapT<GlobalOrdinal, Node>(map));
  }

  void doImport(const DistObject<int, int, GlobalOrdinal, Node> &source,
                const Import<int, GlobalOrdinal, Node> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraIntMultiVectorT<GlobalOrdinal XPETRA_COMMA Node>, source, tSource, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraIntMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal XPETRA_COMMA Node>, importer, tImporter, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntMultiVector &v = *tSource.getEpetra_IntMultiVector();
    int err                        = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void doExport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &dest,
                const Import<int, GlobalOrdinal, Node> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraIntMultiVectorT<GlobalOrdinal XPETRA_COMMA Node>, dest, tDest, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraIntMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal XPETRA_COMMA Node>, importer, tImporter, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntMultiVector &v = *tDest.getEpetra_IntMultiVector();
    int err                        = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void doImport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &source,
                const Export<int, GlobalOrdinal, Node> &exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraIntMultiVectorT<GlobalOrdinal XPETRA_COMMA Node>, source, tSource, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraIntMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal XPETRA_COMMA Node>, exporter, tExporter, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntMultiVector &v = *tSource.getEpetra_IntMultiVector();
    int err                        = vec_->Import(v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void doExport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &dest,
                const Export<int, GlobalOrdinal, Node> &exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraIntMultiVectorT<GlobalOrdinal XPETRA_COMMA Node>, dest, tDest, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraIntMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal XPETRA_COMMA Node>, exporter, tExporter, "Xpetra::EpetraIntMultiVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntMultiVector &v = *tDest.getEpetra_IntMultiVector();
    int err                        = vec_->Export(v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void replaceMap(const RCP<const Map<int, GlobalOrdinal, Node> > &map) {
    XPETRA_MONITOR("EpetraIntMultiVectorT::replaceMap");
    int err = 0;
    if (!map.is_null()) {
      err = this->getEpetra_IntMultiVector()->ReplaceMap(toEpetra<GlobalOrdinal, Node>(map));

    } else {
      // Replace map with a dummy map to avoid potential hangs later
      Epetra_SerialComm SComm;
      Epetra_Map NewMap((GlobalOrdinal)vec_->MyLength(), (GlobalOrdinal)vec_->Map().IndexBase64(), SComm);
      err = this->getEpetra_IntMultiVector()->ReplaceMap(NewMap);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

 protected:
  /// \brief Implementation of the assignment operator (operator=);
  ///   does a deep copy.
  virtual void
  assign(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &rhs) {
    typedef EpetraIntMultiVectorT<GlobalOrdinal, Node> this_type;
    const this_type *rhsPtr = dynamic_cast<const this_type *>(&rhs);
    TEUCHOS_TEST_FOR_EXCEPTION(
        rhsPtr == NULL, std::invalid_argument,
        "Xpetra::MultiVector::operator=: "
        "The left-hand side (LHS) of the assignment has a different type than "
        "the right-hand side (RHS).  The LHS has type Xpetra::EpetraIntMultiVectorT "
        "(which means it wraps an Epetra_IntMultiVector), but the RHS has some "
        "other type.  This probably means that the RHS wraps either an "
        "Tpetra::MultiVector, or an Epetra_MultiVector.  Xpetra::MultiVector "
        "does not currently implement assignment from a Tpetra object to an "
        "Epetra object, though this could be added with sufficient interest.");

    RCP<const Epetra_IntMultiVector> rhsImpl = rhsPtr->getEpetra_IntMultiVector();
    RCP<Epetra_IntMultiVector> lhsImpl       = this->getEpetra_IntMultiVector();

    TEUCHOS_TEST_FOR_EXCEPTION(
        rhsImpl.is_null(), std::logic_error,
        "Xpetra::MultiVector::operator= "
        "(in Xpetra::EpetraIntMultiVectorT::assign): *this (the right-hand side of "
        "the assignment) has a null RCP<Epetra_IntMultiVector> inside.  Please "
        "report this bug to the Xpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        lhsImpl.is_null(), std::logic_error,
        "Xpetra::MultiVector::operator= "
        "(in Xpetra::EpetraIntMultiVectorT::assign): The left-hand side of the "
        "assignment has a null RCP<Epetra_IntMultiVector> inside.  Please report "
        "this bug to the Xpetra developers.");

    // Epetra_IntMultiVector's assignment operator does a deep copy.
    *lhsImpl = *rhsImpl;
  }

 private:
  //! The Epetra_IntMultiVector which this class wraps.
  RCP<Epetra_IntMultiVector> vec_;
};
#endif

}  // namespace Xpetra

#endif  // XPETRA_EPETRAINTMULTIVECTOR_HPP
