// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_EPETRAINTVECTOR_HPP
#define XPETRA_EPETRAINTVECTOR_HPP

#include "Xpetra_EpetraConfigDefs.hpp"

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_Vector.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_EpetraMap.hpp"
#include "Xpetra_EpetraMultiVector.hpp"
#include "Epetra_IntVector.h"

#if defined(XPETRA_ENABLE_DEPRECATED_CODE)
#ifdef __GNUC__
#if defined(Xpetra_SHOW_DEPRECATED_WARNINGS)
#warning "The header file Trilinos/packages/xpetra/src/Vector/Xpetra_EpetraIntVector.hpp is deprecated."
#endif
#endif
#else
#error "The header file Trilinos/packages/xpetra/src/Vector/Xpetra_EpetraIntVector.hpp is deprecated."
#endif

namespace Xpetra {

// TODO: move that elsewhere
template <class GlobalOrdinal, class Node>
XPETRA_DEPRECATED Epetra_IntVector &toEpetra(Vector<int, int, GlobalOrdinal, Node> &);

template <class GlobalOrdinal, class Node>
const XPETRA_DEPRECATED Epetra_IntVector &toEpetra(const Vector<int, int, GlobalOrdinal, Node> &);
//

// stub implementation for EpetraIntVectorT
template <class EpetraGlobalOrdinal, class Node>
class XPETRA_DEPRECATED EpetraIntVectorT
  : public Vector<int, int, EpetraGlobalOrdinal, Node> {
  typedef int Scalar;
  typedef int LocalOrdinal;
  typedef EpetraGlobalOrdinal GlobalOrdinal;

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Sets all vector entries to zero.
  explicit EpetraIntVectorT(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, bool zeroOut = true) {}

  //! Destructor.
  ~EpetraIntVectorT(){};

  //@}

  //! @name Mathematical methods
  //@{

  //! TODO missing comment
  int dot(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &a) const {
    XPETRA_MONITOR("EpetraIntVectorT::dot");
    TEUCHOS_TEST_FOR_EXCEPTION(-1, Xpetra::Exceptions::NotImplemented, "TODO");
    TEUCHOS_UNREACHABLE_RETURN(-1);
  }

  //! Return 1-norm of this Vector.
  Teuchos::ScalarTraits<int>::magnitudeType norm1() const {
    XPETRA_MONITOR("EpetraIntVectorT::norm1");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
    TEUCHOS_UNREACHABLE_RETURN(-1);
  }

  //! Compute 2-norm of this Vector.
  Teuchos::ScalarTraits<int>::magnitudeType norm2() const {
    XPETRA_MONITOR("EpetraIntVectorT::norm2");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
    TEUCHOS_UNREACHABLE_RETURN(-1);
  }

  //! Compute Inf-norm of this Vector.
  Teuchos::ScalarTraits<int>::magnitudeType normInf() const {
    XPETRA_MONITOR("EpetraIntVectorT::normInf");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
    TEUCHOS_UNREACHABLE_RETURN(-1);
  }

  //! Compute mean (average) value of this Vector.
  int meanValue() const {
    XPETRA_MONITOR("EpetraIntVectorT::meanValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
    TEUCHOS_UNREACHABLE_RETURN(-1);
  }

  //! Compute max value of this Vector.
  int maxValue() const {
    XPETRA_MONITOR("EpetraIntVectorT::maxValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
    TEUCHOS_UNREACHABLE_RETURN(-1);
  }

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Replace current value at the specified location with specified value.
  void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntVectorT::replaceGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Adds specified value to existing value at the specified location.
  void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntVectorT::sumIntoGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Replace current value at the specified location with specified values.
  void replaceLocalValue(LocalOrdinal myRow, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntVectorT::replaceLocalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Adds specified value to existing value at the specified location.
  void sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntVectorT::sumIntoLocalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Initialize all values in a multi-vector with specified value.
  void putScalar(const int &value) {}

  //! Set multi-vector values to random numbers.
  void randomize(bool bUseXpetraImplementation = true) {
    XPETRA_MONITOR("EpetraIntVectorT::randomize");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVectorT::randomize(): Functionnality not available in Epetra");
  }

  //! Set multi-vector values to random numbers.
  void randomize(const Scalar & /*minVal*/, const Scalar & /*maxVal*/, bool bUseXpetraImplementation = true) {
    XPETRA_MONITOR("EpetraIntVectorT::randomize");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVectorT::randomize(): Functionnality not available in Epetra");
  }

  //! Set seed for Random function.
  /** Note: this method does not exist in Tpetra interface. Added for MueLu. */
  void setSeed(unsigned int seed) {
    XPETRA_MONITOR("EpetraIntVectorT::setSeed");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVectorT::setSeed(): Functionnality not available in Epetra");
  }

  //@}

  //! @name Data Copy and View get methods
  //@{

  //! Return a Vector which is a const view of column j.
  Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getVector(size_t j) const {
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Return a Vector which is a nonconst view of column j.
  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getVectorNonConst(size_t j) {
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Const Local vector access function.
  //! View of the local values in a particular vector of this multi-vector.
  Teuchos::ArrayRCP<const int> getData(size_t j) const { return Teuchos::ArrayRCP<const int>(); }

  //! Local vector access function.
  //! View of the local values in a particular vector of this multi-vector.
  Teuchos::ArrayRCP<int> getDataNonConst(size_t j) { return Teuchos::ArrayRCP<int>(); }

  //@}

  //! @name Mathematical methods
  //@{
  //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
  void dot(const MultiVector<int, int, GlobalOrdinal, Node> &A, const Teuchos::ArrayView<int> &dots) const {}

  //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
  void abs(const MultiVector<int, int, GlobalOrdinal, Node> &A) {}

  //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
  void reciprocal(const MultiVector<int, int, GlobalOrdinal, Node> &A) {}

  //! Scale the current values of a multi-vector, this = alpha*this.
  void scale(const int &alpha) {}

  //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
  void scale(Teuchos::ArrayView<const int> alpha) {
    XPETRA_MONITOR("EpetraIntVectorT::scale");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
  void update(const int &alpha, const MultiVector<int, int, GlobalOrdinal, Node> &A, const int &beta) {
    XPETRA_MONITOR("EpetraIntVectorT::update");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
  void update(const int &alpha, const MultiVector<int, int, GlobalOrdinal, Node> &A, const int &beta, const MultiVector<int, int, GlobalOrdinal, Node> &B, const int &gamma) {
    XPETRA_MONITOR("EpetraIntVectorT::update");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, B, eB, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute 1-norm of each vector in multi-vector.
  void norm1(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const {
    XPETRA_MONITOR("EpetraIntVectorT::norm1");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute 2-norm of each vector in multi-vector.
  void norm2(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const {
    XPETRA_MONITOR("EpetraIntVectorT::norm2");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute Inf-norm of each vector in multi-vector.
  void normInf(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const {
    XPETRA_MONITOR("EpetraIntVectorT::normInf");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute mean (average) value of each vector in multi-vector.
  void meanValue(const Teuchos::ArrayView<int> &means) const {
    XPETRA_MONITOR("EpetraIntVectorT::meanValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute max value of each vector in multi-vector.
  void maxValue(const Teuchos::ArrayView<int> &maxs) const {
    XPETRA_MONITOR("EpetraIntVectorT::maxValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
  void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const int &alpha, const MultiVector<int, int, GlobalOrdinal, Node> &A, const MultiVector<int, int, GlobalOrdinal, Node> &B, const int &beta) {
    XPETRA_MONITOR("EpetraIntVectorT::multiply");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Not available in Epetra");
  }

  //! Element-wise multiply of a Vector A with a EpetraMultiVector B.
  void elementWiseMultiply(int scalarAB, const Vector<int, int, GlobalOrdinal, Node> &A, const MultiVector<int, int, GlobalOrdinal, Node> &B, int scalarThis) {
    XPETRA_MONITOR("EpetraIntVectorT::elementWiseMultiply");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra_EpetraIntVector: elementWiseMultiply not implemented because Epetra_IntVector does not support this operation");
  }

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Replace value, using global (row) index.
  void replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntVectorT::replaceGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Add value to existing value, using global (row) index.
  void sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntVectorT::sumIntoGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Replace value, using local (row) index.
  void replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntVectorT::replaceLocalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Add value to existing value, using local (row) index.
  void sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntVectorT::sumIntoLocalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //@}

  //! @name Attribute access functions
  //@{

  //! Returns the number of vectors in the multi-vector.
  size_t getNumVectors() const {
    XPETRA_MONITOR("EpetraIntVectorT::getNumVectors");
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

  RCP<Epetra_IntVector> getEpetra_IntVector() const { return Teuchos::null; }

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
  //! The Epetra_IntVector which this class wraps.
  // RCP< Epetra_IntVector > vec_;

};  // class EpetraIntVectorT

// specialization on GO=int and Node=Serial
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
template <>
class EpetraIntVectorT<int, EpetraNode>
  : public virtual Vector<int, int, int, EpetraNode> {
  typedef int Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef EpetraNode Node;

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Sets all vector entries to zero.
  explicit EpetraIntVectorT(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, bool zeroOut = true) {
    vec_ = rcp(new Epetra_IntVector(toEpetra<GlobalOrdinal, Node>(map), zeroOut));
  }

  //! Destructor.
  ~EpetraIntVectorT(){};

  //@}

  //! @name Mathematical methods
  //@{

  //! TODO missing comment
  int dot(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & /* a */) const {
    XPETRA_MONITOR("EpetraIntVectorT::dot");
    TEUCHOS_TEST_FOR_EXCEPTION(-1, Xpetra::Exceptions::NotImplemented, "TODO");
    TEUCHOS_UNREACHABLE_RETURN(-1);
  }

  //! Return 1-norm of this Vector.
  Teuchos::ScalarTraits<int>::magnitudeType norm1() const {
    XPETRA_MONITOR("EpetraIntVectorT::norm1");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
    TEUCHOS_UNREACHABLE_RETURN(-1);
  }

  //! Compute 2-norm of this Vector.
  Teuchos::ScalarTraits<int>::magnitudeType norm2() const {
    XPETRA_MONITOR("EpetraIntVectorT::norm2");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
    TEUCHOS_UNREACHABLE_RETURN(-1);
  }

  //! Compute Inf-norm of this Vector.
  Teuchos::ScalarTraits<int>::magnitudeType normInf() const {
    XPETRA_MONITOR("EpetraIntVectorT::normInf");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
    TEUCHOS_UNREACHABLE_RETURN(-1);
  }

  //! Compute mean (average) value of this Vector.
  int meanValue() const {
    XPETRA_MONITOR("EpetraIntVectorT::meanValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
    TEUCHOS_UNREACHABLE_RETURN(-1);
  }

  //! Compute max value of this Vector.
  int maxValue() const {
    XPETRA_MONITOR("EpetraIntVectorT::maxValue");
    return vec_->MaxValue();
  }

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Replace current value at the specified location with specified value.
  void replaceGlobalValue(GlobalOrdinal /* globalRow */, const Scalar & /* value */) {
    XPETRA_MONITOR("EpetraIntVectorT::replaceGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Adds specified value to existing value at the specified location.
  void sumIntoGlobalValue(GlobalOrdinal /* globalRow */, const Scalar & /* value */) {
    XPETRA_MONITOR("EpetraIntVectorT::sumIntoGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Replace current value at the specified location with specified values.
  void replaceLocalValue(LocalOrdinal myRow, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntVectorT::replaceLocalValue");
    (*vec_)[myRow] = value;
  }

  //! Adds specified value to existing value at the specified location.
  void sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntVectorT::sumIntoLocalValue");
    (*vec_)[myRow] += value;
  }

  //! Initialize all values in a multi-vector with specified value.
  void putScalar(const int &value) { vec_->PutValue(value); }

  //! Set multi-vector values to random numbers.
  void randomize(bool /* bUseXpetraImplementation */ = true) {
    XPETRA_MONITOR("EpetraIntVectorT::randomize");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVectorT::randomize(): Functionnality not available in Epetra");
  }

  //! Set multi-vector values to random numbers.
  void randomize(const Scalar & /*minVal*/, const Scalar & /*maxVal*/, bool /* bUseXpetraImplementation */ = true) {
    XPETRA_MONITOR("EpetraIntVectorT::randomize");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVectorT::randomize(): Functionnality not available in Epetra");
  }

  //! Set seed for Random function.
  /** Note: this method does not exist in Tpetra interface. Added for MueLu. */
  void setSeed(unsigned int /* seed */) {
    XPETRA_MONITOR("EpetraIntVectorT::setSeed");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVectorT::setSeed(): Functionnality not available in Epetra");
  }

  typedef typename Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type dual_view_type;

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
    // access Epetra vector data
    Scalar *data = NULL;
    vec_->ExtractView(&data);
    int localLength = vec_->MyLength();

    // create view
    epetra_view_type test                  = epetra_view_type(data, localLength, 1);
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
  Teuchos::ArrayRCP<const int> getData(size_t /* j */) const {
    XPETRA_MONITOR("EpetraIntVectorT::getData");

    int *data       = vec_->Values();
    int localLength = vec_->MyLength();

    return ArrayRCP<int>(data, 0, localLength, false);  // not ownership
  }

  //! Local vector access function.
  //! View of the local values in a particular vector of this multi-vector.
  Teuchos::ArrayRCP<int> getDataNonConst(size_t /* j */) {
    XPETRA_MONITOR("EpetraIntVectorT::getDataNonConst");

    int *data       = vec_->Values();
    int localLength = vec_->MyLength();

    return ArrayRCP<int>(data, 0, localLength, false);  // not ownership
  }

  //@}

  //! @name Mathematical methods
  //@{
  //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
  void dot(const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const Teuchos::ArrayView<int> & /* dots */) const {
    XPETRA_MONITOR("EpetraIntVectorT::dot");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
  void abs(const MultiVector<int, int, GlobalOrdinal, Node> & /* A */) {
    XPETRA_MONITOR("EpetraIntVectorT::abs");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
  void reciprocal(const MultiVector<int, int, GlobalOrdinal, Node> & /* A */) {
    XPETRA_MONITOR("EpetraIntVectorT::reciprocal");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Scale the current values of a multi-vector, this = alpha*this.
  void scale(const int & /* alpha */) {
    XPETRA_MONITOR("EpetraIntVectorT::scale");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
  void scale(Teuchos::ArrayView<const int> /* alpha */) {
    XPETRA_MONITOR("EpetraIntVectorT::scale");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
  void update(const int & /* alpha */, const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const int & /* beta */) {
    XPETRA_MONITOR("EpetraIntVectorT::update");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
  void update(const int & /* alpha */, const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const int & /* beta */, const MultiVector<int, int, GlobalOrdinal, Node> & /* B */, const int & /* gamma */) {
    XPETRA_MONITOR("EpetraIntVectorT::update");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, B, eB, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute 1-norm of each vector in multi-vector.
  void norm1(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> & /* norms */) const {
    XPETRA_MONITOR("EpetraIntVectorT::norm1");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute 2-norm of each vector in multi-vector.
  void norm2(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> & /* norms */) const {
    XPETRA_MONITOR("EpetraIntVectorT::norm2");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute Inf-norm of each vector in multi-vector.
  void normInf(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> & /* norms */) const {
    XPETRA_MONITOR("EpetraIntVectorT::normInf");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute mean (average) value of each vector in multi-vector.
  void meanValue(const Teuchos::ArrayView<int> & /* means */) const {
    XPETRA_MONITOR("EpetraIntVectorT::meanValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute max value of each vector in multi-vector.
  void maxValue(const Teuchos::ArrayView<int> & /* maxs */) const {
    XPETRA_MONITOR("EpetraIntVectorT::maxValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
  void multiply(Teuchos::ETransp /* transA */, Teuchos::ETransp /* transB */, const int & /* alpha */, const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const MultiVector<int, int, GlobalOrdinal, Node> & /* B */, const int & /* beta */) {
    XPETRA_MONITOR("EpetraIntVectorT::multiply");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Not available in Epetra");
  }

  //! Element-wise multiply of a Vector A with a EpetraMultiVector B.
  void elementWiseMultiply(int /* scalarAB */, const Vector<int, int, GlobalOrdinal, Node> & /* A */, const MultiVector<int, int, GlobalOrdinal, Node> & /* B */, int /* scalarThis */) {
    XPETRA_MONITOR("EpetraIntVectorT::elementWiseMultiply");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra_EpetraIntVector: elementWiseMultiply not implemented because Epetra_IntVector does not support this operation");
  }

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Replace value, using global (row) index.
  void replaceGlobalValue(GlobalOrdinal /* globalRow */, size_t /* vectorIndex */, const Scalar & /* value */) {
    XPETRA_MONITOR("EpetraIntVectorT::replaceGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Add value to existing value, using global (row) index.
  void sumIntoGlobalValue(GlobalOrdinal /* globalRow */, size_t /* vectorIndex */, const Scalar & /* value */) {
    XPETRA_MONITOR("EpetraIntVectorT::sumIntoGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Replace value, using local (row) index.
  void replaceLocalValue(LocalOrdinal /* myRow */, size_t /* vectorIndex */, const Scalar & /* value */) {
    XPETRA_MONITOR("EpetraIntVectorT::replaceLocalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Add value to existing value, using local (row) index.
  void sumIntoLocalValue(LocalOrdinal /* myRow */, size_t /* vectorIndex */, const Scalar & /* value */) {
    XPETRA_MONITOR("EpetraIntVectorT::sumIntoLocalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //@}

  //! @name Attribute access functions
  //@{

  //! Returns the number of vectors in the multi-vector.
  size_t getNumVectors() const {
    XPETRA_MONITOR("EpetraIntVectorT::getNumVectors");
    return 1;
  }

  //! Returns the local vector length on the calling processor of vectors in the multi-vector.
  size_t getLocalLength() const { return vec_->MyLength(); }

  //! Returns the global vector length of vectors in the multi-vector.
  global_size_t getGlobalLength() const { return vec_->GlobalLength64(); }

  //! Checks to see if the local length, number of vectors and size of Scalar type match
  bool isSameSize(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &vec) const {
    XPETRA_MONITOR("EpetraIntVectorT::isSameSize");
    const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> *asvec = dynamic_cast<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> *>(&vec);
    if (!asvec) return false;
    auto vv = toEpetra(*asvec);
    return ((vec_->MyLength() == vv.MyLength()) && (getNumVectors() == vec.getNumVectors()));
  }

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const {
    XPETRA_MONITOR("EpetraIntVectorT::description");

    // This implementation come from Epetra_Vector_def.hpp (without modification)
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{length=" << this->getGlobalLength()
        << "}";
    return oss.str();
  }

  //! Print the object with some verbosity level to an FancyOStream object.
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {
    XPETRA_MONITOR("EpetraIntVectorT::describe");

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

  RCP<Epetra_IntVector> getEpetra_IntVector() const { return vec_; }

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
    XPETRA_MONITOR("EpetraIntVectorT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVectorT<GlobalOrdinal XPETRA_COMMA Node>, source, tSource, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraIntVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal XPETRA_COMMA Node>, importer, tImporter, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntVector &v = *tSource.getEpetra_IntVector();
    int err                   = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void doExport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &dest,
                const Import<int, GlobalOrdinal, Node> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVectorT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVectorT<GlobalOrdinal XPETRA_COMMA Node>, dest, tDest, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraIntVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal XPETRA_COMMA Node>, importer, tImporter, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntVector &v = *tDest.getEpetra_IntVector();
    int err                   = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void doImport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &source,
                const Export<int, GlobalOrdinal, Node> &exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVectorT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVectorT<GlobalOrdinal XPETRA_COMMA Node>, source, tSource, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraIntVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal XPETRA_COMMA Node>, exporter, tExporter, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntVector &v = *tSource.getEpetra_IntVector();
    int err                   = vec_->Import(v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void doExport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &dest,
                const Export<int, GlobalOrdinal, Node> &exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVectorT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVectorT<GlobalOrdinal XPETRA_COMMA Node>, dest, tDest, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraIntVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal XPETRA_COMMA Node>, exporter, tExporter, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntVector &v = *tDest.getEpetra_IntVector();
    int err                   = vec_->Export(v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void replaceMap(const RCP<const Map<int, GlobalOrdinal, Node> > & /* map */) {
    // do nothing
  }

 protected:
  /// \brief Implementation of the assignment operator (operator=);
  ///   does a deep copy.
  virtual void
  assign(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &rhs) {
    typedef EpetraIntVectorT<GlobalOrdinal, Node> this_type;
    const this_type *rhsPtr = dynamic_cast<const this_type *>(&rhs);
    TEUCHOS_TEST_FOR_EXCEPTION(
        rhsPtr == NULL, std::invalid_argument,
        "Xpetra::MultiVector::operator=: "
        "The left-hand side (LHS) of the assignment has a different type than "
        "the right-hand side (RHS).  The LHS has type Xpetra::EpetraIntVectorT "
        "(which means it wraps an Epetra_IntVector), but the RHS has some "
        "other type.  This probably means that the RHS wraps either an "
        "Tpetra::MultiVector, or an Epetra_MultiVector.  Xpetra::MultiVector "
        "does not currently implement assignment from a Tpetra object to an "
        "Epetra object, though this could be added with sufficient interest.");

    RCP<const Epetra_IntVector> rhsImpl = rhsPtr->getEpetra_IntVector();
    RCP<Epetra_IntVector> lhsImpl       = this->getEpetra_IntVector();

    TEUCHOS_TEST_FOR_EXCEPTION(
        rhsImpl.is_null(), std::logic_error,
        "Xpetra::MultiVector::operator= "
        "(in Xpetra::EpetraIntVectorT::assign): *this (the right-hand side of "
        "the assignment) has a null RCP<Epetra_IntVector> inside.  Please "
        "report this bug to the Xpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        lhsImpl.is_null(), std::logic_error,
        "Xpetra::MultiVector::operator= "
        "(in Xpetra::EpetraIntVectorT::assign): The left-hand side of the "
        "assignment has a null RCP<Epetra_IntVector> inside.  Please report "
        "this bug to the Xpetra developers.");

    // Epetra_IntVector's assignment operator does a deep copy.
    *lhsImpl = *rhsImpl;
  }

 private:
  //! The Epetra_IntVector which this class wraps.
  RCP<Epetra_IntVector> vec_;
};
#endif

// specialization on GO=long long and Node=Serial
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
template <>
class EpetraIntVectorT<long long, EpetraNode>
  : public virtual Vector<int, int, long long, EpetraNode> {
  typedef int Scalar;
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Sets all vector entries to zero.
  explicit EpetraIntVectorT(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, bool zeroOut = true) {
    vec_ = rcp(new Epetra_IntVector(toEpetra<GlobalOrdinal, Node>(map), zeroOut));
  }

  //! Destructor.
  ~EpetraIntVectorT(){};

  //@}

  //! @name Mathematical methods
  //@{

  //! TODO missing comment
  int dot(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & /* a */) const {
    XPETRA_MONITOR("EpetraIntVectorT::dot");
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::NotImplemented, "TODO"); /* return -1; */
  }

  //! Return 1-norm of this Vector.
  Teuchos::ScalarTraits<int>::magnitudeType norm1() const {
    XPETRA_MONITOR("EpetraIntVectorT::norm1");
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::NotImplemented, "TODO"); /* return -1; */
  }

  //! Compute 2-norm of this Vector.
  Teuchos::ScalarTraits<int>::magnitudeType norm2() const {
    XPETRA_MONITOR("EpetraIntVectorT::norm2");
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::NotImplemented, "TODO"); /* return -1; */
  }

  //! Compute Inf-norm of this Vector.
  Teuchos::ScalarTraits<int>::magnitudeType normInf() const {
    XPETRA_MONITOR("EpetraIntVectorT::normInf");
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::NotImplemented, "TODO"); /* return -1; */
  }

  //! Compute mean (average) value of this Vector.
  int meanValue() const {
    XPETRA_MONITOR("EpetraIntVectorT::meanValue");
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::NotImplemented, "TODO"); /* return -1; */
  }

  //! Compute max value of this Vector.
  int maxValue() const {
    XPETRA_MONITOR("EpetraIntVectorT::maxValue");
    return Teuchos::as<int>(vec_->MaxValue());
  }

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Replace current value at the specified location with specified value.
  void replaceGlobalValue(GlobalOrdinal /* globalRow */, const Scalar & /* value */) {
    XPETRA_MONITOR("EpetraIntVectorT::replaceGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Adds specified value to existing value at the specified location.
  void sumIntoGlobalValue(GlobalOrdinal /* globalRow */, const Scalar & /* value */) {
    XPETRA_MONITOR("EpetraIntVectorT::sumIntoGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Replace current value at the specified location with specified values.
  void replaceLocalValue(LocalOrdinal myRow, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntVectorT::replaceLocalValue");
    (*vec_)[myRow] = value;
  }

  //! Adds specified value to existing value at the specified location.
  void sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value) {
    XPETRA_MONITOR("EpetraIntVectorT::sumIntoLocalValue");
    (*vec_)[myRow] += value;
  }

  //! Initialize all values in a multi-vector with specified value.
  void putScalar(const int &value) { vec_->PutValue(value); }

  //! Set multi-vector values to random numbers.
  void randomize(bool /* bUseXpetraImplementation */ = true) {
    XPETRA_MONITOR("EpetraIntVectorT::randomize");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVectorT::randomize(): Functionnality not available in Epetra");
  }

  //! Set multi-vector values to random numbers.
  void randomize(const Scalar & /*minVal*/, const Scalar & /*maxVal*/, bool /* bUseXpetraImplementation */ = true) {
    XPETRA_MONITOR("EpetraIntVectorT::randomize");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVectorT::randomize(): Functionnality not available in Epetra");
  }

  //! Set seed for Random function.
  /** Note: this method does not exist in Tpetra interface. Added for MueLu. */
  void setSeed(unsigned int /* seed */) {
    XPETRA_MONITOR("EpetraIntVectorT::setSeed");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVectorT::setSeed(): Functionnality not available in Epetra");
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
  Teuchos::ArrayRCP<const int> getData(size_t /* j */) const {
    XPETRA_MONITOR("EpetraIntVectorT::getData");

    int *data       = vec_->Values();
    int localLength = vec_->MyLength();

    return ArrayRCP<int>(data, 0, localLength, false);  // not ownership
  }

  //! Local vector access function.
  //! View of the local values in a particular vector of this multi-vector.
  Teuchos::ArrayRCP<int> getDataNonConst(size_t /* j */) {
    XPETRA_MONITOR("EpetraIntVectorT::getDataNonConst");

    int *data       = vec_->Values();
    int localLength = vec_->MyLength();

    return ArrayRCP<int>(data, 0, localLength, false);  // not ownership
  }

  //@}

  //! @name Mathematical methods
  //@{
  //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
  void dot(const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const Teuchos::ArrayView<int> & /* dots */) const {
    XPETRA_MONITOR("EpetraIntVectorT::dot");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
  void abs(const MultiVector<int, int, GlobalOrdinal, Node> & /* A */) {
    XPETRA_MONITOR("EpetraIntVectorT::abs");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
  void reciprocal(const MultiVector<int, int, GlobalOrdinal, Node> & /* A */) {
    XPETRA_MONITOR("EpetraIntVectorT::reciprocal");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Scale the current values of a multi-vector, this = alpha*this.
  void scale(const int & /* alpha */) {
    XPETRA_MONITOR("EpetraIntVectorT::scale");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
  void scale(Teuchos::ArrayView<const int> /* alpha */) {
    XPETRA_MONITOR("EpetraIntVectorT::scale");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
  void update(const int & /* alpha */, const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const int & /* beta */) {
    XPETRA_MONITOR("EpetraIntVectorT::update");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
  void update(const int & /* alpha */, const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const int & /* beta */, const MultiVector<int, int, GlobalOrdinal, Node> & /* B */, const int & /* gamma */) {
    XPETRA_MONITOR("EpetraIntVectorT::update");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, B, eB, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute 1-norm of each vector in multi-vector.
  void norm1(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> & /* norms */) const {
    XPETRA_MONITOR("EpetraIntVectorT::norm1");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute 2-norm of each vector in multi-vector.
  void norm2(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> & /* norms */) const {
    XPETRA_MONITOR("EpetraIntVectorT::norm2");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute Inf-norm of each vector in multi-vector.
  void normInf(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> & /* norms */) const {
    XPETRA_MONITOR("EpetraIntVectorT::normInf");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute mean (average) value of each vector in multi-vector.
  void meanValue(const Teuchos::ArrayView<int> & /* means */) const {
    XPETRA_MONITOR("EpetraIntVectorT::meanValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Compute max value of each vector in multi-vector.
  void maxValue(const Teuchos::ArrayView<int> & /* maxs */) const {
    XPETRA_MONITOR("EpetraIntVectorT::maxValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
  void multiply(Teuchos::ETransp /* transA */, Teuchos::ETransp /* transB */, const int & /* alpha */, const MultiVector<int, int, GlobalOrdinal, Node> & /* A */, const MultiVector<int, int, GlobalOrdinal, Node> & /* B */, const int & /* beta */) {
    XPETRA_MONITOR("EpetraIntVectorT::multiply");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Not available in Epetra");
  }

  //! Element-wise multiply of a Vector A with a EpetraMultiVector B.
  void elementWiseMultiply(int /* scalarAB */, const Vector<int, int, GlobalOrdinal, Node> & /* A */, const MultiVector<int, int, GlobalOrdinal, Node> & /* B */, int /* scalarThis */) {
    XPETRA_MONITOR("EpetraIntVectorT::elementWiseMultiply");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra_EpetraIntVector: elementWiseMultiply not implemented because Epetra_IntVector does not support this operation");
  }

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Replace value, using global (row) index.
  void replaceGlobalValue(GlobalOrdinal /* globalRow */, size_t /* vectorIndex */, const Scalar & /* value */) {
    XPETRA_MONITOR("EpetraIntVectorT::replaceGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Add value to existing value, using global (row) index.
  void sumIntoGlobalValue(GlobalOrdinal /* globalRow */, size_t /* vectorIndex */, const Scalar & /* value */) {
    XPETRA_MONITOR("EpetraIntVectorT::sumIntoGlobalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Replace value, using local (row) index.
  void replaceLocalValue(LocalOrdinal /* myRow */, size_t /* vectorIndex */, const Scalar & /* value */) {
    XPETRA_MONITOR("EpetraIntVectorT::replaceLocalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //! Add value to existing value, using local (row) index.
  void sumIntoLocalValue(LocalOrdinal /* myRow */, size_t /* vectorIndex */, const Scalar & /* value */) {
    XPETRA_MONITOR("EpetraIntVectorT::sumIntoLocalValue");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  //@}

  //! @name Attribute access functions
  //@{

  //! Returns the number of vectors in the multi-vector.
  size_t getNumVectors() const {
    XPETRA_MONITOR("EpetraIntVectorT::getNumVectors");
    return 1;
  }

  //! Returns the local vector length on the calling processor of vectors in the multi-vector.
  size_t getLocalLength() const { return vec_->MyLength(); }

  //! Returns the global vector length of vectors in the multi-vector.
  global_size_t getGlobalLength() const { return vec_->GlobalLength64(); }

  //! Checks to see if the local length, number of vectors and size of Scalar type match
  bool isSameSize(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &vec) const {
    XPETRA_MONITOR("EpetraIntVectorT::isSameSize");
    const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> *asvec = dynamic_cast<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> *>(&vec);
    if (!asvec) return false;
    auto vv = toEpetra(*asvec);
    return ((vec_->MyLength() == vv.MyLength()) && (getNumVectors() == vec.getNumVectors()));
  }
  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const {
    XPETRA_MONITOR("EpetraIntVectorT::description");

    // This implementation come from Epetra_Vector_def.hpp (without modification)
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{length=" << this->getGlobalLength()
        << "}";
    return oss.str();
  }

  //! Print the object with some verbosity level to an FancyOStream object.
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {
    XPETRA_MONITOR("EpetraIntVectorT::describe");

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

  RCP<Epetra_IntVector> getEpetra_IntVector() const { return vec_; }

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
    XPETRA_MONITOR("EpetraIntVectorT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVectorT<GlobalOrdinal XPETRA_COMMA Node>, source, tSource, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraIntVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal XPETRA_COMMA Node>, importer, tImporter, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntVector &v = *tSource.getEpetra_IntVector();
    int err                   = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void doExport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &dest,
                const Import<int, GlobalOrdinal, Node> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVectorT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVectorT<GlobalOrdinal XPETRA_COMMA Node>, dest, tDest, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraIntVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal XPETRA_COMMA Node>, importer, tImporter, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntVector &v = *tDest.getEpetra_IntVector();
    int err                   = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void doImport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &source,
                const Export<int, GlobalOrdinal, Node> &exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVectorT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVectorT<GlobalOrdinal XPETRA_COMMA Node>, source, tSource, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraIntVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal XPETRA_COMMA Node>, exporter, tExporter, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntVector &v = *tSource.getEpetra_IntVector();
    int err                   = vec_->Import(v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void doExport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &dest,
                const Export<int, GlobalOrdinal, Node> &exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVectorT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVectorT<GlobalOrdinal XPETRA_COMMA Node>, dest, tDest, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraIntVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal XPETRA_COMMA Node>, exporter, tExporter, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntVector &v = *tDest.getEpetra_IntVector();
    int err                   = vec_->Export(v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void replaceMap(const RCP<const Map<int, GlobalOrdinal, Node> > & /* map */) {
    // do nothing
  }

 protected:
  /// \brief Implementation of the assignment operator (operator=);
  ///   does a deep copy.
  virtual void
  assign(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &rhs) {
    typedef EpetraIntVectorT<GlobalOrdinal, Node> this_type;
    const this_type *rhsPtr = dynamic_cast<const this_type *>(&rhs);
    TEUCHOS_TEST_FOR_EXCEPTION(
        rhsPtr == NULL, std::invalid_argument,
        "Xpetra::MultiVector::operator=: "
        "The left-hand side (LHS) of the assignment has a different type than "
        "the right-hand side (RHS).  The LHS has type Xpetra::EpetraIntVectorT "
        "(which means it wraps an Epetra_IntVector), but the RHS has some "
        "other type.  This probably means that the RHS wraps either an "
        "Tpetra::MultiVector, or an Epetra_MultiVector.  Xpetra::MultiVector "
        "does not currently implement assignment from a Tpetra object to an "
        "Epetra object, though this could be added with sufficient interest.");

    RCP<const Epetra_IntVector> rhsImpl = rhsPtr->getEpetra_IntVector();
    RCP<Epetra_IntVector> lhsImpl       = this->getEpetra_IntVector();

    TEUCHOS_TEST_FOR_EXCEPTION(
        rhsImpl.is_null(), std::logic_error,
        "Xpetra::MultiVector::operator= "
        "(in Xpetra::EpetraIntVectorT::assign): *this (the right-hand side of "
        "the assignment) has a null RCP<Epetra_IntVector> inside.  Please "
        "report this bug to the Xpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
        lhsImpl.is_null(), std::logic_error,
        "Xpetra::MultiVector::operator= "
        "(in Xpetra::EpetraIntVectorT::assign): The left-hand side of the "
        "assignment has a null RCP<Epetra_IntVector> inside.  Please report "
        "this bug to the Xpetra developers.");

    // Epetra_IntVector's assignment operator does a deep copy.
    *lhsImpl = *rhsImpl;
  }

 private:
  //! The Epetra_IntVector which this class wraps.
  RCP<Epetra_IntVector> vec_;
};
#endif

}  // namespace Xpetra

#endif  // XPETRA_EPETRAINTVECTOR_HPP
