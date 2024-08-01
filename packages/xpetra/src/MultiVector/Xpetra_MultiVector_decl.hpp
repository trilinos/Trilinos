// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_MULTIVECTOR_DECL_HPP
#define XPETRA_MULTIVECTOR_DECL_HPP

#include <Teuchos_LabeledObject.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_Range1D.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DistObject.hpp"
#include "Xpetra_Map_decl.hpp"

#include "Xpetra_Access.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_ArithTraits.hpp>

namespace Xpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// forward declaration of Vector, needed to prevent circular inclusions
template <class S, class LO, class GO, class N>
class Vector;
#endif

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class MultiVector
  : public DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  typedef Scalar scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

  //! @name Constructor/Destructor Methods
  //@{

  //! Destructor.
  virtual ~MultiVector();

  /// \brief Assignment operator: Does a deep copy.
  ///
  /// The assignment operator does a deep copy, just like
  /// subclasses' copy constructors.
  ///
  /// \note This currently only works if both <tt>*this</tt> and the
  ///   input argument are instances of the same subclass.  We do
  ///   not currently allow assignment between an
  ///   Xpetra::TpetraMultiVector and an Xpetra::EpetraMultiVector,
  ///   or vice versa, for example.
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
  operator=(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs);

  //@}
  //! @name Post-construction modification routines
  //@{

  //! Replace value, using global (row) index.
  virtual void replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar& value) = 0;

  //! Add value to existing value, using global (row) index.
  virtual void sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar& value) = 0;

  //! Replace value, using local (row) index.
  virtual void replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar& value) = 0;

  //! Add value to existing value, using local (row) index.
  virtual void sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar& value) = 0;

  //! Set all values in the multivector with the given value.
  virtual void putScalar(const Scalar& value) = 0;

  //@}

  //! @name Data Copy and View get methods
  //@{

  //! Return a Vector which is a const view of column j.
  virtual Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> getVector(size_t j) const = 0;

  //! Return a Vector which is a nonconst view of column j.
  virtual Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> getVectorNonConst(size_t j) = 0;

  //! Const view of the local values in a particular vector of this multivector.
  virtual Teuchos::ArrayRCP<const Scalar> getData(size_t j) const = 0;

  //! View of the local values in a particular vector of this multivector.
  virtual Teuchos::ArrayRCP<Scalar> getDataNonConst(size_t j) = 0;

  //@}

  //! @name Mathematical methods
  //@{

  //! Compute dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i]).
  virtual void dot(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const Teuchos::ArrayView<Scalar>& dots) const = 0;

  //! Put element-wise absolute values of input Multi-vector in target: A = abs(this).
  virtual void abs(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A) = 0;

  //! Put element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
  virtual void reciprocal(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A) = 0;

  //! Scale the current values of a multi-vector, this = alpha*this.
  virtual void scale(const Scalar& alpha) = 0;

  //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
  virtual void scale(Teuchos::ArrayView<const Scalar> alpha) = 0;

  //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
  virtual void update(const Scalar& alpha, const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const Scalar& beta) = 0;

  //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
  virtual void update(const Scalar& alpha, const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                      const Scalar& beta, const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                      const Scalar& gamma) = 0;

  //! Compute 1-norm of each vector in multi-vector.
  virtual void norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const = 0;

  //!
  virtual void norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const = 0;

  //! Compute Inf-norm of each vector in multi-vector.
  virtual void normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const = 0;

  //! Compute mean (average) value of each vector in multi-vector. The outcome of this routine is undefined for non-floating point scalar types
  //! (e.g., int).
  virtual void meanValue(const Teuchos::ArrayView<Scalar>& means) const = 0;

  //! Matrix-matrix multiplication: this = beta*this + alpha*op(A)*op(B).
  virtual void multiply(Teuchos::ETransp transA,
                        Teuchos::ETransp transB,
                        const Scalar& alpha,
                        const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                        const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                        const Scalar& beta) = 0;

  /*!
  \brief Multiply a Vector A elementwise by a MultiVector B.

  Compute <tt>this = scalarThis * this + scalarAB * B @ A</tt>
  where <tt>@</tt> denotes element-wise multiplication.  In
  pseudocode, if C denotes <tt>*this</tt> MultiVector:
  \code
  C(i,j) = scalarThis * C(i,j) + scalarAB * B(i,j) * A(i,1);
  \endcode
  for all rows i and columns j of C.

  B must have the same dimensions as <tt>*this</tt>, while A
  must have the same number of rows but a single column.

  We do not require that A, B, and <tt>*this</tt> have
  compatible Maps, as long as the number of rows in A, B, and
  <tt>*this</tt> on each process is the same.  For example, one
  or more of these vectors might have a locally replicated Map,
  or a Map with a local communicator (<tt>MPI_COMM_SELF</tt>).
  This case may occur in block relaxation algorithms when
  applying a diagonal scaling.

  \param[in] scalarAB Scaling factor applied to result of elementwise multiplication of A and B
  \param[in] A Vector A to be element-wise multiplied with B
  \param[in] B MultiVector B to be element-wise multiplied with A
  \param[in] scalarThis Scaling factor for existing values in <tt>*this<\tt> MultiVector C
  */
  virtual void elementWiseMultiply(Scalar scalarAB,
                                   const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                                   const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                                   Scalar scalarThis) = 0;

  //@}

  //! @name Attribute access functions
  //@{

  //! Number of columns in the multivector.
  virtual size_t getNumVectors() const = 0;

  //! Local number of rows on the calling process.
  virtual size_t getLocalLength() const = 0;

  //! Global number of rows in the multivector.
  virtual global_size_t getGlobalLength() const = 0;

  // \brief Checks to see if the local length, number of vectors and size of Scalar type match
  virtual bool isSameSize(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec) const = 0;

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! A simple one-line description of this object.
  virtual std::string description() const = 0;

  //! Print the object with the given verbosity level to a FancyOStream.
  virtual void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const = 0;

  virtual void replaceMap(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map) = 0;

  //@}

  //! @name Xpetra specific
  //@{

  //! Set seed for Random function.
  virtual void setSeed(unsigned int seed) = 0;

  //! Set multi-vector values to random numbers.
  virtual void randomize(bool bUseXpetraImplementation = false) = 0;

  //! Set multi-vector values to random numbers.
  virtual void randomize(const Scalar& minVal, const Scalar& maxVal, bool bUseXpetraImplementation = false) = 0;

  //! Set multi-vector values to random numbers. XPetra implementation
  virtual void Xpetra_randomize();

  //! Set multi-vector values to random numbers. XPetra implementation
  virtual void Xpetra_randomize(const Scalar& minVal, const Scalar& maxVal);

  using impl_scalar_type     = typename Kokkos::ArithTraits<Scalar>::val_type;
  using dual_view_type       = Kokkos::DualView<impl_scalar_type**, Kokkos::LayoutStride, typename node_type::device_type, Kokkos::MemoryUnmanaged>;
  using dual_view_type_const = Kokkos::DualView<const impl_scalar_type**, Kokkos::LayoutStride, typename node_type::device_type, Kokkos::MemoryUnmanaged>;
  using host_execution_space = typename dual_view_type::host_mirror_space;
  using dev_execution_space  = typename dual_view_type::t_dev::execution_space;

  virtual typename dual_view_type::t_host_const_um getHostLocalView(Access::ReadOnlyStruct) const {
    throw std::runtime_error("Dummy function getHostLocalView(Access::ReadOnlyStruct), should be overwritten at" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
#ifndef __NVCC__
    typename dual_view_type::t_host_um test;
#endif
    TEUCHOS_UNREACHABLE_RETURN(test);
  }

  virtual typename dual_view_type::t_dev_const_um getDeviceLocalView(Access::ReadOnlyStruct) const {
    throw std::runtime_error("Dummy function getDeviceLocalView(Access::ReadOnlyStruct), should be overwritten at" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
#ifndef __NVCC__
    typename dual_view_type::t_dev_um test;
#endif
    TEUCHOS_UNREACHABLE_RETURN(test);
  }

  virtual typename dual_view_type::t_host_um getHostLocalView(Access::OverwriteAllStruct) const {
    throw std::runtime_error("Dummy function getHostLocalView(Access::OverwriteAllStruct), should be overwritten at" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
#ifndef __NVCC__
    typename dual_view_type::t_host_um test;
#endif
    TEUCHOS_UNREACHABLE_RETURN(test);
  }

  virtual typename dual_view_type::t_dev_um getDeviceLocalView(Access::OverwriteAllStruct) const {
    throw std::runtime_error("Dummy function getDeviceLocalView(Access::OverwriteAllStruct), should be overwritten at" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
#ifndef __NVCC__
    typename dual_view_type::t_dev_um test;
#endif
    TEUCHOS_UNREACHABLE_RETURN(test);
  }

  virtual typename dual_view_type::t_host_um getHostLocalView(Access::ReadWriteStruct) const {
    throw std::runtime_error("Dummy function getHostLocalView(Access::ReadWriteStruct), should be overwritten at" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
#ifndef __NVCC__
    typename dual_view_type::t_host_um test;
#endif
    TEUCHOS_UNREACHABLE_RETURN(test);
  }

  virtual typename dual_view_type::t_dev_um getDeviceLocalView(Access::ReadWriteStruct) const {
    throw std::runtime_error("Dummy function getDeviceLocalView(Access::ReadWriteStruct), should be overwritten at" + std::string(__FILE__) + ":" + std::to_string(__LINE__));
#ifndef __NVCC__
    typename dual_view_type::t_dev_um test;
#endif
    TEUCHOS_UNREACHABLE_RETURN(test);
  }

  //@}

 protected:
  /// \brief Implementation of the assignment operator (operator=);
  ///   does a deep copy.
  ///
  /// Each subclass must implement this.  This includes
  /// Xpetra::EpetraMultiVector and Xpetra::TpetraMultiVector.
  virtual void assign(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs) = 0;

};  // MultiVector class

}  // namespace Xpetra

#define XPETRA_MULTIVECTOR_SHORT

#endif  // XPETRA_MULTIVECTOR_DECL_HPP
