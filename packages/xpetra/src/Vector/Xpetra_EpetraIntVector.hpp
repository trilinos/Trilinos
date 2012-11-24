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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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

namespace Xpetra {

  class EpetraIntVector
    : public Vector<int,int,int>
  {

    typedef int Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;

  public:

    //! @name Constructor/Destructor Methods
    //@{

    //! Sets all vector entries to zero.
    explicit EpetraIntVector(const Teuchos::RCP<const Map<int,int> > &map, bool zeroOut=true)
    {
      XPETRA_RCP_DYNAMIC_CAST(const EpetraMap, map, eMap, "Xpetra::EpetraCrsMatrix constructors only accept Xpetra::EpetraMap as input arguments.");
      vec_ = rcp(new Epetra_IntVector(eMap->getEpetra_BlockMap(), zeroOut));
    }

    //! Destructor.
    ~EpetraIntVector() {  };

    //@}

    //! @name Mathematical methods
    //@{

    //! TODO missing comment
    int dot(const Vector<int,int,int> &a) const;

    //! Return 1-norm of this Vector.
    Teuchos::ScalarTraits<int>::magnitudeType norm1() const;

    //! Compute 2-norm of this Vector.
    Teuchos::ScalarTraits<int>::magnitudeType norm2() const;

    //! Compute Inf-norm of this Vector.
    Teuchos::ScalarTraits<int>::magnitudeType normInf() const;

    //! Compute Weighted 2-norm (RMS Norm) of this Vector.
    Teuchos::ScalarTraits<int>::magnitudeType normWeighted(const Vector<int,int,int> &weights) const;

    //! Compute mean (average) value of this Vector.
    int meanValue() const;

    //! Compute max value of this Vector.
    int maxValue() const;

    //@}

    //! @name Post-construction modification routines
    //@{

    //! Replace current value at the specified location with specified value.
    void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value);

    //! Adds specified value to existing value at the specified location.
    void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value);

    //! Replace current value at the specified location with specified values.
    void replaceLocalValue(LocalOrdinal myRow, const Scalar &value);

    //! Adds specified value to existing value at the specified location.
    void sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value);

    //! Initialize all values in a multi-vector with specified value.
    void putScalar(const int &value) {  vec_->PutValue(value); }

    //! Set multi-vector values to random numbers.
    void randomize(bool bUseXpetraImplementation = true);

    //! Set seed for Random function.
    /** Note: this method does not exist in Tpetra interface. Added for MueLu. */
    void setSeed(unsigned int seed);

    //@}

    //! @name Data Copy and View get methods
    //@{

    //! Return a Vector which is a const view of column j.
    Teuchos::RCP< const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > > getVector(size_t j) const;

    //! Return a Vector which is a nonconst view of column j.
    Teuchos::RCP< Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > > getVectorNonConst(size_t j);

    //! Const Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    Teuchos::ArrayRCP<const int> getData(size_t j) const;

    //! Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    Teuchos::ArrayRCP<int> getDataNonConst(size_t j);

    //@}

    //! @name Mathematical methods
    //@{
    //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    void dot(const MultiVector<int,int,int,Kokkos::DefaultNode::DefaultNodeType> &A, const Teuchos::ArrayView<int> &dots) const;

    //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
    void abs(const MultiVector<int,int,int,Kokkos::DefaultNode::DefaultNodeType> &A);

    //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    void reciprocal(const MultiVector<int,int,int,Kokkos::DefaultNode::DefaultNodeType> &A);

    //! Scale the current values of a multi-vector, this = alpha*this.
    void scale(const int &alpha);

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    void update(const int &alpha, const MultiVector<int,int,int,Kokkos::DefaultNode::DefaultNodeType> &A, const int &beta);

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    void update(const int &alpha, const MultiVector<int,int,int,Kokkos::DefaultNode::DefaultNodeType> &A, const int &beta, const MultiVector<int,int,int,Kokkos::DefaultNode::DefaultNodeType> &B, const int &gamma);

    //! Compute 1-norm of each vector in multi-vector.
    void norm1(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const;

    //! Compute 2-norm of each vector in multi-vector.
    void norm2(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const;

    //! Compute Inf-norm of each vector in multi-vector.
    void normInf(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const;

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    void normWeighted(const MultiVector<int,int,int,Kokkos::DefaultNode::DefaultNodeType> &weights, const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const;

    //! Compute mean (average) value of each vector in multi-vector.
    void meanValue(const Teuchos::ArrayView<int> &means) const;

    //! Compute max value of each vector in multi-vector.
    void maxValue(const Teuchos::ArrayView<int> &maxs) const;

    //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
    void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const int &alpha, const MultiVector<int,int,int,Kokkos::DefaultNode::DefaultNodeType> &A, const MultiVector<int,int,int,Kokkos::DefaultNode::DefaultNodeType> &B, const int &beta);

    //! Element-wise multiply of a Vector A with a EpetraMultiVector B.
    void elementWiseMultiply(int scalarAB, const Vector<int,int,int,Kokkos::DefaultNode::DefaultNodeType> &A, const MultiVector<int,int,int,Kokkos::DefaultNode::DefaultNodeType> &B, int scalarThis);
    //@}

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

    //@}

    //! @name Attribute access functions
    //@{

    //! Returns the number of vectors in the multi-vector.
    size_t getNumVectors() const;

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    size_t getLocalLength() const {  return vec_->MyLength(); }

    //! Returns the global vector length of vectors in the multi-vector.
    global_size_t getGlobalLength() const {  return vec_->GlobalLength(); }

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

    //@}

    RCP< Epetra_IntVector > getEpetra_IntVector() const {  return vec_; }

    const RCP<const Comm<int> > getComm() const {
      TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO getComm Epetra MultiVector not implemented");
    }

    // Implementing DistObject
    const Teuchos::RCP<const Map<int,int> > getMap() const {
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(vec_->Map()));
      return rcp ( new Xpetra::EpetraMap(map) );
    }

    void doImport(const DistObject<int, int, int> &source, const Import<int, int> &importer, CombineMode CM);

    void doExport(const DistObject<int, int, int> &dest, const Import<int, int>& importer, CombineMode CM);

    void doImport(const DistObject<int, int, int> &source, const Export<int, int>& exporter, CombineMode CM);

    void doExport(const DistObject<int, int, int> &dest, const Export<int, int>& exporter, CombineMode CM);

  private:
    RCP< Epetra_IntVector > vec_;

  }; // class EpetraIntVector

} // namespace Xpetra

#endif // XPETRA_EPETRAINTVECTOR_HPP
