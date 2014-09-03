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
#include "Xpetra_EpetraVector.hpp"

//TODO: replace double -> Scalar etc.

namespace Xpetra {

  EpetraVector::EpetraVector(const Teuchos::RCP<const Map<int,int> > &map, bool zeroOut) : EpetraMultiVector(map,1,zeroOut) { }

  void EpetraVector::replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { XPETRA_MONITOR("EpetraVector::replaceGlobalValue"); getEpetra_MultiVector()->ReplaceGlobalValue(globalRow, 0, value); }

  void EpetraVector::sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { XPETRA_MONITOR("EpetraVector::sumIntoGlobalValue");getEpetra_MultiVector()->SumIntoGlobalValue(globalRow, 0, value); }

  void EpetraVector::replaceLocalValue(LocalOrdinal myRow, const Scalar &value) { XPETRA_MONITOR("EpetraVector::replaceLocalValue");getEpetra_MultiVector()->ReplaceMyValue(myRow, 0, value); }

  void EpetraVector::sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value) { XPETRA_MONITOR("EpetraVector::sumIntoLocalValue");getEpetra_MultiVector()->SumIntoMyValue(myRow, 0, value); }

  double EpetraVector::dot(const Vector<double,int,int>& a) const {
    XPETRA_MONITOR("EpetraVector::dot");

    XPETRA_DYNAMIC_CAST(const EpetraVector, a, tA, "This Xpetra::EpetraVector method only accept Xpetra::EpetraVector as input arguments.");
    //      return getEpetra_Vector()->Dot(*tA.getEpetra_Vector());

    // other way: use the MultiVector Dot instead of VectorDot:
    double r;
    getEpetra_MultiVector()->Epetra_MultiVector::Dot(*tA.getEpetra_MultiVector(), &r);
    return r;
  }

  Teuchos::ScalarTraits<double>::magnitudeType EpetraVector::norm1() const { XPETRA_MONITOR("EpetraVector::norm1"); double r; getEpetra_MultiVector()->Norm1(&r); return r; }

  Teuchos::ScalarTraits<double>::magnitudeType EpetraVector::norm2() const { XPETRA_MONITOR("EpetraVector::norm2"); double r; getEpetra_MultiVector()->Norm2(&r); return r; }

  Teuchos::ScalarTraits<double>::magnitudeType EpetraVector::normInf() const { XPETRA_MONITOR("EpetraVector::normInf"); double r; getEpetra_MultiVector()->NormInf(&r); return r; }

  Teuchos::ScalarTraits<double>::magnitudeType
  EpetraVector::normWeighted (const Vector<double,int,int>& weights) const
  {
    XPETRA_MONITOR("EpetraVector::normWeighted");
    XPETRA_DYNAMIC_CAST(const EpetraVector, weights, tWeights, "This Xpetra::EpetraVector method only accept Xpetra::EpetraVector as input arguments.");
    double r;
    getEpetra_MultiVector()->NormWeighted(*tWeights.getEpetra_MultiVector(), &r); return r;
  }

  double EpetraVector::meanValue() const {
    XPETRA_MONITOR("EpetraVector::meanValue");
    double r;
    getEpetra_MultiVector()->MeanValue(&r);
    return r;
  }

  std::string EpetraVector::description() const {
    XPETRA_MONITOR("EpetraVector::description");
    // This implementation come from Epetra_Vector_def.hpp (without modification)
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{length="<<this->getGlobalLength()
        << "}";
    return oss.str();
  }

    void EpetraVector::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
      XPETRA_MONITOR("EpetraVector::describe");

      if (verbLevel > Teuchos::VERB_NONE) {
        getEpetra_Vector()->Print (out);
      }
    }

  EpetraVector::EpetraVector(const RCP<Epetra_MultiVector> &mv, size_t j)
    : EpetraMultiVector(rcp((*mv)(j), false)), // view of the vector number j. false == I do not own the data.
      internalRefToBaseMV_(mv)                 // keep an internal reference to the initial MultiVector to avoid desallocation of the view.
  {
    // The view of the internal data of 'mv' is only valid until the destruction of 'mv'.
    // The new vector hold an internal reference to 'mv' in order to keep the view valid after disappearance of 'mv' references in user code.
    // This implements the logic of subArray rcp (as required by the Tpetra interface).
  }

  // TODO: move that elsewhere
  Epetra_Vector & toEpetra(Vector<double, int, int> &x) {
    XPETRA_DYNAMIC_CAST(      EpetraVector, x, tX, "toEpetra");
    return *tX.getEpetra_Vector();
  }

  const Epetra_Vector & toEpetra(const Vector<double, int, int> &x) {
    XPETRA_DYNAMIC_CAST(const EpetraVector, x, tX, "toEpetra");
    return *tX.getEpetra_Vector();
  }
  //

}
