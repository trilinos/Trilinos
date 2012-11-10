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
#include "Xpetra_EpetraVector.hpp"

//TODO: replace double -> Scalar etc.

namespace Xpetra {

  EpetraVector::EpetraVector(const Teuchos::RCP<const Map<int,int> > &map, bool zeroOut) : EpetraMultiVector(map,1,zeroOut) { }

  void EpetraVector::replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { XPETRA_MONITOR("EpetraVector::replaceGlobalValue"); getEpetra_MultiVector()->ReplaceGlobalValue(globalRow, 0, value); }

  void EpetraVector::sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { XPETRA_MONITOR("EpetraVector::sumIntoGlobalValue");getEpetra_MultiVector()->SumIntoGlobalValue(globalRow, 0, value); }

  void EpetraVector::replaceLocalValue(LocalOrdinal myRow, const Scalar &value) { XPETRA_MONITOR("EpetraVector::replaceLocalValue");getEpetra_MultiVector()->ReplaceMyValue(myRow, 0, value); }

  void EpetraVector::sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value) { XPETRA_MONITOR("EpetraVector::sumIntoLocalValue");getEpetra_MultiVector()->SumIntoMyValue(myRow, 0, value); }

  double EpetraVector::dot(const Vector<double,int,int,Kokkos::DefaultNode::DefaultNodeType> &a) const {
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

  Teuchos::ScalarTraits<double>::magnitudeType EpetraVector::normWeighted(const Vector<double,int,int,Kokkos::DefaultNode::DefaultNodeType> &weights) const {
      XPETRA_MONITOR("EpetraVector::normWeighted");
      XPETRA_DYNAMIC_CAST(const EpetraVector, weights, tWeights, "This Xpetra::EpetraVector method only accept Xpetra::EpetraVector as input arguments.");
      double r;
      getEpetra_MultiVector()->NormWeighted(*tWeights.getEpetra_MultiVector(), &r); return r;
    }

    double EpetraVector::meanValue() const { XPETRA_MONITOR("EpetraVector::meanValue"); double r; getEpetra_MultiVector()->MeanValue(&r); return r; }

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

      if (verbLevel > Teuchos::VERB_NONE)
        getEpetra_Vector()->Print(out);

//       typedef Kokkos::MultiVector<double> KMV;
//       typedef Kokkos::DefaultArithmetic<KMV>   MVT;

//       // This implementation come from Tpetra_Vector_def.hpp (without modification)
//       using std::endl;
//       using std::setw;
//       using Teuchos::VERB_DEFAULT;
//       using Teuchos::VERB_NONE;
//       using Teuchos::VERB_LOW;
//       using Teuchos::VERB_MEDIUM;
//       using Teuchos::VERB_HIGH;
//       using Teuchos::VERB_EXTREME;
//       Teuchos::EVerbosityLevel vl = verbLevel;

//       TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");

//       if (vl == VERB_DEFAULT) vl = VERB_LOW;
//       Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getMap()->getComm();
//       const int myImageID = comm->getRank(),
//       numImages = comm->getSize();
//       size_t width = 1;
//       for (size_t dec=10; dec<this->getGlobalLength(); dec *= 10) {
//         ++width;
//       }
//       Teuchos::OSTab tab(out);
//       if (vl != VERB_NONE) {
//         // VERB_LOW and higher prints description()
//         if (myImageID == 0) out << this->description() << std::endl;
//         for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
//           if (myImageID == imageCtr) {
//             if (vl != VERB_LOW) {
//               // VERB_MEDIUM and higher prints getLocalLength()
//               out << "node " << setw(width) << myImageID << ": local length=" << this->getLocalLength() << endl;
//               if (vl != VERB_MEDIUM) {
//                 // VERB_HIGH and higher prints isConstantStride() and stride()
//                 if (vl == VERB_EXTREME && this->getLocalLength() > 0) {
//                   Teuchos::RCP<Node> node = this->lclMV_.getNode();
//                   KOKKOS_NODE_TRACE("Vector::describe()")
//                     Teuchos::ArrayRCP<const double> myview = node->template viewBuffer<double>(
//                                                                                                this->getLocalLength(),
//                                                                                                MVT::getValues(this->lclMV_) );
//                   // VERB_EXTREME prints values
//                   for (size_t i=0; i<this->getLocalLength(); ++i) {
//                     out << setw(width) << this->getMap()->getGlobalElement(i)
//                         << ": "
//                         << myview[i] << endl;
//                   }
//                   myview = Teuchos::null;
//                 }
//               }
//               else {
//                 out << endl;
//               }
//             }
//           }
//         }
//       }
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
