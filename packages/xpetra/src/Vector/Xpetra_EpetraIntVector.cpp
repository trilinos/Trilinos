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
#include "Xpetra_EpetraIntVector.hpp"
#include "Xpetra_EpetraImport.hpp"
#include "Xpetra_EpetraExport.hpp"

namespace Xpetra {

  void EpetraIntVector::replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { XPETRA_MONITOR("EpetraIntVector::replaceGlobalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { XPETRA_MONITOR("EpetraIntVector::sumIntoGlobalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::replaceLocalValue(LocalOrdinal myRow, const Scalar &value) { XPETRA_MONITOR("EpetraIntVector::replaceLocalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value) { XPETRA_MONITOR("EpetraIntVector::sumIntoLocalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  int EpetraIntVector::dot(const Vector<int,int,int> &a) const { XPETRA_MONITOR("EpetraIntVector::dot"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  Teuchos::ScalarTraits<int>::magnitudeType EpetraIntVector::norm1() const { XPETRA_MONITOR("EpetraIntVector::norm1"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  Teuchos::ScalarTraits<int>::magnitudeType EpetraIntVector::norm2() const { XPETRA_MONITOR("EpetraIntVector::norm2"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  Teuchos::ScalarTraits<int>::magnitudeType EpetraIntVector::normInf() const { XPETRA_MONITOR("EpetraIntVector::normInf"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  Teuchos::ScalarTraits<int>::magnitudeType EpetraIntVector::normWeighted(const Vector<int,int,int> &weights) const { XPETRA_MONITOR("EpetraIntVector::normWeighted"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  int EpetraIntVector::meanValue() const { XPETRA_MONITOR("EpetraIntVector::meanValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  int EpetraIntVector::maxValue() const { XPETRA_MONITOR("EpetraIntVector::maxValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  void EpetraIntVector::randomize(bool bUseXpetraImplementation) { XPETRA_MONITOR("EpetraIntVector::randomize"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVector::randomize(): Functionnality not available in Epetra"); }

  void EpetraIntVector::setSeed(unsigned int seed) { XPETRA_MONITOR("EpetraIntVector::setSeed"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVector::setSeed(): Functionnality not available in Epetra"); }

  Teuchos::ArrayRCP<const int> EpetraIntVector::getData(size_t j) const {
    XPETRA_MONITOR("EpetraIntVector::getData");

    int * data = vec_->Values();
    int localLength = vec_->MyLength();

    return ArrayRCP<int>(data, 0, localLength, false); // not ownership
  }

  Teuchos::ArrayRCP<int> EpetraIntVector::getDataNonConst(size_t j) {
    XPETRA_MONITOR("EpetraIntVector::getDataNonConst");

    int * data = vec_->Values();
    int localLength = vec_->MyLength();

    return ArrayRCP<int>(data, 0, localLength, false); // not ownership
  }

  void EpetraIntVector::dot(const MultiVector<int,int,int> &A, const Teuchos::ArrayView<int> &dots) const {
    XPETRA_MONITOR("EpetraIntVector::dot");

    //XPETRA_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  void EpetraIntVector::abs(const MultiVector<int,int,int> &A) {
    XPETRA_MONITOR("EpetraIntVector::abs");

    //XPETRA_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  void EpetraIntVector::reciprocal(const MultiVector<int,int,int> &A) {
    XPETRA_MONITOR("EpetraIntVector::reciprocal");

    //XPETRA_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  void EpetraIntVector::scale(const int &alpha) { XPETRA_MONITOR("EpetraIntVector::scale"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::update(const int &alpha, const MultiVector<int,int,int> &A, const int &beta) {
    XPETRA_MONITOR("EpetraIntVector::update");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  void EpetraIntVector::update(const int &alpha, const MultiVector<int,int,int> &A, const int &beta, const MultiVector<int,int,int> &B, const int &gamma) {
    XPETRA_MONITOR("EpetraIntVector::update");

    //XPETRA_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    //XPETRA_DYNAMIC_CAST(const EpetraMultiVector, B, eB, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  void EpetraIntVector::norm1(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const { XPETRA_MONITOR("EpetraIntVector::norm1"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::norm2(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const { XPETRA_MONITOR("EpetraIntVector::norm2"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::normInf(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const { XPETRA_MONITOR("EpetraIntVector::normInf"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::normWeighted(const MultiVector<int,int,int> &weights, const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const { XPETRA_MONITOR("EpetraIntVector::normWeighted"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::meanValue(const Teuchos::ArrayView<int> &means) const { XPETRA_MONITOR("EpetraIntVector::meanValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::maxValue(const Teuchos::ArrayView<int> &maxs) const { XPETRA_MONITOR("EpetraIntVector::maxValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const int &alpha, const MultiVector<int,int,int> &A, const MultiVector<int,int,int> &B, const int &beta) { XPETRA_MONITOR("EpetraIntVector::multiply"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Not available in Epetra"); }

  void EpetraIntVector::elementWiseMultiply(int scalarAB, const Vector<int,int,int> &A, const MultiVector<int,int,int> &B, int scalarThis) {
    XPETRA_MONITOR("EpetraIntVector::elementWiseMultiply");

    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra_EpetraIntVector: elementWiseMultiply not implemented because Epetra_IntVector does not support this operation");
  }

  void EpetraIntVector::replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) { XPETRA_MONITOR("EpetraIntVector::replaceGlobalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) { XPETRA_MONITOR("EpetraIntVector::sumIntoGlobalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) { XPETRA_MONITOR("EpetraIntVector::replaceLocalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) { XPETRA_MONITOR("EpetraIntVector::sumIntoLocalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  size_t EpetraIntVector::getNumVectors() const { XPETRA_MONITOR("EpetraIntVector::getNumVectors"); return 1; }

  std::string EpetraIntVector::description() const {
    XPETRA_MONITOR("EpetraIntVector::description");

    // This implementation come from Epetra_Vector_def.hpp (without modification)
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{length="<<this->getGlobalLength()
        << "}";
    return oss.str();
  }

  void EpetraIntVector::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    XPETRA_MONITOR("EpetraIntVector::describe");

    //   typedef Kokkos::MultiVector<double> KMV;
    //   typedef Kokkos::DefaultArithmetic<KMV> MVT;

    // This implementation come from Tpetra_Vector_def.hpp (without modification) // JG: true?
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;

    if (verbLevel > Teuchos::VERB_NONE)
      vec_->Print(out);
  }

  void EpetraIntVector::doImport(const DistObject<int, int, int> &source,
                                 const Import<int, int> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVector::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVector, source, tSource, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraIntVector as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraImport as input arguments.");

    const Epetra_IntVector & v = *tSource.getEpetra_IntVector();
    int err = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraIntVector::doExport(const DistObject<int, int, int> &dest,
                                 const Import<int, int>& importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVector::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVector, dest, tDest, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraIntVector as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraImport as input arguments.");

    const Epetra_IntVector & v = *tDest.getEpetra_IntVector();
    int err = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraIntVector::doImport(const DistObject<int, int, int> &source,
                                 const Export<int, int>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVector::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVector, source, tSource, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraIntVector as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraImport as input arguments.");

    const Epetra_IntVector & v = *tSource.getEpetra_IntVector();
    int err = vec_->Import(v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraIntVector::doExport(const DistObject<int, int, int> &dest,
                                 const Export<int, int>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVector::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVector, dest, tDest, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraIntVector as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraImport as input arguments.");

    const Epetra_IntVector & v = *tDest.getEpetra_IntVector();
    int err = vec_->Export(v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

} // namespace Xpetra
