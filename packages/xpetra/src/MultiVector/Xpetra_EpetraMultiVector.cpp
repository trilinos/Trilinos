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
#include "Xpetra_EpetraMultiVector.hpp"

#include "Xpetra_EpetraImport.hpp"
#include "Xpetra_EpetraExport.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  EpetraMultiVector::EpetraMultiVector(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map, const Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > &ArrayOfPtrs, size_t NumVectors) {
    //TODO: input argument 'NumVectors' is not necessary in both Xpetra and Tpetra interface. Should it be removed?

    const std::string tfecfFuncName("MultiVector(ArrayOfPtrs)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(NumVectors < 1 || NumVectors != Teuchos::as<size_t>(ArrayOfPtrs.size()), std::runtime_error,
                                          ": ArrayOfPtrs.size() must be strictly positive and as large as ArrayOfPtrs.");

#ifdef HAVE_XPETRA_DEBUG
    // This cannot be tested by Epetra itself
    {
      size_t localLength = map->getNodeNumElements();
      for(int j=0; j<ArrayOfPtrs.size(); j++) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Teuchos::as<size_t>(ArrayOfPtrs[j].size()) != localLength, std::runtime_error,
                                              ": ArrayOfPtrs[" << j << "].size() (== " << ArrayOfPtrs[j].size() <<
                                              ") is not equal to getLocalLength() (== " << localLength);

      }
    }
#endif

    // Convert Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > to double**
    Array<const double*> arrayOfRawPtrs(ArrayOfPtrs.size());
    for(int i=0; i<ArrayOfPtrs.size(); i++) {
      arrayOfRawPtrs[i] = ArrayOfPtrs[i].getRawPtr();
    }
    double** rawArrayOfRawPtrs = const_cast<double**>(arrayOfRawPtrs.getRawPtr()); // This const_cast should be fine, because Epetra_DataAccess=Copy.

    vec_ = Teuchos::rcp(new Epetra_MultiVector(Copy, toEpetra(map), rawArrayOfRawPtrs, NumVectors));
  }

  Teuchos::ArrayRCP<const double> EpetraMultiVector::getData(size_t j) const {
    XPETRA_MONITOR("EpetraMultiVector::getData");

    double ** arrayOfPointers;

    vec_->ExtractView(&arrayOfPointers);

    double * data = arrayOfPointers[j];
    int localLength = vec_->MyLength();

    return ArrayRCP<double>(data, 0, localLength, false); // not ownership
  }

  Teuchos::ArrayRCP<double> EpetraMultiVector::getDataNonConst(size_t j) {
    XPETRA_MONITOR("EpetraMultiVector::getDataNonConst");

    double ** arrayOfPointers;

    vec_->ExtractView(&arrayOfPointers);

    double * data = arrayOfPointers[j];
    int localLength = vec_->MyLength();

    return ArrayRCP<double>(data, 0, localLength, false); // not ownership
  }

  void EpetraMultiVector::dot(const MultiVector<double,int,int> &A, const Teuchos::ArrayView<double> &dots) const {
    XPETRA_MONITOR("EpetraMultiVector::dot");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    vec_->Dot(*eA.getEpetra_MultiVector(), dots.getRawPtr());
  }

  void EpetraMultiVector::norm1(const Teuchos::ArrayView< Teuchos::ScalarTraits< Scalar >::magnitudeType > &norms) const { XPETRA_MONITOR("EpetraMultiVector::norm1"); vec_->Norm1(norms.getRawPtr()); }

  void EpetraMultiVector::norm2(const Teuchos::ArrayView< Teuchos::ScalarTraits< Scalar >::magnitudeType > &norms) const { XPETRA_MONITOR("EpetraMultiVector::norm2"); vec_->Norm2(norms.getRawPtr()); }

  void EpetraMultiVector::normInf(const Teuchos::ArrayView< Teuchos::ScalarTraits< Scalar >::magnitudeType > &norms) const { XPETRA_MONITOR("EpetraMultiVector::normInf"); vec_->NormInf(norms.getRawPtr()); }

  void EpetraMultiVector::normWeighted(const MultiVector<double,int,int> &weights, const Teuchos::ArrayView<Teuchos::ScalarTraits<double>::magnitudeType> &norms) const {
    XPETRA_MONITOR("EpetraMultiVector::normWeighted");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVector, weights, eWeights, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    vec_->NormWeighted(*eWeights.getEpetra_MultiVector(), norms.getRawPtr());
  }

  void EpetraMultiVector::meanValue(const Teuchos::ArrayView<double> &means) const { XPETRA_MONITOR("EpetraMultiVector::meanValue"); vec_->MeanValue(means.getRawPtr()); } //TODO: modify ArrayView size ??

  std::string EpetraMultiVector::description() const {
    XPETRA_MONITOR("EpetraMultiVector::description");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
    return "TODO";
  }

  void EpetraMultiVector::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    XPETRA_MONITOR("EpetraMultiVector::describe");
    vec_->Print(out);
  }

  void EpetraMultiVector::doImport(const DistObject<double, int, int> &source, const Import<int, int> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraMultiVector::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVector, source, tSource, "Xpetra::EpetraMultiVector::doImport only accept Xpetra::EpetraMultiVector as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Xpetra::EpetraMultiVector::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<Epetra_MultiVector> v = tSource.getEpetra_MultiVector();
    int err = this->getEpetra_MultiVector()->Import(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraMultiVector::doExport(const DistObject<double, int, int> &dest, const Import<int, int>& importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraMultiVector::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVector, dest, tDest, "Xpetra::EpetraMultiVector::doImport only accept Xpetra::EpetraMultiVector as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Xpetra::EpetraMultiVector::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<Epetra_MultiVector> v = tDest.getEpetra_MultiVector();
    int err = this->getEpetra_MultiVector()->Export(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraMultiVector::doImport(const DistObject<double,int,int> &source, const Export<int, int>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraMultiVector::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVector, source, tSource, "Xpetra::EpetraMultiVector::doImport only accept Xpetra::EpetraMultiVector as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Xpetra::EpetraMultiVector::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<Epetra_MultiVector> v = tSource.getEpetra_MultiVector();
    int err = this->getEpetra_MultiVector()->Import(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraMultiVector::doExport(const DistObject<double, int, int> &dest, const Export<int, int>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraMultiVector::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVector, dest, tDest, "Xpetra::EpetraMultiVector::doImport only accept Xpetra::EpetraMultiVector as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Xpetra::EpetraMultiVector::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<Epetra_MultiVector> v = tDest.getEpetra_MultiVector();
    int err = this->getEpetra_MultiVector()->Export(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  // TODO: move that elsewhere
  const Epetra_MultiVector & toEpetra(const MultiVector<double, int, int> & x) {
    XPETRA_DYNAMIC_CAST(const EpetraMultiVector, x, tX, "toEpetra");
    return *tX.getEpetra_MultiVector();
  }

  Epetra_MultiVector & toEpetra(MultiVector<double, int, int> & x) {
    XPETRA_DYNAMIC_CAST(      EpetraMultiVector, x, tX, "toEpetra");
    return *tX.getEpetra_MultiVector();
  }
  //


} // namespace Xpetra
