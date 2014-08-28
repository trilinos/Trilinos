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
#include "Xpetra_EpetraMultiVector.hpp"

#include "Xpetra_EpetraImport.hpp"
#include "Xpetra_EpetraExport.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_EpetraVector.hpp"

#include "Epetra_SerialComm.h"

namespace Xpetra {

  template<class EpetraGlobalOrdinal>
  EpetraMultiVectorT<EpetraGlobalOrdinal>::EpetraMultiVectorT(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map, const Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > &ArrayOfPtrs, size_t NumVectors) {
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


  template<class EpetraGlobalOrdinal>
  Teuchos::RCP< const Vector< double, int, EpetraGlobalOrdinal > > EpetraMultiVectorT<EpetraGlobalOrdinal>::getVector(size_t j) const {
    XPETRA_MONITOR("EpetraMultiVectorT::getVector");
    return rcp(new EpetraVectorT<GlobalOrdinal>(vec_, j)); // See constructor EpetraVectorT(const RCP<EpetraMultiVectorT> &mv, size_t j) for more info
  }

  template<class EpetraGlobalOrdinal>
  Teuchos::RCP< Vector< double, int, EpetraGlobalOrdinal > > EpetraMultiVectorT<EpetraGlobalOrdinal>::getVectorNonConst(size_t j) {
    XPETRA_MONITOR("EpetraMultiVectorT::getVector");
    return rcp(new EpetraVectorT<GlobalOrdinal>(vec_, j)); // See constructor EpetraVectorT(const RCP<EpetraMultiVectorT> &mv, size_t j) for more info
  }

  template<class EpetraGlobalOrdinal>
  Teuchos::ArrayRCP<const double> EpetraMultiVectorT<EpetraGlobalOrdinal>::getData(size_t j) const {
    XPETRA_MONITOR("EpetraMultiVectorT::getData");

    double ** arrayOfPointers;

    vec_->ExtractView(&arrayOfPointers);

    double * data = arrayOfPointers[j];
    int localLength = vec_->MyLength();

    return ArrayRCP<double>(data, 0, localLength, false); // no ownership
  }

  template<class EpetraGlobalOrdinal>
  Teuchos::ArrayRCP<double> EpetraMultiVectorT<EpetraGlobalOrdinal>::getDataNonConst(size_t j) {
    XPETRA_MONITOR("EpetraMultiVectorT::getDataNonConst");

    double ** arrayOfPointers;

    vec_->ExtractView(&arrayOfPointers);

    double * data = arrayOfPointers[j];
    int localLength = vec_->MyLength();

    return ArrayRCP<double>(data, 0, localLength, false); // no ownership
  }

  template<class EpetraGlobalOrdinal>
  void EpetraMultiVectorT<EpetraGlobalOrdinal>::dot(const MultiVector<double,int,GlobalOrdinal> &A, const Teuchos::ArrayView<double> &dots) const {
    XPETRA_MONITOR("EpetraMultiVectorT::dot");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    vec_->Dot(*eA.getEpetra_MultiVector(), dots.getRawPtr());
  }

  template<class EpetraGlobalOrdinal>
  void EpetraMultiVectorT<EpetraGlobalOrdinal>::norm1(const Teuchos::ArrayView< Teuchos::ScalarTraits< Scalar >::magnitudeType > &norms) const { XPETRA_MONITOR("EpetraMultiVectorT::norm1"); vec_->Norm1(norms.getRawPtr()); }

  template<class EpetraGlobalOrdinal>
  void EpetraMultiVectorT<EpetraGlobalOrdinal>::norm2(const Teuchos::ArrayView< Teuchos::ScalarTraits< Scalar >::magnitudeType > &norms) const { XPETRA_MONITOR("EpetraMultiVectorT::norm2"); vec_->Norm2(norms.getRawPtr()); }

  template<class EpetraGlobalOrdinal>
  void EpetraMultiVectorT<EpetraGlobalOrdinal>::normInf(const Teuchos::ArrayView< Teuchos::ScalarTraits< Scalar >::magnitudeType > &norms) const { XPETRA_MONITOR("EpetraMultiVectorT::normInf"); vec_->NormInf(norms.getRawPtr()); }

  template<class EpetraGlobalOrdinal>
  void EpetraMultiVectorT<EpetraGlobalOrdinal>::normWeighted(const MultiVector<double,int,GlobalOrdinal> &weights, const Teuchos::ArrayView<Teuchos::ScalarTraits<double>::magnitudeType> &norms) const {
    XPETRA_MONITOR("EpetraMultiVectorT::normWeighted");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, weights, eWeights, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    vec_->NormWeighted(*eWeights.getEpetra_MultiVector(), norms.getRawPtr());
  }

  template<class EpetraGlobalOrdinal>
  void EpetraMultiVectorT<EpetraGlobalOrdinal>::meanValue(const Teuchos::ArrayView<double> &means) const { XPETRA_MONITOR("EpetraMultiVectorT::meanValue"); vec_->MeanValue(means.getRawPtr()); } //TODO: modify ArrayView size ??

  template<class EpetraGlobalOrdinal>
  std::string EpetraMultiVectorT<EpetraGlobalOrdinal>::description() const {
    XPETRA_MONITOR("EpetraMultiVectorT::description");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
    return "TODO";
  }

  template<class EpetraGlobalOrdinal>
  void EpetraMultiVectorT<EpetraGlobalOrdinal>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    XPETRA_MONITOR("EpetraMultiVectorT::describe");
    vec_->Print(out);
  }

  template<class EpetraGlobalOrdinal>
  void EpetraMultiVectorT<EpetraGlobalOrdinal>::doImport(const DistObject<double, int, GlobalOrdinal> &source, const Import<int, GlobalOrdinal> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraMultiVectorT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT<GlobalOrdinal>, source, tSource, "Xpetra::EpetraMultiVectorT::doImport only accept Xpetra::EpetraMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal>, importer, tImporter, "Xpetra::EpetraMultiVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<Epetra_MultiVector> v = tSource.getEpetra_MultiVector();
    int err = this->getEpetra_MultiVector()->Import(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraMultiVectorT<EpetraGlobalOrdinal>::doExport(const DistObject<double, int, GlobalOrdinal> &dest, const Import<int, GlobalOrdinal>& importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraMultiVectorT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT<GlobalOrdinal>, dest, tDest, "Xpetra::EpetraMultiVectorT::doImport only accept Xpetra::EpetraMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal>, importer, tImporter, "Xpetra::EpetraMultiVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<Epetra_MultiVector> v = tDest.getEpetra_MultiVector();
    int err = this->getEpetra_MultiVector()->Export(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraMultiVectorT<EpetraGlobalOrdinal>::doImport(const DistObject<double,int,GlobalOrdinal> &source, const Export<int, GlobalOrdinal>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraMultiVectorT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT<GlobalOrdinal>, source, tSource, "Xpetra::EpetraMultiVectorT::doImport only accept Xpetra::EpetraMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal>, exporter, tExporter, "Xpetra::EpetraMultiVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<Epetra_MultiVector> v = tSource.getEpetra_MultiVector();
    int err = this->getEpetra_MultiVector()->Import(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraMultiVectorT<EpetraGlobalOrdinal>::doExport(const DistObject<double, int, GlobalOrdinal> &dest, const Export<int, GlobalOrdinal>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraMultiVectorT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT<GlobalOrdinal>, dest, tDest, "Xpetra::EpetraMultiVectorT::doImport only accept Xpetra::EpetraMultiVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal>, exporter, tExporter, "Xpetra::EpetraMultiVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    RCP<Epetra_MultiVector> v = tDest.getEpetra_MultiVector();
    int err = this->getEpetra_MultiVector()->Export(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraMultiVectorT<EpetraGlobalOrdinal>::replaceMap(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map) {
    int err = 0;
    if (!map.is_null()) {
      err = this->getEpetra_MultiVector()->ReplaceMap(toEpetra(map));

    } else {
      // Replace map with a dummy map to avoid potential hangs later
      Epetra_SerialComm SComm;
      Epetra_Map NewMap((EpetraGlobalOrdinal) vec_->MyLength(), (EpetraGlobalOrdinal) vec_->Map().IndexBase64(), SComm);
      err = this->getEpetra_MultiVector()->ReplaceMap(NewMap);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraMultiVectorT<EpetraGlobalOrdinal>::
  assign (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs)
  {
    typedef EpetraMultiVectorT this_type;
    const this_type* rhsPtr = dynamic_cast<const this_type*> (&rhs);
    TEUCHOS_TEST_FOR_EXCEPTION(
      rhsPtr == NULL, std::invalid_argument, "Xpetra::MultiVector::operator=: "
      "The left-hand side (LHS) of the assignment has a different type than "
      "the right-hand side (RHS).  The LHS has type Xpetra::EpetraMultiVectorT "
      "(which means it wraps an Epetra_MultiVector), but the RHS has some "
      "other type.  This probably means that the RHS wraps a Tpetra::Multi"
      "Vector.  Xpetra::MultiVector does not currently implement assignment "
      "from a Tpetra object to an Epetra object, though this could be added "
      "with sufficient interest.");

    RCP<const Epetra_MultiVector> rhsImpl = rhsPtr->getEpetra_MultiVector ();
    RCP<Epetra_MultiVector> lhsImpl = this->getEpetra_MultiVector ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      rhsImpl.is_null (), std::logic_error, "Xpetra::MultiVector::operator= "
      "(in Xpetra::EpetraMultiVectorT::assign): *this (the right-hand side of "
      "the assignment) has a null RCP<Epetra_MultiVector> inside.  Please "
      "report this bug to the Xpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      lhsImpl.is_null (), std::logic_error, "Xpetra::MultiVector::operator= "
      "(in Xpetra::EpetraMultiVectorT::assign): The left-hand side of the "
      "assignment has a null RCP<Epetra_MultiVector> inside.  Please report "
      "this bug to the Xpetra developers.");

    // Epetra_MultiVector's assignment operator does a deep copy.
    *lhsImpl = *rhsImpl;
  }

  // TODO: move that elsewhere
  template<class GlobalOrdinal>
  const Epetra_MultiVector & toEpetra(const MultiVector<double, int, GlobalOrdinal> & x) {
    XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT<GlobalOrdinal>, x, tX, "toEpetra");
    return *tX.getEpetra_MultiVector();
  }

  template<class GlobalOrdinal>
  Epetra_MultiVector & toEpetra(MultiVector<double, int, GlobalOrdinal> & x) {
    XPETRA_DYNAMIC_CAST(      EpetraMultiVectorT<GlobalOrdinal>, x, tX, "toEpetra");
    return *tX.getEpetra_MultiVector();
  }
  //

  template<class GlobalOrdinal>
  RCP<MultiVector<double, int, GlobalOrdinal> > toXpetra(RCP<Epetra_MultiVector> vec) {
    if (!vec.is_null())
      return rcp(new EpetraMultiVectorT<GlobalOrdinal>(vec));

    return Teuchos::null;
  }

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
template class EpetraMultiVectorT<int>;
template RCP<MultiVector<double, int, int> > toXpetra<int>(RCP<Epetra_MultiVector>);
template const Epetra_MultiVector & toEpetra(const MultiVector<double, int, int> &);
template Epetra_MultiVector & toEpetra(MultiVector<double, int, int> &);
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
template class EpetraMultiVectorT<long long>;
template RCP<MultiVector<double, int, long long> > toXpetra<long long>(RCP<Epetra_MultiVector>);
template const Epetra_MultiVector & toEpetra(const MultiVector<double, int, long long> &);
template Epetra_MultiVector & toEpetra(MultiVector<double, int, long long> &);
#endif

} // namespace Xpetra
