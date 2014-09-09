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
#include "Xpetra_EpetraIntVector.hpp"
#include "Xpetra_EpetraImport.hpp"
#include "Xpetra_EpetraExport.hpp"

namespace Xpetra {

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { XPETRA_MONITOR("EpetraIntVectorT::replaceGlobalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { XPETRA_MONITOR("EpetraIntVectorT::sumIntoGlobalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::replaceLocalValue(LocalOrdinal myRow, const Scalar &value) { XPETRA_MONITOR("EpetraIntVectorT::replaceLocalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value) { XPETRA_MONITOR("EpetraIntVectorT::sumIntoLocalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  int EpetraIntVectorT<EpetraGlobalOrdinal>::dot(const Vector<int,int,GlobalOrdinal,Node> &a) const { XPETRA_MONITOR("EpetraIntVectorT::dot"); TEUCHOS_TEST_FOR_EXCEPTION(-1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  template<class EpetraGlobalOrdinal>
  Teuchos::ScalarTraits<int>::magnitudeType EpetraIntVectorT<EpetraGlobalOrdinal>::norm1() const { XPETRA_MONITOR("EpetraIntVectorT::norm1"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  template<class EpetraGlobalOrdinal>
  Teuchos::ScalarTraits<int>::magnitudeType EpetraIntVectorT<EpetraGlobalOrdinal>::norm2() const { XPETRA_MONITOR("EpetraIntVectorT::norm2"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  template<class EpetraGlobalOrdinal>
  Teuchos::ScalarTraits<int>::magnitudeType EpetraIntVectorT<EpetraGlobalOrdinal>::normInf() const { XPETRA_MONITOR("EpetraIntVectorT::normInf"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  template<class EpetraGlobalOrdinal>
  Teuchos::ScalarTraits<int>::magnitudeType EpetraIntVectorT<EpetraGlobalOrdinal>::normWeighted(const Vector<int,int,GlobalOrdinal,Node> &weights) const { XPETRA_MONITOR("EpetraIntVectorT::normWeighted"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  template<class EpetraGlobalOrdinal>
  int EpetraIntVectorT<EpetraGlobalOrdinal>::meanValue() const { XPETRA_MONITOR("EpetraIntVectorT::meanValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  template<class EpetraGlobalOrdinal>
  int EpetraIntVectorT<EpetraGlobalOrdinal>::maxValue() const { XPETRA_MONITOR("EpetraIntVectorT::maxValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::randomize(bool bUseXpetraImplementation) { XPETRA_MONITOR("EpetraIntVectorT::randomize"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVectorT::randomize(): Functionnality not available in Epetra"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::setSeed(unsigned int seed) { XPETRA_MONITOR("EpetraIntVectorT::setSeed"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVectorT::setSeed(): Functionnality not available in Epetra"); }

  template<class EpetraGlobalOrdinal>
  Teuchos::RCP< const Vector< int,int,EpetraGlobalOrdinal > > EpetraIntVectorT<EpetraGlobalOrdinal>::getVector(size_t j) const {
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  template<class EpetraGlobalOrdinal>
  Teuchos::RCP< Vector< int,int,EpetraGlobalOrdinal > > EpetraIntVectorT<EpetraGlobalOrdinal>::getVectorNonConst(size_t j) {
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  template<class EpetraGlobalOrdinal>
  Teuchos::ArrayRCP<const int> EpetraIntVectorT<EpetraGlobalOrdinal>::getData(size_t j) const {
    XPETRA_MONITOR("EpetraIntVectorT::getData");

    int * data = vec_->Values();
    int localLength = vec_->MyLength();

    return ArrayRCP<int>(data, 0, localLength, false); // not ownership
  }

  template<class EpetraGlobalOrdinal>
  Teuchos::ArrayRCP<int> EpetraIntVectorT<EpetraGlobalOrdinal>::getDataNonConst(size_t j) {
    XPETRA_MONITOR("EpetraIntVectorT::getDataNonConst");

    int * data = vec_->Values();
    int localLength = vec_->MyLength();

    return ArrayRCP<int>(data, 0, localLength, false); // not ownership
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::dot(const MultiVector<int,int,GlobalOrdinal,Node> &A, const Teuchos::ArrayView<int> &dots) const {
    XPETRA_MONITOR("EpetraIntVectorT::dot");

    //XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::abs(const MultiVector<int,int,GlobalOrdinal,Node> &A) {
    XPETRA_MONITOR("EpetraIntVectorT::abs");

    //XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::reciprocal(const MultiVector<int,int,GlobalOrdinal,Node> &A) {
    XPETRA_MONITOR("EpetraIntVectorT::reciprocal");

    //XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::scale(const int &alpha) {
    XPETRA_MONITOR("EpetraIntVectorT::scale");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::scale (Teuchos::ArrayView< const int > alpha) {
    XPETRA_MONITOR("EpetraIntVectorT::scale");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::update(const int &alpha, const MultiVector<int,int,GlobalOrdinal,Node> &A, const int &beta) {
    XPETRA_MONITOR("EpetraIntVectorT::update");

    // XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::update(const int &alpha, const MultiVector<int,int,GlobalOrdinal,Node> &A, const int &beta, const MultiVector<int,int,GlobalOrdinal,Node> &B, const int &gamma) {
    XPETRA_MONITOR("EpetraIntVectorT::update");

    //XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, A, eA, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    //XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT, B, eB, "This Xpetra::EpetraMultiVectorT method only accept Xpetra::EpetraMultiVectorT as input arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::norm1(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const { XPETRA_MONITOR("EpetraIntVectorT::norm1"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::norm2(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const { XPETRA_MONITOR("EpetraIntVectorT::norm2"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::normInf(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const { XPETRA_MONITOR("EpetraIntVectorT::normInf"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::normWeighted(const MultiVector<int,int,GlobalOrdinal,Node> &weights, const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const { XPETRA_MONITOR("EpetraIntVectorT::normWeighted"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::meanValue(const Teuchos::ArrayView<int> &means) const { XPETRA_MONITOR("EpetraIntVectorT::meanValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::maxValue(const Teuchos::ArrayView<int> &maxs) const { XPETRA_MONITOR("EpetraIntVectorT::maxValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const int &alpha, const MultiVector<int,int,GlobalOrdinal,Node> &A, const MultiVector<int,int,GlobalOrdinal,Node> &B, const int &beta) { XPETRA_MONITOR("EpetraIntVectorT::multiply"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Not available in Epetra"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::elementWiseMultiply(int scalarAB, const Vector<int,int,GlobalOrdinal,Node> &A, const MultiVector<int,int,GlobalOrdinal,Node> &B, int scalarThis) {
    XPETRA_MONITOR("EpetraIntVectorT::elementWiseMultiply");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra_EpetraIntVector: elementWiseMultiply not implemented because Epetra_IntVector does not support this operation");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) { XPETRA_MONITOR("EpetraIntVectorT::replaceGlobalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) { XPETRA_MONITOR("EpetraIntVectorT::sumIntoGlobalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) { XPETRA_MONITOR("EpetraIntVectorT::replaceLocalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) { XPETRA_MONITOR("EpetraIntVectorT::sumIntoLocalValue"); TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  template<class EpetraGlobalOrdinal>
  size_t EpetraIntVectorT<EpetraGlobalOrdinal>::getNumVectors() const { XPETRA_MONITOR("EpetraIntVectorT::getNumVectors"); return 1; }

  template<class EpetraGlobalOrdinal>
  std::string EpetraIntVectorT<EpetraGlobalOrdinal>::description() const {
    XPETRA_MONITOR("EpetraIntVectorT::description");

    // This implementation come from Epetra_Vector_def.hpp (without modification)
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{length="<<this->getGlobalLength()
        << "}";
    return oss.str();
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    XPETRA_MONITOR("EpetraIntVectorT::describe");

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

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::doImport(const DistObject<int, int, GlobalOrdinal, Node> &source,
                                 const Import<int, GlobalOrdinal, Node> &importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVectorT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVectorT<GlobalOrdinal>, source, tSource, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraIntVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal>, importer, tImporter, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntVector & v = *tSource.getEpetra_IntVector();
    int err = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::doExport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &dest,
                                 const Import<int, GlobalOrdinal, Node>& importer, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVectorT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVectorT<GlobalOrdinal>, dest, tDest, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraIntVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImportT<GlobalOrdinal>, importer, tImporter, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntVector & v = *tDest.getEpetra_IntVector();
    int err = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::doImport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &source,
                                 const Export<int, GlobalOrdinal, Node>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVectorT::doImport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVectorT<GlobalOrdinal>, source, tSource, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraIntVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal>, exporter, tExporter, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntVector & v = *tSource.getEpetra_IntVector();
    int err = vec_->Import(v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::doExport(const DistObject<int, LocalOrdinal, GlobalOrdinal, Node> &dest,
                                 const Export<int, GlobalOrdinal, Node>& exporter, CombineMode CM) {
    XPETRA_MONITOR("EpetraIntVectorT::doExport");

    XPETRA_DYNAMIC_CAST(const EpetraIntVectorT<GlobalOrdinal>, dest, tDest, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraIntVectorT as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExportT<GlobalOrdinal>, exporter, tExporter, "Xpetra::EpetraIntVectorT::doImport only accept Xpetra::EpetraImportT as input arguments.");

    const Epetra_IntVector & v = *tDest.getEpetra_IntVector();
    int err = vec_->Export(v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  template<class EpetraGlobalOrdinal>
  void EpetraIntVectorT<EpetraGlobalOrdinal>::
  assign (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal>& rhs)
  {
    typedef EpetraIntVectorT<GlobalOrdinal> this_type;
    const this_type* rhsPtr = dynamic_cast<const this_type*> (&rhs);
    TEUCHOS_TEST_FOR_EXCEPTION(
      rhsPtr == NULL, std::invalid_argument, "Xpetra::MultiVector::operator=: "
      "The left-hand side (LHS) of the assignment has a different type than "
      "the right-hand side (RHS).  The LHS has type Xpetra::EpetraIntVectorT "
      "(which means it wraps an Epetra_IntVector), but the RHS has some "
      "other type.  This probably means that the RHS wraps either an "
      "Tpetra::MultiVector, or an Epetra_MultiVector.  Xpetra::MultiVector "
      "does not currently implement assignment from a Tpetra object to an "
      "Epetra object, though this could be added with sufficient interest.");

    RCP<const Epetra_IntVector> rhsImpl = rhsPtr->getEpetra_IntVector ();
    RCP<Epetra_IntVector> lhsImpl = this->getEpetra_IntVector ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      rhsImpl.is_null (), std::logic_error, "Xpetra::MultiVector::operator= "
      "(in Xpetra::EpetraIntVectorT::assign): *this (the right-hand side of "
      "the assignment) has a null RCP<Epetra_IntVector> inside.  Please "
      "report this bug to the Xpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      lhsImpl.is_null (), std::logic_error, "Xpetra::MultiVector::operator= "
      "(in Xpetra::EpetraIntVectorT::assign): The left-hand side of the "
      "assignment has a null RCP<Epetra_IntVector> inside.  Please report "
      "this bug to the Xpetra developers.");

    // Epetra_IntVector's assignment operator does a deep copy.
    *lhsImpl = *rhsImpl;
  }

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
template class EpetraIntVectorT<int>;
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
template class EpetraIntVectorT<long long>;
#endif

} // namespace Xpetra
