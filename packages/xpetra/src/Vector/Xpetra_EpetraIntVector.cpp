#include "Xpetra_EpetraIntVector.hpp"
#include "Xpetra_EpetraImport.hpp"
#include "Xpetra_EpetraExport.hpp"

namespace Xpetra {

  int EpetraIntVector::dot(const Vector<int,int,int> &a) const { 
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
    return -1;
  }

  Teuchos::ScalarTraits<int>::magnitudeType EpetraIntVector::norm1() const {  TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  Teuchos::ScalarTraits<int>::magnitudeType EpetraIntVector::norm2() const {  TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  Teuchos::ScalarTraits<int>::magnitudeType EpetraIntVector::normInf() const {  TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  Teuchos::ScalarTraits<int>::magnitudeType EpetraIntVector::normWeighted(const Vector<int,int,int> &weights) const { TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  int EpetraIntVector::meanValue() const {  TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  int EpetraIntVector::maxValue() const {  TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); return -1; }

  void EpetraIntVector::randomize() {   TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVector::randomize(): Functionnality not available in Epetra"); }

  void EpetraIntVector::setSeed(unsigned int seed) {
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraIntVector::setSeed(): Functionnality not available in Epetra");
  }

  Teuchos::ArrayRCP<const int> EpetraIntVector::getData(size_t j) const {  
    int * data = vec_->Values();
    int localLength = vec_->MyLength();
      
    return ArrayRCP<int>(data, 0, localLength, false); // not ownership
  }

  Teuchos::ArrayRCP<int> EpetraIntVector::getDataNonConst(size_t j) { 
    int * data = vec_->Values();
    int localLength = vec_->MyLength();
      
    return ArrayRCP<int>(data, 0, localLength, false); // not ownership
  }

  void EpetraIntVector::dot(const MultiVector<int,int,int> &A, const Teuchos::ArrayView<int> &dots) const { 
    //XPETRA_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  void EpetraIntVector::abs(const MultiVector<int,int,int> &A) { 
    //XPETRA_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  void EpetraIntVector::reciprocal(const MultiVector<int,int,int> &A) { 
    //XPETRA_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  void EpetraIntVector::scale(const int &alpha) {  TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::update(const int &alpha, const MultiVector<int,int,int> &A, const int &beta) { 
    // XPETRA_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  void EpetraIntVector::update(const int &alpha, const MultiVector<int,int,int> &A, const int &beta, const MultiVector<int,int,int> &B, const int &gamma) {
    //XPETRA_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    //XPETRA_DYNAMIC_CAST(const EpetraMultiVector, B, eB, "This Xpetra::EpetraMultiVector method only accept Xpetra::EpetraMultiVector as input arguments.");
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO");
  }

  void EpetraIntVector::norm1(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const {  TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::norm2(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const {  TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::normInf(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const {  TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::normWeighted(const MultiVector<int,int,int> &weights, const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const { TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::meanValue(const Teuchos::ArrayView<int> &means) const {        TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::maxValue(const Teuchos::ArrayView<int> &maxs) const {        TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO"); }

  void EpetraIntVector::multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const int &alpha, const MultiVector<int,int,int> &A, const MultiVector<int,int,int> &B, const int &beta) { TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Not available in Epetra"); }

  void EpetraIntVector::elementWiseMultiply(int scalarAB, const Vector<int,int,int> &A, const MultiVector<int,int,int> &B, int scalarThis) {   
    TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "Xpetra_EpetraIntVector: elementWiseMultiply not implemented because Epetra_IntVector does not support this operation");
  }

  size_t EpetraIntVector::getNumVectors() const { return 1; }

  std::string EpetraIntVector::description() const {  
    // This implementation come from Epetra_Vector_def.hpp (without modification)
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{length="<<this->getGlobalLength()
        << "}";
    return oss.str();
  }

  void EpetraIntVector::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const { 
//     typedef Kokkos::MultiVector<double> KMV;
//     typedef Kokkos::DefaultArithmetic<KMV>   MVT;

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

    XPETRA_DYNAMIC_CAST(const EpetraIntVector, source, tSource, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraIntVector as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraImport as input arguments.");

    const Epetra_IntVector & v = *tSource.getEpetra_IntVector();
    int err = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM)); 
    TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraIntVector::doExport(const DistObject<int, int, int> &dest,
                                 const Import<int, int>& importer, CombineMode CM) {

    XPETRA_DYNAMIC_CAST(const EpetraIntVector, dest, tDest, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraIntVector as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraImport as input arguments.");

    const Epetra_IntVector & v = *tDest.getEpetra_IntVector();
    int err = vec_->Import(v, *tImporter.getEpetra_Import(), toEpetra(CM)); 
    TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraIntVector::doImport(const DistObject<int, int, int> &source,
                                 const Export<int, int>& exporter, CombineMode CM) {
    
    XPETRA_DYNAMIC_CAST(const EpetraIntVector, source, tSource, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraIntVector as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraImport as input arguments.");
    
    const Epetra_IntVector & v = *tSource.getEpetra_IntVector();
    int err = vec_->Import(v, *tExporter.getEpetra_Export(), toEpetra(CM)); 
    TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraIntVector::doExport(const DistObject<int, int, int> &dest,
                                 const Export<int, int>& exporter, CombineMode CM) {

    XPETRA_DYNAMIC_CAST(const EpetraIntVector, dest, tDest, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraIntVector as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Xpetra::EpetraIntVector::doImport only accept Xpetra::EpetraImport as input arguments.");

    const Epetra_IntVector & v = *tDest.getEpetra_IntVector();
    int err = vec_->Export(v, *tExporter.getEpetra_Export(), toEpetra(CM)); 
    TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

} // namespace Xpetra
