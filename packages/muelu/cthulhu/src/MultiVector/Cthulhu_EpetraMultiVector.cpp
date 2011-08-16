#include "Cthulhu_EpetraMultiVector.hpp"

#include "Cthulhu_EpetraImport.hpp" 
#include "Cthulhu_EpetraExport.hpp" 

namespace Cthulhu {

  Teuchos::ArrayRCP<const double> EpetraMultiVector::getData(size_t j) const { 
    double ** arrayOfPointers;
      
    vec_->ExtractView(&arrayOfPointers);
     
    double * data = arrayOfPointers[j];
    int localLength = vec_->MyLength();
      
    return ArrayRCP<double>(data, 0, localLength, false); // not ownership
  }

  Teuchos::ArrayRCP<double> EpetraMultiVector::getDataNonConst(size_t j) { 
    double ** arrayOfPointers;
      
    vec_->ExtractView(&arrayOfPointers);
     
    double * data = arrayOfPointers[j];
    int localLength = vec_->MyLength();
      
    return ArrayRCP<double>(data, 0, localLength, false); // not ownership
  }

  void EpetraMultiVector::dot(const MultiVector<double,int,int> &A, const Teuchos::ArrayView<double> &dots) const { 
    CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
    vec_->Dot(*eA.getEpetra_MultiVector(), dots.getRawPtr());
  }

  void EpetraMultiVector::normWeighted(const MultiVector<double,int,int> &weights, const Teuchos::ArrayView<Teuchos::ScalarTraits<double>::magnitudeType> &norms) const { 
       
    CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, weights, eWeights, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
    vec_->NormWeighted(*eWeights.getEpetra_MultiVector(), norms.getRawPtr()); 
  }

  void EpetraMultiVector::meanValue(const Teuchos::ArrayView<double> &means) const {  vec_->MeanValue(means.getRawPtr()); } //TODO: modify ArrayView size ??

  std::string EpetraMultiVector::description() const {  
    TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO");
    return "TODO"; 
  }

  void EpetraMultiVector::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const { 
    vec_->Print(out);
  }

  void EpetraMultiVector::doImport(const DistObject<double, int, int> &source, 
                const Import<int, int> &importer, CombineMode CM) {
      
    CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, source, tSource, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraMultiVector as input arguments.");
    CTHULHU_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

    RCP<Epetra_MultiVector> v = tSource.getEpetra_MultiVector();
    int err = this->getEpetra_MultiVector()->Import(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraMultiVector::doExport(const DistObject<double, int, int> &dest,
                const Import<int, int>& importer, CombineMode CM) {
      
    CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, dest, tDest, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraMultiVector as input arguments.");
    CTHULHU_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

    RCP<Epetra_MultiVector> v = tDest.getEpetra_MultiVector();
    int err = this->getEpetra_MultiVector()->Export(*v, *tImporter.getEpetra_Import(), toEpetra(CM)); 
    TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraMultiVector::doImport(const DistObject<double,int,int> &source,
                const Export<int, int>& exporter, CombineMode CM) {
      
    CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, source, tSource, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraMultiVector as input arguments.");
    CTHULHU_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

    RCP<Epetra_MultiVector> v = tSource.getEpetra_MultiVector();
    int err = this->getEpetra_MultiVector()->Import(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraMultiVector::doExport(const DistObject<double, int, int> &dest,
                const Export<int, int>& exporter, CombineMode CM) {
      
    CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, dest, tDest, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraMultiVector as input arguments.");
    CTHULHU_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

    RCP<Epetra_MultiVector> v = tDest.getEpetra_MultiVector();
    int err = this->getEpetra_MultiVector()->Export(*v, *tExporter.getEpetra_Export(), toEpetra(CM)); 
    TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  // TODO: move that elsewhere
  const Epetra_MultiVector & toEpetra(const MultiVector<double, int, int> & x) {
    CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, x, tX, "toEpetra");
    return *tX.getEpetra_MultiVector();
  }

  Epetra_MultiVector & toEpetra(MultiVector<double, int, int> & x) {
    CTHULHU_DYNAMIC_CAST(      EpetraMultiVector, x, tX, "toEpetra");
    return *tX.getEpetra_MultiVector();
  }
  //


} // namespace Cthulhu
