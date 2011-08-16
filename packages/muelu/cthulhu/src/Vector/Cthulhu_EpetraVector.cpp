#include "Cthulhu_EpetraVector.hpp"

namespace Cthulhu {

  double EpetraVector::dot(const Vector<double,int,int> &a) const { 
       
      CTHULHU_DYNAMIC_CAST(const EpetraVector, a, tA, "This Cthulhu::EpetraVector method only accept Cthulhu::EpetraVector as input arguments.");
      //      return getEpetra_Vector()->Dot(*tA.getEpetra_Vector()); 

      // other way: use the MultiVector Dot instead of VectorDot:
      double r;
      getEpetra_MultiVector()->Epetra_MultiVector::Dot(*tA.getEpetra_MultiVector(), &r); 
      return r;
    }

    Teuchos::ScalarTraits<double>::magnitudeType EpetraVector::norm1() const {  double r; getEpetra_MultiVector()->Norm1(&r); return r; }

    Teuchos::ScalarTraits<double>::magnitudeType EpetraVector::norm2() const {  double r; getEpetra_MultiVector()->Norm2(&r); return r; }

    Teuchos::ScalarTraits<double>::magnitudeType EpetraVector::normInf() const {  double r; getEpetra_MultiVector()->NormInf(&r); return r; }

    Teuchos::ScalarTraits<double>::magnitudeType EpetraVector::normWeighted(const Vector<double,int,int> &weights) const { 
      CTHULHU_DYNAMIC_CAST(const EpetraVector, weights, tWeights, "This Cthulhu::EpetraVector method only accept Cthulhu::EpetraVector as input arguments.");
      double r; 
      getEpetra_MultiVector()->NormWeighted(*tWeights.getEpetra_MultiVector(), &r); return r; 
    }

    double EpetraVector::meanValue() const {  double r; getEpetra_MultiVector()->MeanValue(&r); return r; }

    std::string EpetraVector::description() const { 
      // This implementation come from Epetra_Vector_def.hpp (without modification)
      std::ostringstream oss;
      oss << Teuchos::Describable::description();
      oss << "{length="<<this->getGlobalLength()
      << "}";
      return oss.str();
    }

    void EpetraVector::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const { 
       
      
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

//       TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO");

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
    CTHULHU_DYNAMIC_CAST(      EpetraVector, x, tX, "toEpetra");
    return *tX.getEpetra_Vector();
  }
  
  const Epetra_Vector & toEpetra(const Vector<double, int, int> &x) {
    CTHULHU_DYNAMIC_CAST(const EpetraVector, x, tX, "toEpetra");
    return *tX.getEpetra_Vector();
  }
  //

}
