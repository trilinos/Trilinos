#ifndef CTHULHU_EPETRAVECTOR_HPP
#define CTHULHU_EPETRAVECTOR_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_EPETRA
#error This file should be included only if HAVE_CTHULHU_EPETRA is defined.
#endif

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_MultiVector.hpp"
#include "Cthulhu_Vector.hpp"
#include "Cthulhu_Exceptions.hpp"

#include "Cthulhu_EpetraMap.hpp"
#include "Cthulhu_EpetraMultiVector.hpp"
#include "Cthulhu_EpetraImport.hpp"
#include "Cthulhu_EpetraExport.hpp"

#include "Cthulhu_CombineMode.hpp"

#include "Epetra_Vector.h"


namespace Cthulhu {

  //! \brief A class for constructing and using dense, distributors vectors.
  class EpetraVector : public Vector<double,int,int>,  public EpetraMultiVector {
    
  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Sets all vector entries to zero.
    explicit EpetraVector(const Teuchos::RCP<const Map<int,int> > &map, bool zeroOut=true) 
      : EpetraMultiVector(map,1,zeroOut)
    {
      CTHULHU_DEBUG_ME;
    }
    
#ifdef CTHULHU_NOT_IMPLEMENTED
    //! EpetraVector copy constructor.
    EpetraVector(const EpetraVector &source);
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! \brief Set multi-vector values from an array using Teuchos memory management classes. (copy)
    EpetraVector(const Teuchos::RCP<const Map<int,int> > &map, const Teuchos::ArrayView<const double> &A);
#endif

    //! Destructor.  
    inline ~EpetraVector() { CTHULHU_DEBUG_ME; };

    //@}

    //! @name Post-construction modification routines
    //@{ 
#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Replace current value at the specified location with specified value.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void replaceGlobalValue(int globalRow, const double &value) { CTHULHU_DEBUG_ME; this->EpetraMultiVector::replaceGlobalValue(globalRow,0,value); };

    //! Adds specified value to existing value at the specified location.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void sumIntoGlobalValue(int globalRow, const double &value) { CTHULHU_DEBUG_ME; this->EpetraMultiVector::sumIntoGlobalValue(globalRow, 0, value); };

    //! Replace current value at the specified location with specified values.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void replaceLocalValue(int myRow, const double &value) { CTHULHU_DEBUG_ME; this->EpetraMultiVector::replaceLocalValue(myRow, 0, value); };

    //! Adds specified value to existing value at the specified location.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void sumIntoLocalValue(int myRow, const double &value) { CTHULHU_DEBUG_ME; this->EpetraMultiVector::sumIntoLocalValue(myRow, 0, value); };
#endif
    //@}

#ifdef CTHULHU_TODO
    //! @name Extraction methods
    //@{
    using MultiVector<double,int,int>::get1dCopy; // overloading, not hiding
    //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
    void get1dCopy(Teuchos::ArrayView<double> A) const;
    //@}
#endif

    //! @name Mathematical methods
    //@{ 
    using EpetraMultiVector::dot; // overloading, not hiding
    //! Computes dot product of this Vector against input Vector x.
    double dot(const Vector<double,int,int> &a) const { 
      CTHULHU_DEBUG_ME; 
      CTHULHU_DYNAMIC_CAST(const EpetraVector, a, tA, "This Cthulhu::EpetraVector method only accept Cthulhu::EpetraVector as input arguments.");
      //      return getEpetra_Vector()->Dot(*tA.getEpetra_Vector()); 

      // other way: use the MultiVector Dot instead of VectorDot:
      double r;
      getEpetra_MultiVector()->Epetra_MultiVector::Dot(*tA.getEpetra_MultiVector(), &r); 
      return r;
    }

    using EpetraMultiVector::norm1; // overloading, not hiding
    //! Return 1-norm of this Vector.
    Teuchos::ScalarTraits<double>::magnitudeType norm1() const { CTHULHU_DEBUG_ME; double r; getEpetra_MultiVector()->Norm1(&r); return r; }

    using EpetraMultiVector::norm2; // overloading, not hiding
    //! Compute 2-norm of this Vector.
    Teuchos::ScalarTraits<double>::magnitudeType norm2() const { CTHULHU_DEBUG_ME; double r; getEpetra_MultiVector()->Norm2(&r); return r; }

    using EpetraMultiVector::normInf; // overloading, not hiding
    //! Compute Inf-norm of this Vector.
    Teuchos::ScalarTraits<double>::magnitudeType normInf() const { CTHULHU_DEBUG_ME; double r; getEpetra_MultiVector()->NormInf(&r); return r; }

    using EpetraMultiVector::normWeighted; // overloading, not hiding
    //! Compute Weighted 2-norm (RMS Norm) of this Vector.
    Teuchos::ScalarTraits<double>::magnitudeType normWeighted(const Vector<double,int,int> &weights) const { CTHULHU_DEBUG_ME;
      CTHULHU_DYNAMIC_CAST(const EpetraVector, weights, tWeights, "This Cthulhu::EpetraVector method only accept Cthulhu::EpetraVector as input arguments.");
      double r; 
      getEpetra_MultiVector()->NormWeighted(*tWeights.getEpetra_MultiVector(), &r); return r; 
    }

    using EpetraMultiVector::meanValue; // overloading, not hiding
    //! Compute mean (average) value of this Vector.
    double meanValue() const { CTHULHU_DEBUG_ME; double r; getEpetra_MultiVector()->MeanValue(&r); return r; }

    // Note: Added to Cthulhu. Not present in Tpetra
    using EpetraMultiVector::maxValue; // overloading, not hiding
    //! Compute max value of this Vector.
    double maxValue() const { CTHULHU_DEBUG_ME; double r; getEpetra_MultiVector()->MaxValue(&r); return r; }
    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    inline std::string description() const { 
      CTHULHU_DEBUG_ME; 

      // This implementation come from Epetra_Vector_def.hpp (without modification)
      std::ostringstream oss;
      oss << Teuchos::Describable::description();
      oss << "{length="<<this->getGlobalLength()
      << "}";
      return oss.str();
    };

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { 
      CTHULHU_DEBUG_ME; 
      
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

    //@}

    // protected:

    //     //! Advanced constructor accepting parallel buffer view.
    //     Vector(const Teuchos::RCP<const Map<int,int> > &map, Teuchos::ArrayRCP<double> data) { CTHULHU_DEBUG_ME; vec_->(); };

    Epetra_Vector * getEpetra_Vector() const { CTHULHU_DEBUG_ME; return (*this->EpetraMultiVector::getEpetra_MultiVector())(0); }

    using EpetraMultiVector::doImport; // overloading, not hiding
    void doImport(const Vector<double, int, int> &source, 
                  const Import<int, int> &importer, CombineMode CM) {
      CTHULHU_DEBUG_ME;

      CTHULHU_DYNAMIC_CAST(const EpetraVector, source, tSource, "Cthulhu::EpetraVector::doImport only accept Cthulhu::EpetraVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Cthulhu::EpetraVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

      const Epetra_Vector * v = tSource.getEpetra_Vector();
      int err = this->EpetraMultiVector::getEpetra_MultiVector()->Import(*v, *tImporter.getEpetra_Import(), Cthulhu2Epetra_CombineMode(CM));
      TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
    }

    using EpetraMultiVector::doExport; // overloading, not hiding
    void doExport(const Vector<double, int, int> &dest,
                  const Import<int, int>& importer, CombineMode CM) {
      CTHULHU_DEBUG_ME;
      
      CTHULHU_DYNAMIC_CAST(const EpetraVector, dest, tDest, "Cthulhu::EpetraVector::doImport only accept Cthulhu::EpetraVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Cthulhu::EpetraVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

      const Epetra_Vector * v = tDest.getEpetra_Vector();
      int err = this->EpetraMultiVector::getEpetra_MultiVector()->Export(*v, *tImporter.getEpetra_Import(), Cthulhu2Epetra_CombineMode(CM)); 
      TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
    }


    //    using EpetraMultiVector::doExport; // overloading, not hiding
    void doExport(const Vector<double, int, int> &dest,
                  const Export<int, int>& exporter, CombineMode CM) {
      CTHULHU_DEBUG_ME;
      
      CTHULHU_DYNAMIC_CAST(const EpetraVector, dest, tDest, "Cthulhu::EpetraVector::doImport only accept Cthulhu::EpetraVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Cthulhu::EpetraVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

      const Epetra_Vector * v = tDest.getEpetra_Vector();
      int err = this->EpetraMultiVector::getEpetra_MultiVector()->Export(*v, *tExporter.getEpetra_Export(), Cthulhu2Epetra_CombineMode(CM)); 
      TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
    }



    
  }; // class EpetraVector

#ifdef CTHULHU_NOT_IMPLEMENTED
  /** \brief Non-member function to create a Vector from a specified Map.
  
  \relates Vector
  */
  template <class double, class int, class int, class Node>
  Teuchos::RCP< Vector<double,int,int> >
  createVector(const Teuchos::RCP< const Map<int,int> > &map);
#endif
} // namespace Cthulhu

#endif // CTHULHU_VECTOR_DECL_HPP
