#ifndef CTHULHU_EPETRAVECTOR_DECL_HPP
#define CTHULHU_EPETRAVECTOR_DECL_HPP

#ifndef HAVE_CTHULHU_EPETRA
#error This file should be included only if HAVE_CTHULHU_EPETRA is defined.
#endif

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_MultiVector.hpp"
#include "Cthulhu_Vector.hpp"
#include "Cthulhu_Exceptions.hpp"

#include "Cthulhu_EpetraMap.hpp"
#include "Epetra_Vector.h"

namespace Cthulhu {

  //! \brief A class for constructing and using dense, distributors vectors.
  /*!
  */
  class EpetraVector : public Vector<double,int,int> { // TODO public::EpetraMultiVector ?
    
  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Sets all vector entries to zero.
    explicit EpetraVector(const Teuchos::RCP<const Map<int,int> > &map, bool zeroOut=true) {
      CTHULHU_DEBUG_ME;
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, map, eMap, "Cthulhu::EpetraVector constructors only accept Cthulhu::EpetraMap as input arguments.");
      vec_ = rcp(new Epetra_Vector(eMap->getEpetra_Map(), zeroOut));
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
    inline void replaceGlobalValue(int globalRow, const double &value) { CTHULHU_DEBUG_ME; vec_->replaceGlobalValue(globalRow, value); };
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Adds specified value to existing value at the specified location.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void sumIntoGlobalValue(int globalRow, const double &value) { CTHULHU_DEBUG_ME; vec_->sumIntoGlobalValue(globalRow, value); };
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Replace current value at the specified location with specified values.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void replaceLocalValue(int myRow, const double &value) { CTHULHU_DEBUG_ME; vec_->replaceLocalValue(myRow, value); };
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Adds specified value to existing value at the specified location.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void sumIntoLocalValue(int myRow, const double &value) { CTHULHU_DEBUG_ME; vec_->sumIntoLocalValue(myRow, value); };
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
#ifdef CTHULHU_TODO
    using MultiVector<double,int,int>::dot; // overloading, not hiding
    //! Computes dot product of this Vector against input Vector x.
    double dot(const Vector<double,int,int> &a) const;

    using MultiVector<double,int,int>::norm1; // overloading, not hiding
    //! Return 1-norm of this Vector.
    typename Teuchos::doubleTraits<double>::magnitudeType norm1() const;

    using MultiVector<double,int,int>::norm2; // overloading, not hiding
    //! Compute 2-norm of this Vector.
    typename Teuchos::doubleTraits<double>::magnitudeType norm2() const;

    using MultiVector<double,int,int>::normInf; // overloading, not hiding
    //! Compute Inf-norm of this Vector.
    typename Teuchos::doubleTraits<double>::magnitudeType normInf() const;

    using MultiVector<double,int,int>::normWeighted; // overloading, not hiding
    //! Compute Weighted 2-norm (RMS Norm) of this Vector.
    typename Teuchos::doubleTraits<double>::magnitudeType normWeighted(const Vector<double,int,int> &weights) const;

    using MultiVector<double,int,int>::meanValue; // overloading, not hiding
    //! Compute mean (average) value of this Vector.
    double meanValue() const;
#endif // CTHULHU_TODO
    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{
#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    /** \brief Return a simple one-line description of this object. */
    inline std::string description() const { CTHULHU_DEBUG_ME; return vec_->description(); };
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { CTHULHU_DEBUG_ME; return vec_->describe(out, verbLevel); };
#endif

    //@}

//   protected:

//     //! Advanced constructor accepting parallel buffer view.
//     Vector(const Teuchos::RCP<const Map<int,int> > &map, Teuchos::ArrayRCP<double> data) { CTHULHU_DEBUG_ME; vec_->(); };

    RCP<Epetra_Vector> getEpetra_Vector() const { CTHULHU_DEBUG_ME; return vec_; }
    
  private:
    RCP<Epetra_Vector> vec_;


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
