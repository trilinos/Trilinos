#ifndef CTHULHU_TPETRAVECTOR_HPP
#define CTHULHU_TPETRAVECTOR_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_TPETRA
#error This file should be included only if HAVE_CTHULHU_TPETRA is defined.
#endif

#include <Teuchos_ScalarTraitsDecl.hpp> //TODO: useless?

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_MultiVector.hpp"
#include "Cthulhu_Vector.hpp"
#include "Cthulhu_Exceptions.hpp"

#include "Cthulhu_TpetraMap.hpp"
#include "Cthulhu_TpetraMultiVector.hpp"
#include "Cthulhu_TpetraImport.hpp"
#include "Cthulhu_TpetraExport.hpp"

#include "Cthulhu_CombineMode.hpp"

#include "Tpetra_Vector.hpp"

namespace Cthulhu {

  //! \brief A class for constructing and using dense, distributors vectors.
  /*!
    This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
    The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
    type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template<class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal, class Node=Kokkos::DefaultNode::DefaultNodeType>
  class TpetraVector : public Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>, public TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
    
    // // need this so that MultiVector::operator() can call Vector's private constructor
    // friend class MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

    // The following typedef is used by the CTHULHU_DYNAMIC_CAST() macro.
    //typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass; //not needed??
    
  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Sets all vector entries to zero.
    TpetraVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, bool zeroOut=true)
      : TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map,1,zeroOut)
    {
      CTHULHU_DEBUG_ME;
    }
    
#ifdef CTHULHU_NOT_IMPLEMENTED
    //! TpetraVector copy constructor.
    TpetraVector(const TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source);
#endif

    //! \brief Set multi-vector values from an array using Teuchos memory management classes. (copy)
    TpetraVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, const Teuchos::ArrayView<const Scalar> &A)
      : TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map,A,map->getNodeNumElements(),1)
    {
      CTHULHU_DEBUG_ME;
    }
    
    /** \brief TpetraVector constructor to wrap a Tpetra::Vector object
     */
    TpetraVector(const Teuchos::RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &vec) : TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(vec) { CTHULHU_DEBUG_ME; }
    //TODO: removed const of Tpetra::Vector

    //! Destructor.  
    inline ~TpetraVector() { CTHULHU_DEBUG_ME; };

    //@}

    //! @name Post-construction modification routines
    //@{ 

    //! Replace current value at the specified location with specified value.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { CTHULHU_DEBUG_ME; this->TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValue(globalRow,0,value); };

    //! Adds specified value to existing value at the specified location.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { CTHULHU_DEBUG_ME; this->TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue(globalRow, 0, value); };

    //! Replace current value at the specified location with specified values.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void replaceLocalValue(LocalOrdinal myRow, const Scalar &value) { CTHULHU_DEBUG_ME; this->TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceLocalValue(myRow, 0, value); };

    //! Adds specified value to existing value at the specified location.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value) { CTHULHU_DEBUG_ME; this->TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoLocalValue(myRow, 0, value); };

    //@}

#ifdef CTHULHU_TODO
    //! @name Extraction methods
    //@{
    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get1dCopy; // overloading, not hiding
    //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
    void get1dCopy(Teuchos::ArrayView<Scalar> A) const;
    //@}
#endif

    //! @name Mathematical methods
    //@{ 
    using TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot; // overloading, not hiding
    //! Computes dot product of this Vector against input Vector x.
    Scalar dot(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &a) const { 
      CTHULHU_DEBUG_ME; 
      CTHULHU_DYNAMIC_CAST(const TpetraVector, a, tA, "This Cthulhu::TpetraVector method only accept Cthulhu::TpetraVector as input arguments.");
      return getTpetra_Vector()->dot(*tA.getTpetra_Vector()); 
    }

    using TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1; // overloading, not hiding
    //! Return 1-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm1() const { CTHULHU_DEBUG_ME; return getTpetra_Vector()->norm1(); }

    using TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2; // overloading, not hiding
    //! Compute 2-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm2() const { CTHULHU_DEBUG_ME; return getTpetra_Vector()->norm2(); }

    using TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf; // overloading, not hiding
    //! Compute Inf-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType normInf() const { CTHULHU_DEBUG_ME; return getTpetra_Vector()->normInf(); }

    using TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normWeighted; // overloading, not hiding
    //! Compute Weighted 2-norm (RMS Norm) of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType normWeighted(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights) const { 
      CTHULHU_DEBUG_ME; 
      CTHULHU_DYNAMIC_CAST(const TpetraVector, weights, tWeights, "This Cthulhu::TpetraVector method only accept Cthulhu::TpetraVector as input arguments.");
      return getTpetra_Vector()->normWeighted(*tWeights.getTpetra_Vector()); 
    }

    using TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::meanValue; // overloading, not hiding
    //! Compute mean (average) value of this Vector.
    Scalar meanValue() const { CTHULHU_DEBUG_ME; return getTpetra_Vector()->meanValue(); }

    // Note: Added to Cthulhu. Not present in Tpetra
    using TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::maxValue; // overloading, not hiding
    //! Compute max value of this Vector.
    Scalar maxValue() const { CTHULHU_DEBUG_ME; TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO");
      //return getTpetra_Vector()->maxValue(); 
    }
    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    inline std::string description() const { 
      CTHULHU_DEBUG_ME; 
      return getTpetra_Vector()->description();       
    };

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { 
      getTpetra_Vector()->describe(out, verbLevel);
    }

    //@}

    // protected:

    //     //! Advanced constructor accepting parallel buffer view.
    //     Vector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, Teuchos::ArrayRCP<Scalar> data) { CTHULHU_DEBUG_ME; vec_->(); };

    RCP< Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getTpetra_Vector() const { CTHULHU_DEBUG_ME; return this->TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getTpetra_MultiVector()->getVectorNonConst(0); }
    
  }; // class TpetraVector

#ifdef CTHULHU_NOT_IMPLEMENTED
  /** \brief Non-member function to create a Vector from a specified Map.
  
  \relates Vector
  */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  createVector(const Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> > &map);
#endif
} // namespace Cthulhu

#define CTHULHU_TPETRAVECTOR_SHORT
#endif // CTHULHU_VECTOR_DECL_HPP


//JG TODO: overload constructor to have a real vector as underlying obj
