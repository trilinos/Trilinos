#ifndef CTHULHU_TPETRAVECTOR_DECL_HPP
#define CTHULHU_TPETRAVECTOR_DECL_HPP

#ifndef HAVE_CTHULHU_TPETRA
#error This file should be included only if HAVE_CTHULHU_TPETRA is defined.
#endif

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_MultiVector.hpp"
#include "Cthulhu_Vector.hpp"
#include "Cthulhu_Exceptions.hpp"

#include "Cthulhu_ConfigDefs.hpp"

#include "Cthulhu_TpetraMap.hpp"
#include "Tpetra_Vector.hpp"

namespace Cthulhu {

  //! \brief A class for constructing and using dense, distributors vectors.
  /*!
    This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
    The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
    type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template<class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal, class Node=Kokkos::DefaultNode::DefaultNodeType>
  class TpetraVector : public Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> { // TODO public::TpetraMultiVector ?
    
    // // need this so that MultiVector::operator() can call Vector's private constructor
    // friend class MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

    // The following typedef are used by the CTHULHU_DYNAMIC_CAST() macro.
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;
    
  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Sets all vector entries to zero.
    explicit TpetraVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, bool zeroOut=true) {
      CTHULHU_DEBUG_ME;
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMap, map, tMap, "Cthulhu::TpetraVector constructors only accept Cthulhu::TpetraMap as input arguments.");
      vec_ = rcp(new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(tMap->getTpetra_Map(), zeroOut));
    }
    
#ifdef CTHULHU_NOT_IMPLEMENTED
    //! TpetraVector copy constructor.
    TpetraVector(const TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source);
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! \brief Set multi-vector values from an array using Teuchos memory management classes. (copy)
    TpetraVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, const Teuchos::ArrayView<const Scalar> &A);
#endif

    //! Destructor.  
    inline ~TpetraVector() { CTHULHU_DEBUG_ME; };

    //@}

    //! @name Post-construction modification routines
    //@{ 

    //! Replace current value at the specified location with specified value.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { CTHULHU_DEBUG_ME; vec_->replaceGlobalValue(globalRow, value); };

    //! Adds specified value to existing value at the specified location.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value) { CTHULHU_DEBUG_ME; vec_->sumIntoGlobalValue(globalRow, value); };

    //! Replace current value at the specified location with specified values.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void replaceLocalValue(LocalOrdinal myRow, const Scalar &value) { CTHULHU_DEBUG_ME; vec_->replaceLocalValue(myRow, value); };

    //! Adds specified value to existing value at the specified location.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value) { CTHULHU_DEBUG_ME; vec_->sumIntoLocalValue(myRow, value); };

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
#ifdef CTHULHU_TODO
    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot; // overloading, not hiding
    //! Computes dot product of this Vector against input Vector x.
    Scalar dot(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &a) const;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1; // overloading, not hiding
    //! Return 1-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm1() const;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2; // overloading, not hiding
    //! Compute 2-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm2() const;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf; // overloading, not hiding
    //! Compute Inf-norm of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType normInf() const;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normWeighted; // overloading, not hiding
    //! Compute Weighted 2-norm (RMS Norm) of this Vector.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType normWeighted(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights) const;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::meanValue; // overloading, not hiding
    //! Compute mean (average) value of this Vector.
    Scalar meanValue() const;
#endif // CTHULHU_TODO
    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    inline std::string description() const { CTHULHU_DEBUG_ME; return vec_->description(); };

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { CTHULHU_DEBUG_ME; return vec_->describe(out, verbLevel); };

    //@}

//   protected:

//     //! Advanced constructor accepting parallel buffer view.
//     Vector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, Teuchos::ArrayRCP<Scalar> data) { CTHULHU_DEBUG_ME; vec_->(); };

    RCP< Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getTpetra_Vector() const { CTHULHU_DEBUG_ME; return vec_; }
    
  private:
    RCP< Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > vec_;


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

#endif // CTHULHU_VECTOR_DECL_HPP
