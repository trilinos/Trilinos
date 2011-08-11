#ifndef CTHULHU_VECTOR_HPP
#define CTHULHU_VECTOR_HPP

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_MultiVector.hpp"
#include "Cthulhu_Import.hpp"
#include "Cthulhu_Export.hpp"

namespace Cthulhu {

  //! \brief A class for constructing and using dense, distributors vectors.
  /*!
    This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
    The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
    type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template<class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal, class Node=Kokkos::DefaultNode::DefaultNodeType>
  class Vector : public virtual MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

    // // need this so that MultiVector::operator() can call Vector's private constructor
    // friend class MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    
  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Destructor.  
    virtual ~Vector() {  }

    //@}

    //! @name Mathematical methods
    //@{ 
    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot; // overloading, not hiding
    //! Computes dot product of this Vector against input Vector x.
    virtual Scalar dot(const Cthulhu::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &a) const = 0;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1; // overloading, not hiding
    //! Return 1-norm of this Vector.
    virtual typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm1() const = 0;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2; // overloading, not hiding
    //! Compute 2-norm of this Vector.
    virtual typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm2() const = 0;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf; // overloading, not hiding
    //! Compute Inf-norm of this Vector.
    virtual typename Teuchos::ScalarTraits<Scalar>::magnitudeType normInf() const = 0;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normWeighted; // overloading, not hiding
    //! Compute Weighted 2-norm (RMS Norm) of this Vector.
    virtual typename Teuchos::ScalarTraits<Scalar>::magnitudeType normWeighted(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights) const = 0;

    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::meanValue; // overloading, not hiding
    //! Compute mean (average) value of this Vector.
    virtual Scalar meanValue() const = 0;

    // Note: Added to Cthulhu. Not present in Tpetra
    using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::maxValue; // overloading, not hiding
    //! Compute max value of this Vector.
    virtual Scalar maxValue() const = 0;

    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    virtual std::string description() const = 0;

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const = 0;

    //@}

  }; // class Vector

}

#define CTHULHU_VECTOR_SHORT
#endif // CTHULHU_VECTOR_DECL_HPP
