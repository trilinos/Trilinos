#ifndef CTHULHU_EPETRAVECTOR_HPP
#define CTHULHU_EPETRAVECTOR_HPP

#include "Cthulhu_EpetraConfigDefs.hpp"

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

  class EpetraVector
    : public Vector<double,int,int>,  public EpetraMultiVector
  {
    
  public:
    using EpetraMultiVector::dot; // overloading, not hiding
    using EpetraMultiVector::norm1; // overloading, not hiding
    using EpetraMultiVector::norm2; // overloading, not hiding
    using EpetraMultiVector::normInf; // overloading, not hiding
    using EpetraMultiVector::normWeighted; // overloading, not hiding
    using EpetraMultiVector::meanValue; // overloading, not hiding

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Sets all vector entries to zero.
    explicit EpetraVector(const Teuchos::RCP<const Map<int,int> > &map, bool zeroOut=true) 
      : EpetraMultiVector(map,1,zeroOut)
    {
      
    }
    
    //! EpetraVector constructor.
    EpetraVector(const RCP<Epetra_Vector> &source)
    : EpetraMultiVector(source)
    {
      
    }

    //! Destructor.  
    ~EpetraVector() {  };

    //@}

    //! @name Mathematical methods
    //@{ 

    //! Computes dot product of this Vector against input Vector x.
    double dot(const Vector<double,int,int> &a) const;

    //! Return 1-norm of this Vector.
    Teuchos::ScalarTraits<double>::magnitudeType norm1() const;

    //! Compute 2-norm of this Vector.
    Teuchos::ScalarTraits<double>::magnitudeType norm2() const;

    //! Compute Inf-norm of this Vector.
    Teuchos::ScalarTraits<double>::magnitudeType normInf() const;

    //! Compute Weighted 2-norm (RMS Norm) of this Vector.
    Teuchos::ScalarTraits<double>::magnitudeType normWeighted(const Vector<double,int,int> &weights) const;


    //! Compute mean (average) value of this Vector.
    double meanValue() const;

    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const;

    //@}

    Epetra_Vector * getEpetra_Vector() const {  return (*this->EpetraMultiVector::getEpetra_MultiVector())(0); }

  }; // class EpetraVector

} // namespace Cthulhu

#endif // CTHULHU_EPETRAVECTOR_HPP
