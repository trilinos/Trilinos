// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_COMPOBJECT_HPP_
#define _TEUCHOS_COMPOBJECT_HPP_

#include "Teuchos_Object.hpp"
#include "Teuchos_Flops.cpp"

namespace Teuchos
{


//! Teuchos_CompObject: Functionality and data that is common to all computational classes.

/*! The Teuchos_CompObject is a base class for all Teuchos computational objects.  It provides the basic
    mechanisms and interface specifications for floating point operations using Teuchos_Flops objects. It
    should be noted that currently, Teuchos_CompObject is an almost exact duplicate of Epetra_CompObject.
*/

class CompObject
{

  public:

  //@{ \name Constructors/Destructor.
  //! Basic Teuchos::CompObject constuctor.
  CompObject();

  //! Teuchos::CompObject copy constructor.
  
  CompObject(const CompObject &source);
  
  
  //! Teuchos::CompObject destructor.  
  virtual ~CompObject();
  //@}

  //@{ \name Set/Get counter method.
  //! Set the internal Teuchos::Flops() pointer.
  void setFlopCounter(const Flops &FlopCounter);
  
  //! Set the internal Teuchos::Flops() pointer to the flop counter of another Teuchos::CompObject.
  void setFlopCounter(const CompObject &CompObject);
  
  //! Set the internal Teuchos::Flops() pointer to 0 (no flops counted).
  void unsetFlopCounter();
  
  //! Get the pointer to the Teuchos::Flops() object associated with this object, returns 0 if none.
  Flops * getFlopCounter() const;
  //@}

  //@{ \name Set flop count methods.
  //! Resets the number of floating point operations to zero for \e this multi-vector.
  void resetFlops() const;

  //! Returns the number of floating point operations with \e this multi-vector.
  double flops() const;
  //@}

  //@{ \name Update flop count methods.
  //! Increment Flop count for \e this object
  void updateFlops(int Flops) const;

  //! Increment Flop count for \e this object
  void updateFlops(long int Flops) const;

  //! Increment Flop count for \e this object
  void updateFlops(double Flops) const;

  //! Increment Flop count for \e this object
  void updateFlops(float Flops) const;
  //@}

 protected:


  Flops *flopCounter_;

};

  // #include "Teuchos_CompObject.cpp"

} // namespace Teuchos

#endif // end of _TEUCHOS_COMPOBJECT_HPP_
