/*Paul
27-May-2002 General cleanup. Checked for newNamingConvention (already done). Moved some code into Tpetra_CompObject.cpp
06-August-2002 Changed to images (nothing changed).
*/

#ifndef _TPETRA_COMPOBJECT_H_
#define _TPETRA_COMPOBJECT_H_

#include "Tpetra_Object.h"
#include "Tpetra_Flops.h"

namespace Tpetra
{


//! Tpetra_CompObject: Functionality and data that is common to all computational classes.

/*! The Tpetra_CompObject is a base class for all Tpetra computational objects.  It provides the basic
    mechanisms and interface specifications for floating point operations using Tpetra_Flops objects. It
    should be noted that currently, Tpetra_CompObject is an almost exact duplicate of Epetra_CompObject.
*/

class CompObject
{

  public:

  //@{ \name Constructors/Destructor.
  //! Basic Tpetra::CompObject constuctor.
  CompObject();

  //! Tpetra::CompObject copy constructor.
  
  CompObject(const CompObject &source);
  
  
  //! Tpetra::CompObject destructor.  
  virtual ~CompObject();
  //@}

  //@{ \name Set/Get counter method.
  //! Set the internal Tpetra::Flops() pointer.
  void setFlopCounter(const Flops &FlopCounter);
  
  //! Set the internal Tpetra::Flops() pointer to the flop counter of another Tpetra::CompObject.
  void setFlopCounter(const CompObject &CompObject);
  
  //! Set the internal Tpetra::Flops() pointer to 0 (no flops counted).
  void unsetFlopCounter();
  
  //! Get the pointer to the Tpetra::Flops() object associated with this object, returns 0 if none.
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


} // namespace Tpetra

#include "Tpetra_CompObject.cpp"

#endif // end of _TPETRA_COMPOBJECT_H_
