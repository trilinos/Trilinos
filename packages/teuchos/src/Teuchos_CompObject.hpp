// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef TEUCHOS_COMPOBJECT_HPP
#define TEUCHOS_COMPOBJECT_HPP

#include "Teuchos_Object.hpp"
#include "Teuchos_Flops.hpp"

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
  void setFlopCounter(const Flops &FlopCounter) {flopCounter_= (Flops *) &FlopCounter; return;};
  
  //! Set the internal Teuchos::Flops() pointer to the flop counter of another Teuchos::CompObject.
  void setFlopCounter(const CompObject &compObject) {flopCounter_= (Flops *) (compObject.getFlopCounter()); return;};
  
  //! Set the internal Teuchos::Flops() pointer to 0 (no flops counted).
  void unsetFlopCounter() {flopCounter_=0; return;};
  
  //! Get the pointer to the Teuchos::Flops() object associated with this object, returns 0 if none.
  Flops * getFlopCounter() const {return(flopCounter_);};
  //@}

  //@{ \name Set flop count methods.
  //! Resets the number of floating point operations to zero for \e this multi-vector.
  void resetFlops() const {if (flopCounter_!=0) flopCounter_->resetFlops(); return;};

  //! Returns the number of floating point operations with \e this multi-vector.
  double flops() const {if (flopCounter_!=0) return(flopCounter_->flops()); else return(0.0);};
  //@}

  //@{ \name Update flop count methods.
  //! Increment Flop count for \e this object
  void updateFlops(int addflops) const { if (flopCounter_!=0) flopCounter_->updateFlops(addflops); return;};

  //! Increment Flop count for \e this object
  void updateFlops(long int addflops) const { if (flopCounter_!=0) flopCounter_->updateFlops(addflops); return;};

  //! Increment Flop count for \e this object
  void updateFlops(double addflops) const { if (flopCounter_!=0) flopCounter_->updateFlops(addflops); return;};

  //! Increment Flop count for \e this object
  void updateFlops(float addflops) const {if (flopCounter_!=0) flopCounter_->updateFlops(addflops); return;};
  //@}

 protected:


  Flops *flopCounter_;

};

  // #include "Teuchos_CompObject.cpp"

} // namespace Teuchos

#endif // end of TEUCHOS_COMPOBJECT_HPP
