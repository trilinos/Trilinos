
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_COMPOBJECT_H_
#define _EPETRA_COMPOBJECT_H_

//! Epetra_CompObject: Functionality and data that is common to all computational classes.

/*! The Epetra_CompObject is a base class for all Epetra computational objects.  It provides the basic
    mechanisms and interface specifications for floating point operations using Epetra_Flops objects.

*/
#include "Epetra_Object.h"
#include "Epetra_Flops.h"
//==========================================================================
class Epetra_CompObject {

  public:

  //@{ \name Constructors/Destructor.
  //! Basic Epetra_CompObject constuctor.
  Epetra_CompObject();

  //! Epetra_CompObject copy constructor.
  
  Epetra_CompObject(const Epetra_CompObject& Source);
  
  
  //! Epetra_CompObject destructor.  
  virtual ~Epetra_CompObject();
  //@}

  //@{ \name Set/Get counter method.
  //! Set the internal Epetra_Flops() pointer.
  void SetFlopCounter(const Epetra_Flops & FlopCounter) {FlopCounter_= (Epetra_Flops *) &FlopCounter; return;};
  //! Set the internal Epetra_Flops() pointer to the flop counter of another Epetra_CompObject.
  void SetFlopCounter(const Epetra_CompObject & CompObject) {FlopCounter_= (Epetra_Flops *) (CompObject.GetFlopCounter()); return;};
  //! Set the internal Epetra_Flops() pointer to 0 (no flops counted).
  void UnsetFlopCounter() {FlopCounter_= 0; return;};
  //! Get the pointer to the  Epetra_Flops() object associated with this object, returns 0 if none.
  Epetra_Flops * GetFlopCounter() const {return(FlopCounter_);};
  //@}

  //@{ \name Set flop count methods.
  //! Resets the number of floating point operations to zero for \e this multi-vector.
  void ResetFlops() const {if (FlopCounter_!=0) FlopCounter_->ResetFlops(); return;};

  //! Returns the number of floating point operations with \e this multi-vector.
  double Flops() const {if (FlopCounter_!=0) return(FlopCounter_->Flops()); else return(0.0);};
  //@}

  //@{ \name Update flop count methods.
  //! Increment Flop count for \e this object
  void UpdateFlops(int Flops) const {if (FlopCounter_!=0) FlopCounter_->UpdateFlops(Flops); return;};

  //! Increment Flop count for \e this object
  void UpdateFlops(long int Flops) const {if (FlopCounter_!=0) FlopCounter_->UpdateFlops(Flops); return;};

  //! Increment Flop count for \e this object
  void UpdateFlops(double Flops) const {if (FlopCounter_!=0) FlopCounter_->UpdateFlops(Flops); return;};

  //! Increment Flop count for \e this object
  void UpdateFlops(float Flops) const {if (FlopCounter_!=0) FlopCounter_->UpdateFlops(Flops); return;};
  //@}

 protected:


  Epetra_Flops * FlopCounter_;

};

#endif /* _EPETRA_COMPOBJECT_H_ */
