/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_COMPOBJECT_H
#define EPETRA_COMPOBJECT_H

//! Epetra_CompObject: Functionality and data that is common to all computational classes.

/*! The Epetra_CompObject is a base class for all Epetra computational objects.  It provides the basic
    mechanisms and interface specifications for floating point operations using Epetra_Flops objects.

*/
#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"
#include "Epetra_Flops.h"
//==========================================================================
class EPETRA_LIB_DLL_EXPORT Epetra_CompObject {

  public:

    //! @name Constructors/Destructor
  //@{ 
  //! Basic Epetra_CompObject constuctor.
  Epetra_CompObject();

  //! Epetra_CompObject copy constructor.
  
  Epetra_CompObject(const Epetra_CompObject& Source);
  
  
  //! Epetra_CompObject destructor.  
  virtual ~Epetra_CompObject();
  //@}

  //! @name Set/Get counter method
  //@{ 
  //! Set the internal Epetra_Flops() pointer.
  void SetFlopCounter(const Epetra_Flops & FlopCounter_in) {FlopCounter_= (Epetra_Flops *) &FlopCounter_in; return;}
  //! Set the internal Epetra_Flops() pointer to the flop counter of another Epetra_CompObject.
  void SetFlopCounter(const Epetra_CompObject & CompObject) {FlopCounter_= (Epetra_Flops *) (CompObject.GetFlopCounter()); return;}
  //! Set the internal Epetra_Flops() pointer to 0 (no flops counted).
  void UnsetFlopCounter() {FlopCounter_= 0; return;}
  //! Get the pointer to the  Epetra_Flops() object associated with this object, returns 0 if none.
  Epetra_Flops * GetFlopCounter() const {return(FlopCounter_);}
  //@}

  //! @name Set flop count methods
  //@{ 
  //! Resets the number of floating point operations to zero for \e this multi-vector.
  void ResetFlops() const {if (FlopCounter_!=0) FlopCounter_->ResetFlops(); return;}

  //! Returns the number of floating point operations with \e this multi-vector.
  double Flops() const {if (FlopCounter_!=0) return(FlopCounter_->Flops()); else return(0.0);}
  //@}

  //! @name Update flop count methods
  //@{ 
  //! Increment Flop count for \e this object
  void UpdateFlops(int Flops_in) const {if (FlopCounter_!=0) FlopCounter_->UpdateFlops(Flops_in); return;}

  //! Increment Flop count for \e this object
  void UpdateFlops(long int Flops_in) const {if (FlopCounter_!=0) FlopCounter_->UpdateFlops(Flops_in); return;}

  //! Increment Flop count for \e this object
  void UpdateFlops(double Flops_in) const {if (FlopCounter_!=0) FlopCounter_->UpdateFlops(Flops_in); return;}

  //! Increment Flop count for \e this object
  void UpdateFlops(float Flops_in) const {if (FlopCounter_!=0) FlopCounter_->UpdateFlops(Flops_in); return;}
  //@}

  Epetra_CompObject& operator=(const Epetra_CompObject& src)
    {
      FlopCounter_ = src.FlopCounter_;
      return(*this);
    }

 protected:


  Epetra_Flops * FlopCounter_;

};

#endif /* EPETRA_COMPOBJECT_H */
