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

#ifndef EPETRA_FLOPS_H
#define EPETRA_FLOPS_H

#include "Epetra_ConfigDefs.h"

//! Epetra_Flops:  The Epetra Floating Point Operations Class.
/*! The Epetra_Flops class provides basic support and consistent interfaces
    for counting and reporting floating point operations performed in 
    the Epetra computational classes.  All classes based on the Epetra_CompObject
    can count flops by the user creating an Epetra_Flops object and calling the SetFlopCounter()
    method for an Epetra_CompObject.
  
*/

class EPETRA_LIB_DLL_EXPORT Epetra_Flops {
    
  public:
  //! Epetra_Flops Constructor.
  /*! Creates a Epetra_Flops instance. This instance can be queried for
      the number of floating point operations performed for the associated
      \e this object.
  */
  Epetra_Flops(void);

  //! Epetra_Flops Copy Constructor.
  /*! Makes an exact copy of an existing Epetra_Flops instance.
  */
  Epetra_Flops(const Epetra_Flops& Flops_in);

  //! Returns the number of floating point operations with \e this object and resets the count.
  double Flops() const {double tmp = Flops_; Flops_ = 0.0; return(tmp);};

  //! Resets the number of floating point operations to zero for \e this multi-vector.
  void ResetFlops() {Flops_=0.0;};

  //! Epetra_Flops Destructor.
  /*! Completely deletes a Epetra_Flops object.  
  */
  virtual ~Epetra_Flops(void);

  Epetra_Flops& operator=(const Epetra_Flops& src)
    {
      Flops_ = src.Flops_;
      return(*this);
    }

  friend class Epetra_CompObject;

 protected:
  mutable double Flops_;
  //! Increment Flop count for \e this object from an int
  void UpdateFlops(int Flops_in) const {Flops_ += (double) Flops_in;};
  //! Increment Flop count for \e this object from a long int
  void UpdateFlops(long int Flops_in) const {Flops_ += (double) Flops_in;};
  //! Increment Flop count for \e this object from a double
  void UpdateFlops(double Flops_in) const {Flops_ += Flops_in;};
  //! Increment Flop count for \e this object from a float
  void UpdateFlops(float Flops_in) const {Flops_ +=(double) Flops_in;};
  

 private:
  
};

#endif /* EPETRA_FLOPS_H */
