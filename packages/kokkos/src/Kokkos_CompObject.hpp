//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef KOKKOS_COMPOBJECT_H
#define KOKKOS_COMPOBJECT_H

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_Flops.hpp"

namespace Kokkos {

//! Kokkos::CompObject: Functionality and data that is common to all computational classes.

/*! The Kokkos::CompObject is a base class for all Kokkos computational objects.  It provides the basic
    mechanisms and interface specifications for floating point operations using Kokkos::Flops objects.

*/
class CompObject {

  public:

  //@{ \name Constructors/Destructor.
  //! Basic CompObject constuctor.
  CompObject() {
    flopCounter_ = 0;
  };

  //! CompObject copy constructor.
  
  CompObject(const CompObject& source) 
    : flopCounter_(source.flopCounter_){};
  
  
  //! CompObject destructor.  
  virtual ~CompObject(){flopCounter_=0;};
  //@}

  //@{ \name Set/Get counter method.

  //! Set the internal Flops() pointer.
  void setFlopCounter(const Flops & flopCounter) {flopCounter_= (Flops *) &flopCounter; return;};

  //! Set the internal Flops() pointer to the flop counter of another CompObject.
  void setFlopCounter(const CompObject & CompObject) {flopCounter_= (Flops *) (CompObject.getFlopCounter()); return;};

  //! Set the internal Flops() pointer to 0 (no flops counted).
  void unsetFlopCounter() {flopCounter_= 0; return;};

  //! Get the pointer to the  Flops() object associated with this object, returns 0 if none.
  Flops * getFlopCounter() const {return(flopCounter_);};

  //@}

  //@{ \name Set flop count methods.

  //! Resets the number of floating point operations to zero for \e this multi-vector.
  void resetFlops() const {if (flopCounter_!=0) flopCounter_->resetFlops(); return;};

  //! Returns the number of floating point operations with \e this multi-vector.
  double getFlops() const {if (flopCounter_!=0) return(flopCounter_->getFlops()); else return(0.0);};
  //@}

  //@{ \name Update flop count methods.
  //! Increment flop count for \e this object
  void updateFlops(int flops) const {if (flopCounter_!=0) flopCounter_->updateFlops(flops); return;};

  //! Increment flop count for \e this object
  void updateFlops(long int flops) const {if (flopCounter_!=0) flopCounter_->updateFlops(flops); return;};

  //! Increment flop count for \e this object
  void updateFlops(double flops) const {if (flopCounter_!=0) flopCounter_->updateFlops(flops); return;};

  //! Increment flop count for \e this object
  void updateFlops(float flops) const {if (flopCounter_!=0) flopCounter_->updateFlops(flops); return;};
  //@}

 protected:


  Flops * flopCounter_;

};

} // namespace Kokkos

#endif /* KOKKOS_COMPOBJECT_H */
