// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
// @HEADER

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef TEUCHOS_COMPOBJECT_HPP
#define TEUCHOS_COMPOBJECT_HPP

/*! \file Teuchos_CompObject.hpp
    \brief Object for storing data and providing functionality that is common to all 
	computational classes.
*/

#include "Teuchos_Object.hpp"
#include "Teuchos_Flops.hpp"

/*! \class Teuchos::CompObject
    \brief Functionality and data that is common to all computational classes.

    The Teuchos::CompObject is a base class for all Teuchos computational objects.  It provides the basic
    mechanisms and interface specifications for floating point operations using Teuchos::Flops objects. 
*/

namespace Teuchos
{
class TEUCHOSNUMERICS_LIB_DLL_EXPORT CompObject
{

  public:

    //! @name Constructors/Destructor.
  //@{ 

  //! Default constructor
  CompObject();

  //! Copy Constructor
  CompObject(const CompObject &source);
    
  //! Destructor
  virtual ~CompObject();
  //@}

  //! @name Set/Get counter method.
  //@{ 
  //! Set the internal Teuchos::Flops() pointer.
  void setFlopCounter(const Flops &FlopCounter) {flopCounter_= (Flops *) &FlopCounter; return;}
  
  //! Set the internal Teuchos::Flops() pointer to the flop counter of another Teuchos::CompObject.
  void setFlopCounter(const CompObject &compObject) {flopCounter_= (Flops *) (compObject.getFlopCounter()); return;}
  
  //! Set the internal Teuchos::Flops() pointer to 0 (no flops counted).
  void unsetFlopCounter() {flopCounter_=0; return;}
  
  //! Get the pointer to the Teuchos::Flops() object associated with this object, returns 0 if none.
  Flops * getFlopCounter() const {return(flopCounter_);}
  //@}

  //! @name Set flop count methods.
  //@{ 
  //! Resets the number of floating point operations to zero for \e this multi-std::vector.
  void resetFlops() const {if (flopCounter_!=0) flopCounter_->resetFlops(); return;}

  //! Returns the number of floating point operations with \e this multi-std::vector.
  double getFlops() const {if (flopCounter_!=0) return(flopCounter_->flops()); else return(0.0);}
  //@}

  //! @name Update flop count methods.
  //@{ 
  //! Increment Flop count for \e this object
  void updateFlops(int addflops) const { if (flopCounter_!=0) flopCounter_->updateFlops(addflops); return;}

  //! Increment Flop count for \e this object
  void updateFlops(long int addflops) const { if (flopCounter_!=0) flopCounter_->updateFlops(addflops); return;}

  //! Increment Flop count for \e this object
  void updateFlops(double addflops) const { if (flopCounter_!=0) flopCounter_->updateFlops(addflops); return;}

  //! Increment Flop count for \e this object
  void updateFlops(float addflops) const {if (flopCounter_!=0) flopCounter_->updateFlops(addflops); return;}
  //@}

 protected:

  Flops *flopCounter_;

};

  // #include "Teuchos_CompObject.cpp"

} // namespace Teuchos

#endif // end of TEUCHOS_COMPOBJECT_HPP
