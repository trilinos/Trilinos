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

#ifndef TEUCHOS_FLOPS_HPP
#define TEUCHOS_FLOPS_HPP

/*! \file Teuchos_Flops.hpp 
    \brief Object for providing basic support and consistent interfaces for 
	counting/reporting floating-point operations performed in Teuchos computational
	classes.
*/

/*! \class Teuchos::Flops
    \brief The Teuchos Floating Point Operations Class.

    The Teuchos_Flops class provides basic support and consistent interfaces
    for counting and reporting floating point operations performed in 
    the Teuchos computational classes.  All classes based on the Teuchos::CompObject
    can count flops by the user creating an Teuchos::Flops object and calling the SetFlopCounter()
    method for an Teuchos_CompObject. 
*/

namespace Teuchos
{
class Flops
{    
  public:

    //! @name Constructor/Destructor.
  //@{ 

  //! Default Constructor.
  /*! Creates a Flops instance. This instance can be queried for
      the number of floating point operations performed for the associated
      \e this object.
  */
  Flops();

  //! Copy Constructor.
  /*! Makes an exact copy of an existing Flops instance.
  */
  Flops(const Flops &flops);

  //! Destructor.
  /*! Completely deletes a Flops object.
  */
  virtual ~Flops();

  //@}

  //! @name Accessor methods.
  //@{ 

  //! Returns the number of floating point operations with \e this object and resets the count.
  double flops() const { return flops_; }

  //@}

  //! @name Reset methods.
  //@{ 

  //! Resets the number of floating point operations to zero for \e this multi-std::vector.
  void resetFlops() {flops_ = 0.0;}

  //@}

  friend class CompObject;

 protected:

  mutable double flops_;

  //! @name Updating methods.
  //@{ 
  //! Increment Flop count for \e this object from an int
  void updateFlops(int addflops) const {flops_ += (double) addflops; }

  //! Increment Flop count for \e this object from a long int
  void updateFlops(long int addflops) const {flops_ += (double) addflops; }

  //! Increment Flop count for \e this object from a double
  void updateFlops(double addflops) const {flops_ += (double) addflops; }

  //! Increment Flop count for \e this object from a float
  void updateFlops(float addflops) const {flops_ += (double) addflops; }

  //@}

 private:
  
};

  // #include "Teuchos_Flops.cpp"

} // namespace Teuchos

#endif // end of TEUCHOS_FLOPS_HPP
