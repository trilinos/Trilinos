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


#ifndef KOKKOS_FLOPS_H
#define KOKKOS_FLOPS_H

#include "Kokkos_ConfigDefs.hpp"

namespace Kokkos {

//! Kokkos::Flops:  The Kokkos Floating Point Operations Class.
/*! The Kokkos::Flops class provides basic support and consistent interfaces
  for counting and reporting floating point operations performed in 
  the Kokkos computational classes.  All classes based on the Kokkos::CompObject
  can count flops by the user creating an Kokkos::Flops object and calling the SetFlopCounter()
  method for an Kokkos::CompObject.
  
*/

  class Flops {
    
  public:
    //@{ \name Constructors/Destructor.
    //! Flops Constructor.
    /*! Creates a Flops instance. This instance can be queried for
      the number of floating point operations performed for the associated
      \e this object.
    */
    Flops(void)
      : flops_(0.0){};

    //! Flops Copy Constructor.
    /*! Makes a copy of an existing Flops instance.
      \warning The counter of the Flops copy is reset to zero.
    */
    Flops(const Flops& flops)
      : flops_(0.0){};

    //! Flops Destructor.
    virtual ~Flops(void){};
    //@}


    //@{ \name Attribute access/reset methods.
    //! Returns the number of floating point operations with \e this object and resets the count.
    double getFlops() const {double tmp = flops_; flops_ = 0.0; return(tmp);};

    //! Resets the number of floating point operations to zero for \e this multi-vector.
    void resetFlops() {flops_=0.0;};
    //@}

    //@{ \name Friend class.
    //! CompObject needs to be able to set the flopCounter_ attribute of Flops
    friend class CompObject;
    //@}

  protected:

    mutable double flops_;

    //! Increment flop count for \e this object from an int
    void updateFlops(int flops) const {flops_ += (double) flops;};

    //! Increment flop count for \e this object from a long int
    void updateFlops(long int flops) const {flops_ += (double) flops;};

    //! Increment flop count for \e this object from a double
    void updateFlops(double flops) const {flops_ += flops;};

    //! Increment flop count for \e this object from a float
    void updateFlops(float flops) const {flops_ +=(double) flops;};
  
  };

} // namespace Kokkos

#endif /* KOKKOS_FLOPS_H */
