
//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#ifndef EPETRA_FLOPS_H
#define EPETRA_FLOPS_H

//! Epetra_Flops:  The Epetra Floating Point Operations Class.
/*! The Epetra_Flops class provides basic support and consistent interfaces
    for counting and reporting floating point operations performed in 
    the Epetra computational classes.  All classes based on the Epetra_CompObject
    can count flops by the user creating an Epetra_Flops object and calling the SetFlopCounter()
    method for an Epetra_CompObject.
  
*/

class Epetra_Flops {
    
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
  Epetra_Flops(const Epetra_Flops& Flops);

  //! Returns the number of floating point operations with \e this object and resets the count.
  double Flops() const {double tmp = Flops_; Flops_ = 0.0; return(tmp);};

  //! Resets the number of floating point operations to zero for \e this multi-vector.
  void ResetFlops() {Flops_=0.0;};

  //! Epetra_Flops Destructor.
  /*! Completely deletes a Epetra_Flops object.  
  */
  virtual ~Epetra_Flops(void);

  friend class Epetra_CompObject;

 protected:
  mutable double Flops_;
  //! Increment Flop count for \e this object from an int
  void UpdateFlops(int Flops) const {Flops_ += (double) Flops;};
  //! Increment Flop count for \e this object from a long int
  void UpdateFlops(long int Flops) const {Flops_ += (double) Flops;};
  //! Increment Flop count for \e this object from a double
  void UpdateFlops(double Flops) const {Flops_ += Flops;};
  //! Increment Flop count for \e this object from a float
  void UpdateFlops(float Flops) const {Flops_ +=(double) Flops;};
  

 private:
  
};

#endif /* EPETRA_FLOPS_H */
