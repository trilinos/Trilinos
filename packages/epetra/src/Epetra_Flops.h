
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

#ifndef _EPETRA_FLOPS_H_
#define _EPETRA_FLOPS_H_

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

#endif /* _EPETRA_FLOPS_H_ */
