// 27-May-2002 General cleanup. Changed method names to fit namingConvention (already done). Moved some code to Tpetra_Flops.cpp

#ifndef _TPETRA_FLOPS_H_
#define _TPETRA_FLOPS_H_

//! Tpetra_Flops:  The Tpetra Floating Point Operations Class.
/*! The Tpetra_Flops class provides basic support and consistent interfaces
    for counting and reporting floating point operations performed in 
    the Tpetra computational classes.  All classes based on the Tpetra_CompObject
    can count flops by the user creating an Tpetra_Flops object and calling the SetFlopCounter()
    method for an Tpetra_CompObject. It should be noted that currently, Tpetra_Flops is an
    almost exact duplicate of Epetra_Flops.
*/

namespace Tpetra
{

class Flops
{
    
  public:
  //! Flops Constructor.
  /*! Creates a Flops instance. This instance can be queried for
      the number of floating point operations performed for the associated
      \e this object.
  */
  Flops(void);

  //! Flops Copy Constructor.
  /*! Makes an exact copy of an existing Flops instance.
  */
  Flops(const Flops &flops);

  //! Returns the number of floating point operations with \e this object and resets the count.
  double flops() const;

  //! Resets the number of floating point operations to zero for \e this multi-vector.
  void resetFlops();

  //! Flops Destructor.
  /*! Completely deletes a Flops object.  Tpetra::
  */
  virtual ~Flops();

  friend class CompObject;

 protected:
  mutable double flopCounter_;
  //! Increment Flop count for \e this object from an int
  void updateFlops(int flops) const;
  //! Increment Flop count for \e this object from a long int
  void updateFlops(long int flops) const;
  //! Increment Flop count for \e this object from a double
  void updateFlops(double flops) const;
  //! Increment Flop count for \e this object from a float
  void updateFlops(float flops) const;

 private:
  
};

#include "Tpetra_Flops.cpp"

} // namespace Tpetra

#endif // end of _TPETRA_FLOPS_H_
