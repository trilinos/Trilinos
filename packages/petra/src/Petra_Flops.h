#ifndef _PETRA_FLOPS_H_
#define _PETRA_FLOPS_H_

//! Petra_Flops:  The Petra Floating Point Operations Class.
/*! The Petra_Flops class provides basic support and consistent interfaces
    for counting and reporting floating point operations performed in 
    the Petra computational classes.
  
*/

#include "Petra_Petra.h"
#include "Petra_Comm.h"

class Petra_Flops {
    
  public:
  //! Petra_Flops Constructor.
  /*! Creates a Petra_Flops instance. This instance can be queried for
      the number of floating point operations performed for the associated
      \e this object.
  */
  Petra_Flops(void);

  //! Petra_Flops Copy Constructor.
  /*! Makes an exact copy of an existing Petra_Flops instance.
  */
  Petra_Flops(const Petra_Flops& Flops);

  //! Returns the number of floating point operations with \e this object and resets the count.
  double Flops() const {double tmp = Flops_; Flops_ = 0.0; return(tmp);};

  //! Resets the number of floating point operations to zero for \e this multi-vector.
  void ResetFlops() {Flops_=0.0;};

  //! Petra_Flops Destructor.
  /*! Completely deletes a Petra_Flops object.  
  */
  virtual ~Petra_Flops(void);

 protected:
  mutable double Flops_;
  //! Increment Flop count for \e this object from an int
  void UpdateFlops(int Flops) const {Flops_ += (double) Flops;};
  //! Increment Flop count for \e this object from a long int
  void UpdateFlops(long int Flops) const {Flops_ += (double) Flops;};
  //! Increment Flop count for \e this object from a double
  void UpdateFlops(double Flops) const {Flops_ += Flops;};
  //! Increment Flop count for \e this object from a double
  void UpdateFlops(float Flops) const {Flops_ +=(double) Flops;};
  

 private:
  
};

#endif /* _PETRA_FLOPS_H_ */
