// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef TEUCHOS_FLOPS_HPP
#define TEUCHOS_FLOPS_HPP

namespace Teuchos
{

//! Teuchos_Flops:  The Teuchos Floating Point Operations Class.
/*! The Teuchos_Flops class provides basic support and consistent interfaces
    for counting and reporting floating point operations performed in 
    the Teuchos computational classes.  All classes based on the Teuchos_CompObject
    can count flops by the user creating an Teuchos_Flops object and calling the SetFlopCounter()
    method for an Teuchos_CompObject. It should be noted that currently, Teuchos_Flops is an
    almost exact duplicate of Epetra_Flops.
*/

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
  double flops() const {double tmp = flops_; flops_= 0.0; return(tmp); };

  //! Resets the number of floating point operations to zero for \e this multi-vector.
  void resetFlops() {flops_ = 0.0;};

  //! Flops Destructor.
  /*! Completely deletes a Flops object.
  */
  virtual ~Flops();

  friend class CompObject;

 protected:
  mutable double flops_;
  //! Increment Flop count for \e this object from an int
  void updateFlops(int addflops) const {flops_ += (double) addflops; };
  //! Increment Flop count for \e this object from a long int
  void updateFlops(long int addflops) const {flops_ += (double) addflops; };
  //! Increment Flop count for \e this object from a double
  void updateFlops(double addflops) const {flops_ += (double) addflops; };
  //! Increment Flop count for \e this object from a float
  void updateFlops(float addflops) const {flops_ += (double) addflops; };

 private:
  
};

  // #include "Teuchos_Flops.cpp"

} // namespace Teuchos

#endif // end of TEUCHOS_FLOPS_HPP
