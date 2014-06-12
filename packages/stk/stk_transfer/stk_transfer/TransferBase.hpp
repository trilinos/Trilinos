/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef  STK_TRANSFERBASE_HPP
#define  STK_TRANSFERBASE_HPP


namespace stk {
namespace transfer {

class TransferBase {
public :
  TransferBase(){};
  virtual ~TransferBase(){};
  void initialize() {
    coarse_search();
    communication();
    local_search();
  }
  virtual void coarse_search() = 0;
  virtual void communication() = 0;
  virtual void local_search()  = 0;
  virtual void apply()         = 0;
};
}
}
#endif

