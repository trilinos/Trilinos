/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef STK_Transfer_hpp
#define STK_Transfer_hpp

#include <Intrepid_FieldContainer.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_mesh/base/Comm.hpp>


namespace stk {
namespace transfer {

class Transfer {
public:
  typedef Intrepid::FieldContainer<double> MDArray;

  struct Hints {
    Hints(stk::diag::Timer *timer=0) : 
      Tolerance      (0.1), 
      ExpansionFactor(1.5),
      Timer          (timer ? stk::diag::Timer("Transfer",*timer) : 
                              stk::diag::createRootTimer("Transfer",stk::diag::TimerSet(0))){}
    double Tolerance;
    double ExpansionFactor;
    stk::diag::Timer Timer;
  };

  // Constructor is trivial but could be expanded in the
  // future to allocate and initialize input data.
  Transfer(const std::string &name,
           const Hints& hints=Hints());

  // Given values defined at points in arrays FromValues and
  // FromPoints which have to be of the same length, and given
  // points to transfer in in ToPoints, find the N+1 nearest
  // points to each ToPoints from the FromPoints and do linear
  // interpolation of the FromValues to the ToValues.  N is the
  // spatical dimension and the second dimension of the FromPoints
  // and ToPoints matrix.

  // Is Intrepid::FieldContainer the standard for large amounts
  // of ordered data? Seems that STK might have it's own class
  // for this. Should PointToPoint depend on Intrepid this way?
  void PointToPoint(      MDArray &ToValues,
                    const MDArray &ToPoints,
                    const MDArray &FromValues,
                    const MDArray &FromPoints,
                    const stk::ParallelMachine  comm);
private :
  const std::string Name;
  const double Tolerance;
  const double ExpansionFactor;
  stk::diag::Timer RootTimer;
  stk::diag::Timer GhostingTimer;

  template <unsigned DIM>
  void PToP(MDArray &ToValues,
            const MDArray &ToPoints,
            const MDArray &FromValues,
            const MDArray &FromPoints,
            const stk::ParallelMachine  comm);

  template <unsigned DIM>
  void PPToP(MDArray &ToValues,
            const MDArray &ToPoints,
            const MDArray &FromValues,
            const MDArray &FromPoints,
            const stk::ParallelMachine  comm);
};

}}
#endif
