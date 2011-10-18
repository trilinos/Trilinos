#ifndef _ZOLTAN2_ALGSCOTCH_HPP_
#define _ZOLTAN2_ALGSCOTCH_HPP_

#include <Zoltan2_Standards.hpp>

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_Scotch.hpp
//! \brief Parallel graph partitioning using Scotch.

using std::cout;
using std::endl;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
//! \class AlgScotch
//! \brief Problem base class from which other classes (PartitioningProblem, 
//!        ColoringProblem, OrderingProblem, MatchingProblem, etc.) derive.
template<typename Adapter>
class AlgScotch {
public:
  typedef Adapter::user_t user_t;

  AlgScotch(const RCP<GraphModel<Adapter<user_t> > &, 
            const RCP<PartitioningSolution<user_t> > &,
            const RCP<Teuchos::ParameterList> &);

  // Destructor
  ~AlgScotch() {};

protected:
private:

};


////////////////////////////////////////////////////////////////////////
template <typename Adapter>
AlgScotch<Adapter>::AlgScotch(
  const RCP<GraphModel<Adapter> > &model, 
  const RCP<PartitioningSolution<user_t> > &solution,
  const RCP<Teuchos::ParameterList> &pl
) 
{
  HELLO;
}

}
#endif
