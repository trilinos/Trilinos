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
template<Z2CLASS_TEMPLATE>
class AlgScotch {
public:
  AlgScotch(RCP<Model<Z2PARAM_TEMPLATE> >, //TODO Should be GraphModel, but don't know how to do the Adapter template
            RCP<PartitioningSolution<Z2PARAM_TEMPLATE> >,
            RCP<Teuchos::ParameterList>);

  // Destructor
  ~AlgScotch() {};

protected:
private:

};


////////////////////////////////////////////////////////////////////////
//! Problem class constructor:  Tpetra matrix input must be converted
//! to XpetraMatrixAdapter.
template <Z2FN_TEMPLATE>
AlgScotch<Z2PARAM_TEMPLATE>::AlgScotch(
  RCP<Model<Z2PARAM_TEMPLATE> > model,  //TODO Should be GraphModel, but don't know how to do the adapter template
  RCP<PartitioningSolution<Z2PARAM_TEMPLATE> > solution,
  RCP<Teuchos::ParameterList> pl
) 
{
  HELLO;
}

}
#endif
