#ifndef _ZOLTAN2_ALGBLOCK_HPP_
#define _ZOLTAN2_ALGBLOCK_HPP_

#include <Zoltan2_Standards.hpp>

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgBlock.hpp
//! \brief Parallel block partitioning.

using std::cout;
using std::endl;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
//! \class AlgBlock
//! \brief Problem base class from which other classes (PartitioningProblem, 
//!        ColoringProblem, OrderingProblem, MatchingProblem, etc.) derive.
template<Z2CLASS_TEMPLATE>
class AlgBlock {
public:
  AlgBlock(RCP<Model<Z2PARAM_TEMPLATE> >, //TODO Should be IdModel, but don't know how to do the Adapter template
            RCP<PartitioningSolution<Z2PARAM_TEMPLATE> >,
            RCP<Teuchos::ParameterList>);

  // Destructor
  ~AlgBlock() {};

protected:
private:

};

