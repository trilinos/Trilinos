// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_ALGFORTESTINGONLY_HPP_
#define _ZOLTAN2_ALGFORTESTINGONLY_HPP_

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgForTestingOnly.hpp
//! \brief NOT a real algorithm; this algorithm is used
//!        to force hard-coded results for testing

////////////////////////////////////////////////////////////////////////

namespace Zoltan2 {

template <typename Adapter>
class AlgForTestingOnly : public Algorithm<Adapter>
{
private:  
  typedef typename Adapter::part_t part_t;
  const RCP<const Environment> env;
  const RCP<const Comm<int> > comm;
  const RCP<const typename Adapter::base_adapter_t> adapter;

public:
  
  AlgForTestingOnly(
     const RCP<const Environment> &env__,
     const RCP<const Comm<int> > &problemComm__,
     const RCP<const typename Adapter::base_adapter_t> &adapter__):
     env(env__), comm(problemComm__), adapter(adapter__)
  {
  }

  /*! \brief Set up validators specific to this algorithm
    */
  static void getValidParameters(ParameterList & pl)
  {
    RCP<Teuchos::EnhancedNumberValidator<int>> forTestingOnlyFlag_Validator =
      Teuchos::rcp( new Teuchos::EnhancedNumberValidator<int>(0, 1000, 1, 0) );
    pl.set("forTestingOnlyFlag", 0, "Used only for testing; look at "
      "Zoltan2_AlgForTestingOnly for interpretations",
      forTestingOnlyFlag_Validator);
  }

  void partition(const RCP<PartitioningSolution<Adapter> > &solution)
  {
    size_t nObj = adapter->getLocalNumIDs();
    ArrayRCP<part_t> partList(new part_t[nObj], 0, nObj, true);
    size_t nGlobalParts = solution->getTargetGlobalNumberOfParts();

    const Teuchos::ParameterEntry *pe = 
          env->getParameters().getEntryPtr("forTestingOnlyFlag");
    int forTestingOnlyFlag = pe->getValue<int>(&forTestingOnlyFlag);
    
    switch (forTestingOnlyFlag) {
    case 0:
      // rank 0 has all objects in part 0
      // all other ranks assign to {0, nGlobalParts-1, 0, nGlobalParts-1, ..}
      if (comm->getRank() == 0) {
        for (size_t i = 0; i < nObj; i++) partList[i] = 0;
      }
      else {
        for (size_t i = 0; i < nObj; i++) 
          if (i % 2) partList[i] = nGlobalParts - 1;
          else partList[i] = 0;
      }
      break;
    case 1:
      // rank 0 has all objects in part 0
      // all other ranks assign to {nGlobalParts-1, 0, nGlobalParts-1, 0, ..}
      if (comm->getRank() == 0) {
        for (size_t i = 0; i < nObj; i++) partList[i] = 0;
      }
      else {
        for (size_t i = 0; i < nObj; i++) 
          if (i % 2) partList[i] = 0;
          else partList[i] = nGlobalParts - 1;
      }
      break;
    default:
      throw std::runtime_error("invalid forTestingOnlyFlag value");
    }
      
    std::cout << comm->getRank() << " forTestingOnly " << forTestingOnlyFlag
              << " partList:  ";
    for (size_t i = 0; i < nObj; i++)
      std::cout << partList[i] << " ";
    std::cout << std::endl;

    solution->setParts(partList);
  }

};
}

#endif
