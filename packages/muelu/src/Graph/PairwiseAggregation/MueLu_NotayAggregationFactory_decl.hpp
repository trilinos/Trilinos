// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_NOTAYAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_NOTAYAGGREGATIONFACTORY_DECL_HPP_


#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_NotayAggregationFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

template<class Scalar = DefaultScalar,
         class LocalOrdinal = DefaultLocalOrdinal,
         class GlobalOrdinal = DefaultGlobalOrdinal,
         class Node = DefaultNode>
class NotayAggregationFactory : public SingleLevelFactoryBase {
#undef MUELU_NOTAYAGGREGATIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  NotayAggregationFactory() { };

  //! Destructor.
  virtual ~NotayAggregationFactory() { }

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! @name Set/get methods.
  //@{

  // Options shared by all aggregation algorithms

  //! Input
  //@{

  void DeclareInput(Level &currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*! @brief Build aggregates. */
  void Build(Level &currentLevel) const;


  void Build_InitialAggregation(const Teuchos::ParameterList& params,
				const RCP<const Matrix>& A,
				Aggregates& aggregates,
				std::vector<unsigned>& aggStat,
				LO& numNonAggregatedNodes) const;

  
  //@}

private:
}; // class NotayAggregationFactory

}

#define MUELU_NOTAYAGGREGATIONFACTORY_SHORT
#endif /* MUELU_NOTAYAGGREGATIONFACTORY_DECL_HPP_ */
