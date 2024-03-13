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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_SEGREGATEDAFACTORY_DECL_HPP
#define MUELU_SEGREGATEDAFACTORY_DECL_HPP

#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLu {

/*!
    @class  SegregatedAFactory class.
    @brief  Factory for building a new "segregated" A operator. Here, "segregated" means that the user
            provides some map(s) (containing a subset of GIDs of the input matrix A) and the factory
            drops entries depending on the dropping scheme.

            ## Idea ##

            The idea is to use the output matrix A as input for the aggregation factory to have control over
            the aggregates and make sure that aggregates do not cross certain areas.

            ## Remarks ##

            This factory supports multiple dropping schemes based on different inputs. They are:

            - blockmap: Based on the user provided "blockmap", the off-diagonal entries (a,b) and (b,a) in A are dropped.
            "a" denotes a GID entry in the provided map and "b" denotes a GID that is not contained in the provided map.
            In this use case the Factory expects a dropMap1 (==blockmap).
            The blockmap scheme also doesn't support the "Call ReduceAll on dropMap1/2" options.

            - map-pair: Based on a "map-pair", the user provides two maps "dropMap1" and "dropMap2",
            which specify global row/column pairs in the operator A to be dropped.
            The Factory drops any possible combination of the dropMaps 1 and 2. To ensure that entry A(a,b) is
            dropped, as well as entry A(b,a), there is an option to create redundant dropMaps on all Procs.
            This ensures that entries aren't overlooked due to the local rowmaps of the operator A.

            Note: we have to drop the entries (i.e. not just set them to zero) as the CoalesceDropFactory
                 does not distinguish between matrix entries which are zero and nonzero.

    ## Input/output of this factory ##

    ### User parameters of SegregatedAFactory ###
    Parameter | type | default | master.xml | validated | requested | description
    ----------|------|---------|:----------:|:---------:|:---------:|------------
    A             | Factory | null  |   | * | * | Generating factory of the matrix A
    droppingScheme| string  | vague |   | * | * | Strategy to drop entries from matrix A based on the input of some map(s) [blockmap, map-pair]
    dropMap1      | Factory | null  |   | * | * | Generating factory for dropMap1
    dropMap2      | Factory | null  |   | * | * | Generating factory for dropMap2
    Call ReduceAll on dropMap1 | bool |   | * |  | Boolean for calling reduceAll on dropMap1
    Call ReduceAll on dropMap2 | bool |   | * |  | Boolean for calling reduceAll on dropMap2

    The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
    The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see @c GetValidParameters() ).<br>
    The * in the @c requested column states that the data is requested as input with all dependencies (see @c DeclareInput() ).

    ### Variables provided by this factory ###

    After SegregatedAFactory::Build the following data is available (if requested)

    Parameter | generated by | description
    ----------|--------------|------------
    A         | SegregatedAFactory | Provides a filtered Matrix, where all possible combinations of the entries of the dropMap(s) have been removed from the input matrix A
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class SegregatedAFactory : public SingleLevelFactoryBase {
#undef MUELU_SEGREGATEDAFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! Input
  //@{

  RCP<const ParameterList> GetValidParameterList() const override;

  void DeclareInput(Level& currentLevel) const override;

  //! @name Build methods.

  /*!
    @brief Build method.

    Builds filtered matrix and returns it in <tt>currentLevel</tt>.
    */
  void Build(Level& currentLevel) const override;

 private:
  void BuildBasedOnBlockmap(Level& currentLevel) const;

  void BuildBasedOnMapPair(Level& currentLevel) const;

  //  RCP<const Map> CreateRedundantMaps(Teuchos::RCP<const Map> localDropMap, Teuchos::RCP<const Matrix> Ain) const;

};  // class SegregatedAFactory

}  // namespace MueLu

#define MUELU_SEGREGATEDAFACTORY_SHORT
#endif  // MUELU_SEGREGATEDAFACTORY_DECL_HPP
