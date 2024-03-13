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
#ifndef MUELU_COMBINEPFACTORY_DECL_HPP
#define MUELU_COMBINEPFACTORY_DECL_HPP

#include <Teuchos_SerialDenseVector.hpp>

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_CombinePFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"

namespace MueLuTests {
// Forward declaration of friend tester class used to UnitTest CombinePFactory
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class CombinePFactoryTester;
}  // namespace MueLuTests

namespace MueLu {

/*!
  @class CombinePFactory
  @ingroup MueLuTransferClasses
  @brief Prolongator factory that replicates 'Psubblock' matrix to create new prolongator suitable for PDE systems


  Takes a set of previously generated prolongators each for a subsystem multiphysics PDE and effectively makes a block diagonal prolongator
  by combining the "combo: npdes" times so that it can be used with a PDE system. A normal use case
  would be to run an existing MueLu multigrid algorithm on a set of scalar PDEs to generate a hierarchy.  Then use
  something like

  MueLu::HierarchyUtils<SC,LO,GO,NO>::CopyBetweenHierarchies(*scalarHierarchy0,*systemHierarchy, "P",  "Psubblock0", "RCP<Matrix>");
                                                                             .                                   .
                                                                             .                                   .
                                                                             .                                   .
  MueLu::HierarchyUtils<SC,LO,GO,NO>::CopyBetweenHierarchies(*scalarHierarchy7,*systemHierarchy, "P",  "Psubblock6", "RCP<Matrix>");

  to copy them to a new hierarchy. This new hierarchy would then have <Parameter name="multigrid algorithm" type="string"  value="combo"/>
  in its parameter list to then invoke ComboPFactory. Currently, this is used in src/Operators/MueLu_Multiphysics_def.hpp with an example
  in test/multiphysics.

  The end result is a block diagonal matrix corresponding to
                   [Psubblock0                                                                  ]
                   [           Psubblock1                                                       ]
                   [                      Psubblock2                                            ]
          P   =    [                                 Psubblock3                                 ]
                   [                                            Psubblock4                      ]
                   [                                                       Psubblock5           ]
                   [                                                                  Psubblock6]

  ## Input/output of CombinePFactory ##

  ### User parameters of SemiCoarsenPFactory ###
  | Parameter   | type    | default   | master.xml | validated | requested | description                                                                      |
  | ------------|---------|-----------|:----------:|:---------:|:---------:|----------------------------------------------------------------------------------|
  |combo: npdes | int     |  1        |     *      |     *     |           | Specifies the number of Psubblock# matrices to expect                            |
  | Psubblock0   | Matrix |           |            |           |           | Matrix defining P(0,0) block                                                     |
  | Psubblock1   | Matrix |           |            |           |           | Matrix defining P(1,1) block                                                     |
  |      .       |    .   |           |            |           |           |            .                                                                     |
  |      .       |    .   |           |            |           |           |            .                                                                     |
  |      .       |    .   |           |            |           |           |            .                                                                     |
  | Psubblockk   | Matrix |           |            |           |           | Matrix defining P(k,k) block                                                     |



  The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see CombineCoarsenPFactory::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see CombineCoarsenPFactory::DeclareInput).

  ### Variables provided by CombinePFactory ###

  After CombinePFactory::Build the following data is available (if requested)

  | Parameter         | generated by             | description
  |-------------------|--------------------------|------------------------------------------------------------------------------------------------------------------|
  | P                 | CombinePFactory          | Prolongator                                                                                                      |

*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class CombinePFactory : public PFactory {
#undef MUELU_COMBINEPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  friend class MueLuTests::CombinePFactoryTester<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  CombinePFactory() {}

  //! Destructor.
  virtual ~CombinePFactory() {}
  //@}

  RCP<const ParameterList> GetValidParameterList() const;

  //! Input
  //@{

  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  void Build(Level& fineLevel, Level& coarseLevel) const;
  void BuildP(Level& fineLevel, Level& coarseLevel) const;

  //@}

 private:
  int numPDEs_;

};  // class CombinePFactory

}  // namespace MueLu

#define MUELU_COMBINEPFACTORY_SHORT
#endif  // MUELU_COMBINEPFACTORY_DECL_HPP
