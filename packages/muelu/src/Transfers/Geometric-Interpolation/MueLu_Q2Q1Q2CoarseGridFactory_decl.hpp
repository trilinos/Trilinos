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

#ifndef MUELU_Q2Q1Q2COARSEGRIDFACTORY_DECL_HPP
#define MUELU_Q2Q1Q2COARSEGRIDFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
//#include "MueLu_Q2Q1Q2CoarseGridFactory_fwd.hpp"

#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_Level_fwd.hpp"
#include "MueLu_PFactory.hpp"

namespace MueLu {

/*!
  @class Q2Q1Q2CoarseGridFactory class.
  @brief Factory for GMG Q2-Q1-Q2 interpolation
  @ingroup MueLuTransferClasses
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class Q2Q1Q2CoarseGridFactory : public TwoLevelFactoryBase {
#undef MUELU_Q2Q1Q2COARSEGRIDFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  Q2Q1Q2CoarseGridFactory();

  //! Destructor.
  virtual ~Q2Q1Q2CoarseGridFactory();

  //@}

  //! @name Input
  //@{

  /*! @brief Specifies the data that this class needs, and the factories that generate that data.

      If the Build method of this class requires some data, but the
      generating factory is not specified in DeclareInput, then this
      class will fall back to the settings in FactoryManager.
  */
  void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  //! Build an object with this factory.
  void Build(Level &fineLevel, Level &coarseLevel) const;

  void BuildCoarseGrid(Level &fineLevel, Level &coarseLevel) const;

  /* For our geometric multigrid needs, we will explicitly build the
     coarse grid here and store it as level data. There are two
     things we care about here: coordinates of the dofs and
     element-dof lists.

     We will assume only that the elements are numbered
     lexicographically, left-to-right, bottom-to-top. The order of
     the degrees of freedom will only be assumed on coarse grids. It
     can be anything on the finest grid (fineLevel.GetLevelID()==1).

     We further assume that the x- and y-components of velocity are
     zippered together: UX1 UY1 UX2 UY2 ETC...
  */

  //@}

 private:
  // No parameters need to be stored... To save "P", use
  //    coarseLevel.Set("P", finalP, this);

};  // class Q2Q1Q2CoarseGridFactory

}  // namespace MueLu

#define MUELU_Q2Q1Q2COARSEGRIDFACTORY_SHORT
#endif  // MUELU_Q2Q1Q2COARSEGRIDFACTORY_DECL_HPP
