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
#ifndef MUELU_AMALGAMATIONFACTORY_DECL_HPP
#define MUELU_AMALGAMATIONFACTORY_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

/*!
  @class AmalgamationFactory
  @brief AmalgamationFactory for subblocks of strided map based amalgamation data

  Class generates unamalgamation information using matrix A with strided maps.
  It stores the output information within an AmalgamationInfo object as "UnAmalgamationInfo".
  This object contains

  \li \c nodegid2dofgids_ a map of all node ids of which the current proc has corresponding DOF gids (used by \c TentativePFactory).
  \li \c gNodeIds vector of all node ids on the current proc (may be less than nodegid2dofgids_.size()). These nodes are stored on the current proc.

*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class AmalgamationFactory : public SingleLevelFactoryBase {
#undef MUELU_AMALGAMATIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  AmalgamationFactory() = default;

  //! Destructor
  virtual ~AmalgamationFactory() = default;

  RCP<const ParameterList> GetValidParameterList() const override;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const override;

  //@}

  void Build(Level& currentLevel) const override;

  /*! @brief Translate global (row/column) id to global amalgamation block id
   *
   * @note Assume that the node map has the same \c indexBase as the dof map
   *
   * @param gid (GlobalOrdinal): input global id (row gid or column gid)
   * @param blockSize (LocalOrdinal): block size (needed for constant block size)
   * @param offset (GlobalOrdinal): global offset for dofs (stored in strided map, default = 0)
   * @param indexBase (GlobalOrdinal): indexBase for DOF map (and node map, default = 0)
   */
  static const GlobalOrdinal DOFGid2NodeId(GlobalOrdinal gid, LocalOrdinal blockSize, const GlobalOrdinal offset /*= 0*/,
                                           const GlobalOrdinal indexBase /* = 0*/);

  /*! @brief Method to create merged  map for systems of PDEs.
   *
   * @param sourceMap (const Map&): source map with dofs which shall be amalgamated to a node map
   * @param A (const Matrix&): operator A (matrix) with striding information (if available)
   * @param amalgamatedMap (const Map&): amalgamated node based map
   * @param translation (Array<LO>&): array storing local node ids given local dof ids (needed in CoalesceDropFactory)
   */
  static void AmalgamateMap(const Map& sourceMap, const Matrix& A, RCP<const Map>& amalgamatedMap, Array<LO>& translation);

};  // class AmalgamationFactory

}  // namespace MueLu

#define MUELU_AMALGAMATIONFACTORY_SHORT
#endif  // MUELU_AMALGAMATIONFACTORY_DECL_HPP
