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
#ifndef MUELU_CLASSICALPFACTORY_DECL_HPP
#define MUELU_CLASSICALPFACTORY_DECL_HPP

#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_ClassicalPFactory_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"
#include "MueLu_CoarseMapFactory_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_LWGraph.hpp"
#include "MueLu_LWGraph_kokkos.hpp"

#include "MueLu_Level_fwd.hpp"

namespace MueLu {

  template <class Scalar = DefaultScalar,
            class LocalOrdinal = DefaultLocalOrdinal,
            class GlobalOrdinal = DefaultGlobalOrdinal,
            class Node = DefaultNode>
  class ClassicalPFactory : public PFactory {
#undef MUELU_CLASSICALPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    // Defining types that require the short names included above
    using local_graph_type = typename LWGraph_kokkos::local_graph_type;
    using colors_view_type = Kokkos::View<typename local_graph_type::entries_type::data_type,
                                          typename local_graph_type::device_type::memory_space>;


    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    ClassicalPFactory() { }

    //! Destructor.
    virtual ~ClassicalPFactory() { }
    //@}

    RCP<const ParameterList> GetValidParameterList() const;

    //! Input
    //@{

    void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{

    void Build (Level& fineLevel, Level& coarseLevel) const;
    void BuildP(Level& fineLevel, Level& coarseLevel) const;
    
  private:
    typedef enum {UNASSIGNED=0, F_PT, C_PT, DIRICHLET_PT} point_type;

    // Utility algorithms
    void GenerateCoarseMap(const Map & fineMap, const Teuchos::Array<point_type> & fc_splitting, RCP<const Map> & coarseMap, Array<LO> & cpoint2pcol, Teuchos::Array<LO> & pcol2cpoint) const;
    void DoGraphColoring(const GraphBase & graph, colors_view_type & myColors, LO & numColors) const;    
    void GenerateStrengthFlags(const Matrix & A,const GraphBase & graph, Teuchos::Array<size_t> & eis_rowptr, Teuchos::Array<bool> & edgeIsStrong) const;

    // Coarsening algorithms
    void Coarsen_Ext_Plus_I(const Matrix & A,const GraphBase & graph, RCP<const Map> & coarseColMap, RCP<const Map> & coarseDomainMap, LO num_c_points, LO num_f_points, Teuchos::Array<point_type> & myPointType, const Teuchos::Array<LO> & cpoint2pcol, const Teuchos::Array<LO> & pcol2cpoint, Teuchos::Array<size_t> & eis_rowptr, Teuchos::Array<bool> & edgeIsStrong, RCP<Matrix> & P) const;

    //@}

  }; //class ClassicalPFactory

} //namespace MueLu

#define MUELU_CLASSICALPFACTORY_SHORT
#endif // MUELU_CLASSICALPFACTORY_DECL_HPP
