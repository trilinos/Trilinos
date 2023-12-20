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
#ifndef MUELU_SMOOVECCOALESCEDROPFACTORY_DECL_HPP
#define MUELU_SMOOVECCOALESCEDROPFACTORY_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_SmooVecCoalesceDropFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Graph_fwd.hpp"
#include "MueLu_LWGraph_fwd.hpp"
#include "MueLu_PreDropFunctionBaseClass_fwd.hpp"

namespace MueLu {

/*!
  @class SmooVecCoalesceDropFactory
  @brief Factory for creating a graph base on a given matrix.

  Factory for creating graphs from matrices with entries selectively dropped.

  ## Code paths ##

  Experimental dropping function based on taking a set of random vectors u, running
  a smoother on A u = 0, and then basing the drop decisions on "how smooth" the vectors
  are local. Neighobring regions where the vectors are smooth can be aggregated
  together and so these are kept in the associated drop matrix. Areas that are
  not smooth should end up in different aggregates and so the A_ij representing
  these should be dropped.  This Factory can address both PDE systems and
  scalar PDEs, always creating a matrix reprsenting nodal connections as opposed
  to dof connections.

   To enter this factor as opposed to the more standard CoalesceDropFactory() one
   must set "aggregation: drop scheme" to "unsupported vector smoothing". In this
   case some of the parameter options associated with CoalesceDropFactory (e.g.,
   "aggregation: drop tol", "aggregation: Dirichlet threshold", "lightweight wrap")
   will cause parameter validator errors.

  ## Input/output of SmooVecCoalesceDropFactory ##

  ### User parameters of SmooVecCoalesceDropFactory ###
  Parameter                  | type      | default   | master.xml | validated | requested | description
  ---------------------------|-----------|-----------|:----------:|:---------:|:---------:|------------
   A                         |Factory    | null      |            | *         | *         | Generating factory of the operator A
   "aggregation: drop scheme"|std::string|"classical"| *          | *         | *         | Must choose "unsupported vector smoothing"
   "aggregation: number of times to pre or post smooth"|int| 10|* |           | *         | Amount of pre or post smoothing invocations
   "aggregation: number of random vectors"|int| 10   |          * | *         | *         | Number of random vectors
   "aggregation: penalty parameters"|Array(double)|{12.0,-.20}| * | *         | *         | Ultimately determines how much dropping is done

  The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see SmooVecCoalesceDropFactory::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see SmooVecCoalesceDropFactory::DeclareInput).

  ### Variables provided by UncoupledAggregationFactory ###

  After SmooVecCoalesceDropFactory::Build the following data is available (if requested)

  Parameter | generated by | description
  ----------|--------------|------------
  Graph   | SmooVecCoalesceDropFactory   | Graph of matrix A
  DofsPerNode | SmooVecCoalesceDropFactory | number of DOFs per node. Note, that we assume a constant number of DOFs per node for all nodes associated with the operator A.

*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class SmooVecCoalesceDropFactory : public SingleLevelFactoryBase {
#undef MUELU_SMOOVECCOALESCEDROPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  SmooVecCoalesceDropFactory();

  //! Destructor
  virtual ~SmooVecCoalesceDropFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  /// set predrop function
  void SetPreDropFunction(const RCP<MueLu::PreDropFunctionBaseClass<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& predrop) { predrop_ = predrop; }

  //@}

  void Build(Level& currentLevel) const;  // Build

 private:
  // pre-drop function
  mutable RCP<PreDropFunctionBaseClass> predrop_;

  //! Methods to support compatible-relaxation style dropping
  void badGuysCoalesceDrop(const Matrix& Amat, Teuchos::ArrayRCP<Scalar>& dropParams, LO nPDEs, const MultiVector& smoothedTVecs, const MultiVector& smoothedNull, RCP<GraphBase>& filteredGraph) const;
  void badGuysDropfunc(LO row, const Teuchos::ArrayView<const LocalOrdinal>& indices, const Teuchos::ArrayView<const Scalar>& vals, const MultiVector& smoothedTVecs, LO nPDEs, Teuchos::ArrayRCP<Scalar>& penalties, const MultiVector& smoothedNull, Teuchos::ArrayRCP<LO>& Bcols, Teuchos::ArrayRCP<bool>& keepOrNot, LO& Nbcols, LO nLoc) const;

};  // class SmooVecCoalesceDropFactory

}  // namespace MueLu

#define MUELU_SMOOVECCOALESCEDROPFACTORY_SHORT
#endif  // MUELU_SMOOVECCOALESCEDROPFACTORY_DECL_HPP
