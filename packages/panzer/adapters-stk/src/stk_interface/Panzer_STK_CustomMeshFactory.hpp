// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_CustomMeshFactory_hpp__
#define __Panzer_STK_CustomMeshFactory_hpp__

#include <Panzer_Traits.hpp> // for panzer::GlobalOrdinal
#include <Panzer_STK_MeshFactory.hpp>
#include <Panzer_STK_Interface.hpp>

namespace panzer_stk {

  class STK_Interface;

  /** This builds a parallel mesh object with given mesh array structures
   * Parameters
   */
  class CustomMeshFactory : public STK_MeshFactory {
  public:
    //! Constructor
    CustomMeshFactory();

    //! Destructor
    virtual ~CustomMeshFactory();

    //! Build the mesh object
    Teuchos::RCP<STK_Interface> buildMesh(stk::ParallelMachine parallelMach) const;

    virtual Teuchos::RCP<STK_Interface> buildUncommitedMesh(stk::ParallelMachine parallelMach) const;
    virtual void completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const;

    //! From ParameterListAcceptor
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList);

    //! From ParameterListAcceptor
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  protected:
    void initializeWithDefaults();

    void buildMetaData(STK_Interface & mesh) const;

    void buildElements(STK_Interface & mesh) const;
    void addSideSets(STK_Interface & mesh) const;

    void fillSolutionFieldData(STK_Interface & mesh) const;

    int Dimension_;

    int NumBlocks_;

    int NumNodesPerProc_;
    int *Nodes_;

    double *Coords_;

    int NumElementsPerProc_;
    int *BlockIDs_;
    int *Element2Nodes_;

    int OffsetToGlobalElementIDs_;

    double *ChargeDensity_;
    double *ElectricPotential_;

    mutable unsigned int machRank_, machSize_;
  };

}

#endif
