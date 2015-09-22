// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_STK_CustomMeshFactory_hpp__
#define __Panzer_STK_CustomMeshFactory_hpp__

#include <Panzer_Traits.hpp> // for Ordinal64
#include <Panzer_STK_MeshFactory.hpp>
#include <Panzer_STK_Interface.hpp>

namespace panzer_stk_classic {

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
    Teuchos::RCP<STK_Interface> buildMesh(stk_classic::ParallelMachine parallelMach) const;

    virtual Teuchos::RCP<STK_Interface> buildUncommitedMesh(stk_classic::ParallelMachine parallelMach) const;
    virtual void completeMeshConstruction(STK_Interface & mesh,stk_classic::ParallelMachine parallelMach) const;

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

    // wrapper for 
    struct FieldContainer {
      std::size_t _dim0, _dim1;
      double *_buffer;

      FieldContainer(const std::size_t dim0, const std::size_t dim1, double *buffer) 
        : _dim0(dim0), _dim1(dim1), _buffer(buffer) { }

      // row-major indexing: Intrepid_FieldContainerDef.hpp
      double operator()(const std::size_t i, 
                        const std::size_t j) const { return _buffer[i*_dim1+j]; }

      double& operator()(const std::size_t i, 
                         const std::size_t j) { return _buffer[i*_dim1+j]; }
      
    };

    mutable unsigned int machRank_, machSize_;
  };

}

#endif
