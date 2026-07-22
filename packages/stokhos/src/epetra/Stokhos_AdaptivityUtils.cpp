// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_AdaptivityUtils.hpp"
#include "Stokhos_BasisInteractionGraph.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"

#include "Epetra_IntVector.h"
#include "Epetra_Import.h"

Teuchos::RCP<Epetra_Map> Stokhos::adapt_utils::buildAdaptedRowMapAndOffsets(
        const Epetra_Comm & Comm,
        const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & per_dof_row_basis,
        std::vector<int> & myRowGidOffsets)
{
   // build row offsets
   int totalSz = 0;
   std::vector<int> myRowOffsets(per_dof_row_basis.size());
   for(std::size_t i=0;i<per_dof_row_basis.size();i++) {
      myRowOffsets[i] = totalSz;
      totalSz += per_dof_row_basis[i]->size();
   }

   // build row map 
   Teuchos::RCP<Epetra_Map> rowMap = Teuchos::rcp(new Epetra_Map(-1,totalSz,0,Comm));

   myRowGidOffsets.resize(per_dof_row_basis.size());
   for(std::size_t i=0;i<myRowGidOffsets.size();i++)
      myRowGidOffsets[i] = rowMap->GID(myRowOffsets[i]);

   return rowMap;
}

Teuchos::RCP<Epetra_Map> Stokhos::adapt_utils::buildAdaptedRowMap(
        const Epetra_Comm & Comm,
        const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & per_dof_row_basis)
{
   // build row offsets
   int totalSz = 0;
   for(std::size_t i=0;i<per_dof_row_basis.size();i++)
      totalSz += per_dof_row_basis[i]->size();

   // build row map 
   return Teuchos::rcp(new Epetra_Map(-1,totalSz,0,Comm));
}

void Stokhos::adapt_utils::buildAdaptedColOffsets(
        const Epetra_CrsGraph & determGraph,
        const std::vector<int> & myRowGidOffsets,
        std::vector<int> & myColGidOffsets)
{
   // build global distributed offsets index
   Epetra_IntVector adaptRowGids(determGraph.RowMap()); // uniquely defined row GIDs that will
                                                        // be used to determine col GIDs
   for(std::size_t i=0;i<myRowGidOffsets.size();i++)
      adaptRowGids[i] = myRowGidOffsets[i];

   // perform global communication to determine offsets
   Epetra_Import importer(determGraph.ColMap(),determGraph.RowMap());
   Epetra_IntVector adaptColGids(determGraph.ColMap());
   adaptColGids.Import(adaptRowGids,importer,Insert);

   // copy offsets into std::vector
   myColGidOffsets.resize(adaptColGids.MyLength());
   for(std::size_t i=0;i<myColGidOffsets.size();i++)
      myColGidOffsets[i] = adaptColGids[i]; 
}

void Stokhos::adapt_utils::buildColBasisFunctions(
        const Epetra_CrsGraph & determGraph,
        const Teuchos::RCP<const Stokhos::ProductBasis<int,double> > & masterBasis,
        const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & per_dof_row_basis,
        std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & per_dof_col_basis)
{
   using Teuchos::RCP;

   const Epetra_BlockMap & rowMap = determGraph.RowMap();
   const Epetra_BlockMap & colMap = determGraph.ColMap();

   // this is block size
   int numStochDim = masterBasis->dimension();

   // build blocked maps 
   Epetra_BlockMap blkRowMap(-1,rowMap.NumMyElements(),rowMap.MyGlobalElements(),numStochDim,0,rowMap.Comm());
   Epetra_BlockMap blkColMap(-1,colMap.NumMyElements(),colMap.MyGlobalElements(),numStochDim,0,colMap.Comm());
   // construct int vectors to determine order
   Epetra_IntVector stochRowOrders(blkRowMap),
                    stochColOrders(blkColMap);  // to build built by Import

   // loop over row degrees of freedom building Row Orders vector
   for(std::size_t dof=0;dof<per_dof_row_basis.size();dof++) {
      RCP<const Stokhos::ProductBasis<int,double> > rowBasis = per_dof_row_basis[dof];
      TEUCHOS_TEST_FOR_EXCEPTION(rowBasis->dimension()!=masterBasis->dimension(),std::invalid_argument,
                      "Stokhos::adapt_utils::buildColBasisFunctions: Row basis must match dimension of master basis!");
   
      Teuchos::Array<RCP<const OneDOrthogPolyBasis<int,double> > > onedBasis 
            = rowBasis->getCoordinateBases();

      TEUCHOS_TEST_FOR_EXCEPTION(onedBasis.size()!=numStochDim,std::logic_error,
                      "Stokhos::adapt_utils::buildColBasisFunctions: Wrong number of dimensions from row basis!");

      // fill stochastic orders vector
      for(int i=0;i<numStochDim;i++) 
         stochRowOrders[i+dof*numStochDim] = onedBasis[i]->order();
   }

   // do communication to determine row maps
   Epetra_Import importer(blkColMap,blkRowMap);
   stochColOrders.Import(stochRowOrders,importer,Insert);

   Teuchos::Array<RCP<const OneDOrthogPolyBasis<int,double> > > oneDMasterBasis
         = masterBasis->getCoordinateBases();

   // now construct column basis functions
   std::vector<int> polyOrder(numStochDim,0);
   per_dof_col_basis.resize(blkColMap.NumMyElements());
   for(std::size_t col=0;col<per_dof_col_basis.size();col++) {
      int colOffset = numStochDim*col;
      
      // extract polynomial order for this column
      for(int o=0;o<numStochDim;o++)
         polyOrder[o] = stochColOrders[colOffset+o];

      Teuchos::Array<Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > newBasisArray(numStochDim);
      for(Teuchos::Ordinal dim=0;dim<oneDMasterBasis.size();dim++) {
         RCP<const OneDOrthogPolyBasis<int,double> > oneDBasis = oneDMasterBasis[dim];

         newBasisArray[dim] = oneDBasis->cloneWithOrder(polyOrder[dim]);
      }

      per_dof_col_basis[col] = Teuchos::rcp(
           new Stokhos::CompletePolynomialBasis<int,double>(newBasisArray)); 
   }
}

Teuchos::RCP<Epetra_CrsGraph> Stokhos::adapt_utils::buildAdaptedGraph(
        const Epetra_CrsGraph & determGraph,
        const Teuchos::RCP<const Stokhos::ProductBasis<int,double> > & masterBasis,
        const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & per_dof_row_basis,
        bool onlyUseLinear,
        int kExpOrder)
{
   std::vector<int> myRowGidOffsets, myColGidOffsets;

   return buildAdaptedGraph(determGraph,masterBasis,per_dof_row_basis,myRowGidOffsets,myColGidOffsets,onlyUseLinear,kExpOrder);
}

Teuchos::RCP<Epetra_CrsGraph> Stokhos::adapt_utils::buildAdaptedGraph(
        const Epetra_CrsGraph & determGraph,
        const Teuchos::RCP<const Stokhos::ProductBasis<int,double> > & masterBasis,
        const std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > & per_dof_row_basis,
        std::vector<int> & myRowGidOffsets,std::vector<int> & myColGidOffsets,
        bool onlyUseLinear,
        int kExpOrder)
{
   TEUCHOS_TEST_FOR_EXCEPTION(int(per_dof_row_basis.size())!=determGraph.NumMyRows(),std::logic_error,
                      "Stokhos::adapt_utils::buildAdaptedGraph: per_dof_row_basis.size()!=determGraph.NumMyRows()");

   myRowGidOffsets.clear();
   myColGidOffsets.clear();

   std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > per_dof_col_basis;
   Teuchos::RCP<Epetra_Map> rowMap;

   // build basis functions associated with the columns
   buildColBasisFunctions(determGraph,masterBasis,per_dof_row_basis,per_dof_col_basis);

   // build row map, and row and column offsets.
   rowMap = buildAdaptedRowMapAndOffsets(determGraph.Comm(),per_dof_row_basis,myRowGidOffsets);
   buildAdaptedColOffsets(determGraph,myRowGidOffsets,myColGidOffsets);

   Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*rowMap,0));

   Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk;
   if(kExpOrder<0)
      Cijk = masterBasis->computeTripleProductTensor();
   else
      Cijk = masterBasis->computeLinearTripleProductTensor();

   // iterate over nonzero structure of graph
   int maxNNZ = determGraph.MaxNumNonzeros();
   std::vector<int> determGraphCols(maxNNZ);
   std::vector<int> graphCols;
   for(int lRID=0;lRID<determGraph.NumMyRows();lRID++) {
      int gRID = determGraph.GRID(lRID);
      int numIndices = -1;
      int rowOffsetIndex = myRowGidOffsets[lRID];

      // grab row nonzero structure
      determGraph.ExtractGlobalRowCopy(gRID,maxNNZ,numIndices,&determGraphCols[0]);
      
      for(int i=0;i<numIndices;i++) {
         int gCID = determGraphCols[i];
         int lCID = determGraph.LCID(gCID);
         int colOffsetIndex = myColGidOffsets[lCID];

         Stokhos::BasisInteractionGraph interactGraph(*masterBasis,*Cijk,*per_dof_row_basis[lRID],
                                                                   *per_dof_col_basis[lCID],
                                                                   onlyUseLinear,kExpOrder);
         for(std::size_t basisRow=0;basisRow<interactGraph.rowCount();basisRow++) {
            const std::vector<std::size_t> & basisCols = interactGraph.activeIndices(basisRow);
            graphCols.resize(basisCols.size());
            for(std::size_t g=0;g<basisCols.size();g++) 
               graphCols[g] = basisCols[g] + colOffsetIndex;

            graph->InsertGlobalIndices(basisRow+rowOffsetIndex,graphCols.size(),&graphCols[0]);            
         }
      }
   }
   
   graph->FillComplete(); 

   return graph;
}
