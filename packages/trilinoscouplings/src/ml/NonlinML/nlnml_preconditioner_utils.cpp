/*
#@HEADER
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file nlnml_preconditioner_utils.cpp
 *
 * \brief ML nonlinear preconditioner and solver utilities
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
// ML-headers
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "Epetra_Comm.h"

#include "nlnml_preconditioner_utils.H"

/*----------------------------------------------------------------------*
 |  make a deep copy of a graph                              m.gee 01/05|
 |  allocate the new graph                                              |
 *----------------------------------------------------------------------*/
Epetra_CrsGraph* NLNML::deepcopy_graph(const Epetra_CrsGraph* oldgraph)
{
   int  i;
   int  nrows = oldgraph->NumMyRows();
   int* nIndicesperRow = new int[nrows];

   for (i=0; i<nrows; i++)
      nIndicesperRow[i] = oldgraph->NumMyIndices(i);
   Epetra_CrsGraph* graph = new Epetra_CrsGraph(Copy,oldgraph->RowMap(),oldgraph->ColMap(),
                                                &(nIndicesperRow[0]));
   delete [] nIndicesperRow;
   nIndicesperRow = 0;

   for (i=0; i<nrows; i++)
   {
      int  numIndices;
      int* Indices=0;
      (void) oldgraph->ExtractMyRowView(i,numIndices,Indices);
      (void) graph->InsertMyIndices(i,numIndices,Indices);
   }

   graph->FillComplete();
   return graph;
}

/*----------------------------------------------------------------------*
 |  color a graph in a collapsed way                         m.gee 05/05|
 |  note:                                                               |
 |  this routine makes 2 assumptions:                                   |
 |  - the blocksize is constant everywhere                              |
 |  - degrees of freedom (rows) on a node(block) a contigous            |
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_MapColoring>  NLNML::Collapsedcoloring(Epetra_CrsGraph* cgraph,
                                                          const int bsize,
                                                          bool diagonalonly,
                                                          int printlevel)
{
  // create a new rangemap for the amalgamated graph
  int new_nummyrows     = cgraph->NumMyRows();
  int new_numglobalrows = cgraph->NumGlobalRows();
  int lok=1;
  int gok=1;
  if (new_nummyrows % bsize != 0 || new_numglobalrows % bsize != 0)
  {
    lok = 0;
    if (printlevel>5)
    std::cout << "**WRN**: NLNML::Collapsedcoloring:\n"
         << "**WRN**: cannot amalgamate cgraph\n"
         << "**WRN**: numlocalrows%bsize= " << (new_nummyrows % bsize) << " numglobalrows%bsize= " << (new_numglobalrows % bsize) << std::endl;
  }
  cgraph->Comm().MinAll(&lok,&gok,1);
  if (!gok)
    return null;

  lok=1;
  gok=1;
  // number of local and global rows
  new_nummyrows /= bsize;
  new_numglobalrows /= bsize;
  //int* myRows = new int [new_nummyrows];
  std::vector<int> myRows(new_nummyrows);
  int  counter=0;
  // calculate which global rows are mine
  lok=1;
  gok=1;
  for (int i=0; i<cgraph->RowMap().NumMyElements(); ++i)
  {
    int old_grow = cgraph->RowMap().GID(i);
    if (old_grow<0)
    {
      if (printlevel>5)
      std::cout << "**WRN**: NLNML::Collapsedcoloring:\n"
           << "**WRN**: cgraph->RowMap().GID() returned " <<  old_grow << std::endl
           << "**WRN**: switching to none-collapsed coloring\n"
           << "**WRN**: file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      lok=0;
    }
    if (lok==0)
      break;
    if (old_grow % bsize != 0)
      continue;
    int new_grow = old_grow/bsize;
    myRows[counter] = new_grow;
    ++counter;
  }
  cgraph->Comm().MinAll(&lok,&gok,1);
  if (!gok)
    return null;

  if (counter != new_nummyrows)
  {
    if (printlevel>5)
    std::cout << "**WRN**: NLNML::Collapsedcoloring:\n"
         << "**WRN**: collapsed coloring problem, switching to standard coloring\n"
         << "**WRN**: file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    lok=0;
  }
  cgraph->Comm().MinAll(&lok,&gok,1);
  if (!gok)
    return null;

  // create the BlockMap
  Epetra_BlockMap newrowmap(new_numglobalrows,new_nummyrows,&myRows[0],1,0,cgraph->Comm());
  myRows.clear();

  // create the nodal cgraph
  Epetra_CrsGraph nodegraph(Copy,newrowmap,27);

  // loop over NumMyElements in old graph and insert every bsize row in
  // new graph
  int  new_length=200;
  std::vector<int> new_gindices(new_length);
  //int* new_gindices = new int[new_length];

  lok=1;
  gok=1;
  for (int i=0; i<cgraph->RowMap().NumMyElements(); ++i)
  {
    int old_grow = cgraph->RowMap().GID(i);
    if (old_grow<0)
    {
      if (printlevel>5)
      std::cout << "**WRN**: NLNML::Collapsedcoloring:\n"
           << "**WRN**: cgraph->RowMap().GID() returned " <<  old_grow << std::endl
           << "**WRN**: switching to none-collapsed coloring\n"
           << "**WRN**: file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      lok=0;
    }
    if (old_grow % bsize != 0)
      continue;
    if (lok==0)
      break;
    int new_grow = old_grow/bsize;

    // extract the row from the old cgraph
    int old_numindices;
    int* old_colindices;
    int err = cgraph->ExtractMyRowView(i,old_numindices,old_colindices);
    if (err)
    {
      if (printlevel>5)
      std::cout << "**WRN**: NLNML::Collapsedcoloring:\n"
           << "**WRN**: cgraph->ExtractMyRowView returned " <<  err << std::endl
           << "**WRN**: switching to none-collapsed coloring\n"
           << "**WRN**: file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      lok=0;
    }
    if (lok==0)
      break;

    // check whether the allocated row for the new cgraph is long enough
    if (old_numindices>new_length)
    {
      new_length = old_numindices;
      new_gindices.resize(new_length);
    }

    // calculate the row in global numbering for the new cgraph
    int new_numindices=0;
    for (int j=0; j<old_numindices; ++j)
    {
      int old_gcol = cgraph->ColMap().GID(old_colindices[j]);
      if (old_gcol<0)
      {
        if (printlevel>5)
        std::cout << "**WRN**: NLNML::Collapsedcoloring:\n"
             << "**WRN**: cgraph->ColMap().GID(old_colindices[j]) returned " <<  old_gcol << std::endl
             << "**WRN**: file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      }
      if (lok==0)
        break;
      if (old_gcol % bsize != 0)
        continue;

      int new_gcol = old_gcol/bsize;
      new_gindices[new_numindices] = new_gcol;
      ++new_numindices;
    }

    // insert the row in the new graph
    err = nodegraph.InsertGlobalIndices(new_grow,new_numindices,&new_gindices[0]);
    if (err!=0 && err !=1)
    {
      if (printlevel>5)
      std::cout << "**WRN**: NLNML::Collapsedcoloring:\n"
           << "**WRN**: nodegraph.InsertGlobalIndices returned " <<  err << std::endl
           << "**WRN**: switching to none-collapsed coloring\n"
           << "**WRN**: file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      lok=0;
    }
    if (lok==0)
      break;
  }
  // catch all errors that may have occurred in this loop
  cgraph->Comm().MinAll(&lok,&gok,1);
  if (!gok)
    return null;

  nodegraph.FillComplete();
  nodegraph.OptimizeStorage();

  // tidy up
  new_gindices.clear();

  EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType =
                                  EpetraExt::CrsGraph_MapColoring::JONES_PLASSMAN;

  RefCountPtr<EpetraExt::CrsGraph_MapColoring> node_MapColoring =
    rcp(new EpetraExt::CrsGraph_MapColoring(algType,0,diagonalonly,0));

  Epetra_MapColoring* node_colorMap = &(*node_MapColoring)(nodegraph);

  // now get the colors out of the node_colorMap and create a colorMap for the original system
  int  node_ncolors = node_colorMap->NumColors();
  int* node_colors  = node_colorMap->ElementColors();
  int* node_loc     = node_colorMap->ListOfColors();
  std::vector<int> col_colors(cgraph->ColMap().NumMyElements());

  // node_ncolors is a local value, we need the global highest color number here
  int lmax = 0;
  int gmax = 0;
  for (int i=0; i<node_ncolors; ++i)
    if (lmax<node_loc[i]) lmax = node_loc[i];
  cgraph->Comm().MaxAll(&lmax,&gmax,1);

  // expand the colors to point wise columns
  int j=0;
  for (int i=0; i<nodegraph.ColMap().NumMyElements(); ++i)
    for (int k=0; k<bsize; ++k)
    {
       col_colors[j] = node_colors[i]+k*gmax;
       ++j;
    }
  lok=1;
  gok=1;
  if (j!=cgraph->ColMap().NumMyElements())
  {
    if (printlevel>5)
    std::cout << "**ERR**: NLNML::Collapsedcoloring:\n"
         << "**ERR**: error expanding node colors to column colors\n"
         << "**WRN**: switching to none-collapsed coloring\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    lok=0;
  }
  cgraph->Comm().MinAll(&lok,&gok,1);
  if (!gok)
    return null;

  RefCountPtr<Epetra_MapColoring> colorMap =
    rcp(new Epetra_MapColoring(cgraph->ColMap(),&col_colors[0]));

  return colorMap;
}

/*----------------------------------------------------------------------*
 |  color a graph in a standard way                          m.gee 05/05|
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_MapColoring> NLNML::Standardcoloring(
                                                    Epetra_CrsGraph* graph,
                                                    bool diagonalonly)
{
  EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType =
                                  EpetraExt::CrsGraph_MapColoring::JONES_PLASSMAN;

  RefCountPtr<EpetraExt::CrsGraph_MapColoring> MapColoring =
    rcp(new EpetraExt::CrsGraph_MapColoring(algType,0,diagonalonly,0));

  RefCountPtr<Epetra_MapColoring> colorMap = rcp(&(*MapColoring)(*graph));

  return colorMap;
}

/*----------------------------------------------------------------------*
 |                                                           m.gee 11/05|
 *----------------------------------------------------------------------*/
#include <cstdio>
bool NLNML::Print_Epetra_CrsMatrix(Epetra_CrsMatrix& matrix)
{
  using namespace std;
  for (int i=0; i<matrix.NumMyRows(); ++i)
  {
    printf("Lrow %5d:  ",i); fflush(stdout);
    int numentries;
    int* indices;
    double* values;
    matrix.ExtractMyRowView(i,numentries,values,indices);
    for (int j=0; j<numentries; ++j)
    printf("%5d %10.3e   ",indices[j],values[j]);
    printf("\n"); fflush(stdout);
  }
  return true;
}

/*----------------------------------------------------------------------*
 |                                                           m.gee 11/05|
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_CrsMatrix> NLNML::StripZeros(
                             RefCountPtr<Epetra_CrsMatrix> in, double eps)
{
  RefCountPtr<Epetra_CrsMatrix> out =
                           rcp(new Epetra_CrsMatrix(Copy,in->RowMap(),200));
  for (int lrow=0; lrow<in->NumMyRows(); ++lrow)
  {
    int grow = in->GRID(lrow);
    if (grow<0) { std::cout << "ERROR: grow<0 \n"; exit(0); }
    int numentries;
    int* lindices;
    double* values;
    int err  = in->ExtractMyRowView(lrow,numentries,values,lindices);
    if (err) { std::cout << "ExtractMyRowView returned " << err << std::endl; exit(0);}
    for (int j=0; j<numentries; ++j)
    {
      int lcol = lindices[j];
      int gcol = in->GCID(lcol);
      if (gcol<0) { std::cout << "ERROR: gcol<0 \n"; exit(0); }
      if (fabs(values[j])<eps && gcol != grow)
        continue;
      int err = out->InsertGlobalValues(grow,1,&values[j],&gcol);
      if (err) { std::cout << "InsertGlobalValues returned " << err << std::endl; exit(0);}
    }
  }
  out->FillComplete();
  out->OptimizeStorage();
  return out;
}

#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
