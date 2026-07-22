// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
// ML-headers
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#ifdef ML_MPI  // FIXME: do I need this?
#include <mpi.h>
#endif

// this class
#include "ml_nox_matrixfreelevel.H"

/*----------------------------------------------------------------------*
 |  the class defining a matrix-free coarse level            m.gee 12/04|
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     m.gee 12/04|
 |  IMPORTANT:                                                          |
 |  No matter on which level we are here, the vector xfine is ALWAYS    |
 |  a fine grid vector here!                                            |
 *----------------------------------------------------------------------*/
ML_NOX::ML_Nox_MatrixfreeLevel::ML_Nox_MatrixfreeLevel(int level, int nlevel, int plevel, 
                                 ML* ml, ML_Aggregate* ag, Epetra_CrsMatrix** P,
                                 ML_NOX::Ml_Nox_Fineinterface& interface, const Epetra_Comm& comm,
                                 const Epetra_Vector& xfine, double fd_alpha, double fd_beta,
                                 bool fd_centered, bool isDiagonalOnly, int bsize) 
: fineinterface_(interface),
comm_(comm)                                                          
{
   level_          = level;
   nlevel_         = nlevel;
   ml_printlevel_  = plevel;
   ml_             = ml;
   ag_             = ag;
   fd_alpha_       = fd_alpha;
   fd_beta_        = fd_beta;
   fd_centered_    = fd_centered;
   isDiagonalOnly_ = isDiagonalOnly;
   A_              = 0;
   coarseinterface_= 0;
   bsize_          = bsize;
   

   // we need the graph of the operator on this level. On the fine grid we can just ask the
   // fineinterface for it, on the coarser levels it has to be extracted from the ML-hierachy
   if (level_==0)
   {
      // the Epetra_CrsGraph-copy-constructor shares data with the original one.
      // We want a really deep copy here so we cannot use it
      // graph_ will be given to the FiniteDifferencing class and will be destroyed by it
      graph_ = ML_NOX::deepcopy_graph(interface.getGraph());
   }
   else
   {
      // Note that ML has no understanding of global indices, so it makes up new GIDs
      // (This also holds for the Prolongators P)
      Epetra_CrsMatrix* tmpMat  = 0;
      int               maxnnz  = 0;
      double            cputime = 0.0;
      ML_Operator2EpetraCrsMatrix(&(ml_->Amat[level_]), tmpMat, maxnnz, false, cputime);
      // copy the graph
      double t0 = GetClock();
      graph_ = ML_NOX::deepcopy_graph(&(tmpMat->Graph()));
      // delete the copy of the Epetra_CrsMatrix
      if (tmpMat) delete tmpMat; tmpMat = 0;
      double t1 = GetClock();
      if (ml_printlevel_ > 0 && 0 == comm_.MyPID())
         cout << "matrixfreeML (level " << level_ << "): extraction/copy of Graph in " << cputime+t1-t0 << " sec\n"
              << "                        max-nonzeros in Graph: " << maxnnz << "\n";
   }
   
   // create this levels coarse interface
   coarseinterface_ = new ML_NOX::Nox_CoarseProblem_Interface(fineinterface_,level_,ml_printlevel_,
                                                              P,&(graph_->RowMap()),nlevel_);

   // restrict the xfine-vector to this level 
   Epetra_Vector* xthis = coarseinterface_->restrict_fine_to_this(xfine);
   if (!xthis)
   {
     cout << "**ERR**: ML_Epetra::ML_Nox_MatrixfreeLevel::ML_Nox_MatrixfreeLevel:\n"
          << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_fine_to_this returned NULL on level " 
          << level_ << "\n"
          << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }

   Epetra_Vector* xc = new Epetra_Vector(graph_->RowMap(),false);
   // FIXME: after intesive testing, this test might be obsolet
#if 0
   bool samemap = xc->Map().PointSameAs(xthis->Map());
   if (samemap)
   {
#endif
      xc->Update(1.0,*xthis,0.0);
#if 0
   }
   else
   {
      cout << "**WRN** Maps are not equal in\n"
           << "**WRN** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      // import the xthis vector in the Map that ML produced for graph_
      Epetra_Import* importer = new Epetra_Import(graph_->RowMap(),xthis->Map());  
      int ierr = xc->Import(*xthis,*importer,Insert);
      if (ierr)
      {
         cout << "**ERR**: ML_Epetra::ML_Nox_MatrixfreeLevel::ML_Nox_MatrixfreeLevel:\n"
              << "**ERR**: export from xthis to xc returned err=" << ierr <<"\n"
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      if (importer) delete importer; importer = 0;
   }
#endif   
   if (xthis) delete xthis; xthis = 0;

   // create the coloring of the graph
   if (ml_printlevel_>0 && comm_.MyPID()==0) 
   {
      cout << "matrixfreeML (level " << level_ << "): Entering Coloring on level " << level_ << "\n";
      fflush(stdout);
   }
   double t0 = GetClock();
   colorMap_ = ML_NOX::ML_Nox_collapsedcoloring(graph_,bsize_,isDiagonalOnly,ml_printlevel_);
   if (!colorMap_) colorMap_ = ML_NOX::ML_Nox_standardcoloring(graph_,isDiagonalOnly);
   colorMapIndex_ = new EpetraExt::CrsGraph_MapColoringIndex(*colorMap_);
   colorcolumns_  = &(*colorMapIndex_)(*graph_);
   double t1 = GetClock();
   if (ml_printlevel_>0 && comm_.MyPID()==0)
   {
      cout << "matrixfreeML (level " << level_ << "): Proc " << comm_.MyPID() <<" Coloring time is " << (t1-t0) << " sec\n";
      fflush(stdout);
   }

   // construct the FiniteDifferenceColoring-Matrix
   if (ml_printlevel_>0 && comm_.MyPID()==0)
   {
      cout << "matrixfreeML (level " << level_ << "): Entering Construction FD-Operator on level " << level_ << "\n";
      fflush(stdout);
   }

   t0 = GetClock();

#if 1 // FD-operator with coloring
   FD_ = new NOX::EpetraNew::FiniteDifferenceColoring(*coarseinterface_,
                                                      *xc,
                                                      *graph_,
                                                      *colorMap_,
                                                      *colorcolumns_,
                                                      true,
                                                      isDiagonalOnly_,
                                                      fd_beta_,fd_alpha_);
#else // FD-operator without coloring
   FD_ = new NOX::EpetraNew::FiniteDifference(*coarseinterface_,
                                              *xc,
                                              *graph_,
                                              fd_beta_,fd_alpha_);
#endif
   // set differencing method
   if (fd_centered_)
      FD_->setDifferenceMethod(NOX::EpetraNew::FiniteDifferenceColoring::Centered);

   bool err = FD_->computeJacobian(*xc); 
   if (err==false)
   {
     cout << "**ERR**: ML_NOX::ML_Nox_MatrixfreeLevel::ML_Nox_MatrixfreeLevel:\n"
          << "**ERR**: NOX::Epetra::FiniteDifferenceColoring returned an error on level " 
          << level_ << "\n"
          << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }

   // print number of calls to the coarse interface
   if (ml_printlevel_>0 && comm_.MyPID()==0)
      cout << "matrixfreeML (level " << level_ << "): Calls to coarse-computeF in FD-Operator: "
           << coarseinterface_->numcallscomputeF() << "\n"; 

   t1 = GetClock();
   if (ml_printlevel_>0 && comm_.MyPID()==0)
   {
      cout << "matrixfreeML (level " << level_ << "): Proc " << comm_.MyPID() <<" colored Finite Differencing time is " << (t1-t0) << " sec\n";
      fflush(stdout);
   }

   // get ref to computed Epetra_CrsMatrix  
   A_ = dynamic_cast<Epetra_CrsMatrix*>(&(FD_->getUnderlyingMatrix()));

   // set counter for number of calls to the coarseinterface_->computeF back to zero
   coarseinterface_->resetnumcallscomputeF();   
   
   // tidy up 
   if (xc) delete xc; xc = 0;
   
   return;
}

/*----------------------------------------------------------------------*
 |  recreate this level (public)                             m.gee 01/05|
 |  this function assumes, that the graph of the fine level problem has |
 |  not changed since call to the constructor and therefore             |
 | the graph and it's coloring do not have to be recomputed             |  
 |  IMPORTANT:                                                          |
 |  No matter on which level we are here, the vector xfine is ALWAYS    |
 |  a fine grid vector here!                                            |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_MatrixfreeLevel::recreateLevel(int level, int nlevel, int plevel, 
                                                   ML* ml, ML_Aggregate* ag, Epetra_CrsMatrix** P,
                                                   ML_NOX::Ml_Nox_Fineinterface& interface, 
                                                   const Epetra_Comm& comm, const Epetra_Vector& xfine) 
{
   // make some tests
   if (level != level_)
   {
     cout << "**ERR**: ML_Epetra::ML_Nox_MatrixfreeLevel::recreateLevel:\n"
          << "**ERR**: level_ " << level_ << " not equal level " << level << "\n"
          << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   if (nlevel != nlevel_)
   {
     cout << "**ERR**: ML_Epetra::ML_Nox_MatrixfreeLevel::recreateLevel:\n"
          << "**ERR**: nlevel_ " << nlevel_ << " not equal nlevel " << nlevel << "\n" 
          << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   // printlevel might have changed
   ml_printlevel_ = plevel;
   ml_            = ml;
   ag_            = ag;
   destroyP();  // safer to use the new Ps
   setP(NULL);

   // we need the graph of the operator on this level. On the fine grid we can just ask the
   // fineinterface for it, on the coarser levels it has to be extracted from the ML-hierachy
   bool same;
   if (level_==0)
   {
      const Epetra_CrsGraph* graph = interface.getGraph();
      // check whether the old graph matches the new one
      same = compare_graphs(graph,graph_);
      destroyFD(); // we are here to recompute the FD-operator (this destroys graph_)
      graph_ = ML_NOX::deepcopy_graph(graph);
   }
   else
   {
      // Note that ML has no understanding of global indices, so it makes up new GIDs
      // (This also holds for the Prolongators P)
      Epetra_CrsMatrix* tmpMat  = 0;
      int               maxnnz  = 0;
      double            cputime = 0.0;
      ML_Operator2EpetraCrsMatrix(&(ml_->Amat[level_]), tmpMat, maxnnz, false, cputime);
      // get a view from the graph
      const Epetra_CrsGraph& graph = tmpMat->Graph();
      // compare the graph to the existing one
      same = compare_graphs(&graph,graph_);
      destroyFD(); // we are here to recompute the FD-operator (this destroys graph_)
      double t0 = GetClock();
      graph_ = ML_NOX::deepcopy_graph(&graph);
      // delete the copy of the Epetra_CrsMatrix
      if (tmpMat) delete tmpMat; tmpMat = 0;
      double t1 = GetClock();
      if (ml_printlevel_ > 0 && 0 == comm_.MyPID())
         cout << "matrixfreeML (level " << level_ << "): extraction/copy of Graph in " << cputime+t1-t0 << " sec\n"
              << "                        max-nonzeros in Graph: " << maxnnz << "\n";
   }
   
   // recreate this levels coarse interface
   if (same)
      coarseinterface_->recreate(ml_printlevel_,P,&(graph_->RowMap()));
   else
   {
      delete coarseinterface_;
      coarseinterface_ = new ML_NOX::Nox_CoarseProblem_Interface(fineinterface_,level_,ml_printlevel_,
                                                                 P,&(graph_->RowMap()),nlevel_);
   }
   
   // restrict the xfine-vector to this level 
   Epetra_Vector* xthis = coarseinterface_->restrict_fine_to_this(xfine);
   if (!xthis)
   {
     cout << "**ERR**: ML_Epetra::ML_Nox_MatrixfreeLevel::ML_Nox_MatrixfreeLevel:\n"
          << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_fine_to_this returned NULL on level " 
          << level_ << "\n"
          << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }

   Epetra_Vector* xc = new Epetra_Vector(graph_->RowMap(),false);
   // FIXME: after intesive testing, this test might be obsolet
#if 0
   bool samemap = xc->Map().PointSameAs(xthis->Map());
   if (samemap)
   {
#endif
      xc->Update(1.0,*xthis,0.0);
#if 0
   }
   else
   {
      cout << "**WRN** Maps are not equal in\n"
           << "**WRN** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      // import the xthis vector in the Map that ML produced for graph_
      Epetra_Import* importer = new Epetra_Import(graph_->RowMap(),xthis->Map());  
      int ierr = xc->Import(*xthis,*importer,Insert);
      if (ierr)
      {
         cout << "**ERR**: ML_Epetra::ML_Nox_MatrixfreeLevel::ML_Nox_MatrixfreeLevel:\n"
              << "**ERR**: export from xthis to xc returned err=" << ierr <<"\n"
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      if (importer) delete importer; importer = 0;
   }
#endif   
   if (xthis) delete xthis; xthis = 0;

   // create the coloring of the graph
   if (ml_printlevel_>0 && comm_.MyPID()==0)
   {
      cout << "matrixfreeML (level " << level_ << "): Entering Recoloring on level " << level_ << "\n";
      fflush(stdout);
   }
   double t0 = GetClock();
   if (!same) // te graph has obviously changed, so we need to recolor
   {
      if (colorMap_) delete colorMap_; colorMap_ = 0;
      if (colorMapIndex_) delete colorMapIndex_; colorMapIndex_ = 0;
      if (colorcolumns_) delete colorcolumns_; colorcolumns_ = 0;

      colorMap_ = ML_NOX::ML_Nox_collapsedcoloring(graph_,bsize_,isDiagonalOnly_,ml_printlevel_);
      if (!colorMap_) colorMap_ = ML_NOX::ML_Nox_standardcoloring(graph_,isDiagonalOnly_);
      colorMapIndex_ = new EpetraExt::CrsGraph_MapColoringIndex(*colorMap_);
      colorcolumns_  = &(*colorMapIndex_)(*graph_);
    }
    else if (ml_printlevel_>0 && comm_.MyPID()==0)
      cout << "matrixfreeML (level " << level_ << "): Reusing existing Coloring on level " << level_ << "\n";
    double t1 = GetClock();
    if (ml_printlevel_>5)
    {
      cout << "matrixfreeML (level " << level_ << "): Proc " << comm_.MyPID() <<" (Re)Coloring time is " << (t1-t0) << " sec\n";
      fflush(stdout);
    }

#if 0
   // print the colorMap_
   if (comm_.MyPID()==0)
   cout << "colorMap_\n";
   cout << *colorMap_;
   for (int i=0; i<colorcolumns_->size(); i++)
   {
      if (comm_.MyPID()==0)
         cout << "the " << i << " th colorcolumn_ - vector\n";
      cout << colorcolumns_->operator[](i);
   }
#endif
   
   // construct the FiniteDifferenceColoring-Matrix
   if (ml_printlevel_>0 && comm_.MyPID()==0)
   {
      cout << "matrixfreeML (level " << level_ << "): Entering Construction FD-Operator on level " << level_ << "\n";
      fflush(stdout);
   }

   t0 = GetClock();
#if 1 // FD-operator with coloring (see the #if 1 in ml_nox_matrixfreelevel.H as well!)
   FD_ = new NOX::EpetraNew::FiniteDifferenceColoring(*coarseinterface_,
                                                      *xc,
                                                      *graph_,
                                                      *colorMap_,
                                                      *colorcolumns_,
                                                      true,
                                                      isDiagonalOnly_,
                                                      fd_beta_,fd_alpha_);
#else // FD-operator without coloring
   FD_ = new NOX::EpetraNew::FiniteDifference(*coarseinterface_,
                                              *xc,
                                              *graph_,
                                              fd_beta_,fd_alpha_);
#endif

   // set differencing method
   if (fd_centered_)
   {
      FD_->setDifferenceMethod(NOX::EpetraNew::FiniteDifferenceColoring::Centered);
   }

   bool err = FD_->computeJacobian(*xc); 
   if (err==false)
   {
     cout << "**ERR**: ML_Epetra::ML_Nox_MatrixfreeLevel::ML_Nox_MatrixfreeLevel:\n"
          << "**ERR**: NOX::Epetra::FiniteDifferenceColoring returned an error on level " 
          << level_ << "\n"
          << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }

   t1 = GetClock();
   if (ml_printlevel_>5)
      cout << "matrixfreeML (level " << level_ << "): Proc " << comm_.MyPID() 
           <<" Finite Differencing operator constr. in " << (t1-t0) << " sec\n";

   // get ref to computed Epetra_CrsMatrix  
   A_ = dynamic_cast<Epetra_CrsMatrix*>(&(FD_->getUnderlyingMatrix()));
   
   // print number of calls to the coarse interface
   if (ml_printlevel_>5 && comm_.MyPID()==0)
      cout << "matrixfreeML (level " << level_ << "): Calls to coarse-computeF in FD-Operator: "
           << coarseinterface_->numcallscomputeF() << "\n"; 
   
   // set counter for number of calls to the coarseinterface_->computeF back to zero
   coarseinterface_->resetnumcallscomputeF();   
   
   // tidy up 
   if (xc)       delete xc;       xc = 0;
   
   return true;
}

/*----------------------------------------------------------------------*
 |  compare 2 graphs pointwise in indices                    m.gee 01/05|
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_MatrixfreeLevel::compare_graphs(const Epetra_CrsGraph* newgraph, 
                                                    const Epetra_CrsGraph* oldgraph)
{
   int i,j;
   if (newgraph->NumMyRows()         != oldgraph->NumMyRows() ||
       newgraph->NumMyNonzeros()     != oldgraph->NumMyNonzeros() ||
       newgraph->NumGlobalNonzeros() != oldgraph->NumGlobalNonzeros())
      return false;
   for (i=0; i<newgraph->NumMyRows(); i++)
   {
      if (newgraph->GRID(i) != oldgraph->GRID(i))
         return false;
      int  new_Numindices;
      int* new_Indices;
      newgraph->ExtractMyRowView(i,new_Numindices,new_Indices);
      int  old_Numindices;
      int* old_Indices;
      oldgraph->ExtractMyRowView(i,old_Numindices,old_Indices);
      if (new_Numindices != old_Numindices)
         return false;
      for (j=0; j<new_Numindices; j++)
         if (new_Indices[j] != old_Indices[j])
            return false;
   }
   return true;
}                                                     
 
/*----------------------------------------------------------------------*
 |  Destructor (public)                                      m.gee 12/04|
 *----------------------------------------------------------------------*/
ML_NOX::ML_Nox_MatrixfreeLevel::~ML_Nox_MatrixfreeLevel()
{

// graph_ and A_ are part of FD and are destroyed there 
  if (FD_)              
     destroyFD();

//  if (MapColoring_)     
//     delete MapColoring_;
//  MapColoring_ = 0;
 
  if (colorMap_)        
     delete colorMap_;
  colorMap_ = 0;
 
  if (colorMapIndex_)   
     delete colorMapIndex_;
  colorMapIndex_ = 0;
 
  if (colorcolumns_)    
     delete colorcolumns_;
  colorcolumns_ = 0;
 
  if (coarseinterface_) 
     delete coarseinterface_;
  coarseinterface_ = 0;
 
  return;
}


#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) & defined(HAVE_ML_EPETRAEXT)
