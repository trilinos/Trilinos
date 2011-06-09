/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Glen Hansen (Glen.Hansen@inl.gov)
#
# ************************************************************************
#@HEADER
*/
/*!
 * \file TwoSquares.cpp
 *
 * \brief Simple serial example showing Moertel usage and solver interfaces
 *
 * \date Last update do Doxygen: 20-March-06
 *
 */
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// MOERTEL headers
#include "mrtr_manager.H"
#include "mrtr_segment_linear1D.H"

// Galeri headers
#include "Galeri_Utils.h"
#include "Galeri_FiniteElements.h"

//#if defined(HAVE_TRIOS_EXODUS) || defined(HAVE_MOERTEL_EXODUS)
#ifdef HAVE_MOERTEL_EXODUS
#include "ExodusInterface.h"
#endif

using namespace Galeri;
using namespace Galeri::FiniteElements;

// ==========================================================
// This file solves the scalar problem
//
//   - \mu \nabla u + \sigma u = f    on \Omega
//                           u = g    on \partial \Omega
// 
// where \Omega is a 2D rectangle, divided into triangles.
// The input grid should be generated using similar 
// conventions of file galeri/data/TwoSquares.m:
// - the bc ID of 10 and 20 are for the mortar interface
// - the bc ID of 0 is for the external domain.
// ==========================================================

double Diffusion(const double& x, const double& y, const double& z)
{
  return (1.0);
}

double Source(const double& x, const double& y, const double& z)
{
  return (0.0);
}

double Force(const double& x, const double& y, const double& z)
{
  if (y < 0.8)
    return(1.0);
  else
    return(0.0);
}

// Specifies the boundary condition.
int BoundaryType(const int& Patch)
{
  if (Patch == 10 || Patch == 20)
    return(GALERI_DO_NOTHING);
  else
    return(GALERI_DIRICHLET);
}

// Specifies the boundary condition.
double BoundaryValue(const double& x, const double& y, 
                     const double& z, const int& Patch)
{
  if (x == -1.0 || x == 1.0)
    return(0.0);
  else if (Patch == 10 || Patch == 20)
    return(1.0);
  else
    return (0.0);
}

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int status = 0; // return status

  try {

    // this example is in serial only
    if (Comm.NumProc()>1) exit(0);

    // read grid from file, see also TwoSquares.m used to generate the grids
    /*
    %  +------+
    %  |  S2  |
    %  +------+ <- edge has tag 20 and 9 elements
    %  +------+ <- edge has tag 10 and 10 elements
    %  |  S1  |
    %  +------+
    %
    % where S1 = (-1,1) x (-1,1) and S2 = (-1, 1) x (1, 3).
    */

    FileGrid Grid(Comm, "TwoSquares.grid");
    
    
    
    
    // create a list of all nodes that are linked to a face
    // with tag 10 and tag 20
    map<int,int> nodes10;
    map<int,int> nodes20;
    for (int i=0; i<Grid.NumMyBoundaryFaces(); ++i)
    {
      int tag;
      int nodeids[2];
      Grid.FaceVertices(i,tag,nodeids);
      if (tag==10)
      {
        nodes10[nodeids[0]] = nodeids[0];
        nodes10[nodeids[1]] = nodeids[1];
      }
      else if (tag==20)
      {
        nodes20[nodeids[0]] = nodeids[0];
        nodes20[nodeids[1]] = nodeids[1];
      }
      else 
        continue;
    }
    
    // ------------------------------------------------------------- //
    // create an empty MOERTEL::Interface, in this example just one
    // ------------------------------------------------------------- //
    int printlevel = 0; // ( moertel takes values 0 - 10 )
    //int printlevel = 9; // ( moertel takes values 0 - 10 )
    MOERTEL::Interface interface(0,true,Comm,printlevel);
    
    // ------------------------------------------------------------- //
    // Add nodes on both sides of interface to interface
    // loop all nodes in the maps nodes10 and nodes20 and add them
    // to the interface with unique ids
    // tag 10 will become interface side 0
    // tag 20 will become interface side 1
    // ------------------------------------------------------------- //
    map<int,int>::iterator curr;
    // do tag==10 or interface side 0
    for (curr = nodes10.begin(); curr != nodes10.end(); ++curr)
    {
      // get unique node id (here it's equal to the degree of freedom on that node)
      int nodeid = curr->second;
      // get node coordinates
      double coord[3];
      Grid.VertexCoord(nodeid,coord);
      // get dirichlet boundary conditions
      double bou = BoundaryValue(coord[0],coord[1],coord[2],10);
      bool dboundary = false;
      if (bou==0.0) 
        dboundary = true;
      // create a moertel node
      MOERTEL::Node node(nodeid,coord,1,&nodeid,dboundary,printlevel);
      // add node to the interface on side 0
      interface.AddNode(node,0);
    }
    nodes10.clear();
    
    // do tag==20 or interface side 1
    for (curr = nodes20.begin(); curr != nodes20.end(); ++curr)
    {
      // get unique node id (here it's equal to the degree of freedom on that node)
      int nodeid = curr->second;
      // get node coordinates
      double coord[3];
      Grid.VertexCoord(nodeid,coord);
      // get dirichlet boundary conditions
      double bou = BoundaryValue(coord[0],coord[1],coord[2],20);
      bool dboundary = false;
      if (bou==0.0) 
        dboundary = true;
      // create a moertel node
      MOERTEL::Node node(nodeid,coord,1,&nodeid,dboundary,printlevel);
      // add node to the interface on side 1
      interface.AddNode(node,1);
    }
    nodes20.clear();
    
    // ------------------------------------------------------------- //
    // add segments on both sides of the interface to the interface
    // ------------------------------------------------------------- //
    for (int i=0; i<Grid.NumMyBoundaryFaces(); ++i)
    {
      int tag;
      int nodeids[2];
      Grid.FaceVertices(i,tag,nodeids);
      if (tag != 10 && tag != 20)
        continue;
      // create a segment (galeri calls it a face)
      MOERTEL::Segment_Linear1D segment(i,2,nodeids,printlevel);
      
      // add it to the interface on side 0
      if (tag==10)
        interface.AddSegment(segment,0);
      // add it to the interface on side 1
      else if (tag==20)
        interface.AddSegment(segment,1);
    }

    // ------------------------------------------------------------- //
    // choose the mortar side of the interface (0 or 1)
    // here: let the package choose it (-2)
    // ------------------------------------------------------------- //
    interface.SetMortarSide(-2);

    // ------------------------------------------------------------- //
    // As we do not know the mortar side yet (we decided to le the
    // package choose it), we can not set a dual trace function (mortar space)
    // as we don't know the side to set it to
    // so we just give orders for the function type
    // ------------------------------------------------------------- //
    interface.SetFunctionTypes(MOERTEL::Function::func_Linear1D,       // primal trace space
                               MOERTEL::Function::func_DualLinear1D);  // dual mortar space (recommended)
                               //MOERTEL::Function::func_Linear1D);    // mortar space (not recommended)

    // ------------------------------------------------------------- //
    // complete the interface
    // ------------------------------------------------------------- //
    if (!interface.Complete())
    {
       cout << "Interface completion returned false\n";
       exit(EXIT_FAILURE);
    }

    // ------------------------------------------------------------- //
    // create an empty MOERTEL::Manager for 2D problems
    // It organizes everything from integration to solution
    // ------------------------------------------------------------- //
    MOERTEL::Manager manager(Comm,printlevel);
    manager.SetDimension(MOERTEL::Manager::manager_2D);
    
    // ------------------------------------------------------------- //
    // Add the interface to the manager
    // ------------------------------------------------------------- //
    manager.AddInterface(interface);

    // ------------------------------------------------------------- //
    // for mortar integration, the mortar manager needs to know about
    // the rowmap of the original (uncoupled) problem because it will
    // create coupling matrices D and M matching that rowmap
    // ------------------------------------------------------------- //
    manager.SetProblemMap(&Grid.RowMap());

    // ============================================================= //
    // choose integration parameters
    // ============================================================= //
    Teuchos::ParameterList& moertelparams = manager.Default_Parameters();
    // this does not affect this 2D case
    moertelparams.set("exact values at gauss points",true);
    // 1D interface possible values are 1,2,3,4,5,6,7,8,10 (2 recommended with linear shape functions)
    moertelparams.set("number gaussian points 1D",2);
    // 2D interface possible values are 3,6,12,13,16,19,27
    moertelparams.set("number gaussian points 2D",27);

    // ============================================================= //
    // Here we are done with the construction phase of the interface
    // so we can integrate the mortar integrals
    // (Note we have not yet evaluated the PDE at all!)
    // ============================================================= //
    manager.Mortar_Integrate();
    
    // print interface information
    // (Manager, Interface, Segment, Node implement the << operator)
    if (printlevel) cout << manager;
        
    // ======================================================== //
    // Prepares the linear system. This requires the definition //
    // of a quadrature formula compatible with the grid, a      //
    // variational formulation, and a problem object which take //
    // care of filling matrix and right-hand side.              //
    // NOTE:
    // we are doing this AFTER we did all the mortar stuff to
    // show that the mortar integration is actually PDE-independent
    // ======================================================== //
    Epetra_CrsMatrix A(Copy, Grid.RowMap(), 0);
    Epetra_Vector    LHS(Grid.RowMap(),true);
    Epetra_Vector    RHS(Grid.RowMap());

    int NumQuadratureNodes = 3;

    GalerkinVariational<TriangleQuadrature>
      Laplace2D(NumQuadratureNodes, Diffusion, Source, Force, 
                BoundaryValue, BoundaryType);

    LinearProblem FiniteElementProblem(Grid, Laplace2D, A, LHS, RHS); 
    FiniteElementProblem.Compute();

    // ============================================================= //
    // this is Galeri's dense solve method if you'd like to see how
    // the uncoupled solution looks like
    // ============================================================= //
    //Solve(&A, &LHS, &RHS);

    // ============================================================= //
    // Since we now have all the pieces together, let's use the 
    // MOERTEL interface to other Trilinos packages to solve the
    // problem
    // ============================================================= //
    
    // ------------------------------------------------------------- //
    // Create a Teuchos::ParameterList to hold solver arguments and also
    // to hold arguments for connected packages AztecOO, ML and Amesos
    // ------------------------------------------------------------- //
    Teuchos::ParameterList list;
    
    // ------------------------------------------------------------- //
    // Choose which type of system of equations to generate
    // Note that only when using DUAL mortar spaces an spd system 
    // can be generated
    // ------------------------------------------------------------- //
    //list.set("System","SaddleSystem");
    list.set("System","SPDSystem");
    
    // ------------------------------------------------------------- //
    // choose solver, currently there is a choice of Amesos and ML/AztecOO
    // Note that if "SaddleSystem" was chosen as system of equations
    // ML/AztecOO doesn't work
    // ------------------------------------------------------------- //
    list.set("Solver","Amesos");
    //list.set("Solver","ML/Aztec");
    
    // ------------------------------------------------------------- //
    // create sublists for packages Amesos, ML, AztecOO. they will be
    // passed on to the individual package that is used
    // ------------------------------------------------------------- //

    // Amesos parameters:
    Teuchos::ParameterList& amesosparams = list.sublist("Amesos");
    amesosparams.set("Solver","Amesos_Klu");
    amesosparams.set("PrintTiming",true);
    amesosparams.set("PrintStatus",true);
    amesosparams.set("UseTranspose",true);
    
    // AztecOO parameters
    Teuchos::ParameterList& aztecparams = list.sublist("Aztec");
    aztecparams.set("AZ_solver","AZ_cg");
    // This will involve ML as preconditioner
    // See the AztecOO manual for other options
    aztecparams.set("AZ_precond","AZ_user_precond");
    aztecparams.set("AZ_max_iter",1200);
    aztecparams.set("AZ_output",100);
    aztecparams.set("AZ_tol",1.0e-7);
    aztecparams.set("AZ_scaling","AZ_none");
        
    // ML parameters
    // As Moertel comes with his own special mortar multigrid hierachy
    // based on ML's smoothed aggregation, not all ML parameters are recognized
    // It basically recognizes everything that recognized by ML's MLAPI
    // (ML Application Programming Interface), see MLAPI documentation
    Teuchos::ParameterList& mlparams = list.sublist("ML");
    ML_Epetra::SetDefaults("SA",mlparams);
    mlparams.set("output",10);
    mlparams.set("print unused",1/*-2*/);
    mlparams.set("increasing or decreasing","increasing");
    mlparams.set("PDE equations",1);
    mlparams.set("max levels",10);
    mlparams.set("coarse: max size",80);
    mlparams.set("aggregation: type","Uncoupled");
    mlparams.set("aggregation: damping factor",1.33);
    
    // original   : The unmodified ML (smoothed) aggregation prolongator
    // mod_simple : ( R * (I-B*W^T) )^T
    // mod_middle : ( (I - R B*W^T*P) * R * (I-B*W^T) )^T
    // mod_full   : ( (I - R B*W^T*P) * R * (I-B*W^T) )^T + ( R B*W^T*P * R * B*W^T )^T
    mlparams.set("prolongator: type","mod_middle"); 
    
    // solvers/smoothers currently recognized by the MLAPI_InverseOperator are
    // Ifpack:
    //         "Jacobi" "Gauss-Seidel" "symmetric Gauss-Seidel"
    //         "ILU" "ILUT" "IC" "ICT" "LU" "Amesos" "Amesos-KLU"
    //         and accompanying parameters
    // ML:
    //         "MLS" "ML MLS" "ML symmetric Gauss-Seidel"
    //         "ML Gauss-Seidel"
    //         and accompanying parameters
    mlparams.set("coarse: type","Amesos-KLU"); 
    mlparams.set("smoother: type","symmetric Gauss-Seidel"); 
    mlparams.set("smoother: MLS polynomial order",3);
    mlparams.set("relaxation: min diagonal value",0.1); 
    mlparams.set("smoother: damping factor",0.67);
    mlparams.set("smoother: sweeps",1);
    mlparams.set("smoother: pre or post","both");
    // the ns for Laplace is the constant
    int dimnullspace = 1;
    int nummyrows = manager.ProblemMap()->NumMyElements();
    int dimnsp    = dimnullspace*nummyrows;
    double* nsp   = new double[dimnsp];
    for (int i=0; i<dimnsp; ++i) nsp[i] = 1.;
    mlparams.set("null space: type","pre-computed");
    mlparams.set("null space: add default vectors",false);
    mlparams.set("null space: dimension",dimnullspace);
    mlparams.set("null space: vectors",nsp);
        
    // ------------------------------------------------------------- //
    // Pass input matrix to Moertel, 
    // Moertel does NOT take ownership of A!
    // ------------------------------------------------------------- //
    manager.SetInputMatrix(&A,false);
    
    // ============================================================= //
    // Solve
    // ============================================================= //
    manager.Solve(list,LHS,RHS);

    // ------------------------------------------------------------- //
    // One can reset the solver, change parameters and/or matrix (with the
    // same rowmap) and solve again if needed.
    // If no ResetSolver() is called, the same matrix and preconditioner
    // will be used to solve for multiple rhs
    // ------------------------------------------------------------- //
    //manager.ResetSolver();
    //LHS.PutScalar(0.0);
    //manager.SetInputMatrix(&A,false);
    //manager.Solve(list,LHS,RHS);
	
//#if defined(HAVE_TRIOS_EXODUS) || defined(HAVE_MOERTEL_EXODUS)
#ifdef HAVE_MOERTEL_EXODUS

    // ==================    //
    // Output using ExodusII //
    // ==================    //
    ExodusInterface exodus(Comm);
    exodus.Write(Grid, "output", LHS);
#else

    // ================== //
    // Output using MEDIT //
    // ================== //
    
    MEDITInterface MEDIT(Comm);
    MEDIT.Write(Grid, "output", LHS);
#endif

  }
  catch (int e) {
    cerr << "Caught exception, value = " << e << endl;
	status = 1;
  }
  catch (Galeri::Exception& rhs) {

			cerr << "Caught Galeri exception: ";
			rhs.Print();
	status = 1;
  }
  catch (...) {
    cerr << "Caught generic exception" << endl;
	status = 1;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

    if (status == 0) 
      std::cout << "\nTest passed!" << endl;
    else
      std::cout << "\nTest Failed!" << endl;


 // Final return value (0 = successfull, non-zero = failure)

  return status;
}
