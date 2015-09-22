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
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/

#include "Moertel_config.h"

#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#include <assert.h>
#include "Epetra_MpiComm.h"
#include "mrtr_manager.H"
#include "mrtr_segment_linear1D.H"
#endif

//
// two dimensional example (1D interface) run on two processors
//

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  bool success = false;
  bool verbose = false;
  try {
    int i, MyPID, NumProc, nedge=0, side=0, nnode, ScanSum;
    MPI_Init(&argc,&argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    MyPID = Comm.MyPID();
    NumProc = Comm.NumProc();
    if (NumProc != 2)
      throw "test1.exe needs to be run on 2 processors\n";
    // ------------------------------------------------------------- //
    // create an empty MOERTEL::Interface, in this example just one
    // ------------------------------------------------------------- //
    int printlevel = 8; // ( moertel takes values 0 - 10 )
    MOERTEL::Interface interface(0, true, Comm, printlevel);
    // ------------------------------------------------------------- //
    // Add nodes on both sides of interface to interface
    // ------------------------------------------------------------- //
    if (MyPID == 0) {
      nedge = 4;
      side = 0;
    }
    if (MyPID == 1) {
      nedge = 5;
      side = 1;
    }
    nnode = nedge + 1;
    Comm.ScanSum(&nnode, &ScanSum, 1);
    int *nodeid = new int[nnode];
    double L = 1/double(nedge);
    double coord[3];
    coord[0] = 1;
    coord[2] = 0;
    //
    // nodes on interface
    //
    for (i=0; i<nnode; i++) {
      coord[1] = i*L;
      nodeid[i] = ScanSum - nnode + i;
      MOERTEL::Node node(nodeid[i], coord, 1, &nodeid[i], false, printlevel);
      if(!interface.AddNode(node, side))
        throw "interface AddNode returned false\n";
    }
    //
    // segments on interface
    //
    Comm.ScanSum(&nedge, &ScanSum, 1);
    int econ[2], seg_id;
    for (i=0; i<nedge; i++) {
      // segments have to be defined in anti-clockwise order for each subdomain!
      if (MyPID==0)
      {
        econ[0] = nodeid[i];
        econ[1] = nodeid[i+1];
      }
      else
      {
        econ[1] = nodeid[i];
        econ[0] = nodeid[i+1];
      }
      seg_id = ScanSum - nedge + i + 1;
      MOERTEL::Segment_Linear1D segment(seg_id, 2, econ, printlevel);
      if(!interface.AddSegment(segment, side))
        throw "interface AddSegment returned false\n";
    }
    //
    // mortar (master) side is 0
    //
    interface.SetMortarSide(0);
    interface.SetFunctionTypes(MOERTEL::Function::func_Linear1D,       // primal trace space
        MOERTEL::Function::func_DualLinear1D);  // dual mortar space (recommended)
    if (!interface.Complete())
      throw "Interface completion returned false\n";
    // ------------------------------------------------------------- //
    // create an empty MOERTEL::Manager for 2D problems
    // It organizes everything from integration to solution
    // ------------------------------------------------------------- //
    MOERTEL::Manager manager(Comm, printlevel);
    manager.SetDimension(MOERTEL::Manager::manager_2D);
    // ------------------------------------------------------------- //
    // Add the interface to the manager
    // ------------------------------------------------------------- //
    manager.AddInterface(interface);

    // ------------------------------------------------------------- //
    // we have to supply a rowmap of the whole problem so Moertel knows
    // how the system of equations looks like. Constraint equations will
    // have rows n to n + nlm-1 where n is the number of degrees of freedom
    // and nlm is the number of lagrange multipliers
    // ------------------------------------------------------------- //
    int gnnode;
    Comm.SumAll(&nnode, &gnnode, 1);
    Epetra_Map map(gnnode,nnode,nodeid,0,Comm);
    manager.SetProblemMap(&map);

    // ============================================================= //
    // choose integration parameters
    // ============================================================= //
    //Teuchos::ParameterList& moertelparams = manager.Default_Parameters();
    // this takes effect in 3D only
    //moertelparams.set("exact values at gauss points",true);
    // 1D interface possible values are 1,2,3,4,5,6,7,8,10 (2 recommended with linear shape functions)
    //moertelparams.set("number gaussian points 1D",2);
    // 2D interface possible values are 3,6,12,13,16,19,27 (12 recommended with linear functions)
    //moertelparams.set("number gaussian points 2D",12);

    // ============================================================= //
    // Here we are done with the construction phase of the interface
    // so we can integrate the mortar integrals
    // (Note we have not yet evaluated the PDE at all!)
    // ============================================================= //
    manager.Mortar_Integrate();

    // print interface information
    // (Manager, Interface, Segment, Node implement the << operator)
    if (printlevel) std::cout << manager;


    // get the 2 pieces of the constraint equation
    const Epetra_CrsMatrix* D = manager.D();
    const Epetra_CrsMatrix* M = manager.M();

    // create an empty Epetra_CrsMatrix and add D and M
    Epetra_CrsMatrix constraints(Copy,D->RowMap(),5);
    MOERTEL::MatrixMatrixAdd(*D,false,1.0,constraints,0.0);
    MOERTEL::MatrixMatrixAdd(*M,false,1.0,constraints,1.0);
    constraints.FillComplete(D->DomainMap(),D->RangeMap());
    constraints.OptimizeStorage();

    std::cout << constraints;

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  MPI_Finalize();

  if (success)
    std::cout << "\nTest passed!" << std::endl;

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );

#else

  return EXIT_FAILURE;

#endif

}
