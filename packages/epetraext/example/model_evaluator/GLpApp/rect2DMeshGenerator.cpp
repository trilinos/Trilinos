/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "rect2DMeshGenerator.hpp"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RefCountPtr.hpp"

void GLpApp::rect2DMeshGenerator(
  const int                      numProc
  ,const int                     procRank
  ,const double                  len_x
  ,const double                  len_y
  ,const int                     local_nx
  ,const int                     local_ny
  ,const int                     bndy_marker
  ,Epetra_IntSerialDenseVector   *ipindx_out
  ,Epetra_SerialDenseMatrix      *ipcoords_out
  ,Epetra_IntSerialDenseVector   *pindx_out
  ,Epetra_SerialDenseMatrix      *pcoords_out
  ,Epetra_IntSerialDenseMatrix   *t_out
  ,Epetra_IntSerialDenseMatrix   *e_out
  ,std::ostream                  *out_arg
  ,const bool                    dumpAll
  )
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::getFancyOStream(Teuchos::rcp(out_arg,false));
  Teuchos::OSTab tab(out);
  if(out.get())
    *out << "\nGenerating rectangular mesh with len_x="<<len_x<<", len_y="<<len_y<<", local_nx="<<local_nx<<", and local_ny="<<local_ny<<" ...\n";
  //
  // Validate input
  //
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(len_x <= 0.0);
  TEST_FOR_EXCEPT(len_y <= 0.0);
  TEST_FOR_EXCEPT(local_nx <= 0);
  TEST_FOR_EXCEPT(local_ny <= 0);
  TEST_FOR_EXCEPT(ipindx_out==NULL);
  TEST_FOR_EXCEPT(ipcoords_out==NULL);
  TEST_FOR_EXCEPT(pindx_out==NULL);
  TEST_FOR_EXCEPT(pcoords_out==NULL);
  TEST_FOR_EXCEPT(t_out==NULL);
  TEST_FOR_EXCEPT(e_out==NULL);
#endif
  //
  // Get local references
  //
  Epetra_IntSerialDenseVector   &ipindx = *ipindx_out;
  Epetra_SerialDenseMatrix      &ipcoords = *ipcoords_out;
  Epetra_IntSerialDenseVector   &pindx = *pindx_out;
  Epetra_SerialDenseMatrix      &pcoords = *pcoords_out;
  Epetra_IntSerialDenseMatrix   &t = *t_out;
  Epetra_IntSerialDenseMatrix   &e = *e_out;
  //
  // Get dimensions
  //
  const int
    numip = ( local_nx + 1 ) * ( local_ny + ( procRank == numProc-1 ? 1 : 0 ) ),
    numcp = ( procRank == numProc-1 ? 0 : (local_nx+1) ),
    nump  = numip + numcp,
    numelems = 2*local_nx*local_ny,
    numedges = 2*local_ny + ( procRank==0 ? local_nx : 0 ) + ( procRank==numProc-1 ? local_nx : 0 );
  const double
    delta_x = len_x / local_nx,
    delta_y = (len_y/numProc) / local_ny;
  if(out.get()) {
    *out
      << "\nnumip = " << numip
      << "\nnumcp = " << numcp
      << "\nnump = " << nump
      << "\nnumelems = " << numelems
      << "\nnumedges = " << numedges
      << "\ndelta_x = " << delta_x
      << "\ndelta_y = " << delta_y
      << "\n";
  }
  //
  // Resize mesh storage
  //
  ipindx.Size(numip);
  ipcoords.Shape(numip,2);
  pindx.Size(nump);
  pcoords.Shape(nump,2);
  t.Shape(numelems,3);
  e.Shape(numedges,3);
  //
  // Get the global offsets for the local nodes IDs and y coordinates
  //
  const int     localNodeOffset   = procRank*(local_nx+1)*local_ny;
  const double  localYCoordOffset = procRank*delta_y*local_ny;
  //
  // Set the local node IDs and their cooridinates
  //
  if(out.get()) *out << "\nGenerating the node list ...\n";
  tab.incrTab();
  // Owned and shared
  int node_i  = 0;
  int global_node_id = localNodeOffset+1;
  for( int i_y = 0; i_y < local_ny + 1; ++i_y ) {
    for( int i_x = 0; i_x < local_nx + 1; ++i_x ) {
      pindx(node_i) = global_node_id;
      pcoords(node_i,0) = i_x * delta_x;
      pcoords(node_i,1) = i_y * delta_y + localYCoordOffset;
      if(out.get() && dumpAll)
        *out << "node_i="<<node_i<<",global_node_id="<<global_node_id<<",("<<pcoords(node_i,0)<<","<<pcoords(node_i,1)<<")\n";
      ++node_i;
      ++global_node_id;
    }
  }
  tab.incrTab(-1);
  TEST_FOR_EXCEPT(node_i != nump);
  // Locally owned only
  for( int i = 0; i < numip; ++i ) {
    ipindx(i) = pindx(i);
    ipcoords(i,0) = pcoords(i,0);
    ipcoords(i,1) = pcoords(i,1);
  }
  //
  // Set the elements
  //
  if(out.get()) *out << "\nGenerating the element list ...\n";
  tab.incrTab();
  int ele_k = 0;
  global_node_id = localNodeOffset+1;
  for( int i_y = 0; i_y < local_ny; ++i_y ) {
    for( int i_x = 0; i_x < local_nx; ++i_x ) {
      t(ele_k,0) = global_node_id;
      t(ele_k,1) = global_node_id + (local_nx + 1);
      t(ele_k,2) = global_node_id + 1;
      if(out.get() && dumpAll)
        *out << "ele_k="<<ele_k<<",("<<t(ele_k,0)<<","<<t(ele_k,1)<<","<<t(ele_k,2)<<")\n";
      ++ele_k;
      t(ele_k,0) = global_node_id + 1;
      t(ele_k,1) = global_node_id + (local_nx + 1);
      t(ele_k,2) = global_node_id + (local_nx + 1) + 1;
      if(out.get() && dumpAll)
        *out << "ele_k="<<ele_k<<",("<<t(ele_k,0)<<","<<t(ele_k,1)<<","<<t(ele_k,2)<<")\n";
      ++ele_k;
      ++global_node_id;
    }
    ++global_node_id;
  }
  tab.incrTab(-1);
  TEST_FOR_EXCEPT(ele_k != numelems);
  //
  // Set the edges
  //
  int edge_j = 0;
  if(procRank==0) {
    // Bottom edges
    if(out.get()) *out << "\nGenerating the bottom edges ...\n";
    tab.incrTab();
    global_node_id = localNodeOffset+1;
    for( int i_x = 0; i_x < local_nx; ++i_x ) {
      e(edge_j,0) = global_node_id;
      e(edge_j,1) = global_node_id + 1;
      e(edge_j,2) = bndy_marker;
      if(out.get() && dumpAll)
        *out << "edge_j="<<edge_j<<",("<<e(edge_j,0)<<","<<e(edge_j,1)<<"),"<<e(edge_j,2)<<"\n";
      ++edge_j;
      global_node_id += 1;
    }
    tab.incrTab(-1);
  }
  // Left edges
  if(out.get()) *out << "\nGenerating the left edges ...\n";
  tab.incrTab();
  global_node_id = localNodeOffset+1;
  for( int i_y = 0; i_y < local_ny; ++i_y ) {
    e(edge_j,0) = global_node_id;
    e(edge_j,1) = global_node_id + (local_nx + 1);
    e(edge_j,2) = bndy_marker;
    if(out.get() && dumpAll)
      *out << "edge_j="<<edge_j<<",("<<e(edge_j,0)<<","<<e(edge_j,1)<<"),"<<e(edge_j,2)<<"\n";
    ++edge_j;
    global_node_id += (local_nx + 1);
  }
  tab.incrTab(-1);
  // Right edges
  if(out.get()) *out << "\nGenerating the right edges ...\n";
  tab.incrTab();
  global_node_id = localNodeOffset + 1 + local_nx;
  for( int i_y = 0; i_y < local_ny; ++i_y ) {
    e(edge_j,0) = global_node_id;
    e(edge_j,1) = global_node_id + (local_nx + 1);
    e(edge_j,2) = bndy_marker;
    if(out.get() && dumpAll)
      *out << "edge_j="<<edge_j<<",("<<e(edge_j,0)<<","<<e(edge_j,1)<<"),"<<e(edge_j,2)<<"\n";
    ++edge_j;
    global_node_id += (local_nx + 1);
  }
  tab.incrTab(-1);
  if(procRank==numProc-1) {
    // Top edges
    if(out.get()) *out << "\nGenerating the top edges ...\n";
    tab.incrTab();
    global_node_id = localNodeOffset+1+(local_nx+1)*local_ny;
    for( int i_x = 0; i_x < local_nx; ++i_x ) {
      e(edge_j,0) = global_node_id;
      e(edge_j,1) = global_node_id + 1;
      e(edge_j,2) = bndy_marker;
      if(out.get() && dumpAll)
        *out << "edge_j="<<edge_j<<",("<<e(edge_j,0)<<","<<e(edge_j,1)<<"),"<<e(edge_j,2)<<"\n";
      ++edge_j;
      global_node_id += 1;
    }
    tab.incrTab(-1);
  }
  TEST_FOR_EXCEPT(edge_j != numedges);
}
