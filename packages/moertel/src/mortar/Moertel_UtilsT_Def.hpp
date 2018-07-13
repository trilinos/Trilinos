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
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include <sstream>
#include <string>

#include "Moertel_UtilsT.hpp"
#include "mrtr_segment.H"
#include "mrtr_segment_linear1D.H"
#include "mrtr_segment_bilineartri.H"
#include "mrtr_segment_bilinearquad.H"
#include "mrtr_node.H"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseSolver.hpp>
#include <TpetraExt_MatrixMatrix_def.hpp>


/*----------------------------------------------------------------------*
  | do 3x3 solve                                              mwgee 10/05|
 *----------------------------------------------------------------------*/
template <class LO, class ST>
bool 
MoertelT::solve33T(const double A[][3], double* x, const double* b)
{
  Teuchos::SerialDenseMatrix<LO, ST> AA(3,3, false);
  Teuchos::SerialDenseMatrix<LO, ST> XX(3,1, false);
  Teuchos::SerialDenseMatrix<LO, ST> BB(3,1, false);
  for (int i=0; i<3; ++i)
  {
    BB(i,0) = b[i];
    for (int j=0; j<3; ++j)
      AA(i,j) = A[i][j];
  }
  Teuchos::SerialDenseSolver<LO, ST> solver;
  solver.setMatrix(AA);
  solver.setVectors(XX,BB);
  solver.factorWithEquilibration(true);
  solver.factor();
  int err = solver.solve();
  if (err)
  {
    std::stringstream oss;
    oss << AA;
    oss << BB;
    oss << XX;
    oss << "***WRN*** MoertelT::solve33T:\n"
      << "***WRN*** Teuchos::SerialDenseSolver::solve returned " << err << "\n"
      << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    throw MOERTEL::ReportError(oss);
  }
  for (int i=0; i<3; ++i)
    x[i] = XX(i,0);

  return true;
}


/*----------------------------------------------------------------------*
  |                                                                 08/05|
  |  modified version of the epetraext matrixmatrixadd                   |
  |  NOTE:                                                               |
  |  - A has to be FillComplete, B must NOT be FillComplete()            |
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
int 
MoertelT::MatrixMatrixAdd(const Tpetra::CrsMatrix<ST, LO, GO, N>& A, bool transposeA, double scalarA,
    Tpetra::CrsMatrix<ST, LO, GO, N>& B, double scalarB )
{
  //
  //This method forms the matrix-matrix sum B = scalarA * op(A) + scalarB * B, where

  Tpetra::MatrixMatrix::Add(A, transposeA, scalarA, B, scalarB);

  return(0);

}


/*----------------------------------------------------------------------*
  | Multiply matrices A*B                                     mwgee 01/06|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> >
MoertelT::MatMatMult(const Tpetra::CrsMatrix<ST, LO, GO, N>& A, bool transA,
    const Tpetra::CrsMatrix<ST, LO, GO, N>& B, bool transB,
    int outlevel)
{


  // create resultmatrix with correct rowmap
  Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> > C;

  if (!transA)
    C = Teuchos::rcp(new Tpetra::CrsMatrix<ST, LO, GO, N>(A.getRangeMap(), 20));
  else
    C = Teuchos::rcp(new Tpetra::CrsMatrix<ST, LO, GO, N>(A.getDomainMap(), 20));

  Tpetra::MatrixMatrix::Multiply(A, transA, B, transB, *C);

  return C;

}

/*----------------------------------------------------------------------*
  | Allocate and return a matrix padded with zeros on the diagonal  06/06|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> > 
MoertelT::PaddedMatrix(const Tpetra::Map<LO, GO, N>& rowmap, double val, const int numentriesperrow)
{
  Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> > tmp 
       = Teuchos::rcp(new Tpetra::CrsMatrix<ST, LO, GO, N>(rowmap, numentriesperrow));
  const int numrows = tmp->getGlobalNumRows();
  for (int i=0; i<numrows; ++i)
  {
    int grid = tmp->getRangeMap()->getGlobalElement(i);
    int err = tmp->insertGlobalValues(grid,1,&val,&grid);
    if (err<0)
    {
      std::cout << "***ERR*** MoertelT::PaddedMatrix:\n"
        << "***ERR*** Cannot insert values into matrix\n"
        << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      tmp = Teuchos::null; // release the memory
      return Teuchos::null;
    }
  }
  return tmp;
}

/*----------------------------------------------------------------------*
  | strip out zeros from a matrix                             m.gee 01/06|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> > 
MoertelT::StripZeros(const Tpetra::CrsMatrix<ST, LO, GO, N>& A, double eps)
{
  Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> > out 
      = Teuchos::rcp(new Tpetra::CrsMatrix<ST, LO, GO, N>(A.getRowMap(),10));
      // Same as Epetra RowMap()
  for (size_t lrow=0; lrow<A.getNodeNumRows(); ++lrow) // Same as Epetra NumMyRows()
  {
    GO grow = A.getRowMap()->getGlobalElement(lrow); // Same as Epetra GRID()
    if (grow<0)
    {
      std::cout << "***ERR*** MoertelT::StripZeros:\n"
        << "***ERR*** Cannot find global row indes from local row index\n"
        << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      out = Teuchos::null; // Free memory
      return Teuchos::null;
    }
    int numentries;
    const int* lindices;
    const double* values;
//    int err  = A.ExtractMyRowView(lrow,numentries,values,lindices);
    LO err  = A.getLocalRowViewRaw(lrow,numentries,lindices,values);
    if (err)
    {
      std::cout << "***ERR*** MoertelT::StripZeros:\n"
        << "***ERR*** A.ExtractMyRowView returned " << err << std::endl
        << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      out = Teuchos::null;
      return Teuchos::null;
    }
    for (int j=0; j<numentries; ++j)
    {
      int lcol = lindices[j];
//      int gcol = A.GCID(lcol);
      GO gcol = A.getColMap()->getGlobalElement(lcol);
      if (gcol<0) {
        std::stringstream oss;
        oss << "ERROR: gcol<0 \n";
        throw MOERTEL::ReportError(oss);
      }
      if (abs(values[j])<eps)
        continue;
      out->insertGlobalValues(grow,1,&values[j],&gcol);
    }
  }
//  out->fillComplete(A.OperatorDomainMap(),A.OperatorRangeMap());
  out->fillComplete(A.getDomainMap(),A.getRangeMap());
  return out;
}

/*----------------------------------------------------------------------*
  | print matrix                                              m.gee 01/06|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::Print_Matrix(std::string name, const Tpetra::CrsMatrix<ST, LO, GO, N>& A, int ibase)
{
  char mypidc[100];
  sprintf(mypidc,"%d",A.Comm().getRank());
  name = name + mypidc + ".mtx";
  const char* nameptr = &name[0];
  FILE* out = fopen(nameptr,"w");
  if (!out)
  {
    std::cout << "***ERR*** MoertelT::Print_Matrix:\n"
      << "***ERR*** Cannot open file " << name << "\n"
      << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  // write global and local dimensions of this operator
  fprintf(out,"%d %d 0\n",A.RangeMap().NumGlobalElements(),A.DomainMap().NumGlobalElements());
  for (int lrow=0; lrow<A.NumMyRows(); ++lrow)
  {
    int grow = A.GRID(lrow);
    if (grow<0)
    {
      std::cout << "***ERR*** MoertelT::Print_Matrix:\n"
        << "***ERR*** Cannot gind global row index from local row index\n"
        << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      fclose(out);
      return false;
    }
    int numentries;
    int* lindices;
    double* values;
    int err  = A.extractMyRowView(lrow,numentries,values,lindices);
    if (err)
    {
      std::cout << "***ERR*** MoertelT::Print_Matrix:\n"
        << "***ERR*** A.ExtractMyRowView returned " << err << std::endl
        << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      fclose(out);
      return false;
    }
    for (int j=0; j<numentries; ++j)
    {
      int lcol = lindices[j];
      int gcol = A.GCID(lcol);
      if (gcol<0) {
        std::stringstream oss;
        oss << "ERROR: gcol<0 \n";
        throw MOERTEL::ReportError(oss);
      }
      fprintf(out," %d   %d   %20.10e\n",grow+ibase,gcol+ibase,values[j]);
    }
  }
  fflush(out);
  fclose(out);
  std::cout << "Tpetra_CrsMatrix is written to file " << name << std::endl;
  fflush(stdout);
  return true;
}


/*----------------------------------------------------------------------*
  | print vector                                              m.gee 02/06|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::Print_Vector(std::string name, const Tpetra::Vector<ST, LO, GO, N>& v, int ibase)
{
  char mypidc[100];
  sprintf(mypidc,"%d",v.Comm().getRank());
  name = name + mypidc + ".vec";
  const char* nameptr = &name[0];
  FILE* out = fopen(nameptr,"w");
  Teuchos::ArrayRCP<const double> vv = v.get1dView();

  if (!out)
  {
    std::cout << "***ERR*** MoertelT::Print_Vector:\n"
      << "***ERR*** Cannot open file " << name << "\n"
      << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  for (int lrow=0; lrow < v->getLocalLength(); ++lrow)
  {
    fprintf(out,"  %20.10e\n",vv[lrow]);
  }
  fflush(out);
  fclose(out);
  std::cout << "Tpetra_Vector is written to file " << name << std::endl;
  fflush(stdout);
  return true;
}


/*----------------------------------------------------------------------*
  | print graph                                               m.gee 04/06|
 *----------------------------------------------------------------------*/
template <class LO,
          class GO,
          class N >
bool 
MoertelT::Print_Graph(std::string name, const Tpetra::CrsGraph<LO, GO, N>& A, int ibase)
{
  char mypidc[100];
  sprintf(mypidc,"%d",A.Comm().getRank());
  name = name + mypidc + ".mtx";
  const char* nameptr = &name[0];
  FILE* out = fopen(nameptr,"w");
  if (!out)
  {
    std::cout << "***ERR*** MoertelT::Print_Graph:\n"
      << "***ERR*** Cannot open file " << name << "\n"
      << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  // write global and local dimensions of this operator
  fprintf(out,"%d %d 0\n",A.RangeMap().NumGlobalElements(),A.DomainMap().NumGlobalElements());
  for (int lrow=0; lrow<A.NumMyRows(); ++lrow)
  {
    int grow = A.GRID(lrow);
    if (grow<0)
    {
      std::cout << "***ERR*** MoertelT::Print_Graph:\n"
        << "***ERR*** Cannot gind global row index from local row index\n"
        << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      fclose(out);
      return false;
    }
    int numentries;
    int* lindices;
    int err  = A.extractMyRowView(lrow,numentries,lindices);
    if (err)
    {
      std::cout << "***ERR*** MoertelT::Print_Graph:\n"
        << "***ERR*** A.ExtractMyRowView returned " << err << std::endl
        << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      fclose(out);
      return false;
    }
    for (int j=0; j<numentries; ++j)
    {
      int lcol = lindices[j];
      int gcol = A.GCID(lcol);
      if (gcol<0) {
        std::stringstream oss;
        oss << "ERROR: gcol<0 \n";
        throw MOERTEL::ReportError(oss);
      }
      fprintf(out," %d   %d   %20.10e\n",grow+ibase,gcol+ibase,1.0);
    }
  }
  fflush(out);
  fclose(out);
  std::cout << "Tpetra_CrsGraph is written to file " << name << std::endl;
  fflush(stdout);
  return true;
}

/*----------------------------------------------------------------------*
  | split matrix into 2x2 block system with given rowmap A22rowmap  06/06|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::SplitMatrix2x2(Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> > A,
    Teuchos::RCP<Tpetra::Map<LO, GO, N> >& A11rowmap,
    Teuchos::RCP<Tpetra::Map<LO, GO, N> >& A22rowmap,
    Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> >& A11,
    Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> >& A12,
    Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> >& A21,
    Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> >& A22)
{
  if (A==Teuchos::null)
  {
    std::stringstream oss;
    oss << "***ERR*** MoertelT::SplitMatrix2x2_A22row_given:\n"
      << "***ERR*** A == null on entry\n"
      << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    throw MOERTEL::ReportError(oss);
  }
  if (A11rowmap==Teuchos::null && A22rowmap != Teuchos::null)
    A11rowmap = Teuchos::rcp(MoertelT::SplitMap(A->RowMap(),*A22rowmap));
  else if (A11rowmap != Teuchos::null && A22rowmap != Teuchos::null);
  else if (A11rowmap != Teuchos::null && A22rowmap == Teuchos::null)
    A22rowmap = Teuchos::rcp(MoertelT::SplitMap(A->RowMap(),*A11rowmap));
  else
  {
    std::stringstream oss;
    oss << "***ERR*** MoertelT::SplitMatrix2x2_A22row_given:\n"
      << "***ERR*** Either A11rowmap OR A22rowmap or both have to be not null"
      << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    throw MOERTEL::ReportError(oss);
  }

  const Teuchos::Comm<LO>& Comm   = A->Comm();
  const Tpetra::Map<LO, GO, N>&  A22map = *(A22rowmap.get());
  const Tpetra::Map<LO, GO, N>&  A11map = *(A11rowmap.get());

  //----------------------------- create a parallel redundant map of A22map
  std::map<int,int> a22gmap;
  {
    std::vector<int> a22global(A22map.NumGlobalElements());
    int count=0;
    for (int proc=0; proc<Comm.getSize(); ++proc)
    {
      int length = 0;
      if (proc==Comm.getRank())
      {
        for (int i=0; i<A22map.NumMyElements(); ++i)
        {
          a22global[count+length] = A22map.GID(i);
          ++length;
        }
      }
      Comm.broadcast(proc, 1, &length);
      Comm.broadcast(proc, length, &a22global[count]);
      count += length;
    }
    if (count != A22map.NumGlobalElements())
    {
      std::stringstream oss;
      oss << "***ERR*** MoertelT::SplitMatrix2x2_A22row_given:\n"
        << "***ERR*** mismatch in dimensions\n"
        << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      throw MOERTEL::ReportError(oss);
    }
    // create the map
    for (int i=0; i<count; ++i)
      a22gmap[a22global[i]] = 1;
    a22global.clear();
  }

  //--------------------------------------------------- create matrix A22
  A22 = Teuchos::rcp(new Tpetra::CrsMatrix<ST, LO, GO, N>(A22map,100));
  {
    std::vector<int>    a22gcindices(100);
    std::vector<double> a22values(100);
    for (int i=0; i<A->NumMyRows(); ++i)
    {
      const int grid = A->GRID(i);
      if (A22map.MyGID(grid)==false)
        continue;
      //std::cout << "Row " << grid << " in A22 Columns ";
      int     numentries;
      double* values;
      int*    cindices;
      int err = A->extractMyRowView(i,numentries,values,cindices);
      if (err)
      {
        std::stringstream oss;
        oss << "***ERR*** MoertelT::SplitMatrix2x2_A22row_given:\n"
          << "***ERR*** A->ExtractMyRowView returned " << err << std::endl
          << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        throw MOERTEL::ReportError(oss);
      }
      if (numentries>(int)a22gcindices.size())
      {
        a22gcindices.resize(numentries);
        a22values.resize(numentries);
      }
      int count=0;
      for (int j=0; j<numentries; ++j)
      {
        const int gcid = A->ColMap().GID(cindices[j]);
        // see whether we have gcid in a22gmap
        std::map<int,int>::iterator curr = a22gmap.find(gcid);
        if (curr==a22gmap.end()) continue;
        //std::cout << gcid << " ";
        a22gcindices[count] = gcid;
        a22values[count]    = values[j];
        ++count;
      }
      //std::cout << std::endl; fflush(stdout);
      // add this filtered row to A22
      err = A22->insertGlobalValues(grid,count,&a22values[0],&a22gcindices[0]);
      if (err<0)
      {
        std::stringstream oss;
        oss << "***ERR*** MoertelT::SplitMatrix2x2_A22row_given:\n"
          << "***ERR*** A22->InsertGlobalValues returned " << err << std::endl
          << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        throw MOERTEL::ReportError(oss);
      }
    } //for (int i=0; i<A->NumMyRows(); ++i)
    a22gcindices.clear();
    a22values.clear();
  }
  A22->FillComplete();
  A22->OptimizeStorage();

  //----------------------------------------------------- create matrix A11
  A11 = Teuchos::rcp(new Tpetra::CrsMatrix<ST, LO, GO, N>(A11map,100));
  {
    std::vector<int>    a11gcindices(100);
    std::vector<double> a11values(100);
    for (int i=0; i<A->NumMyRows(); ++i)
    {
      const int grid = A->GRID(i);
      if (A11map.MyGID(grid)==false) continue;
      int     numentries;
      double* values;
      int*    cindices;
      int err = A->extractMyRowView(i,numentries,values,cindices);
      if (err)
      {
        std::stringstream oss;
        oss << "***ERR*** MoertelT::SplitMatrix2x2_A22row_given:\n"
          << "***ERR*** A->ExtractMyRowView returned " << err << std::endl
          << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        throw MOERTEL::ReportError(oss);
      }
      if (numentries>(int)a11gcindices.size())
      {
        a11gcindices.resize(numentries);
        a11values.resize(numentries);
      }
      int count=0;
      for (int j=0; j<numentries; ++j)
      {
        const int gcid = A->ColMap().GID(cindices[j]);
        // see whether we have gcid as part of a22gmap
        std::map<int,int>::iterator curr = a22gmap.find(gcid);
        if (curr!=a22gmap.end()) continue;
        a11gcindices[count] = gcid;
        a11values[count] = values[j];
        ++count;
      }
      err = A11->insertGlobalValues(grid,count,&a11values[0],&a11gcindices[0]);
      if (err<0)
      {
        std::stringstream oss;
        oss << "***ERR*** MoertelT::SplitMatrix2x2_A22row_given:\n"
          << "***ERR*** A11->InsertGlobalValues returned " << err << std::endl
          << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        throw MOERTEL::ReportError(oss);
      }
    } // for (int i=0; i<A->NumMyRows(); ++i)
    a11gcindices.clear();
    a11values.clear();
  }
  A11->FillComplete();
  A11->OptimizeStorage();

  //---------------------------------------------------- create matrix A12
  A12 = Teuchos::rcp(new Tpetra::CrsMatrix<ST, LO, GO, N>(A11map,100));
  {
    std::vector<int>    a12gcindices(100);
    std::vector<double> a12values(100);
    for (int i=0; i<A->NumMyRows(); ++i)
    {
      const int grid = A->GRID(i);
      if (A11map.MyGID(grid)==false) continue;
      int     numentries;
      double* values;
      int*    cindices;
      int err = A->extractMyRowView(i,numentries,values,cindices);
      if (err)
      {
        std::stringstream oss;
        oss << "***ERR*** MoertelT::SplitMatrix2x2_A22row_given:\n"
          << "***ERR*** A->ExtractMyRowView returned " << err << std::endl
          << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        throw MOERTEL::ReportError(oss);
      }
      if (numentries>(int)a12gcindices.size())
      {
        a12gcindices.resize(numentries);
        a12values.resize(numentries);
      }
      int count=0;
      for (int j=0; j<numentries; ++j)
      {
        const int gcid = A->ColMap().GID(cindices[j]);
        // see whether we have gcid as part of a22gmap
        std::map<int,int>::iterator curr = a22gmap.find(gcid);
        if (curr==a22gmap.end()) continue;
        a12gcindices[count] = gcid;
        a12values[count] = values[j];
        ++count;
      }
      err = A12->insertGlobalValues(grid,count,&a12values[0],&a12gcindices[0]);
      if (err<0)
      {
        std::stringstream oss;
        oss << "***ERR*** MoertelT::SplitMatrix2x2_A22row_given:\n"
          << "***ERR*** A12->InsertGlobalValues returned " << err << std::endl
          << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        throw MOERTEL::ReportError(oss);
      }
    } // for (int i=0; i<A->NumMyRows(); ++i)
    a12values.clear();
    a12gcindices.clear();
  }
  A12->FillComplete(A22map,A11map);
  A12->OptimizeStorage();

  //----------------------------------------------------------- create A21
  A21 = Teuchos::rcp(new Tpetra::CrsMatrix<ST, LO, GO, N>(A22map,100));
  {
    std::vector<int>    a21gcindices(100);
    std::vector<double> a21values(100);
    for (int i=0; i<A->NumMyRows(); ++i)
    {
      const int grid = A->GRID(i);
      if (A22map.MyGID(grid)==false) continue;
      int     numentries;
      double* values;
      int*    cindices;
      int err = A->extractMyRowView(i,numentries,values,cindices);
      if (err)
      {
        std::stringstream oss;
        oss << "***ERR*** MoertelT::SplitMatrix2x2_A22row_given:\n"
          << "***ERR*** A->ExtractMyRowView returned " << err << std::endl
          << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        throw MOERTEL::ReportError(oss);
      }
      if (numentries>(int)a21gcindices.size())
      {
        a21gcindices.resize(numentries);
        a21values.resize(numentries);
      }
      int count=0;
      for (int j=0; j<numentries; ++j)
      {
        const int gcid = A->ColMap().GID(cindices[j]);
        // see whether we have gcid as part of a22gmap
        std::map<int,int>::iterator curr = a22gmap.find(gcid);
        if (curr!=a22gmap.end()) continue;
        a21gcindices[count] = gcid;
        a21values[count] = values[j];
        ++count;
      }
      err = A21->insertGlobalValues(grid,count,&a21values[0],&a21gcindices[0]);
      if (err<0)
      {
        std::stringstream oss;
        oss << "***ERR*** MoertelT::SplitMatrix2x2_A22row_given:\n"
          << "***ERR*** A12->InsertGlobalValues returned " << err << std::endl
          << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        throw MOERTEL::ReportError(oss);
      }
    } // for (int i=0; i<A->NumMyRows(); ++i)
    a21values.clear();
    a21gcindices.clear();
  }
  A21->FillComplete(A11map,A22map);
  A21->OptimizeStorage();

  //-------------------------------------------------------------- tidy up
  a22gmap.clear();
  return true;
}


/*----------------------------------------------------------------------*
  | split a map into 2 pieces with given Agiven                     06/06|
 *----------------------------------------------------------------------*/
template <class LO,
          class GO,
          class N >
Teuchos::RCP<Tpetra::Map<LO, GO, N> > MoertelT::SplitMap(const Tpetra::Map<LO, GO, N>& Amap,
    const Tpetra::Map<LO, GO, N>& Agiven)
{
  const Teuchos::Comm<LO>& Comm = Amap.Comm();
  const Tpetra::Map<LO, GO, N>&  Ag = Agiven;

  int count=0;
  std::vector<int> myaugids(Amap.NumMyElements());
  for (int i=0; i<Amap.NumMyElements(); ++i)
  {
    const int gid = Amap.GID(i);
    if (Ag.MyGID(gid)) continue;
    myaugids[count] = gid;
    ++count;
  }
  myaugids.resize(count);
  int gcount;
  Teuchos::reduceAll(Comm, Teuchos::REDUCE_SUM, 1, &count, &gcount);
  Teuchos::RCP<Tpetra::Map<LO, GO, N> > Aunknown 
        = Teuchos::rcp(new Tpetra::Map<LO, GO, N>(gcount,count,&myaugids[0],0,Comm));
  myaugids.clear();
  return Aunknown;
}

/*----------------------------------------------------------------------*
  | split a vector into 2 pieces with given submaps                 06/06|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::SplitVector(const Tpetra::Vector<ST, LO, GO, N>& x,
    const Tpetra::Map<LO, GO, N>& x1map,
    const Teuchos::RCP<Tpetra::Vector<ST, LO, GO, N> >&   x1,
    const Tpetra::Map<LO, GO, N>& x2map,
    const Teuchos::RCP<Tpetra::Vector<ST, LO, GO, N> >&   x2)
{
  x1 = Teuchos::rcp(new Tpetra::Vector<ST, LO, GO, N>(x1map,false));
  x2 = Teuchos::rcp(new Tpetra::Vector<ST, LO, GO, N>(x2map,false));

  //use an exporter or importer object
  Tpetra::Export<LO, GO, N> exporter_x1(x.Map(),x1map);
  Tpetra::Export<LO, GO, N> exporter_x2(x.Map(),x2map);

  int err = x1->Export(x,exporter_x1,Insert);
  if (err)
  {
    std::stringstream oss;
    oss << "***ERR*** MoertelT::SplitVector:\n"
      << "***ERR*** Export returned " << err << std::endl
      << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    throw MOERTEL::ReportError(oss);
  }

  err = x2->Export(x,exporter_x2,Insert);
  if (err)
  {
    std::stringstream oss;
    oss << "***ERR*** MoertelT::SplitVector:\n"
      << "***ERR*** Export returned " << err << std::endl
      << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    throw MOERTEL::ReportError(oss);
  }

  return true;
}

/*----------------------------------------------------------------------*
  | merge content of 2 vectors into one (assumes matching submaps)  06/06|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::MergeVector(const Tpetra::Vector<ST, LO, GO, N>& x1,
    const Tpetra::Vector<ST, LO, GO, N>& x2,
    Tpetra::Vector<ST, LO, GO, N>& xresult)
{
  //use an exporter or importer object
  Tpetra::Export<LO, GO, N> exporter_x1(x1.Map(),xresult.Map());
  Tpetra::Export<LO, GO, N> exporter_x2(x2.Map(),xresult.Map());

  int err = xresult.Export(x1,exporter_x1,Insert);
  if (err)
  {
    std::stringstream oss;
    oss << "***ERR*** MoertelT::SplitVector:\n"
      << "***ERR*** Export returned " << err << std::endl
      << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    throw MOERTEL::ReportError(oss);
  }

  err = xresult.Export(x2,exporter_x2,Insert);
  if (err)
  {
    std::stringstream oss;
    oss << "***ERR*** MoertelT::SplitVector:\n"
      << "***ERR*** Export returned " << err << std::endl
      << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    throw MOERTEL::ReportError(oss);
  }

  return true;
}
