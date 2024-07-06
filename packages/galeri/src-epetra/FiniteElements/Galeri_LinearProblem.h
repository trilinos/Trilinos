// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_SCALARPROBLEM_H
#define GALERI_SCALARPROBLEM_H

/*!
 * \file Galeri_LinearProblem.h
 */

#include "Epetra_Vector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Export.h"
#include "Galeri_Utils.h"
#include "Galeri_AbstractGrid.h"
#include "Galeri_AbstractVariational.h"
#include "Galeri_AbstractProblem.h"
#include <limits>

namespace Galeri {
namespace FiniteElements {

/*!
 * \class LinearProblem
 *
 * \brief Basic implementation of scalar finite element problem.
 *
 * This class fill the linea system matrix (defined as an Epetra_CrsMatrix),
 * the right-hand side (defined as an Epetra_MultiVector) and the
 * starting solution (defined as a zero Epetra_MultiVector). In the current
 * implementation, only one rhs is created.
 *
 * \note Neumann boundary conditions are still to be fixed.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last updated on Apr-05.
 */

class LinearProblem : public AbstractProblem
{

public:
      
  //! Constructor.
  /*!
   * \param Grid - (In) Reference to an AbstractGrid object
   *
   * \param Variational - (In) Reference to an AbstractVariational object
   *
   * \param A - (In/Out) Epetra_CrsMatrix, whose Map is Grid().RowMap(),
   *                     that will contain the linear system matrix.
   *
   * \param LHS - (In/Out) Epetra_MultiVector, whose Map is Grid().RowMap(),
   *                       that will contain the starting solution 
   *                       (zero vector).
   * 
   * \param RHS - (In/Out) Epetra_MultiVector, whose Map is Grid().RowMap(),
   *                       that will contain the right-hand side.
   */
  LinearProblem(const AbstractGrid& Grid,
                const AbstractVariational& Variational,
                Epetra_CrsMatrix& A, Epetra_MultiVector& LHS, 
                Epetra_MultiVector& RHS) :
    Grid_(Grid),
    Variational_(Variational),
    A_(A),
    LHS_(LHS),
    RHS_(RHS)
  {}

  //! Destructor.  
  virtual ~LinearProblem() {};

  //! Fills the linear system matrix and the right-hand side, zeros out the solution
  void Compute()
  {
    LHS().PutScalar(0.0);
    RHS().PutScalar(0.0);

    const Epetra_Map& VertexMap = Grid().VertexMap();

    Epetra_CrsMatrix LocalA(Copy, VertexMap, 0);
    Epetra_Vector LocalRHS(VertexMap);

    // get maximum number of unknowns per element
    int size = Grid().NumVerticesPerElement();

    // allocate elemental matrices and RHS
    std::vector<double> ElementMatrix(size * size);
    std::vector<double> ElementRHS(size);

    std::vector<double> x(size);
    std::vector<double> y(size);
    std::vector<double> z(size);
    std::vector<int>    LVID(size);
    std::vector<int>    GVID(size);

    // ==================== //
    // Fill matrix elements //
    // ==================== //

    for (int ie = 0 ; ie < Grid().NumMyElements() ; ++ie) 
    {
      double h = Grid().ElementMaxLength(ie);

      Grid().ElementVertices(ie, &LVID[0]);
      Grid().VertexCoord(size, &LVID[0], &x[0], &y[0], &z[0]);

      // form elemental matrix and rhs for element `ie'
      Variational().IntegrateOverElement(Variational(), &x[0], &y[0], &z[0], &h,
                                         &ElementMatrix[0], &ElementRHS[0]);

      for (int i = 0 ; i < size ; ++i)
      {
        long long LLrow = VertexMap.GID64(LVID[i]);
        if(LLrow > std::numeric_limits<int>::max())
        {
          cerr << "LLrow out of int bound" << endl;
          cerr << "File " << __FILE__ << ", line " << __LINE__ << endl;
          throw(-1);
        }

        int row = (int) LLrow;
        assert (row != -1);
        for (int j = 0 ; j < size ; ++j)
        {
          long long LLcol = VertexMap.GID64(LVID[j]);
          if(LLcol > std::numeric_limits<int>::max())
          {
            cerr << "LLcol out of int bound" << endl;
            cerr << "File " << __FILE__ << ", line " << __LINE__ << endl;
            throw(-1);
          }

          int col = (int) LLcol;

          double mat_value = ElementMatrix[i + j * size];
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
          if (LocalA.SumIntoGlobalValues(row, 1, &mat_value, &LLcol) > 0)
            LocalA.InsertGlobalValues(row, 1, &mat_value, &LLcol);
#else
          if (LocalA.SumIntoGlobalValues(row, 1, &mat_value, &col) > 0)
            LocalA.InsertGlobalValues(row, 1, &mat_value, &col);
#endif
        }
        LocalRHS[LVID[i]] += ElementRHS[i];
      }
    }

    LocalA.FillComplete(VertexMap, VertexMap);

    // ========================== //
    // impose boundary conditions //
    // ========================== //

    size = Grid().NumVerticesPerFace();
    LVID.resize(size);
    GVID.resize(size);

    for (int face = 0 ; face < Grid().NumMyBoundaryFaces() ; ++face) 
    {
      int Patch = Grid().FacePatch(face);
      int BCType = Variational().BC(Patch);

      Grid().FaceVertices(face, Patch, &LVID[0]);
      Grid().VertexCoord(size, &LVID[0], &x[0], &y[0], &z[0]);

      if (BCType == GALERI_DIRICHLET)
      {
        for (int j = 0 ; j< size ;  ++j) 
        {
          int MyRow = LVID[j];
          double value = Variational().BC(x[j], y[j], z[j], Patch);
          LocalRHS[MyRow] = value;
          int NumEntries = 0;
          int* Indices;
          double* Values;
          LocalA.ExtractMyRowView(MyRow, NumEntries, Values, Indices);

          for (int i = 0 ; i < NumEntries ; ++i)
            if (Indices[i] == MyRow) {
              Values[i] = 1.0;
            }
            else
              Values[i] = 0.0;
        }
      }
      else if (BCType == GALERI_NEUMANN)
      {
#if 0
        double Area = 0.0; // Grid.FaceArea(face);
        for (int j = 0 ; j< size ;  ++j) 
        {
          double value = Variational().BC(x[j], y[j], z[j], Patch);
          value /= size;
          value *= Area;
          RHS(Vertices[j]) = value;
        }
#endif
        cerr << "Still to check..." << endl;
        throw(-1);

      }
      else if (BCType == GALERI_DO_NOTHING)
      {
        // do nothing here..
      }
      else
      {
        cerr << "Type of boundary condition not recognized" << endl;
        throw(-1);
      }
    }

    const Epetra_Map& RowMap = Grid().RowMap();
    Epetra_Export Exporter(VertexMap, RowMap);

    CrsA().Export(LocalA, Exporter, Add);
    CrsA().FillComplete();

    RHS().Export(LocalRHS, Exporter, Add);
  }

  //! Computes L2, semi-H1 and H1 norms.
  void ComputeNorms(Epetra_MultiVector& RowMatrixField,
                    int (*ExactSolution)(double, double, double, double *),
                    const bool verbose = true,
                    double* Solution = 0, double* Exact = 0, double* Diff = 0)
  {

    const Epetra_Map& VertexMap = Grid().VertexMap(); 
    //const Epetra_Map& RowMap = Grid().RowMap(); 
    Epetra_MultiVector VertexField(VertexMap, RowMatrixField.NumVectors());

    Grid().ExportToVertexMap(RowMatrixField, VertexField);

    double NormSol[3],    NormSolGlobal[3];
    double NormExact[3],  NormExactGlobal[3];
    double NormDiff[3],   NormDiffGlobal[3];

    for (int i = 0 ; i < 3 ; ++i) 
    {
      NormSol[i]   = 0.0; NormSolGlobal[i]   = 0.0;
      NormExact[i] = 0.0; NormExactGlobal[i] = 0.0;
      NormDiff[i]  = 0.0; NormDiffGlobal[i]  = 0.0;
    }

    int size = Grid().NumVerticesPerElement();
    std::vector<double> x(size);
    std::vector<double> y(size);
    std::vector<double> z(size);
    std::vector<double> LocalSol(size);
    std::vector<int>    Vertices(size);

    for (int i = 0 ; i < size ; ++i) 
    {
      x[i] = 0.0;
      y[i] = 0.0; 
      z[i] = 0.0;
    }

    ///double xq, yq, zq;

    // =========================== //
    // roll over all grid elements //
    // =========================== //

    for (int ie = 0 ; ie < Grid().NumMyElements() ; ++ie) 
    {
      Grid().ElementVertices(ie, &Vertices[0]);
      Grid().VertexCoord(size, &Vertices[0], &x[0], &y[0], &z[0]);
      for (int i = 0 ; i < size ; ++i)
        LocalSol[i] = VertexField[0][Vertices[i]];

      Variational().ElementNorm(&LocalSol[0], &x[0], &y[0], &z[0], NormSol);
      Variational().ElementNorm(ExactSolution, &x[0], &y[0], &z[0], &NormExact[0]);
      Variational().ElementNorm(&LocalSol[0] ,ExactSolution,
                             &x[0], &y[0], &z[0], NormDiff);
    }

    Grid().Comm().SumAll((double*)&NormSol,  NormSolGlobal, 3);
    Grid().Comm().SumAll((double*)&NormExact,NormExactGlobal, 3);
    Grid().Comm().SumAll((double*)&NormDiff, NormDiffGlobal, 3);

    NormSolGlobal[2]   = NormSolGlobal[0]   + NormSolGlobal[1];
    NormExactGlobal[2] = NormExactGlobal[0] + NormExactGlobal[1];
    NormDiffGlobal[2]  = NormDiffGlobal[0]  + NormDiffGlobal[1];

    for (int i = 0 ; i < 3 ; ++i) 
    {
      NormSolGlobal[i]   = sqrt(NormSolGlobal[i]);
      NormExactGlobal[i] = sqrt(NormExactGlobal[i]);
      NormDiffGlobal[i]  = sqrt(NormDiffGlobal[i]);
    }

    if (verbose && Grid().Comm().MyPID() == 0) 
    {
      cout << "|| vector ||_L2              = " << NormSolGlobal[0] << endl;
      cout << "|| vector ||_semi-H1         = " << NormSolGlobal[1] << endl;
      cout << "|| vector ||_H1              = " << NormSolGlobal[2] << endl;
      cout << endl;
      cout << "|| exact solution ||_L2      = " << NormExactGlobal[0] << endl;
      cout << "|| exact solution ||_semi-H1 = " << NormExactGlobal[1] << endl;
      cout << "|| exact solution ||_H1      = " << NormExactGlobal[2] << endl;
      cout << endl;
      cout << "|| error ||_L2               = " << NormDiffGlobal[0] << endl;
      cout << "|| error ||_semi-H1          = " << NormDiffGlobal[1] << endl;
      cout << "|| error ||_H1               = " << NormDiffGlobal[2] << endl;
      cout << endl;
    }

    if (Solution)
      for (int i = 0 ; i < 3 ; ++i)
        Solution[i] = NormSolGlobal[i];

    if (Exact)
      for (int i = 0 ; i < 3 ; ++i)
        Exact[i] = NormExactGlobal[i];

    if (Diff)
      for (int i = 0 ; i < 3 ; ++i)
        Diff[i] = NormDiffGlobal[i];
  }

  //! Returns a reference to the linear system matrix.
  virtual Epetra_RowMatrix& A()
  {
    return(A_);
  }

  //! Returns a reference to the linear system matrix as Epetra_CrsMatrix.
  virtual Epetra_CrsMatrix& CrsA()
  {
    return(A_);
  }

  virtual Epetra_MultiVector& RHS()
  {
    return(RHS_);
  }

  virtual Epetra_MultiVector& LHS()
  {
    return(LHS_);
  }

  virtual const AbstractGrid& Grid() const
  {
    return(Grid_);
  }

  virtual const AbstractVariational& Variational() const
  {
    return(Variational_);
  }

public:

  const AbstractGrid& Grid_;
  const AbstractVariational& Variational_;
  Epetra_CrsMatrix& A_;
  Epetra_MultiVector& LHS_;
  Epetra_MultiVector& RHS_;
};

} // namespace FiniteElements
} // namespace Galeri
#endif
