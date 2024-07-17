// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_PROBLEM_VECTOR_LAPLACIAN_H
#define GALERI_PROBLEM_VECTOR_LAPLACIAN_H

#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"

#include "Galeri_core_Workspace.h"
#include "Galeri_problem_Base.h"
#include "Galeri_quadrature_Segment.h"
#include "Galeri_quadrature_Triangle.h"
#include "Galeri_quadrature_Quad.h"
#include "Galeri_quadrature_Tet.h"
#include "Galeri_quadrature_Hex.h"

#include "Teuchos_RefCountPtr.hpp"

using namespace Teuchos;

namespace Galeri {
namespace problem {
  
template<class T>
class VectorLaplacian // FIXME????? : public Base
{
  public:
  VectorLaplacian(const int numPDEs,
                  const std::string& elementType,
                  const int integrationDegree = Galeri::core::Workspace::MIN, 
                  const int normDegree = Galeri::core::Workspace::MAX) :
    equation_(0),
    numPDEs_(numPDEs)
  {
    if (elementType == "Segment")
    {
      IE_ = rcp(new Galeri::quadrature::Segment(integrationDegree));
      NE_ = rcp(new Galeri::quadrature::Segment(normDegree));
    }
    else if (elementType == "Triangle")
    {
      IE_ = rcp(new Galeri::quadrature::Triangle(integrationDegree));
      NE_ = rcp(new Galeri::quadrature::Triangle(normDegree));
    }
    else if (elementType == "Quad")
    {
      IE_ = rcp(new Galeri::quadrature::Quad(integrationDegree));
      NE_ = rcp(new Galeri::quadrature::Quad(normDegree));
    }
    else if (elementType == "Tet")
    {
      IE_ = rcp(new Galeri::quadrature::Tet(integrationDegree));
      NE_ = rcp(new Galeri::quadrature::Tet(normDegree));
    }
    else if (elementType == "Hex")
    {
      IE_ = rcp(new Galeri::quadrature::Hex(integrationDegree));
      NE_ = rcp(new Galeri::quadrature::Hex(normDegree));
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                         "input elementType not recognized, " << elementType);
  }

  ~VectorLaplacian() {}

  virtual void createGraph(Galeri::grid::Loadable& domain,
                           Epetra_FECrsGraph& A)
  {
    int ierr;

    int numVerticesPerElement = domain.getNumVerticesPerElement();
    int localSize = numVerticesPerElement * numPDEs_;

    Epetra_IntSerialDenseVector vertexList(numVerticesPerElement);

    for (int i = 0; i < domain.getNumMyElements(); ++i)
    {
      int GEID = domain.getGEID(i);
      // load the element vertex IDs
      for (int j = 0; j < numVerticesPerElement; ++j)
        vertexList[j] = domain.getGlobalConnectivity(GEID, j);

      ierr = A.InsertGlobalIndices(numVerticesPerElement, vertexList.Values(),
                                   numVerticesPerElement, vertexList.Values());

      TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::logic_error,
                         "InsertGlobalIndices() returned a negative value, "
                         << ierr);
    }

    // FIXME: :SHOULD BE A PARAMETER...
    A.GlobalAssemble();
  }

  virtual void integrate(Galeri::grid::Loadable& domain,
                         Epetra_RowMatrix& RowA,
                         Epetra_FEVector& RHS)
  {
    int ierr;

    Epetra_FEVbrMatrix& A = dynamic_cast<Epetra_FEVbrMatrix&>(RowA);

    // FIXME: build graph first??
    int numVerticesPerElement = domain.getNumVerticesPerElement();
    int localSize = numVerticesPerElement * numPDEs_;

    Epetra_IntSerialDenseVector vertexList(numVerticesPerElement);
    Epetra_SerialDenseMatrix elementLHS(localSize, localSize);
    Epetra_SerialDenseVector elementRHS(localSize);

    Epetra_SerialDenseMatrix entry(numPDEs_, numPDEs_);

    for (int i = 0; i < domain.getNumMyElements(); ++i)
    {
      int GEID = domain.getGEID(i);
      // load the element vertex IDs
      for (int j = 0; j < numVerticesPerElement; ++j)
        vertexList[j] = domain.getGlobalConnectivity(GEID, j);

      // load the element coordinates
      for (int j = 0; j < numVerticesPerElement; ++j)
        for (int k = 0; k < Galeri::core::Workspace::getNumDimensions(); ++k) 
          (*IE_)(j, k) = domain.getGlobalCoordinates(vertexList[j], k);

      integrateOverElement(*IE_, elementLHS, elementRHS);

      for (int j = 0; j < numVerticesPerElement; ++j)
      {
        ierr = A.BeginSumIntoGlobalValues(vertexList[j], vertexList.Length(),
                                          vertexList.Values());
        TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::logic_error,
                           "SubmitBlockEntry() returned a negative value, "
                           << ierr);

        for (int k = 0; k < numVerticesPerElement; ++k)
        {
          for (int ieq = 0; ieq < numPDEs_; ++ieq)
            for (int jeq = 0; jeq < numPDEs_; ++jeq)
              entry(ieq, jeq) = elementLHS(j * numPDEs_ + ieq, k * numPDEs_ + jeq);

          ierr = A.SubmitBlockEntry(entry.A(), numPDEs_, numPDEs_, numPDEs_);
          TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::logic_error,
                             "SubmitBlockEntry() returned a negative value, "
                             << ierr);
        }

        ierr = A.EndSubmitEntries();
        TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::logic_error,
                           "SubmitBlockEntry() returned a negative value, "
                           << ierr);
      }

      RHS.SumIntoGlobalValues(vertexList, elementRHS);
    }

    A.GlobalAssemble();
    RHS.GlobalAssemble();
  }

  virtual void integrateOverElement(Galeri::quadrature::Element& QE,
                         Epetra_SerialDenseMatrix& ElementLHS, 
                         Epetra_SerialDenseMatrix& ElementRHS)
  {
    ElementLHS.Scale(0.0);
    ElementRHS.Scale(0.0);
#if 0
    for (int i = 0; i < ElementLHS.M(); ++i)
      for (int j = 0; j < ElementLHS.N(); ++j)
        ElementLHS(i, j) = 0.0;

    for (int i = 0; i < ElementRHS.M(); ++i)
      for (int j = 0; j < ElementRHS.N(); ++j)
        ElementRHS(i, j) = 0.0;
#endif

    // cycle over all quadrature nodes

    for (int ii = 0 ; ii < QE.getNumQuadrNodes() ; ii++) 
    {
      double xq, yq, zq;

      QE.computeQuadrNodes(ii, xq, yq, zq);
      QE.computeJacobian(ii);
      QE.computeDerivatives(ii);

      const double& weight = QE.getQuadrWeight(ii);
      const double& det    = QE.getDetJacobian(ii);

      for (int i = 0 ; i < QE.getNumBasisFunctions() ; ++i) 
      {
        const double& phi_i   = QE.getPhi(i);
        const double& phi_x_i = QE.getPhiX(i);
        const double& phi_y_i = QE.getPhiY(i);
        const double& phi_z_i = QE.getPhiZ(i);

        for (int j = 0 ; j < QE.getNumBasisFunctions() ; ++j) 
        {
          const double& phi_j   = QE.getPhi(j);
          const double& phi_x_j = QE.getPhiX(j);
          const double& phi_y_j = QE.getPhiY(j);
          const double& phi_z_j = QE.getPhiZ(j);

          for (int ieq = 0; ieq < numPDEs_; ++ieq)
          {
            for (int jeq = 0; jeq < numPDEs_; ++jeq)
            {
              double contrib = T::getElementLHS(xq, yq, zq, ieq, jeq,
                                                phi_i, phi_x_i, phi_y_i, phi_z_i,
                                                phi_j, phi_x_j, phi_y_j, phi_z_j);
              contrib *= weight * det;

              ElementLHS(i * numPDEs_ + ieq, j * numPDEs_ + jeq) += contrib;
            }

            ElementRHS(i * numPDEs_ + ieq, 0) += 
              weight * det * T::getElementRHS(xq, yq, zq, ieq, phi_i);
          }
        }
      }
    }
  }

  void imposeDirichletBoundaryConditions(Galeri::grid::Loadable& boundary,
                                         Epetra_RowMatrix& RowA,
                                         Epetra_FEVector& RHS,
                                         Epetra_FEVector& LHS)
  {
    Epetra_FEVbrMatrix& A = dynamic_cast<Epetra_FEVbrMatrix&>(RowA);

    TEUCHOS_TEST_FOR_EXCEPTION(A.Filled() == false, std::logic_error,
                       "input matrix must be filled");

    const Epetra_BlockMap& matrixMap = A.RowMatrixRowMap();

    double min = numeric_limits<double>::min();

    Epetra_FEVector rowValues(matrixMap);
    rowValues.PutScalar(min);
    double coord[3];
    coord[0] = 0.0; coord[1] = 0.0; coord[2] = 0.0;

    for (int LEID = 0; LEID < boundary.getNumMyElements(); ++LEID)
    {
      for (int i = 0; i < boundary.getNumVerticesPerElement(); ++i)
      {
        int GVID = boundary.getMyConnectivity(LEID, i);
        for (int j = 0; j < Galeri::core::Workspace::getNumDimensions(); ++j)
          coord[j] = boundary.getGlobalCoordinates(GVID, j);

        for (int ieq = 0; ieq < numPDEs_; ++ieq)
        {
          int j = GVID * numPDEs_ + ieq;
          if (T::getBoundaryType(boundary.getID(), coord[0], coord[1], coord[2], ieq) == 'd')
          {
            double value = T::getBoundaryValue(coord[0], coord[1], coord[2], ieq);
            rowValues.ReplaceGlobalValues(1, &j, &value);
          }
        }
      }
    }

    rowValues.GlobalAssemble(Insert);

    Epetra_Vector colValues(A.RowMatrixColMap());
    colValues.PutScalar(min);
    Epetra_Import importer(A.RowMatrixColMap(), matrixMap);
    colValues.Import(rowValues, importer, Insert); 

    for (int blockRow = 0; blockRow < A.NumMyBlockRows(); ++blockRow)
    {
      int rowDim, numBlockEntries;
      int* blockIndices;
      Epetra_SerialDenseMatrix** values;

      A.ExtractMyBlockRowView(blockRow, rowDim, numBlockEntries, blockIndices,
                              values);

      // cycle over all equations in this block row
      for (int ieq = 0; ieq < numPDEs_; ++ieq)
      {
        // this is the global row ID
        int row = blockRow * numPDEs_ + ieq;

        if (rowValues[0][row] != min)
        {
          for (int j = 0; j < numBlockEntries; ++j)
          {
            for (int jeq = 0; jeq < numPDEs_; ++jeq)
              (*values[j])(ieq, jeq) = 0.0;

            if (blockIndices[j] == blockRow) 
              (*values[j])(ieq, ieq) = 1.0;
          }
          RHS[0][row] = rowValues[0][row]; 
        }
        else
        {
          for (int j = 0; j < numBlockEntries; ++j)
          {
            for (int jeq = 0; jeq < numPDEs_; ++jeq)
            {
              int col = blockIndices[j] * numPDEs_ + jeq;
              if (col == row) continue;
              if (colValues[col] != min) 
              {
                RHS[0][row] -= (*values[j])(ieq, jeq) * colValues[col];
                (*values[j])(ieq, jeq) = 0.0;
              }
            }
          }
        }
      }
    }
  }

  virtual void 
  computeNorms(Galeri::grid::Loadable& domain,
               const Epetra_MultiVector& solution,
               double& solNormL2, double& solSemiNormH1,
               double& exaNormL2, double& exaSemiNormH1,
               double& errNormL2, double& errSemiNormH1)
  {
    // FIXME: should be done only if NumProc > 1.
    const Epetra_Map& vertexMap = domain.getVertexMap();
    Epetra_MultiVector vertexSolution(vertexMap, solution.NumVectors());
    Epetra_Import importer(vertexMap, solution.Map());
    vertexSolution.Import(solution, importer, Insert);

    int numVerticesPerElement = domain.getNumVerticesPerElement();

    Epetra_SerialDenseVector elementSol(numVerticesPerElement);
    Epetra_SerialDenseVector elementNorm(Galeri::core::Workspace::getNumDimensions());

    Epetra_IntSerialDenseVector vertexList(numVerticesPerElement);
    Epetra_SerialDenseMatrix elementLHS(numVerticesPerElement, numVerticesPerElement);
    Epetra_SerialDenseVector elementRHS(numVerticesPerElement);

    exaNormL2 = 0.0, exaSemiNormH1 = 0.0;
    solNormL2 = 0.0, solSemiNormH1 = 0.0;
    errNormL2 = 0.0, errSemiNormH1 = 0.0;

    for (int i = 0; i < domain.getNumMyElements(); ++i)
    {
      for (int j = 0; j < numVerticesPerElement; ++j)
      {
        vertexList[j] = domain.getMyConnectivity(i, j);
        elementSol[j] = vertexSolution[0][vertexMap.LID(vertexList[j])];
      }

      // load the element coordinates
      for (int j = 0; j < numVerticesPerElement; ++j)
        for (int k = 0; k < Galeri::core::Workspace::getNumDimensions(); ++k) 
          (*NE_)(j, k) = domain.getGlobalCoordinates(vertexList[j], k);

      computeNormOverElement((*NE_), elementNorm);
      exaNormL2 += elementNorm[0]; exaSemiNormH1 += elementNorm[1];

      computeErrorOverElement((*NE_), elementSol, elementNorm);
      errNormL2 += elementNorm[0]; errSemiNormH1 += elementNorm[1];

      computeNormOverElement((*NE_), elementSol, elementNorm);
      solNormL2 += elementNorm[0]; solSemiNormH1 += elementNorm[1];
    }

    exaNormL2 = sqrt(exaNormL2); exaSemiNormH1 = sqrt(exaSemiNormH1);
    errNormL2 = sqrt(errNormL2); errSemiNormH1 = sqrt(errSemiNormH1);
    solNormL2 = sqrt(solNormL2); solSemiNormH1 = sqrt(solSemiNormH1);
  }

  virtual void computeNorms(Galeri::grid::Loadable& domain,
                            const Epetra_MultiVector& solution)
  {
    double exaNormL2 = 0.0, exaSemiNormH1 = 0.0;
    double solNormL2 = 0.0, solSemiNormH1 = 0.0;
    double errNormL2 = 0.0, errSemiNormH1 = 0.0;

    computeNorms(domain, solution, solNormL2, solSemiNormH1,
                 exaNormL2, exaSemiNormH1, errNormL2, errSemiNormH1);

    if (solution.Comm().MyPID() == 0)
    {
      cout << endl;
      cout << "||x_h||_2        = " << solNormL2 << endl;
      cout << "||x_ex||_2       = " << exaNormL2 << endl;
      cout << "||x_ex - x_h||_2 = " << errNormL2 << endl;
      cout << endl;
      cout << "|x_h|_1        = " << solSemiNormH1 << endl;
      cout << "|x_ex|_1       = " << exaSemiNormH1 << endl;
      cout << "|x_ex - x_h|_1 = " << errSemiNormH1 << endl;
      cout << endl;
    }
  }

  virtual void computeNormOverElement(Galeri::quadrature::Element& QE,
                                      Epetra_SerialDenseMatrix& elementSol,
                                      Epetra_SerialDenseMatrix& elementNorm)
  {
    elementNorm(0, 0) = 0.0;
    elementNorm(1, 0) = 0.0;

    for (int ii = 0 ; ii < QE.getNumQuadrNodes() ; ii++) 
    {
      double xq, yq, zq;

      QE.computeQuadrNodes(ii, xq, yq, zq);
      QE.computeJacobian(ii);
      QE.computeDerivatives(ii);

      const double& weight = QE.getQuadrWeight(ii);
      const double& det    = QE.getDetJacobian(ii);

      double sol      = 0.0, sol_derx = 0.0;
      double sol_dery = 0.0, sol_derz = 0.0;

      for (int k = 0 ; k < QE.getNumBasisFunctions() ; ++k)
      {
        sol      += QE.getPhi(k)  * elementSol(k, 0);
        sol_derx += QE.getPhiX(k) * elementSol(k, 0);
        sol_dery += QE.getPhiY(k) * elementSol(k, 0);
        sol_derz += QE.getPhiZ(k) * elementSol(k, 0);
      }

      elementNorm(0, 0) += weight * det * sol * sol;
      elementNorm(1, 0) += weight * det * (sol_derx * sol_derx +
                                           sol_dery * sol_dery +
                                           sol_derz * sol_derz);
    }
  }

  virtual void computeErrorOverElement(Galeri::quadrature::Element& QE,
                                       Epetra_SerialDenseMatrix& elementSol,
                                       Epetra_SerialDenseMatrix& elementNorm)
  {
    elementNorm(0, 0) = 0.0;
    elementNorm(1, 0) = 0.0;

    for (int ii = 0; ii < QE.getNumQuadrNodes(); ii++) 
    {
      double xq, yq, zq;

      QE.computeQuadrNodes(ii, xq, yq, zq);
      QE.computeJacobian(ii);
      QE.computeDerivatives(ii);

      const double& weight = QE.getQuadrWeight(ii);
      const double& det    = QE.getDetJacobian(ii);

      double sol      = 0.0, sol_derx = 0.0;
      double sol_dery = 0.0, sol_derz = 0.0;

      for (int k = 0; k < QE.getNumBasisFunctions(); ++k)
      {
        sol      += QE.getPhi(k)  * elementSol(k, 0);
        sol_derx += QE.getPhiX(k) * elementSol(k, 0);
        sol_dery += QE.getPhiY(k) * elementSol(k, 0);
        sol_derz += QE.getPhiZ(k) * elementSol(k, 0);
      }

      sol      -= T::getExactSolution('f', xq, yq, zq, equation_);
      sol_derx -= T::getExactSolution('x', xq, yq, zq, equation_);
      sol_dery -= T::getExactSolution('y', xq, yq, zq, equation_);
      sol_derz -= T::getExactSolution('z', xq, yq, zq, equation_);

      elementNorm(0, 0) += weight * det * sol * sol;
      elementNorm(1, 0) += weight * det * (sol_derx * sol_derx +
                                           sol_dery * sol_dery +
                                           sol_derz * sol_derz);
    }
  }

  virtual void 
  computeNormOverElement(Galeri::quadrature::Element& QE,
                         Epetra_SerialDenseMatrix& elementNorm)
  {
    elementNorm(0, 0) = 0.0;
    elementNorm(1, 0) = 0.0;

    for (int ii = 0; ii < QE.getNumQuadrNodes(); ii++) 
    {
      double xq, yq, zq;

      QE.computeQuadrNodes(ii, xq, yq, zq);
      QE.computeJacobian(ii);
      QE.computeDerivatives(ii);

      const double& weight = QE.getQuadrWeight(ii);
      const double& det    = QE.getDetJacobian(ii);

      double sol      = T::getExactSolution('f', xq, yq, zq, equation_);
      double sol_derx = T::getExactSolution('x', xq, yq, zq, equation_);
      double sol_dery = T::getExactSolution('y', xq, yq, zq, equation_);
      double sol_derz = T::getExactSolution('z', xq, yq, zq, equation_);

      elementNorm(0, 0) += weight * det * sol * sol;
      elementNorm(1, 0) += weight * det * (sol_derx * sol_derx +
                                           sol_dery * sol_dery +
                                           sol_derz * sol_derz);
    }
  }

  void setEquation(const int& equation)
  {
    equation_ = equation;
  }

private:
  int equation_;
  int numPDEs_;
  RefCountPtr<Galeri::quadrature::Element> IE_, NE_;
}; // class Base

} // namespace problem
} // namespace Galeri
#endif

