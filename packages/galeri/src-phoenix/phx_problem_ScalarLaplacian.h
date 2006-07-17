// @HEADER
// ************************************************************************
//
//                  Galeri Matrix Generation Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef PHX_PROBLEM_SCALAR_LAPLACIAN_H
#define PHX_PROBLEM_SCALAR_LAPLACIAN_H

#include "Epetra_Import.h"

#include "phx_core_Constants.h"
#include "phx_problem_Base.h"
#include "phx_quadrature_Segment.h"
#include "phx_quadrature_Triangle.h"
#include "phx_quadrature_Quad.h"
#include "phx_quadrature_Tet.h"
#include "phx_quadrature_Hex.h"

#include "Teuchos_RefCountPtr.hpp"

using namespace Teuchos;

namespace phx {
namespace problem {
  
template<class T>
class ScalarLaplacian : public Base
{
  public:
  ScalarLaplacian(const string& elementType,
                  const int integrationDegree = phx::core::Constants::MIN, 
                  const int normDegree = phx::core::Constants::MAX) 
  {
    if (elementType == "Segment")
    {
      IE_ = rcp(new phx::quadrature::Segment(integrationDegree));
      NE_ = rcp(new phx::quadrature::Segment(normDegree));
    }
    else if (elementType == "Triangle")
    {
      IE_ = rcp(new phx::quadrature::Triangle(integrationDegree));
      NE_ = rcp(new phx::quadrature::Triangle(normDegree));
    }
    else if (elementType == "Quad")
    {
      IE_ = rcp(new phx::quadrature::Quad(integrationDegree));
      NE_ = rcp(new phx::quadrature::Quad(normDegree));
    }
    else if (elementType == "Tet")
    {
      IE_ = rcp(new phx::quadrature::Tet(integrationDegree));
      NE_ = rcp(new phx::quadrature::Tet(normDegree));
    }
    else if (elementType == "Hex")
    {
      IE_ = rcp(new phx::quadrature::Hex(integrationDegree));
      NE_ = rcp(new phx::quadrature::Hex(normDegree));
    }
    else
      TEST_FOR_EXCEPTION(true, std::logic_error,
                         "input elementType not recognized, " << elementType);
  }

  ~ScalarLaplacian() {}

  virtual void integrate(phx::grid::Loadable& domain,
                         Epetra_FECrsMatrix& A,
                         Epetra_FEVector& RHS)
  {
    // FIXME: BUILD GRAPH FIRST???
    int numVerticesPerElement = domain.getNumVerticesPerElement();
    Epetra_IntSerialDenseVector vertexList(numVerticesPerElement);
    Epetra_SerialDenseMatrix elementLHS(numVerticesPerElement, numVerticesPerElement);
    Epetra_SerialDenseVector elementRHS(numVerticesPerElement);

    for (int i = 0; i < domain.getNumMyElements(); ++i)
    {
      int GEID = domain.getGEID(i);
      // load the element vertex IDs
      for (int j = 0; j < numVerticesPerElement; ++j)
        vertexList[j] = domain.getGlobalConnectivity(GEID, j);

      // load the element coordinates
      for (int j = 0; j < numVerticesPerElement; ++j)
        for (int k = 0; k < phx::core::Utils::getNumDimensions(); ++k) 
          (*IE_)(j, k) = domain.getGlobalCoordinates(vertexList[j], k);

      integrateOverElement(*IE_, elementLHS, elementRHS);

      A.InsertGlobalValues(vertexList, elementLHS);
      RHS.SumIntoGlobalValues(vertexList, elementRHS);
    }

    A.GlobalAssemble();
    RHS.GlobalAssemble();
  }

  void imposeDirichletBoundaryConditions(phx::grid::Loadable& boundary,
                                         Epetra_FECrsMatrix& A,
                                         Epetra_FEVector& RHS,
                                         Epetra_FEVector& LHS)
  {
    TEST_FOR_EXCEPTION(A.Filled() == false, std::logic_error,
                       "input matrix must be filled");

    const Epetra_Map& matrixMap = A.RowMatrixRowMap();

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
        for (int j = 0; j < phx::core::Utils::getNumDimensions(); ++j)
          coord[j] = boundary.getGlobalCoordinates(GVID, j);

        if (T::getBoundaryType(boundary.getID(), coord[0], coord[1], coord[2]) == 'd')
        {
          double value = T::getBoundaryValue(coord[0], coord[1], coord[2]);
          rowValues.ReplaceGlobalValues(1, &GVID, &value);
        }
      }
    }

    rowValues.GlobalAssemble(Insert);

    Epetra_Vector colValues(A.RowMatrixColMap());
    colValues.PutScalar(min);
    Epetra_Import importer(A.RowMatrixColMap(), matrixMap);
    colValues.Import(rowValues, importer, Insert); 

    for (int i = 0; i < A.NumMyRows(); ++i)
    {
      int* indices;
      double* values;
      int numEntries;
      A.ExtractMyRowView(i, numEntries, values, indices);

      if (rowValues[0][i] != min)
      {
        for (int j = 0; j < numEntries; ++j)
        {
          if (indices[j] != i) values[j] = 0.0;
          else values[j] = 1.0;
        }
        RHS[0][i] = rowValues[0][i];  // FIXME???
      }
      else
      {
        for (int j = 0; j < numEntries; ++j)
        {
          if (indices[j] == i) continue;
          if (colValues[indices[j]] != min) 
          {
            RHS[0][i] -= values[j] * colValues[indices[j]];
            values[j] = 0.0;
          }
        }
      }
    }
  }

  virtual void 
  computeNorms(phx::grid::Loadable& domain,
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
    Epetra_SerialDenseVector elementNorm(phx::core::Utils::getNumDimensions());

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
        for (int k = 0; k < phx::core::Utils::getNumDimensions(); ++k) 
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

  virtual void computeNorms(phx::grid::Loadable& domain,
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

  virtual void integrateOverElement(phx::quadrature::Element& QE,
                         Epetra_SerialDenseMatrix& ElementLHS, 
                         Epetra_SerialDenseMatrix& ElementRHS)
  {
    for (int i = 0; i < ElementLHS.M(); ++i)
      for (int j = 0; j < ElementLHS.N(); ++j)
        ElementLHS(i, j) = 0.0;

    for (int i = 0; i < ElementRHS.M(); ++i)
      for (int j = 0; j < ElementRHS.N(); ++j)
        ElementRHS(i, j) = 0.0;

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

          double contrib = T::getElementLHS(xq, yq, zq,
                                            phi_i, phi_x_i, phi_y_i, phi_z_i,
                                            phi_j, phi_x_j, phi_y_j, phi_z_j);
          contrib *= weight * det;

          ElementLHS(i, j) += contrib;
        }

        ElementRHS(i, 0) += weight * det * T::getElementRHS(xq, yq, zq, phi_i);
      }
    }
  }

  virtual void computeNormOverElement(phx::quadrature::Element& QE,
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

  virtual void computeErrorOverElement(phx::quadrature::Element& QE,
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

      sol      -= T::getExactSolution('f', xq, yq, zq);
      sol_derx -= T::getExactSolution('x', xq, yq, zq);
      sol_dery -= T::getExactSolution('y', xq, yq, zq);
      sol_derz -= T::getExactSolution('z', xq, yq, zq);

      elementNorm(0, 0) += weight * det * sol * sol;
      elementNorm(1, 0) += weight * det * (sol_derx * sol_derx +
                                           sol_dery * sol_dery +
                                           sol_derz * sol_derz);
    }
  }

  virtual void 
  computeNormOverElement(phx::quadrature::Element& QE,
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

      double sol      = T::getExactSolution('f', xq, yq, zq);
      double sol_derx = T::getExactSolution('x', xq, yq, zq);
      double sol_dery = T::getExactSolution('y', xq, yq, zq);
      double sol_derz = T::getExactSolution('z', xq, yq, zq);

      elementNorm(0, 0) += weight * det * sol * sol;
      elementNorm(1, 0) += weight * det * (sol_derx * sol_derx +
                                           sol_dery * sol_dery +
                                           sol_derz * sol_derz);
    }
  }

private:
  RefCountPtr<phx::quadrature::Element> IE_, NE_;
}; // class Base

} // namespace problem
} // namespace phx
#endif

