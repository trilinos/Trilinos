// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_ABSTRACTPROBLEM_H
#define GALERI_ABSTRACTPROBLEM_H

/*!
 * \file Galeri_AbstractProblem.h
 */

class Epetra_RowMatrix;
class Epetra_MultiVector;

namespace Galeri {
namespace FiniteElements {

class Variational;
class Grid;
/*!
 * \class AbstractProblem
 *
 * \brief Abstract interface to define linear problems.
 *
AbstractProblem defines a set of abstract interfaces, used to construct
the linear system corresponding to the finite element discretization 
of a scalar PDE problem. A concrete implementation will require an
AbstractGrid and an AbstractVariational object; the former is used
to query for the grid elements, the latter to integrate the variational
form over such elements. The role of AbstractProblem is to take the
elemental matrices, given by AbstractVariational, and insert them into
the global, distributed matrix (whose RowMatrixRowMap() is given by
Grid().RowMap()).

 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last updated on Apr-05.
 */

class AbstractProblem 
{

public:
      
  //! Destructor.
  virtual ~AbstractProblem() {};

  //! Returns a reference to the linear system matrix.
  virtual Epetra_RowMatrix& A() = 0;

  //! Returns a reference to the multi-vector of right-hand side.
  virtual Epetra_MultiVector& RHS() = 0;

  //! Returns a reference to the multi-vector of starting solution.
  virtual Epetra_MultiVector& LHS() = 0;

  //! Returns a reference to the grid object.
  virtual const AbstractGrid& Grid() const = 0;

  //! Returns a reference to the variational object.
  virtual const AbstractVariational& Variational() const = 0;

  //! Computes the linear system matrix, LHS and RHS.
  virtual void Compute() = 0;

  //! Computes the norm of computed solution, exact solution, and error.
  /*!
   * \param RowMatrixField - (In) Multi-vector defined on Grid().RowMap()
   *                              which contains the numerical solution.
   *
   * \param ExactSolution - (In) Function defined as in the following example:
   *                             ExactSolution(double x, double y, double x,
   *                             double* sol) will contain the value of the
   *                             solution in sol[0], the x-derivative in
   *                             sol[1], the y-derivative in sol[2], and the
   *                             z-derivative in sol[2].
   *
   * \param verbose - (In) If \c true, prints out the results.
   *
   * \param SolutionNorm - (Out) a double array of size 3, which will contain
   *                             the L2 norm, the semi-H1 norm and the H1-norm
   *                             of the numerical solution.
   *
   * \param ExactNorm - (Out) a double array of size 3, which will contain
   *                          the L2 norm, the semi-H1 norm and the H1-norm
   *                          of the exact solution.
   *
   * \param DiffNorm - (Out) a double array of size 3, which will contain
   *                         the L2 norm, the semi-H1 norm and the H1-norm
   *                         of the error.
   */
  virtual void ComputeNorms(Epetra_MultiVector& RowMatrixField,
                            int (*ExactSolution)(double, double, double, double *),
                            const bool verbose = true,
                            double* SolutionNorm = 0,
                            double* ExactNorm = 0,
                            double* DiffNorm = 0) = 0;

};

} // namespace FiniteElements
} // namespace Galeri
#endif
