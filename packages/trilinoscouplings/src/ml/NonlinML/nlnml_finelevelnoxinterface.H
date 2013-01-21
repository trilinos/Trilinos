/*
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file nlnml_finelevelnoxinterface.H
 *
 * \class NLNML_finelevelnoxinterface
 *
 * \brief virtual class NLNML_finelevelnoxinterface used by the nonlinear preconditioner
 *
 * \date Last update do Doxygen: 31-Mar-06
 *
 */
// ML headers
#include "ml_common.h"
#include "ml_include.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)
#ifndef NLNML_FINELEVELNOXINTERFACE_H
#define NLNML_FINELEVELNOXINTERFACE_H

#include <iostream>
// Epetra headers
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
// NOX interface
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Required.H"

//! NLNML: default namespace for all ML/NOX nonlinear methods
namespace NLNML
{

/*!
 
   \brief Ml_Nox_Fineinterface: a virtual class used by NLNML_Preconditioner

   NLNML::NLNML_FineLevelNoxInterface is a pure virtual class that serves as
   an interface between the NLNML_Preconditioner and any application.
   It implements NOX interfaces NOX::Epetra::Interface::Required and
   NOX::EpetraNew::Interface::Jacobian.
   Besides serving as an interface for the nonlinear iteration, it also serves as an 
   interface to do finite differencing of Jacobians and as the fine level interface all
   coarser level interfaces refer to. It therefore plays a major role in using
   the NLNML_Preconditioner.
   
    To make use of the full functionality of this class it 
    requires ML to be configured with the following options:
    - \c --with-ml_nox
    - \c --enable-epetra
    - \c --enable-epetraext
    - \c --enable-nox
    - \c --enable-nox-epetra
    - \c --enable-prerelease
    - \c --enable-aztecoo
    - \c --enable-amesos
    
    \author Michael Gee, SNL 1414
*/  
class  NLNML_FineLevelNoxInterface : public virtual NOX::Epetra::Interface::Required,
                                     public virtual NOX::Epetra::Interface::Jacobian
{
public:

  //! Constructor
  NLNML_FineLevelNoxInterface()
  { isnewJacobian_=false; numJacobian_=0; printlevel_=0; t_ = 0.; ncalls_computeF_=0; return;}

  //! virtual Destructor
  virtual ~NLNML_FineLevelNoxInterface() {}
  
  //! Compute and return F (derived from NOX::Epetra::Interface::Required)
  /*! Compute and return F (derived from NOX::Epetra::Interface::Required). This method is called
  by the nonlinear preconditioner. It supplies a solution vector x and needs a residual vector F computed
  by this funtion on output.
   \param x        (In)  : current solution vector
   \param F        (Out) : residual vector matching x
   \param fillFlag (In)  : flag indicating for what purpose this routine was called, see NOX documentation of
                           NOX::Epetra::Interface::Required
   */
  virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, 
			const FillType fillFlag) = 0;
  
  //! Compute the current Jacobian (derived from NOX::Epetra::Interface::Jacobian)
  /*! Compute the current Jacobian (derived from NOX::Epetra::Interface::Jacobian). 
      This method is called by the nonlinear preconditioner. It supplies a solution vector x and wants
      a Jacobian computed matching this solution vector on output. 
      getJacobian will be used to obtain the Jacobian
      \sa getJacobian
   \param x        (In)  : current solution vector
   \param Jac      (In)  : the Jacobian to be computed
   */
  virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) = 0;
  
  //! Return the Jacobian from the application.
  /*! Return the Jacobian from the application. The NLNML_Preconditioner uses this
      method to get the jacobian from the application. It will never be used in the matrixfree case.
      \sa NLNML_Preconditioner
   */
  virtual Epetra_CrsMatrix* getJacobian()  = 0;
  
  //! Return the Graph of the fine level problem from the application.
  /*! Return the Graph of the fine level problem from the application. 
      The NLNML_Preconditioner uses this method to get the graph of
      the fine level problem from the application to do finite differencing
   */
  virtual const Epetra_CrsGraph* getGraph() = 0;

  //! Return the modified Graph of the fine level problem from the application.
  /*! Return a graph that has sufficient modifications to reflect 
      constraints that will be applied to the solution and residual.
      See Sandia report for more details
   */
  virtual const Epetra_CrsGraph* getModifiedGraph() = 0;
  
  //! Return the current solution from the application.
  /*! Return the current solution from the application.
      The NLNML_Preconditioner uses this method to get the current solution vector of
      the fine level problem from the application.
   */
  virtual const Epetra_Vector* getSolution() = 0;
  
  //! Return the Map of the current solution from the application.
  /*! Return the Map of the current solution from the application.
      The NLNML_Preconditioner uses this method to get the map of the current solution vector of
      the fine level problem from the application.
   */
  virtual const Epetra_Map& getMap() = 0;

  //! Compute and return the Nullspace of the problem from the application.
  /*! Compute and return the Nullspace of the problem from the application. This method is called
  by the nonlinear preconditioner to obtain the Nullspace of the problem from the application.
  If the application returns NULL on exit, the NLNML_Preconditioner will use a default Nullspace of
  correct size (supplied in parameter list) and will assume constant number of dofs per node (constant block size).
  This will lead to errors in the aggregation process if the problem is NOT of constant block size!
   \param nummyrows     (In)  : dimension of the nullspace vectors
   \param numpde        (In) : number of PDEs (normally dofs per node)
   \param dim_nullsp    (In)  : dimension of the nullspace
   */
  virtual double* Get_Nullspace(const int nummyrows, const int numpde, 
                                const int dim_nullsp) = 0;

  // get block information
  //! Compute and return the block node information from the application.
  /*! Compute and return the block node information from the application. This method is called
  by the nonlinear preconditioner to obtain information about non-constant block sizes. It is used
  only in the case of the VBMETIS aggregation scheme. If on exit blocks and block_pde are NULL, the
  NLNML_Preconditioner will assume a constant vblock size and will use METIS instead of the VBMETIS
  aggregation scheme.
   \param nblocks     (Out)  : on exit, number of variable sized nodal blocks
   \param blocks      (Out)  : on exit, allocated vector of length NumMyRows containing the block indize
                               of each row in C-style numbering (starting with 0)
   \param block_pde   (Out)  : on exit, allocated vector of length NumMyRows containing the number of
                               PDE a certain row belongs to.
   */
  virtual bool getBlockInfo(int *nblocks, std::vector<int>& blocks, std::vector<int>& block_pde) = 0;

  //! Dummy routine to apply constraints
  /*! If the underlying application needs to apply constraints to the
      gradient computed by the preconditioner this method has to be implemented.
   */
  virtual void ApplyAllConstraints(Epetra_Vector& gradient, int level) 
  { return;}

  //! Query whether the current Jacobian matches the current solution vector.
  virtual bool isnewJacobian() { return isnewJacobian_; }

  //! Query how often the Jacobian has been recomputed.
  virtual int getnumJacobian() { return numJacobian_; }

  //! Query the summed time spent in this interface.
  double getsumtime() {return t_;}

  //! Reset the summed time spent in this interface
  void resetsumtime() { t_ = 0.; return; }

  //! Query the number of calls to the computeF method of this class.
  int getnumcallscomputeF() { return ncalls_computeF_; }

  //! Reset the number of calls to the computeF method of this class to a certain value
  bool setnumcallscomputeF(int ncalls) { ncalls_computeF_=ncalls; return true; }

  //! Returns the level of output to be generated [ 0 - 10 ]
  inline int OutLevel() { return printlevel_; }

private:  
  
  NLNML_FineLevelNoxInterface(NLNML_FineLevelNoxInterface& old);

protected:  
  
  bool isnewJacobian_;
  int numJacobian_;
  int printlevel_;
  double t_;
  int ncalls_computeF_;
  
}; // class NLNML_FineLevelNoxInterface
} // namespace NLNML

#endif
#endif
