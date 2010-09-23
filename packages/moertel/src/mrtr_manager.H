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
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file mrtr_manager.H
 *
 * \class MOERTEL::Manager
 *
 * \brief Top level user interface to the mortar package 'Moertel'
 *
 * \date Last update do Doxygen: 09-July-2010
 *
 */
#ifndef MOERTEL_MANAGER_H
#define MOERTEL_MANAGER_H

#include <ctime>
#include <iostream>
#include <vector>

// Trilinos includes
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

// Moertel includes
#include "mrtr_segment.H"
#include "mrtr_node.H"
#include "mrtr_interface.H"
#include "mrtr_lm_selector.H"
#include "mrtr_solver.H"
#include "mrtr_utils.H"




/*!
\brief MOERTEL: namespace of the Moertel package

The Moertel package depends on \ref Epetra, \ref EpetraExt, \ref Teuchos,
\ref Amesos, \ref ML and \ref AztecOO:<br>
Use at least the following lines in the configure of Trilinos:<br>
\code
--enable-moertel 
--enable-epetra 
--enable-epetraext
--enable-teuchos 
--enable-ml
--enable-aztecoo --enable-aztecoo-teuchos 
--enable-amesos
\endcode

*/

namespace MOERTEL
{
// forward declarations
class Solver;

/*!
\class Manager

\brief <b> Top level user interface to mortar package 'Moertel' </b>

This class supplies capabilities for nonconforming mesh tying in 2 and 3 
dimensions using Mortar methods. <br>
It also supports use of mortar methods for contact formulations and 
domain decomposition as well as some limited support to solving 
mortar-coupled systems of equations.

\b Proceeding:
- The user constructs an arbitrary number of conforming or non-conforming 
 interfaces using the \ref MOERTEL::Interface class and passes them to this
 \ref MOERTEL::Manager.
- This class will then choose appropriate mortar and slave sides on all interfaces
 (if not already done by the user) such that conflicts are avoided on 
 corner and edge nodes that are shared among interfaces.
 It will detect such corner/edge nodes based on the user input and will perform
 necessary modifications to the user chosen Lagrange multiplier shape functions
 close to the boundary of an interface.
- It will then perform the construction of a C0-continuous field of normals on
 each slave side of an interface and will use this field to construct the
 Mortar projection operator that projects nodal values from the mortar to the slave side
- It will perform the integration of the Mortar integrals and assemble the
 coupling matrix \b D that couples the interface degrees of freedom of 
 the slave side to the Lagrange multipliers (LMs) and the coupling matrix \b M that
 couples the degrees of freedom of the master side to the LMs.
- A saddle point system of equations can be constructed and returned to the user
 to solve if the matrix of the original (uncoupled) problem is supplied.
- Some limited support in solving the coupled system of equations using direct
 and iterative, multigrid preconditioned solvers is supplied. (not yet implemented)
- The user has access to coupling matrices \b D and \b M to use in contact formulations (not yet impl.).
- It supplies a projection operator that projects arbitrary nodal field values
  from the slave side to the master side and vice versa to use in domain decomposition methods.
  (not yet impl.)    
  
The \ref MOERTEL::Manager class supports the ostream& operator <<

The package and this class make use of the other Trilinos packages 
\b Epetra , \b EpetraExt , \b Teuchos , \b Amesos , \b ML .

For more detailed information on the usage, see "The Moertel User's Guide", Sandia report XXXX

For background information on Mortar methods, we refer to <br>
Wohlmuth, B.:"Discretization Methods and Iterative Solvers Based on Domain Decomposition", 
Springer, New York, 2001, ISBN 3-540-41083-X

\author Glen Hansen (Glen.Hansen@inl.gov)

*/
class Manager 
{
public:

  friend class MOERTEL::Solver;
  
  /*!
  \brief Specification of the problem dimension
  
  */
  enum DimensionType
  {
    manager_none,
    manager_2D,
    manager_3D
  };
  
  // @{ \name Constructors and destructor
  
  /*!
  \brief Creates and empty instance of this class
  
  Constructs an empty instance of this class to be filled by the user with
  information about the non-conforming interface(s). <br>
  This is a collective call for all processors associated with the
  Epetra_Comm 
  
  \param comm : An Epetra_Comm object
  \param dimension : Dimension of the problem
  \param outlevel : The level of output to stdout to be generated ( 0 - 10 )
  */
  explicit Manager(Epetra_Comm& comm, 
                   MOERTEL::Manager::DimensionType dimension, 
                   int outlevel);
                   
  
  /*!
  \brief Creates and empty instance of this class
  
  Constructs an empty instance of this class that is filled in by the user with
  information about the non-conforming interface(s). <br>
  This is a collective call for all processors associated with the
  Epetra_Comm 
  
  \param comm : An Epetra_Comm object
  \param outlevel : The level of output to stdout to be generated ( 0 - 10 )
  */
  explicit Manager(Epetra_Comm& comm, int outlevel);

  /*!
  \brief Destroys an instance of this class
  
  Destructor
  
  */
  virtual ~Manager();
  
  //@}
  
  // @{ \name Query methods
  
  /*!
  \brief Returns the Epetra_Comm object associated with this class
  */
  Epetra_Comm& Comm() const {return comm_;}

  /*!
  \brief Returns the Level of output (0 - 10) the user specified in the constructor
  */
  int OutLevel() { return outlevel_; }

  /*!
  \brief Query the problem dimension

  Query the problem dimension
  
  */
  MOERTEL::Manager::DimensionType Dimension() { return dimensiontype_; }
  
  /*!
  \brief Query the number of interfaces passed to this class via AddInterface

  Returns the number of non-conforming interfaces passed to this class via AddInterface

  */
  int Ninterfaces() { return interface_.size(); }

  /*!
  \brief Print all information stored in this class to stdout

  */
  bool Print() const;
  
  //@}

  // @{ \name Construction methods


  /*!
  \brief Set problem dimension
  
  This class can handle 2D and 3D problems but not both at the same time.
  It is necessary to specify the dimension of the problem
  
  \param type : Dimension of the problem
  */
  void SetDimension(MOERTEL::Manager::DimensionType type) { dimensiontype_ = type; return; }
  
  /*!
  \brief Add an interface class to this class
  
  Add a previously constructed interface to this class. </br>
  Before adding an interface, the interface's public method 
  \ref MOERTEL::Interface::Complete() has to be called.</br>
  This class will not accept interface that are not completed.
  Any number of interfaces can be added. This class does not take ownership over the
  added class. The interface added can be destroyed immediately after the call to this method
  as this class stores a deep copy of it. 
    
  \param interface : A completed interface
  
  \return true when adding the interface was successful and false otherwise
  */
  bool AddInterface(MOERTEL::Interface& interface);

  /*!
  \brief Obtain a parameter list with default values
  
  Obtain a parameter list containing default values for 
  integration of mortar interfaces.
  Default values set are:
  \code
    // When true, values of shape functions are exact at
    // every single integration point.
    // When false, linear shape functions of an overlap 
    // discretization are used to linearly interpolate 
    // shape functions in the overlap region.
    // The latter is significantly cheaper but only 
    // recommended for truly linear functions, e.g. with 
    // 1D interfaces discretized with linear elements
    // 2D interfaces discretized with linear triangle elements
    // It is NOT recommended for 2D interfaces discretized with
    // bilinear quads as shape functions are not truly linear.
    integrationparams_->set("exact values at gauss points",true);

    // 1D interface possible values are 1,2,3,4,5,6,7,8,10
    integrationparams_->set("number Gaussian points 1D",3);

    // 2D interface possible values are 3,6,12,13,16,19,27
    integrationparams_->set("number Gaussian points 2D",12);
  \endcode

  \sa SetProblemMap , AddInterface, Mortar_integrate
    
  */
  Teuchos::ParameterList& Default_Parameters();

  /*!
  \brief Perform integration of the mortar integral on all interfaces
  
  Once all interfaces and a row map of the original system are passed to 
  this class this method evaluates the mortar integrals on all interfaces.
  
  Process:
  - Chooses a slave and a mortar side for each interface of not already
    chosen by the user in the interface class
  - Build a C0-continuous field of normals on the slave side
  - Makes special boundary modifications to the Lagrange multiplier shape functions
  - Performs integration of the mortar integrals on individual interfaces
  - Chooses degrees of freedom for the Lagrange multipliers based on 
    the row map of the original uncoupled problem.
  - Builds a row map of the saddle point system
  - Assembles interface contributions to coupling matrices D and M
  
  \return true if successful, false otherwise
  
  \sa SetProblemMap , AddInterface , Integrate_Interfaces
    
  */
  bool Mortar_Integrate();

  /*!
  \brief Perform the integration of the mortar integral on all interfaces.
  
  Once all interfaces and a row map of the original system are passed to 
  this class this method evaluates the mortar integrals on all interfaces.
  
  Process:
  - Chooses a slave and a mortar side for each interface of not already
    chosen by the user in the interface class
  - Build a C0-continuous field of normals on the slave side
  - Makes special boundary modifications to the Lagrange multiplier shape functions
  - Performs integration of the mortar integrals on individual interfaces
  - Chooses degrees of freedom for the Lagrange multipliers based on 
    the row map of the original uncoupled problem.
  - Builds a row map of the saddle point system
  
  \return true if successful, false otherwise
  
  \sa SetProblemMap , AddInterface , Mortar_Integrate
    
  */

  bool Integrate_Interfaces();

  /*!
  \brief Assembles contributions from the D and M matrices into the JFNK residual vector.
  
  Once all interfaces and a row map of the original system are passed to 
  this class this method evaluates the mortar integrals on all interfaces.
  
  \param sel (in): A user-derived object inheriting from Lmselector. 
  \param rhs (out) : JFNK rhs vector containing accumulated D and M values
  \param soln (in) : Current JFNK solution vector, used to matvec with D and M
  
  \return true if successful, false otherwise
  
  */

  bool AssembleInterfacesIntoResidual(Lmselector *sel);
  
  
  //@}

  // @{ \name Solution methods

  /*!
  \brief Set the RowMap of the original uncoupled problem
  
  Passes in the RowMap of the original (uncoupled) problem.
  In a saddle point system, this would be the row map of the (1,1) block.
  This class does not take ownership of the passed in map. The map can
  be destroyed immediately after the call to this method as this class
  stores a deep copy of it.
  A row map has to be passed to this class before a call to \ref Mortar_Integrate() 
    
  \param map : Row map of the original (uncoupled) problem
  
  \sa Mortar_Integrate()
  
  */
  bool SetProblemMap(const Epetra_Map* map);
  
  /*!
  \brief Return view of row map of the uncoupled problem
  
  return a view of the row map of the uncoupled problem the user provided using
  SetProblemMap. Returns NULL if no row map was previously set.
  
  \sa SetProblemMap

  */
  Epetra_Map* ProblemMap() const { return problemmap_.get(); }
  
  /*!
  \brief Set the Matrix of the original uncoupled problem
  
  Passes in the Matrix of the original (uncoupled) problem.
  In a saddle point system, this would be the (1,1) block.
  This class takes ownership of the matrix if DeepCopy is set to true. 

  The matrix of the original, uncoupled problem is used to
  construct a saddle point system via \ref MakeSaddleProblem() 
  or solve the coupled problem
    
  \param map : Matrix of the original (uncoupled) problem
  
  \sa MakeSaddleProblem() 
  
  */
  bool SetInputMatrix(Epetra_CrsMatrix* inputmatrix, bool DeepCopy = false);
  
  /*!
  \brief Construct a saddle point system of equations
  
  After a call to \ref Mortar_Integrate() a saddle point system of equations
  can be constructed and returned to the user. Prerequisite is that the user has 
  supplied the original uncoupled matrix via \ref SetInputMatrix and
  the integration was performed with \ref Mortar_Integrate().
  This class holds ownership of the saddle point system so the
  user must not destroy it.
  
  */
  Epetra_CrsMatrix* MakeSaddleProblem();

  /*!
  \brief Construct a spd system of equations if possible
  
  If dual Lagrange multipliers are used, a symm. positive definite 
  system of equations can be constructed from the saddle point problem.
  The final solution of the problem is then obtained from<br>
  Atilde x = RHSMatrix * b
  
  \sa GetSPDRHSMatrix()
  */
  Epetra_CrsMatrix* MakeSPDProblem();

  /*!
  \brief returns the righ hand side matrix of the spd system of equations
  
  If dual Lagrange multipliers are used, a symm. positive definite 
  system of equations can be constructed from the saddle point problem.
  The final solution of the problem is then obtained from<br>
  Atilde x = RHSMatrix * b
  
  \sa MakeSPDProblem()
  */
  Epetra_CrsMatrix* GetRHSMatrix();

  /*!
  \brief Set a parameter list holding solver parameters
  
  Set a ParameterList containing solver options. This class does not take
  ownership over params but instead uses a view of it.<br>
  Currently Moertel has interfaces to the Amesos, the AztecOO and the ML package.<br>
  All packages are controlled by this parameter list. The parameter list
  Contains a general part where the package and the type of linear system
  to be generated is chosen and contains one or more sub-parameter lists which
  are passed on to the appropriate package. This way all parameters described
  in the documentation of Amesos, ML and Aztec can be passed to the appropriate
  package.<br>

  Below one example how to choose the parameter list is given:<br>
  \code
  Teuchos::ParameterList params;
  params.set("System","SaddleSystem"); // use a saddle point system of equations in solve
  // or
  params.set("System","SPDSystem");    // use a spd system with condensed Lagrange multipliers in solve

  // choose solver package
  params.set("Solver","Amesos");       // use solver package Amesos
  // or
  params.set("Solver","ML/Aztec");     // use solver packages ML and AztecOO
  
  // argument sublist passed to and used for Amesos
  // see Amesos documentation, all configured Amesos solvers can be used
  // This is an example of possible parameters
  Teuchos::ParameterList& amesosparams = params.sublist("Amesos");
  amesosparams.set("Solver","Amesos_Klu"); // name of Amesos solver to use
  amesosparams.set("PrintTiming",false);   
  amesosparams.set("PrintStatus",false);   
  amesosparams.set("UseTranspose",true);

  // argument sublist passed to and used for ML
  // see ML documentation, all configured parameters recognized by the
  // ML_Epetra::ml_MultiLevelPreconditioner class can be used
  // Moertel sets the ML default parameters first which are then overridden by
  // this list
  // This is an example configuration:
  Teuchos::ParameterList& mlparams = params.sublist("ML");
  ML_Epetra::SetDefaults("SA",mlparams);
  mlparams.set("output",6);
  mlparams.set("print unused",-2);
  mlparams.set("increasing or decreasing","increasing");
  mlparams.set("PDE equations",3);
  mlparams.set("max levels",20);
  mlparams.set("aggregation: type","Uncoupled");
  mlparams.set("coarse: max size",2500);
  mlparams.set("coarse: type","Amesos-KLU");
  mlparams.set("smoother: type","MLS");
  mlparams.set("smoother: MLS polynomial order",3);
  mlparams.set("smoother: sweeps",1);
  mlparams.set("smoother: pre or post","both");
  mlparams.set("null space: type","pre-computed");
  mlparams.set("null space: add default vectors",false);
  int dimnullspace  = 6;
  int dimnsp        = dimnullspace*nummyrows;
  double* nsp       = new double[dimnsp];
  application_compute_nullspace(nsp,dimnsp);
  mlparams.set("null space: dimension",dimnullspace);
  // the user is responsible for freeing nsp after solve
  mlparams.set("null space: vectors",nsp); 

  // argument sublist passed to and used for AztecOO
  // see AztecOO documentation, all configured parameters recognized by the
  // AztecOO class can be used
  // Moertel sets the Aztec default parameters first which are then overridden by
  // this list
  // This is an example configuration:
  Teuchos::ParameterList& aztecparams = params.sublist("Aztec");
  aztecparams.set("AZ_solver","AZ_cg");
  // if "AZ_user_precond" is specified, ML will be used as a preconditioner
  // Any other aztec-buildin preconditioner can be used as well 
  aztecparams.set("AZ_precond","AZ_user_precond");
  aztecparams.set("AZ_max_iter",500);
  aztecparams.set("AZ_output",100);
  aztecparams.set("AZ_tol",1.0e-08);
  aztecparams.set("AZ_scaling","AZ_none");

  \endcode
  
  \param params : Parameter list containing solver options
  
  \warning When using ML and/or Aztec, one has to use "SPDSystem" as system matrix
           as the amg-preconditioned iterative method will currently fail on the 
           saddle point system.
  */
  void SetSolverParameters(Teuchos::ParameterList& params);

  /*!
  \brief Solve the problem
  
  Solve the coupled problem using the solver options supplied in params.
  Succesive calls to this method will reuse the system matrix together with
  the supplied solution and rhs vectors until
  \ref ResetSolver() is called. Then, the linear problem is recreated from
  scratch.
  
  \param params (in) : Solution parameters
  \param sol (out) : The solution
  \param rhs (in) : The right hand side vector
  */
  bool Solve(Teuchos::ParameterList& params, Epetra_Vector& sol, const Epetra_Vector& rhs);

  /*!
  \brief Solve the problem
  
  Solve the coupled problem using the solver options previously 
  supplied in using \ref SetSolverParameters
  Succesive calls to this method will resue the system matrix together with
  the supplied solution and rhs vectors until
  \ref ResetSolver() is called. Then, the linear problem is recreated from
  scratch.
  
  \param sol (out) : The solution
  \param rhs (in) : The right hand side vector
  */
  bool Solve(Epetra_Vector& sol, const Epetra_Vector& rhs);
  
  /*!
  \brief Reset the internal solver
  
  The internal solver (no matter which one) will go on using the
  once constructed (and maybe factored) matrix until this method is
  called. After a call to this method, the linear problem as well as
  the solver are build from scratch.
  
  */
  void ResetSolver();
  
  /*!
  \brief Get pointer to constraints matrix D
  
  */
  const Epetra_CrsMatrix* D() const { return D_.get();}

  /*!
  \brief Get pointer to constraints matrix M
  
  */
  const Epetra_CrsMatrix* M() const { return M_.get();}

  /*!
  \brief Query a managed interface added using AddInterface
  
  Accessor function to the interfaces stored within the manager by
  previous AddInterface calls.

  \param idx : index of the interface to be returned (0 up to Ninterfaces)

  \return the requested interface
  
  \sa AddInterface, Ninterfaces
  */
  Teuchos::RCP<MOERTEL::Interface> GetInterface(int idx) { return interface_[idx]; }

  /*!
   */
  Teuchos::RCP<Epetra_Map> GetSaddleMap() { return saddlemap_; }

  //@}

private:  
  // don't want = operator and copy-ctor
  Manager operator = (const Manager& old);
  Manager(MOERTEL::Manager& old);

  // Build the M and D matrices 2D case
  bool Build_MD_2D();
  
  // just do integration of all interfaces 2D case
  bool Integrate_Interfaces_2D();
  
  // Build the M and D matrices 3D case
  bool Build_MD_3D();
  
  // just do integration of all interfaces 3D case
  bool Integrate_Interfaces_3D();
  
  // build the map for the saddle point problem
  bool BuildSaddleMap();

  // choose the dofs for Lagrange multipliers
  Teuchos::RCP<Epetra_Map> LagrangeMultiplierDofs();
  
  // automatically choose mortar side (called when mortarside==-2 on any interface)
  bool ChooseMortarSide();
  bool ChooseMortarSide_2D(std::vector<Teuchos::RCP<MOERTEL::Interface> >& inter);
  bool ChooseMortarSide_3D(std::vector<Teuchos::RCP<MOERTEL::Interface> >& inter);
  
protected:

  int                           outlevel_;            // output level (0-10)
  Epetra_Comm&                  comm_;                // communicator (global, contains ALL procs)
  DimensionType                 dimensiontype_;       // problem dimension

  std::map<int,Teuchos::RCP<MOERTEL::Interface> >  interface_; // the interfaces

  Teuchos::RCP<Epetra_Map>       problemmap_;          // the rowmap of the input matrix
  Teuchos::RCP<Epetra_CrsMatrix> inputmatrix_;         // the uncoupled matrix on input
  Teuchos::RCP<Epetra_Map>       constraintsmap_;      // the rowmap of M and D (both of them use comm_)
  Teuchos::RCP<Epetra_CrsMatrix> D_;                   // the coupling matrix D
  Teuchos::RCP<Epetra_CrsMatrix> M_;                   // the coupling matrix M
  Teuchos::RCP<Epetra_Map>       saddlemap_;           // the rowmap of the saddlepointproblem
  Teuchos::RCP<Epetra_CrsMatrix> saddlematrix_;        // the matrix of the saddle problem;
  Teuchos::RCP<Epetra_CrsMatrix> spdmatrix_;           // the spd matrix of the problem;
  Teuchos::RCP<Epetra_CrsMatrix> spdrhs_;              // the matrix to left-multiply the rhs vector with for the spd problem
  
  Teuchos::RCP<std::map<int,int> >    lm_to_dof_;           // maps lagrange multiplier dofs to primary dofs (slave side))
  
  Teuchos::RCP<Epetra_CrsMatrix> I_;
  Teuchos::RCP<Epetra_CrsMatrix> WT_;
  Teuchos::RCP<Epetra_CrsMatrix> B_;
  Teuchos::RCP<Epetra_Map>       annmap_;

  Teuchos::RCP<Teuchos::ParameterList> integrationparams_; // parameter list with integration parameters
  Teuchos::RCP<Teuchos::ParameterList> solverparams_;  // solution parameters passed from the user    

  Teuchos::RCP<MOERTEL::Solver>        solver_;        // the mortar solver


};

} // namespace MOERTEL

/*!
\brief << operator

outputs all information stored in the \ref MOERTEL::Manager class to ostream

*/
std::ostream& operator << (ostream& os, const MOERTEL::Manager& node);

#endif // MOERTEL_MANAGER_H
