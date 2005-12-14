/*
#@HEADER
# ************************************************************************
#
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
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
 * \date Last update do Doxygen: 14-Dec-05
 *
 */
#ifdef TRILINOS_PACKAGE

#ifndef MOERTEL_MANAGER_H
#define MOERTEL_MANAGER_H

#include <ctime>
#include <iostream>
#include <vector>

// Trilinos includes
#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include "Teuchos_RefCountPtr.hpp"
// mrtr includes
#include "mrtr_segment.H"
#include "mrtr_node.H"
#include "mrtr_interface.H"



using namespace std;
using namespace Teuchos;

//! MOERTEL: namespace of the Moertel package
namespace MOERTEL
{

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
- This class will then choose approbiate mortar and slave sides on all interfaces
 (if not already done by the user) such that conflicts are avoided on 
 corner and edge nodes that are shared among interfaces.
 It will detect such corner/edge nodes based on the user input and will perform
 neccessary modifications to the user chosen lagrange mutliplier shape functions
 close to the boundary of an interface.
- It will then perform the construction of a C0-continous field of normals on
 each slave side of an interface and will use this field to construct the
 Mortar projection operator that projects nodal values from the mortar to the slave side
- It will perform the integration of the Mortar integrals and assemble the
 coupling matrix \b D that couples the interface degrees of freedom of 
 the slave side to the Lagrange multipliers (LMs) and the coupling matrix \b M that
 couples the degrees of freedom of the master side to the LMs.
- A saddle point system of equations can be constructed and returned to the user
 to solve if the matrix of the original (uncoupled) problem is supplied.
- Some limited support in solving the coupled system of equations using direct
 and iterative, mutligrid preconditioned solvers is supplied. (not yet implemented)
- The user has access to coupling matrices \b D and \b M to use in contact formulations (not yet impl.).
- It supplies a projection operator that projects arbitrary nodal field values
  from the slave side to the master side and vice versa to use in domain decomposition methods.
  (not yet impl.)    
  
The \ref MOERTEL::Manager class supports the << operator  

The package and this class make use of the other Trilinos packages 
\b Epetra , \b EpetraExt , \b Teuchos , \b Amesos , \b ML .

For more detailed information on the usage, see "The Moertel User's Guide", Sandia report XXXX

For background information on Mortar methods, we refer to <br>
Wohlmuth, B.:"Discretization Methods and Iterative Solvers Based on Domain Decomposition", 
Springer, New York, 2001, ISBN 3-540-41083-X

\author Michael W. Gee (mwgee@sandia.gov)

*/
class Manager 
{
public:

  /*!
  \brief Specification of the problem dimension
  
  */
  enum DimensionType
  {
    manager_none,
    manager_2D,
    manager_3D
  };
  
  // @{ \name Constructors and destructors
  
  /*!
  \brief Creates and empty instance of this class
  
  Constructs an empty instance of this class that is then filled buy the user with
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
  Before adding an interface,the interfaces's public method 
  \ref MOERTEL::Interface::Complete() has to be called.</br>
  This class will not accept interface that are not completed.
  Any number of interfaces can be added. This class does not take ownership over the
  added class. The interface added can be destroyed immediately after the call to this method
  as this class stores a deep copy of it. 
    
  \param interface : A completed interface
  
  \return true when adding the interface was successfull and false otherwise
  */
  bool AddInterface(MOERTEL::Interface& interface);

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
  bool SetInputMap(Epetra_Map* map);
  
  /*!
  \brief Set the Matrix of the original uncoupled problem
  
  Passes in the Matrix of the original (uncoupled) problem.
  In a saddle point system, this would be the the (1,1) block.
  This class takes ownership of the matrix if DeepCopy is set to true. 

  The matrix of the original, uncoupled problem is used to
  construct a saddle point system via \ref MakeSaddleProblem() 
  or solve the coupled problem
    
  \param map : Matrix of the original (uncoupled) problem
  
  \sa MakeSaddleProblem() 
  
  */
  bool SetInputMatrix(Epetra_CrsMatrix* inputmatrix, bool DeepCopy = false);
  
  /*!
  \brief Perform integration of the mortar integral on all interfaces
  
  Once all interfaces and a row map of the original system are passed to 
  this class this method evaluates the mortar integrals on all interfaces.
  
  Proceding:
  - Chooses a slave and a mortar side for each interface of not already
    chosen by the user in the interface class
  - Build a C0-continous field of normals on the slave side
  - Makes special boundary modifications to the lagrange mutliplier shape functions
  - Performs integration of the mortar integrals on individual interfaces
  - Chooses degrees of freedom for the lagrange multipliers based on 
    the row map of the original uncoupled problem.
  - Builds a row map of the saddle point system
  - Assembles interface contributions to coupling matrices D and M
  
  \return true if successfull, false otherwise
  
  \sa SetInputMap , AddInterface
    
  */
  bool Mortar_Integrate();
  
  /*!
  \brief Construct a saddle point system of equations
  
  After a call to \ref Mortar_Integrate() a saddle point system of equations
  can be constructed and returned to the user. Prerequisite is that the user has 
  supplied the original uncoupled matrix via \ref SetInputMatrix and
  the integration was performed with \ref Mortar_Integrate().
  This class holds ownership of the saddle point system so the
  user must not destroy it.
  
  \return true if succesfull, false otherwise
  
  */
  Epetra_CrsMatrix* MakeSaddleProblem();
  
  //@}


private:  
  // don't want = operator and copy-ctor
  Manager operator = (const Manager& old);
  Manager(MOERTEL::Manager& old);

  // do integration of all interfaces 2D case
  bool Mortar_Integrate_2D();
  
  // do integration of all interfaces 3D case
  bool Mortar_Integrate_3D();
  
  // build the map for the saddle point problem
  bool BuildSaddleMap();

  // choose the dofs for lagrange multipliers
  RefCountPtr<Epetra_Map> LagrangeMultiplierDofs();
  
  // modified version of the epetraext matrixmatrixadd
  int MatrixMatrixAdd(const Epetra_CrsMatrix& A, bool transposeA,double scalarA,
                      Epetra_CrsMatrix& B,double scalarB );

  // automatically choose mortar side (called when mortarside==-2 on any interface)
  bool ChooseMortarSide();
  bool ChooseMortarSide_2D(vector<RefCountPtr<MOERTEL::Interface> >& inter);
  bool ChooseMortarSide_3D(vector<RefCountPtr<MOERTEL::Interface> >& inter);
  
private:

  int                           outlevel_;            // output level (0-10)
  Epetra_Comm&                  comm_;                // communicator (global, contains ALL procs)
  DimensionType                 dimensiontype_;       // problem dimension

  RefCountPtr<Epetra_Map>       inputmap_;            // the rowmap of the input matrix
  RefCountPtr<Epetra_CrsMatrix> inputmatrix_;         // the uncoupled matrix on input
  RefCountPtr<Epetra_Map>       constraintsmap_;      // the rowmap of M and D (both of them use comm_)
  RefCountPtr<Epetra_CrsMatrix> D_;                   // the coupling matrix D
  RefCountPtr<Epetra_CrsMatrix> M_;                   // the coupling matrix M
  RefCountPtr<Epetra_Map>       saddlemap_;           // the rowmap of the saddlepointproblem
  RefCountPtr<Epetra_CrsMatrix> saddlematrix_;        // the matrix of the saddle problem;
  
  map<int,RefCountPtr<MOERTEL::Interface> >  interface_; // the interfaces
};

} // namespace MOERTEL

/*!
\brief << operator

outputs all information stored in the \ref MOERTEL::Manager class to ostream

*/
std::ostream& operator << (ostream& os, const MOERTEL::Manager& node);

#endif // MOERTEL_MANAGER_H
#endif // TRILINOS_PACKAGE
