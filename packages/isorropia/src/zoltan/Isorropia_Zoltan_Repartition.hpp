//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

#ifndef _Isorropia_Zoltan_Repartition_hpp_
#define _Isorropia_Zoltan_Repartition_hpp_

#include <Isorropia_configdefs.hpp>

#ifdef HAVE_ISORROPIA_ZOLTAN

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_EPETRA
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_LinearProblem;
#endif

#include <Isorropia_ZoltanQuery.h>

/** Isorropia_Zoltan is the namespace that contains isorropia's
  Zoltan-specific classes and functions.
*/
namespace Isorropia_Zoltan {

#ifdef HAVE_EPETRA

/** Calculate a new partitioning, and fill lists with new elements for
    the local partition, as well as export and import lists.
*/
int
repartition(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
	    Teuchos::ParameterList& paramlist,
            std::vector<int>& myNewElements,
            std::map<int,int>& exports,
            std::map<int,int>& imports);

/** Calculate a new partitioning, and fill lists with new elements for
    the local partition, as well as export and import lists.
*/
int
repartition(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
	    Teuchos::ParameterList& paramlist,
            std::vector<int>& myNewElements,
            std::map<int,int>& exports,
            std::map<int,int>& imports);

int
load_balance(MPI_Comm comm,
	     Teuchos::ParameterList& paramlist,
	     Zoltan::QueryObject& queryObject,
	     std::vector<int>& myNewElements,
	     std::map<int,int>& exports,
	     std::map<int,int>& imports);

#endif //HAVE_EPETRA

}//namespace Isorropia_Zoltan

//the following endif closes the '#ifdef HAVE_ISORROPIA_ZOLTAN' block.
#endif

#endif

