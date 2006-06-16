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

#ifndef _Isorropia_Epetra_utils_hpp_
#define _Isorropia_Epetra_utils_hpp_

#include <Isorropia_configdefs.hpp>
#include <Teuchos_RefCountPtr.hpp>

#ifdef HAVE_EPETRA
class Epetra_Comm;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_RowMatrix;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
#endif

/** Isorropia is the namespace that contains isorropia's declarations
  for classes and functions.
*/
namespace Isorropia {
  class Partitioner;

namespace Epetra_Utils {

#ifdef HAVE_EPETRA
/** Given a Partitioner object, create a target map representing the
   new partitioning.
   This function calls partitioner.compute_partitioning() if it has not
   already been called.
*/
Teuchos::RefCountPtr<Epetra_Map>
create_target_map(const Epetra_Comm& comm, Partitioner& partitioner);

/** Return a vector containing weights that are equal to the number of
  nonzeros per row in the input_matrix. The returned vector will have
  the same size and distribution as input_matrix's row-map.
*/
Epetra_Vector* create_row_weights_nnz(const Epetra_RowMatrix& input_matrix);

/** Return a vector containing weights that are equal to the number of
  nonzeros per row in the input_graph. The returned vector will have
  the same size and distribution as input_graph's row-map.
*/
Epetra_Vector* create_row_weights_nnz(const Epetra_CrsGraph& input_graph);

/** Calculate a new partitioning, and fill lists with new elements for
    the local partition, as well as export and import lists.
*/
int
repartition(const Epetra_BlockMap& input_map,
	    const Epetra_Vector& weights,
            std::vector<int>& myNewElements,
            std::map<int,int>& exports,
            std::map<int,int>& imports);

/** Given an Epetra_BlockMap object, fill a vector of length numprocs+1
  with each processor's starting offset into the Epetra_BlockMap's global
  set of elements (the last position will contain num-global-elements).
  Gather the vector of offsets onto all processors.
*/
void gather_all_proc_global_offsets(const Epetra_BlockMap& blkmap,
                                    std::vector<int>& all_proc_offsets);

#endif //HAVE_EPETRA

}//namespace Epetra_Utils
}//namespace Isorropia

#endif

