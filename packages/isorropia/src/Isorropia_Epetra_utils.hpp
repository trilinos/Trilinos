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

#ifdef HAVE_EPETRA
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Vector;
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
#endif

/** Isorropia is the namespace that contains isorropia's declarations
  for classes and functions.
*/
namespace Isorropia {
namespace Epetra_Utils {

#ifdef HAVE_EPETRA
/** Return a vector containing weights that are equal to the number of
  nonzeros per row in the input_matrix. The returned vector will have
  the same size and distribution as input_matrix's row-map.
*/
Epetra_Vector* create_row_weights_nnz(const Epetra_RowMatrix& input_matrix);

/** Return a Epetra_Map that describes a row-distribution that would
  place roughly the same 'weight' on each processor if the
  input-matrix were imported into that distribution. The 'weight' quantity
  is the sum of corresponding entries in the 'weights' vector, which must
  have the same size and distribution as the row-map of 'input_matrix'.
  The returned Epetra_Map object is passed by value, utilizing the
  fact that Epetra_Map is a light-weight reference-counted "front-end"
  object with an underlying data-object.
*/
Epetra_Map create_rowmap_balanced(const Epetra_RowMatrix& input_matrix,
                                  const Epetra_Vector& weights);

/** Given an Epetra_BlockMap object, fill a vector of length num-procs+1
  with each processor's starting offset into the Epetra_BlockMap's global
  set of elements (the last position will contain num-global-elements).
  Gather the vector of offsets onto all processors.
*/
void gather_all_proc_global_offsets(const Epetra_BlockMap& blkmap,
                                    std::vector<int>& all_proc_offsets);

/** Import input_matrix into target_matrix.
  On entry to this function, target_matrix is assumed to be constructed
  with the row-map for the desired distribution.
*/
void import_matrix(const Epetra_CrsMatrix& input_matrix,
                   Epetra_CrsMatrix& target_matrix);

#endif //HAVE_EPETRA

}//namespace Epetra_Utils
}//namespace Isorropia

#endif

