// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#ifndef EpetraExt_PERMUTATION_H
#define EpetraExt_PERMUTATION_H

#include <Epetra_IntVector.h>

#include <EpetraExt_Transform.h>

namespace EpetraExt {

/** Permutation stores and describes a permutation matrix P.
    As described in "Matrix Computations" (Golub and Van Loan), a permutation
    matrix is the identity matrix with its rows re-ordered. The permutation is
    internally stored as an integer vector p, where p[i] is the column-index
    of the "1" in P's i-th row.
    Consider the example of permuting a matrix A by applying the permutation
    matrix P to form the result B. i.e., B = PA.
    If p[i] = j, then row j of A becomes row i of B.

    This Permutation class is templated on the type of the object to be
    permuted. However, not all objects are eligible to be template parameters.
    Currently the following objects may be used: Epetra_CrsMatrix, Epetra_CrsGraph
    and Epetra_MultiVector.

    A test program which exercises this Permutation class is located in the
    directory packages/epetraext/test/Permutation.

<pre>
    Implementation Notes:
    Permutation currently inherits StructuralSameTypeTransform, which in turn
    inherits Transform through SameTypeTransform. Permutation, and its base classes,
    are templates. A couple of noteworthy consequences result from this:

       1. A separate instantiation of Permutation must be created for each type
       of object to be permuted. Example:
         Epetra_CrsGraph& graph = ...
         Epetra_CrsMatrix& A = ...
         Permutation<Epetra_CrsGraph> graph_perm(...);
         Permutation<Epetra_CrsMatrix> matrix_perm(...);

         Epetra_CrsMatrix& PA = matrix_perm(A);
	 Epetra_CrsGraph& Pgraph = graph_perm(graph);

       2. Following the semantics of Transform, when the Permutation class is used
       to create a new permuted copy of an object, ownership of the new copy is
       retained by Permutation. Permutation will destroy the new object. This means
       that only one object should be permuted by a Permutation instance.


    It is not clear that these are desirable behaviors for permutations. It is
    possible that Permutation will be altered to remove these limitations, as
    follows:

       1. If Permutation doesn't inherit Transform, then Permutation need not be
       a template and instead we could either overload or template-ize the
       operator() method member. This would allow a single instantiation of
       Permutation to be used for permuting all of the eligible target types.

       2. Allowing the caller (user) to take ownership of the newly-produced
       permuted objects would allow a single Permutation instance to be used
       repeatedly since it would no longer need to hold a pointer to the new object
       for later deletion.

       Then, example usage could look like this:
         Epetra_CrsMatrix& A = ...
	 Epetra_MultiVector& v = ...
         Permutation P(...);

	 Epetra_CrsMatrix PA = P(A);
	 Epetra_MultiVector Pv = P(v);
</pre>

Questions and comments about this class may be directed to Alan Williams.
*/
template<typename T>
class Permutation : public Epetra_IntVector,
                    public EpetraExt::StructuralSameTypeTransform<T> {
 public:
  /** Constructor

      @param CV Set to either Copy or View.

      @param map Defines the index space to be permuted.

      @param permutation Array defining the permutation.
         The length of this array must be 'map.NumMyElements()'. This array
	 is the local portion of the 'p' vector described in the
	 'Detailed Description' section.
   */
  Permutation(Epetra_DataAccess CV,
	      const Epetra_BlockMap& map,
	      int* permutation);

  /** Constructor. This constructor creates an empty permutation object.
      The contents must then be set using regular Epetra_IntVector methods.

      @param map Defines the index space to be permuted.
  */
  Permutation(const Epetra_BlockMap& map);

  /** Copy Constructor */
  Permutation(const Permutation<T>& src);

  /** Destructor */
  virtual ~Permutation();

  typedef typename EpetraExt::SameTypeTransform<T>::TransformTypeRef OutputRef;
  typedef typename EpetraExt::SameTypeTransform<T>::TransformTypeRef InputRef;
  typedef typename EpetraExt::SameTypeTransform<T>::TransformTypePtr OutputPtr;
  typedef typename EpetraExt::SameTypeTransform<T>::TransformTypePtr InputPtr;

  /** This method creates a new object which is a permuted copy of
      the input argument.

      Notes:
      <ul>
      <li> This is a collective function, so in a parallel setting it must be
      called by all processors before any will complete it. (This is because
      map objects are being created, and import/export operations are being
      performed.)
      <li> The new object that is created, if it is a graph or matrix, has
      already had FillComplete() called before it is returned to the user.
      <li> The new object will be destroyed by this permutation object,
      so the caller should not delete it.
      </ul>

      @param orig Input Object to be permuted.
  */
  OutputRef operator()( InputRef orig );

  /** This method creates a new object which is a permuted copy of
      the input argument.

      The same notes apply to this method (regarding collective communications
      etc.) as to the row-permutation operator method above.

      @param orig Input Object to be permuted.

      @param column_permutation Optional Input, defaults to false if not
      provided. A value of false means that a row-permutation will be 
      performed (result = P*orig), a value of true means that a
      column-permutation will be performed (result = orig*P).
  */
  OutputRef operator()( InputRef orig,
		      bool column_permutation );

 private:
  bool isTypeSupported();

  OutputPtr newObj_;
  InputPtr origObj_;
};

}//namespace EpetraExt

#include <EpetraExt_Permutation.cpp>

#endif

