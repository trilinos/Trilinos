/* Copyright (2001) Sandia Corportation. Under the terms of Contract
 *  * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 *   * work by or on behalf of the U.S. Government.  Export of this program
 *    * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 *
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef EDT_PERMUTATION_H
#define EDT_PERMUTATION_H

#include <Epetra_IntVector.h>

#include <Epetra_Transform.h>

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

Questions and comments about this class may be directed to /dev/null.
Just kidding. Direct them to Alan Williams.
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

  typedef typename EpetraExt::Permutation<T>::NewTypeRef NewTypeRef;
  typedef typename EpetraExt::Permutation<T>::OriginalTypeRef OriginalTypeRef;
  typedef typename EpetraExt::Permutation<T>::NewTypePtr NewTypePtr;
  typedef typename EpetraExt::Permutation<T>::OriginalTypePtr OriginalTypePtr;

  /** This method creates a new object which is a permuted copy of
      the input argument. Note that the new object will be destroyed by this
      permutation object.

      @param orig Input Matrix to be permuted.
  */
  NewTypeRef operator()( OriginalTypeRef orig );

  /** This method creates a new object which is a permuted copy of
      the input argument. Note: Column permutations are not yet implemented.
      Note that the new object will be destroyed by this
      permutation object.

      @param orig Input Matrix to be permuted.

      @param column_permutation Optional Input, defaults to false if not
      provided. A value of false means that a row-permutation will be 
      performed (result = P*orig), a value of true means that a
      column-permutation will be performed (result = orig*P).
  */
  NewTypeRef operator()( OriginalTypeRef orig,
			 bool column_permutation );

 private:
  bool isTypeSupported();

  NewTypePtr newObj_;
  OriginalTypePtr origObj_;
};

}//namespace EpetraExt

#endif

