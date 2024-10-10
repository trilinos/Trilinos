//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef EpetraExt_PERMUTATION_H
#define EpetraExt_PERMUTATION_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include <Epetra_ConfigDefs.h>
#include <Epetra_GIDTypeVector.h>

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
template<typename T, typename int_type>
class TPermutation : public Epetra_GIDTypeVector<int_type>::impl,
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
  TPermutation(Epetra_DataAccess CV,
	      const Epetra_BlockMap& map,
	      int_type* permutation);

  /** Constructor. This constructor creates an empty permutation object.
      The contents must then be set using regular Epetra_IntVector methods.

      @param map Defines the index space to be permuted.
  */
  TPermutation(const Epetra_BlockMap& map);

  /** Copy Constructor */
  TPermutation(const TPermutation<T, int_type>& src);

  /** Destructor */
  virtual ~TPermutation();

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
  TPermutation<T, int_type>& operator=(const TPermutation<T, int_type>& src)
    {
      //not currently supported
      abort();
      return(*this);
    }

  bool isTypeSupported();

  OutputPtr newObj_;
  InputPtr origObj_;
};

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
template<typename T>
class Permutation : public TPermutation<T, int> {
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
	      int* permutation)
		  : TPermutation<T, int>(CV, map, permutation)
  {
  }

  /** Constructor. This constructor creates an empty permutation object.
      The contents must then be set using regular Epetra_IntVector methods.

      @param map Defines the index space to be permuted.
  */
  Permutation(const Epetra_BlockMap& map)
	  : TPermutation<T, int>(map)
  {
  }

  /** Copy Constructor */
  Permutation(const Permutation<T>& src)
	  : TPermutation<T, int>(src)
  {
  }
};
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

template<typename T>
class Permutation64 : public TPermutation<T, long long> {
public:
  /** Constructor

      @param CV Set to either Copy or View.

      @param map Defines the index space to be permuted.

      @param permutation Array defining the permutation.
         The length of this array must be 'map.NumMyElements()'. This array
	 is the local portion of the 'p' vector described in the
	 'Detailed Description' section.
   */
  Permutation64(Epetra_DataAccess CV,
	      const Epetra_BlockMap& map,
	      long long* permutation)
		  : TPermutation<T, long long>(CV, map, permutation)
  {
  }

  /** Constructor. This constructor creates an empty permutation object.
      The contents must then be set using regular Epetra_IntVector methods.

      @param map Defines the index space to be permuted.
  */
  Permutation64(const Epetra_BlockMap& map)
	  : TPermutation<T, long long>(map)
  {
  }

  /** Copy Constructor */
  Permutation64(const Permutation64<T>& src)
	  : TPermutation<T, long long>(src)
  {
  }
};
#endif

}//namespace EpetraExt

#include <EpetraExt_Permutation_impl.h>

#endif
