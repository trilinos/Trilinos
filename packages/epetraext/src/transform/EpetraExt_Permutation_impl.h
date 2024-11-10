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
#ifndef EpetraExt_PERMUTATION_IMPL_H
#define EpetraExt_PERMUTATION_IMPL_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include <EpetraExt_ConfigDefs.h>

#include <EpetraExt_Permutation.h>

#include <Epetra_Export.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_GIDTypeVector.h>

namespace EpetraExt {

/** Define some traits to make it easier to deal with template-parameters which
   are objects to be permuted. Given a template parameter, we'll want to
   have the following operations available:
   <ul>
   <li>determine the type
   <li>construct an instance of it
   <li>replace its row-map
   <li>produce a column-permutation of it
   </ul>

   First the default definition, which catches all types "T", followed by some
   specializations for anticipated types. Any type other than the types
   specifically anticipated will be handled by this default definition,
   allowing the Permutation class to abort or return NULL where appropriate.

   We define these trait structs in this file rather than in a separate file
   in an attempt to avoid some template-instantiation complications...
*/
template<class T>
struct Perm_traits {
  /** return a std::string name for the object type */
  static const char* typeName()
  { static const char name[] = "unknown"; return( name ); }

  /** clone function accepts an example of the object being cloned, and enough
      constructor arguments to be able to create any of these: CrsMatrix,
      CrsGraph, Vector, MultiVector.   And probably more later...

      Why is an example object needed? For instance, if a MultiVector is
      created, we may want to know how many vectors it should contain...
  */
  static T* clone(T* example,
		  Epetra_DataAccess CV,
		  const Epetra_BlockMap& map,
		  int int_argument)
  {  return( NULL ); }

  /** replace the object's row-map (or if it's not a matrix, replace its only
      map)
   */
  static void replaceMap(T* obj, const Epetra_BlockMap& map)
  { std::cerr << "not implemented for unknown type"<<std::endl; }

  /** return new object, which is a column-permutation of srcObj */
  template<typename int_type>
  static T*
  produceColumnPermutation(TPermutation<T, int_type>* perm,
			   T* srcObj)
  { std::cerr << "not implemented for unknown type"<<std::endl; }

};//struct Perm_traits



/** A specialization of Perm_traits for the specific type
    Epetra_CrsMatrix.
 */
template<>
struct Perm_traits<Epetra_CrsMatrix> {

  /** typeName implementation */
  static const char* typeName()
  { static const char name[] = "Epetra_CrsMatrix"; return( name ); }


  /** clone implementation */
  static Epetra_CrsMatrix* clone(Epetra_CrsMatrix* example,
				 Epetra_DataAccess CV,
				 const Epetra_BlockMap& map,
				 int rowLength)
  {
    //don't need the example object currently...
    (void)example;

    //we need a Epetra_Map, rather than a Epetra_BlockMap, to create a
    //Epetra_CrsMatrix.

    const Epetra_Map* pointmap =
      dynamic_cast<const Epetra_Map*>(&map);
    if (pointmap == NULL) {
      std::cerr << "dynamic_cast<const Epetra_Map*> failed."<<std::endl;
      return(NULL);
    }

    return( new Epetra_CrsMatrix(CV, *pointmap, rowLength) );
  }


  /** replaceMap implementation */
  static void replaceMap(Epetra_CrsMatrix* mat, const Epetra_BlockMap& map)
  { mat->ReplaceRowMap(map); }

  /** return new object, which is a column-permutation of srcObj */
  template<typename int_type>
  static Epetra_CrsMatrix*
  TproduceColumnPermutation(TPermutation<Epetra_CrsMatrix, int_type>* perm,
			   Epetra_CrsMatrix* srcObj)
  {
    //First we need to export this permutation to match the column-map of the
    //object being column-permuted. (We need to have locally available all
    //elements of the permutation corresponding to the local columns of the
    //object being permuted.)

    const Epetra_Map& origColMap = srcObj->ColMap();

    TPermutation<Epetra_CrsMatrix, int_type>* colperm =
      new TPermutation<Epetra_CrsMatrix, int_type>(origColMap);
    colperm->PutValue(0);

    Epetra_Export p_exporter(perm->Map(), origColMap);
    colperm->Export(*perm, p_exporter, Add);

    const Epetra_Map& origRowMap = srcObj->RowMap();
    int numMyRows = origRowMap.NumMyElements();
    int_type* myGlobalRows = 0;
    origRowMap.MyGlobalElementsPtr(myGlobalRows);

    //Create the new object, giving it the same map as the original object.

    Epetra_CrsMatrix* result = new Epetra_CrsMatrix(Copy, origRowMap, 1);

    for(int i=0; i<numMyRows; ++i) {
      int_type globalRow = myGlobalRows[i];
      int len = srcObj->NumGlobalEntries(globalRow);

      int numIndices;
      double* src_values = new double[len];
      int_type* src_indices = new int_type[len];
      int err = srcObj->ExtractGlobalRowCopy(globalRow, len, numIndices,
					     src_values, src_indices);
      if (err < 0 || numIndices != len) {
	std::cerr<<"Perm_traits<CrsMatrix>::produceColumnPermutation err("<<err<<") row "
	    <<globalRow<<", len "<<len<<", numIndices "<<numIndices<<std::endl;
      }

      int_type* pindices = new int_type[len];

      const Epetra_BlockMap& pmap = colperm->Map();
      int_type* p = colperm->Values();

      for(int j=0; j<len; ++j) {
	int_type old_col = src_indices[j];

	int lid = pmap.LID(old_col);
	if (lid<0) {
	  std::cerr << "Perm_traits<CrsMatrix>::permuteColumnIndices GID("<<old_col
	       <<") not found"<<std::endl;
	  break;
	}

	pindices[j] = p[lid];
      }

      err = result->InsertGlobalValues(globalRow, len, src_values, pindices);
      if (err < 0) {
	std::cerr << "Perm_traits<CrsMatrix>::permuteColumnIndices err("<<err
	     <<") row "<<globalRow<<std::endl;
      }

      delete [] pindices;
      delete [] src_indices;
      delete [] src_values;
    }

    result->FillComplete();

    delete colperm;

    return(result);
  }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  /** return new object, which is a column-permutation of srcObj */
  static Epetra_CrsMatrix*
  produceColumnPermutation(TPermutation<Epetra_CrsMatrix, int>* perm,
			   Epetra_CrsMatrix* srcObj)
  {
    return TproduceColumnPermutation<int>(perm, srcObj);
  }
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  /** return new object, which is a column-permutation of srcObj */
  static Epetra_CrsMatrix*
  produceColumnPermutation(TPermutation<Epetra_CrsMatrix, long long>* perm,
			   Epetra_CrsMatrix* srcObj)
  {
    return TproduceColumnPermutation<long long>(perm, srcObj);
  }
#endif
};//struct Perm_traits<Epetra_CrsMatrix>



/** A specialization of Perm_traits for the specific type
    Epetra_CrsGraph.
 */
template<>
struct Perm_traits<Epetra_CrsGraph> {

  /** typeName implementation */
  static const char* typeName()
  { static const char name[] = "Epetra_CrsGraph"; return( name ); }


  /** clone implementation */
  static Epetra_CrsGraph* clone(Epetra_CrsGraph* example,
				Epetra_DataAccess CV,
				const Epetra_BlockMap& map,
				int rowLength)
  {
    //don't need the example object currently...
    (void)example;

    return( new Epetra_CrsGraph(CV, map, rowLength) );
  }


  /** replaceMap implementation */
  static void replaceMap(Epetra_CrsGraph* graph, const Epetra_BlockMap& map)
  { graph->ReplaceRowMap(map); }

  /** return new object which is a column-permutation of srcObj */
  template<typename int_type>
  static Epetra_CrsGraph*
  TproduceColumnPermutation(TPermutation<Epetra_CrsGraph, int_type>* perm,
			   Epetra_CrsGraph* srcObj)
  {
    //First we need to export this permutation to match the column-map of the
    //object being column-permuted. (We need to have locally available all
    //elements of the permutation corresponding to the local columns of the
    //object being permuted.)

    const Epetra_BlockMap& origColMap = srcObj->ColMap();

    TPermutation<Epetra_CrsGraph, int_type>* colperm =
      new TPermutation<Epetra_CrsGraph, int_type>(origColMap);
    colperm->PutValue(0);

    Epetra_Export p_exporter(perm->Map(), origColMap);
    colperm->Export(*perm, p_exporter, Add);

    const Epetra_BlockMap& origRowMap = srcObj->RowMap();
    int numMyRows = origRowMap.NumMyElements();
    int_type* myGlobalRows = 0;
    origRowMap.MyGlobalElementsPtr(myGlobalRows);

    //Create the new object, giving it the same map as the original object.

    Epetra_CrsGraph* result = new Epetra_CrsGraph(Copy, origRowMap, 1);

    for(int i=0; i<numMyRows; ++i) {
      int_type globalRow = myGlobalRows[i];
      int len = srcObj->NumGlobalIndices(globalRow);

      int numIndices;
      int_type* src_indices = new int_type[len];
      int err = srcObj->ExtractGlobalRowCopy(globalRow, len, numIndices, src_indices);
      if (err < 0 || numIndices != len) {
	std::cerr<<"Perm_traits<CrsGraph>::produceColumnPermutation err("<<err<<") row "
	  <<globalRow<<", len "<<len<<", numIndices "<<numIndices<<std::endl;
      }

      int_type* pindices = new int_type[len];

      const Epetra_BlockMap& pmap = colperm->Map();
      int_type* p = colperm->Values();

      for(int j=0; j<len; ++j) {
	int_type old_col = src_indices[j];

	int lid = pmap.LID(old_col);
	if (lid<0) {
	  std::cerr << "Perm_traits<CrsGraph>::permuteColumnIndices GID("<<old_col
	       <<") not found"<<std::endl;
	  break;
	}

	pindices[j] = p[lid];
      }

      err = result->InsertGlobalIndices(globalRow, len, pindices);
      if (err < 0) {
	std::cerr << "Perm_traits<CrsGraph>::produceColumnPermutation err("<<err
	     <<") row "<<globalRow<<std::endl;
      }

      delete [] pindices;
      delete [] src_indices;
    }

    result->FillComplete();

    delete colperm;

    return(result);
  }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  /** return new object which is a column-permutation of srcObj */
  static Epetra_CrsGraph*
  produceColumnPermutation(TPermutation<Epetra_CrsGraph, int>* perm,
			   Epetra_CrsGraph* srcObj)
  {
    return TproduceColumnPermutation<int>(perm, srcObj);
  }
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  /** return new object which is a column-permutation of srcObj */
  static Epetra_CrsGraph*
  produceColumnPermutation(TPermutation<Epetra_CrsGraph, long long>* perm,
			   Epetra_CrsGraph* srcObj)
  {
    return TproduceColumnPermutation<long long>(perm, srcObj);
  }
#endif
};//struct Perm_traits<Epetra_CrsGraph>


/** A specialization of Perm_traits for the specific type
    Epetra_MultiVector.
 */
template<>
struct Perm_traits<Epetra_MultiVector> {

  /** typeName implementation */
  static const char* typeName()
  { static const char name[] = "Epetra_MultiVector"; return( name ); }


  /** clone implementation */
  static Epetra_MultiVector* clone(Epetra_MultiVector* example,
				   Epetra_DataAccess /* CV */,
				   const Epetra_BlockMap& map,
				   int /* numVectors */)
  {
    return( new Epetra_MultiVector(map, example->NumVectors()) );
  }


  /** replaceMap implementation */
  static void replaceMap(Epetra_MultiVector* mvec, const Epetra_BlockMap& map)
  { mvec->ReplaceMap(map); }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  /** permute column-indices within a specified row, if applicable*/
  static Epetra_MultiVector*
  produceColumnPermutation(Permutation<Epetra_MultiVector>* /* perm */,
			   Epetra_MultiVector* /* srcObj */)
  {
    std::cerr << "col-permutation not implemented for Epetra_MultiVector"<<std::endl;
    return(NULL);
  }
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  /** permute column-indices within a specified row, if applicable*/
  static Epetra_MultiVector*
  produceColumnPermutation(Permutation64<Epetra_MultiVector>* /* perm */,
			   Epetra_MultiVector* /* srcObj */)
  {
    std::cerr << "col-permutation not implemented for Epetra_MultiVector"<<std::endl;
    return(NULL);
  }
#endif
};//struct Perm_traits<Epetra_CrsGraph>


//-------------------------------------------------------------------------
//Now the method definitions for the EpetraExt::Permutation class.
//-------------------------------------------------------------------------

template<typename T, typename int_type>
TPermutation<T, int_type>::TPermutation(Epetra_DataAccess CV,
                         const Epetra_BlockMap& map,
                         int_type* permutation)
  : Epetra_GIDTypeVector<int_type>::impl(CV, map, permutation),
    newObj_(NULL),
    origObj_(NULL)
{
  if (!isTypeSupported()) {
    std::cerr << "unsupported type for permutation, aborting" << std::endl;
    abort();
  }
}

template<typename T, typename int_type>
TPermutation<T, int_type>::TPermutation(const Epetra_BlockMap& map)
  : Epetra_GIDTypeVector<int_type>::impl(map),
    newObj_(NULL),
    origObj_(NULL)
{
  if (!isTypeSupported()) {
    std::cerr << "unsupported type for permutation, aborting" << std::endl;
    abort();
  }
}

template<typename T, typename int_type>
TPermutation<T, int_type>::TPermutation(const TPermutation& src)
  : Epetra_GIDTypeVector<int_type>::impl((const typename Epetra_GIDTypeVector<int_type>::impl&)src),
    newObj_(NULL),
    origObj_(NULL)
{
  if (!isTypeSupported()) {
    std::cerr << "unsupported type for permutation, aborting" << std::endl;
    abort();
  }
}

template<typename T, typename int_type>
TPermutation<T, int_type>::~TPermutation()
{
  if (newObj_ != NULL) delete newObj_;
}

template<typename T, typename int_type>
bool TPermutation<T, int_type>::isTypeSupported()
{
  const char* type_name = Perm_traits<T>::typeName();
  if (!strcmp(type_name, "unknown")) {
    return(false);
  }

  return( true );
}

template<typename T, typename int_type>
typename TPermutation<T, int_type>::OutputRef
TPermutation<T, int_type>::operator()( typename TPermutation<T, int_type>::InputRef orig )
{
  //In this function we're going to produce a new object which is a
  //row-permutation of the input object (orig).
  //
  //Our permutation inherits IntVector, and the permutation is defined by the
  //contents of the integer vector 'p', such that if p[i] = j then row i of
  //the input object becomes row j of the permuted object.
  //
  //The permutation is accomplished by creating a map defined by the
  //permutation, then using an Epetra_Export operation to move data from the
  //input object into the permuted object.
  //
  //The permutation may be global. In other words, the rows of the object may
  //be arbitrarily rearranged, including across processors.
  //

  origObj_ = &orig;

  //The 'Map()' accessor returns Epetra_DistObject::Map() for CrsGraph and
  //CrsMatrix, which turns out to be the RowMap() for those objects. For
  //MultiVector it returns the correct object because MultiVectors only have
  //one map.

  const Epetra_BlockMap& origMap = orig.Map();

  //Create an Epetra_Map representing the permutation.

  Epetra_Map* pmap = new Epetra_Map((int_type) Epetra_DistObject::Map().NumGlobalPoints64(),
				    Epetra_DistObject::Map().NumMyPoints(),
				    Epetra_GIDTypeVector<int_type>::impl::Values(),
				    (int_type) Epetra_DistObject::Map().IndexBase64(),
				    Epetra_DistObject::Map().Comm());

  TPermutation* p = this;

  //Next check that the maps are compatible. If they aren't, we'll redistribute
  //the permutation to match the distribution of the input object.

  if (!pmap->PointSameAs(origMap)) {
    Epetra_Export p_exporter(Epetra_DistObject::Map(), origMap);
    TPermutation* newp = new TPermutation(origMap);
    newp->Export(*p, p_exporter, Add);
    p = newp;

    delete pmap;
    pmap = new Epetra_Map((int_type) p->Map().NumGlobalPoints64(),
			  p->Map().NumMyPoints(),
			  p->Values(),
			  (int_type) p->Map().IndexBase64(),
			  p->Map().Comm());
  }

  //Create the new object, initially giving it the map defined by the
  //permutation.

  newObj_ = Perm_traits<T>::clone(origObj_, Copy, *pmap, 1);

  //Create an exporter which will export data from the original object to the
  //permuted object.

  Epetra_Export exporter(origMap, *pmap);

  //Now export the original object to the permuted object.

  newObj_->Export(orig, exporter, Add);

  //Now, since the export operation moved not only row-contents but also
  //row-numbering, we need to replace the permuted row-numbering with the
  //original row-numbering. We do this by replacing the permuted map with
  //the original row-map.

  Perm_traits<T>::replaceMap(newObj_, origMap);

  delete pmap;

  if (p != this) {
    delete p; //delete "newp" created if the PointSameAs test failed above
  }

  return( *newObj_ );
}

template<typename T, typename int_type>
typename TPermutation<T, int_type>::OutputRef
TPermutation<T, int_type>::operator()( typename TPermutation<T, int_type>::InputRef orig,
			    bool column_permutation )
{
  origObj_ = &orig;
  newObj_ = NULL;

  if (!column_permutation) {
    return( operator()(orig) );
  }

  if (strcmp("Epetra_CrsMatrix", Perm_traits<T>::typeName()) &&
      strcmp("Epetra_CrsGraph", Perm_traits<T>::typeName())) {
    std::cerr << "Permutation: column-permutation only implemented for"
	 << "CrsMatrix and CrsGraph." << std::endl;
    assert(0);
  }

  newObj_ = Perm_traits<T>::produceColumnPermutation(this, &orig);

  return( *newObj_ );
}

} // namespace EpetraExt

#endif //EpetraExt_PERMUTATION_IMPL_H
