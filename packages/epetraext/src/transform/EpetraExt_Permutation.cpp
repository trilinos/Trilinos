//@HEADER
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
//@HEADER
#ifndef EpetraExt_PERMUTATION_CPP
#define EpetraExt_PERMUTATION_CPP

#include <EpetraExt_ConfigDefs.h>

#include <EpetraExt_Permutation.h>

#include <Epetra_Export.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>

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
  /** return a string name for the object type */
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
  { return( NULL ); }

  /** replace the object's row-map (or if it's not a matrix, replace its only
      map)
   */
  static void replaceMap(T* obj, const Epetra_BlockMap& map)
  { cerr << "not implemented for unknown type"<<endl; }

  /** return new object, which is a column-permutation of srcObj */
  static T*
  produceColumnPermutation(Permutation<T>* perm,
			   T* srcObj)
  { cerr << "not implemented for unknown type"<<endl; }

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
      cerr << "dynamic_cast<const Epetra_Map*> failed."<<endl;
      return(NULL);
    }

    return( new Epetra_CrsMatrix(CV, *pointmap, rowLength) );
  }


  /** replaceMap implementation */
  static void replaceMap(Epetra_CrsMatrix* mat, const Epetra_BlockMap& map)
  { mat->ReplaceRowMap(map); }

  /** return new object, which is a column-permutation of srcObj */
  static Epetra_CrsMatrix*
  produceColumnPermutation(Permutation<Epetra_CrsMatrix>* perm,
			   Epetra_CrsMatrix* srcObj)
  {
    //First we need to export this permutation to match the column-map of the
    //object being column-permuted. (We need to have locally available all
    //elements of the permutation corresponding to the local columns of the
    //object being permuted.)

    const Epetra_Map& origColMap = srcObj->ColMap();

    Permutation<Epetra_CrsMatrix>* colperm =
      new Permutation<Epetra_CrsMatrix>(origColMap);
    colperm->PutValue(0);

    Epetra_Export p_exporter(perm->Map(), origColMap);
    colperm->Export(*perm, p_exporter, Add);

    const Epetra_Map& origRowMap = srcObj->RowMap();
    int numMyRows = origRowMap.NumMyElements();
    const int* myGlobalRows = origRowMap.MyGlobalElements();

    //Create the new object, giving it the same map as the original object.

    Epetra_CrsMatrix* result = new Epetra_CrsMatrix(Copy, origRowMap, 1);

    for(int i=0; i<numMyRows; ++i) {
      int globalRow = myGlobalRows[i];
      int len = srcObj->NumGlobalEntries(globalRow);

      int numIndices;
      double* src_values = new double[len];
      int* src_indices = new int[len];
      int err = srcObj->ExtractGlobalRowCopy(globalRow, len, numIndices,
					     src_values, src_indices);
      if (err < 0 || numIndices != len) {
	cerr<<"Perm_traits<CrsMatrix>::produceColumnPermutation err("<<err<<") row "
	    <<globalRow<<", len "<<len<<", numIndices "<<numIndices<<endl;
      }

      int* pindices = new int[len];

      const Epetra_BlockMap& pmap = colperm->Map();
      int* p = colperm->Values();

      for(int j=0; j<len; ++j) {
	int old_col = src_indices[j];

	int lid = pmap.LID(old_col);
	if (lid<0) {
	  cerr << "Perm_traits<CrsMatrix>::permuteColumnIndices GID("<<old_col
	       <<") not found"<<endl;
	  break;
	}

	pindices[j] = p[lid];
      }

      err = result->InsertGlobalValues(globalRow, len, src_values, pindices);
      if (err < 0) {
	cerr << "Perm_traits<CrsMatrix>::permuteColumnIndices err("<<err
	     <<") row "<<globalRow<<endl;
      }

      delete [] pindices;
      delete [] src_indices;
      delete [] src_values;
    }

    result->FillComplete();

    delete colperm;

    return(result);
  }

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
  static Epetra_CrsGraph*
  produceColumnPermutation(Permutation<Epetra_CrsGraph>* perm,
			   Epetra_CrsGraph* srcObj)
  {
    //First we need to export this permutation to match the column-map of the
    //object being column-permuted. (We need to have locally available all
    //elements of the permutation corresponding to the local columns of the
    //object being permuted.)

    const Epetra_BlockMap& origColMap = srcObj->ColMap();

    Permutation<Epetra_CrsGraph>* colperm =
      new Permutation<Epetra_CrsGraph>(origColMap);
    colperm->PutValue(0);

    Epetra_Export p_exporter(perm->Map(), origColMap);
    colperm->Export(*perm, p_exporter, Add);

    const Epetra_BlockMap& origRowMap = srcObj->RowMap();
    int numMyRows = origRowMap.NumMyElements();
    const int* myGlobalRows = origRowMap.MyGlobalElements();

    //Create the new object, giving it the same map as the original object.

    Epetra_CrsGraph* result = new Epetra_CrsGraph(Copy, origRowMap, 1);

    for(int i=0; i<numMyRows; ++i) {
      int globalRow = myGlobalRows[i];
      int len = srcObj->NumGlobalIndices(globalRow);

      int numIndices;
      int* src_indices = new int[len];
      int err = srcObj->ExtractGlobalRowCopy(globalRow, len, numIndices, src_indices);
      if (err < 0 || numIndices != len) {
	cerr<<"Perm_traits<CrsGraph>::produceColumnPermutation err("<<err<<") row "
	  <<globalRow<<", len "<<len<<", numIndices "<<numIndices<<endl;
      }

      int* pindices = new int[len];

      const Epetra_BlockMap& pmap = colperm->Map();
      int* p = colperm->Values();

      for(int j=0; j<len; ++j) {
	int old_col = src_indices[j];

	int lid = pmap.LID(old_col);
	if (lid<0) {
	  cerr << "Perm_traits<CrsGraph>::permuteColumnIndices GID("<<old_col
	       <<") not found"<<endl;
	  break;
	}

	pindices[j] = p[lid];
      }

      err = result->InsertGlobalIndices(globalRow, len, pindices);
      if (err < 0) {
	cerr << "Perm_traits<CrsGraph>::produceColumnPermutation err("<<err
	     <<") row "<<globalRow<<endl;
      }

      delete [] pindices;
      delete [] src_indices;
    }

    result->FillComplete();

    delete colperm;

    return(result);
  }

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
				   Epetra_DataAccess CV,
				   const Epetra_BlockMap& map,
				   int numVectors)
  {
    return( new Epetra_MultiVector(map, example->NumVectors()) );
  }


  /** replaceMap implementation */
  static void replaceMap(Epetra_MultiVector* mvec, const Epetra_BlockMap& map)
  { mvec->ReplaceMap(map); }

  /** permute column-indices within a specified row, if applicable*/
  static Epetra_MultiVector*
  produceColumnPermutation(Permutation<Epetra_MultiVector>* perm,
			   Epetra_MultiVector* srcObj)
  {
    cerr << "col-permutation not implemented for Epetra_MultiVector"<<endl;
    return(NULL);
  }

};//struct Perm_traits<Epetra_CrsGraph>


//-------------------------------------------------------------------------
//Now the method definitions for the EpetraExt::Permutation class.
//-------------------------------------------------------------------------

template<typename T>
Permutation<T>::Permutation(Epetra_DataAccess CV,
                         const Epetra_BlockMap& map,
                         int* permutation)
  : Epetra_IntVector(CV, map, permutation),
    newObj_(NULL)
{
  if (!isTypeSupported()) {
    cerr << "unsupported type for permutation, aborting" << endl;
    abort();
  }
}

template<typename T>
Permutation<T>::Permutation(const Epetra_BlockMap& map)
  : Epetra_IntVector(map),
    newObj_(NULL)
{
  if (!isTypeSupported()) {
    cerr << "unsupported type for permutation, aborting" << endl;
    abort();
  }
}

template<typename T>
Permutation<T>::Permutation(const Permutation& src)
  : Epetra_IntVector((const Epetra_IntVector&)src),
    newObj_(NULL)
{
  if (!isTypeSupported()) {
    cerr << "unsupported type for permutation, aborting" << endl;
    abort();
  }
}

template<typename T>
Permutation<T>::~Permutation()
{
  if (newObj_ != NULL) delete newObj_;
}

template<typename T>
bool Permutation<T>::isTypeSupported()
{
  const char* type_name = Perm_traits<T>::typeName();
  if (!strcmp(type_name, "unknown")) {
    return(false);
  }

  return( true );
}

template<typename T>
typename Permutation<T>::OutputRef
Permutation<T>::operator()( typename Permutation<T>::InputRef orig )
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

  Epetra_Map* pmap = new Epetra_Map(Map().NumGlobalPoints(),
				    Map().NumMyPoints(),
				    Values(),
				    Map().IndexBase(),
				    Map().Comm());

  Permutation* p = this;

  //Next check that the maps are compatible. If they aren't, we'll redistribute
  //the permutation to match the distribution of the input object.

  if (!pmap->PointSameAs(origMap)) {
    Epetra_Export p_exporter(Map(), origMap);
    Permutation* newp = new Permutation(origMap);
    newp->Export(*p, p_exporter, Add);
    p = newp;

    delete pmap;
    pmap = new Epetra_Map(p->Map().NumGlobalPoints(),
			  p->Map().NumMyPoints(),
			  p->Values(),
			  p->Map().IndexBase(),
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

template<typename T>
typename Permutation<T>::OutputRef
Permutation<T>::operator()( typename Permutation<T>::InputRef orig,
			    bool column_permutation )
{
  origObj_ = &orig;
  newObj_ = NULL;

  if (!column_permutation) {
    return( operator()(orig) );
  }

  if (strcmp("Epetra_CrsMatrix", Perm_traits<T>::typeName()) &&
      strcmp("Epetra_CrsGraph", Perm_traits<T>::typeName())) {
    cerr << "Permutation: column-permutation only implemented for"
	 << "CrsMatrix and CrsGraph." << endl;
    assert(0);
  }

  newObj_ = Perm_traits<T>::produceColumnPermutation(this, &orig);

  return( *newObj_ );
}

} // namespace EpetraExt

//
//Explicit template instantiation requests
//
template class EpetraExt::Permutation<Epetra_MultiVector>;
template class EpetraExt::Permutation<Epetra_CrsMatrix>;
template class EpetraExt::Permutation<Epetra_CrsGraph>;

#endif //EpetraExt_PERMUTATION_CPP
