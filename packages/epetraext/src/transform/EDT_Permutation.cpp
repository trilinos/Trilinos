//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#include <Epetra_ConfigDefs.h>
#include <EDT_Permutation.h>
#include <Epetra_Export.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>

//
//Explicit template instantiation requests
//
template class EpetraExt::Permutation<Epetra_MultiVector>;
template class EpetraExt::Permutation<Epetra_CrsMatrix>;
template class EpetraExt::Permutation<Epetra_CrsGraph>;

namespace EpetraExt {

/** Define some traits to make it easier to deal with template-parameters which
   are objects to be permuted. Given a template parameter, we'll want to know
   what type it is, be able to construct an instance of it, and be able to
   replace its row-map.

   First the default definition, which catches all types "T", followed by some
   specializations for anticipated types. Any type other than the types specifically
   anticipated will be handled by this default definition, allowing the Permutation
   class to abort or return NULL where appropriate.

   We define these trait structs in this file rather than in a separate file
   in an attempt to avoid some template-instantiation complications...
*/
template<class T>
struct Perm_traits {
  /** return a string name for the object type */
  static const char* typeName()
  { static const char name[] = "unknown"; return( name ); }

  /** clone function accepts an example of the object being cloned, and enough
      constructor arguments to be able to create any of these: CrsMatrix, CrsGraph,
      Vector, MultiVector.   And probably more later...

      Why is an example object needed? For instance, if a MultiVector is created, we
      may want to know how many vectors it should contain...
  */
  static T* clone(T* example,
		  Epetra_DataAccess CV,
		  const Epetra_BlockMap& map,
		  int int_argument)
  { return( NULL ); }

  /** replace the object's row-map (or if it's not a matrix, replace its only map)
   */
  static void replaceMap(T* obj, const Epetra_BlockMap& map)
  { cerr << "not implemented for unknown type"<<endl; }
};//struct Perm_traits



/** A specialization of Perm_traits for the specific anticipated type Epetra_CrsMatrix.
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

};//struct Perm_traits<Epetra_CrsMatrix>



/** A specialization of Perm_traits for the specific anticipated type Epetra_CrsGraph.
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

};//struct Perm_traits<Epetra_CrsGraph>


/** A specialization of Perm_traits for the specific anticipated type Epetra_MultiVector.
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
typename Permutation<T>::NewTypeRef Permutation<T>::operator()( OriginalTypeRef orig )
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

  //The 'Map()' accessor returns Epetra_DistObject::Map() for CrsGraph and CrsMatrix,
  //which turns out to be the RowMap() for those objects. For MultiVector it returns
  //the correct object because MultiVectors only have one map.

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

  return( *newObj_ );
}

template<typename T>
typename Permutation<T>::NewTypeRef Permutation<T>::operator()( OriginalTypeRef orig,
						 bool column_permutation )
{
  origObj_ = &orig;
  newObj_ = NULL;

  if (!column_permutation) {
    return( operator()(orig) );
  }
  else {
    cerr <<"Permutation: column-permutations aren't implemented yet"<<endl;
    assert(0);
  }

  return( *newObj_ );
}

} // namespace EpetraExt

