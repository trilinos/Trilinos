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
Foundation, Inc, 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

************************************************************************
*/
//@HEADER

#ifndef _Isorropia_EpetraRedistributor_hpp_
#define _Isorropia_EpetraRedistributor_hpp_

#include <Isorropia_Redistributor.hpp>
#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

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
class Epetra_SrcDistObject;
class Epetra_DistObject;

namespace Isorropia {

namespace Epetra {
  class Partitioner;

/** Class which is constructed with a Partitioner instance, and
     provides several methods for redistributing Epetra objects
     given the partitioning computed by the Partitioner object.
*/

class Redistributor : public Isorropia::Redistributor {
public:
  /** This constructor calls the Isorropia::Epetra::Partitioner::partition
      method on the @c partitioner if it has not already been called.
 
      \param[in] partitioner this input partitioner determines the new partitioning
            to be created when Isorropia::Epetra::Redistributor::redistribute is called
   */
  Redistributor(Teuchos::RefCountPtr<Isorropia::Epetra::Partitioner> partitioner);

  /** Destructor
   */
  virtual ~Redistributor();

  /** Method to redistribute a Epetra_SrcDistObject into a
      Epetra_DistObject. The caller is required to have constructed
      the target object using the correct target map.
  */
  void redistribute(const Epetra_SrcDistObject& src,
		    Epetra_DistObject& target);

  /** Method to accept a Epetra_CrsGraph object, and
      return a redistributed Epetra_CrsGraph object.

      \param[in] input_graph the graph for which we want a new graph that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param[in] callFillComplete The new graph is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_graph has different domain and range maps, in which case
      the range map is preserved.  By default callFillComplete is @c true.

      \return a reference counted pointer to the new redistributed graph 
  */
  Teuchos::RefCountPtr<Epetra_CrsGraph>
     redistribute(const Epetra_CrsGraph& input_graph, bool callFillComplete= true);

  /** Method to accept a Epetra_CrsMatrix object, and
      return a redistributed Epetra_CrsMatrix object.

      \param[in] input_matrix the matrix for which we want a new matrix that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param[in] callFillComplete The new matrix is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_matrix has different domain and range maps, in which case
      the range map is preserved.  By default callFillComplete is @c true.

      \return a reference counted pointer to the new redistributed matrix
  */
  Teuchos::RefCountPtr<Epetra_CrsMatrix>
     redistribute(const Epetra_CrsMatrix& input_matrix, bool callFillComplete= true);

  /** Method to accept a Epetra_RowMatrix object, and
      return a redistributed Epetra_CrsMatrix object.

      \param[in] input_matrix the row matrix for which we want a new matrix that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param[in] callFillComplete The new matrix is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_matrix has different domain and range maps, in which case
      the range map is preserved.  By default callFillComplete is @c true.

      \return a reference counted pointer to the new redistributed matrix
  */
  Teuchos::RefCountPtr<Epetra_CrsMatrix>
     redistribute(const Epetra_RowMatrix& input_matrix, bool callFillComplete= true);

  /** Method to accept a Epetra_Vector object, and
      return a redistributed Epetra_Vector object.

      \param[in] input_vector the vector for which we want a new vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \return a reference counted pointer to the new redistributed vector
  */
  Teuchos::RefCountPtr<Epetra_Vector>
     redistribute(const Epetra_Vector& input_vector);

  /** Method to accept a Epetra_MultiVector object, and
      return a redistributed Epetra_MultiVector object.

      \param[in] input_vector the multi vector for which we want a new multi vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \return a reference counted pointer to the new redistributed multi vector
  */
  Teuchos::RefCountPtr<Epetra_MultiVector>  
     redistribute(const Epetra_MultiVector& input_vector);

  /** Reverse redistribute an Epetra_Vector.

      \param[in] input_vector a vector that is distributed according to the partitioner that was used to create this Redistributor

      \param[out] output_vector a copy of the @c input_vector which has been redistributed according
                    to the reverse of the partitioner that was used to create this Redistributor
  */
  void
     redistribute_reverse(const Epetra_Vector& input_vector, Epetra_Vector& output_vector);

  /** Reverse redistribute an Epetra_MultiVector.

      \param[in] input_vector a multi vector that is distributed according to the partitioner that was used to create this Redistributor

      \param[out] output_vector a copy of the @c input_vector which has been redistributed according
                    to the reverse of the partitioner that was used to create this Redistributor

  */
  void
     redistribute_reverse(const Epetra_MultiVector& input_vector, Epetra_MultiVector& output_vector);
private:
  /** Create an importer object to be used in the redistribution

      \param[in] src_map the map describing the pattern of the import operation
   */
  void create_importer(const Epetra_BlockMap& src_map);

  Teuchos::RefCountPtr<Isorropia::Epetra::Partitioner> partitioner_;
  Teuchos::RefCountPtr<Epetra_Import> importer_;
  Teuchos::RefCountPtr<Epetra_Map> target_map_;

  bool created_importer_;

}; //class Redistributor

}//namespace Epetra

}//namespace Isorropia

#endif //HAVE_EPETRA

#endif

