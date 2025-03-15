//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
//@HEADER

#ifndef _Isorropia_EpetraRedistributor_hpp_
#define _Isorropia_EpetraRedistributor_hpp_

#include <Isorropia_Redistributor.hpp>
#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
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

/** @ingroup partitioning_grp partitioning_rcp_grp partitioning_ptr_grp
     Class which is constructed with a Partitioner instance, and
     provides several methods for redistributing Epetra objects
     given the partitioning computed by the Partitioner object.
*/

class Redistributor : public Isorropia::Redistributor {
public:

  /** @ingroup partitioning_rcp_grp
      This constructor calls the Isorropia::Epetra::Partitioner::partition
      method on the @c partitioner if it has not already been called.
 
      \param partitioner (in) this input partitioner determines the new partitioning
            to be created when Isorropia::Epetra::Redistributor::redistribute is called
   */
  Redistributor(Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner);

  /** @ingroup partitioning_rcp_grp
      This constructor sets the target map for the redistribution.

      \param target_map (in) this input map determines the new matrices/vectors
      to be created when Isorropia::Epetra::Redistributor::redistribute is
      called
   */
  Redistributor(Teuchos::RCP<Epetra_Map> target_map);

  /** @ingroup partitioning_ptr_grp
      This constructor calls the Isorropia::Epetra::Partitioner::partition
      method on the @c partitioner if it has not already been called.
 
      \param partitioner (in) this input partitioner determines the new partitioning
            to be created when Isorropia::Epetra::Redistributor::redistribute is called
   */
  Redistributor(Isorropia::Epetra::Partitioner *partitioner);

  /** @ingroup partitioning_ptr_grp
      This constructor sets the target map for the redistribution.

      \param target_map (in) this input map determines the new matrices/vectors
      to be created when Isorropia::Epetra::Redistributor::redistribute is
      called
   */
  Redistributor(Epetra_Map *target_map);

  /** 
       Destructor
   */
  virtual ~Redistributor();

  /** @ingroup partitioning_grp
      Method to redistribute a Epetra_SrcDistObject into a
      Epetra_DistObject. The caller is required to have constructed
      the target object using the correct target map.
  */
  void redistribute(const Epetra_SrcDistObject& src,
		    Epetra_DistObject& target);

  /** @ingroup partitioning_rcp_grp
      Method to accept a Epetra_CrsGraph object, and
      return a redistributed Epetra_CrsGraph object.

      \param input_graph (in) the graph for which we want a new graph that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param callFillComplete (in) The new graph is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_graph has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.

      \return a reference counted pointer to the new redistributed graph 
  */
  Teuchos::RCP<Epetra_CrsGraph>
     redistribute(const Epetra_CrsGraph& input_graph, bool callFillComplete= true);

  /** @ingroup partitioning_ptr_grp
      Method to accept a Epetra_CrsGraph object, and
      return a redistributed Epetra_CrsGraph object.

      \param input_graph (in) the graph for which we want a new graph that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param outputGraphPtr (out) pointer to the new redistributed graph 

      \param callFillComplete (in) The new graph is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_graph has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.

  */
  void redistribute(const Epetra_CrsGraph& input_graph, Epetra_CrsGraph * &outputGraphPtr, bool callFillComplete= true);

  /** @ingroup partitioning_rcp_grp
      Method to accept a Epetra_CrsMatrix object, and
      return a redistributed Epetra_CrsMatrix object.

      \param input_matrix (in) the matrix for which we want a new matrix that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param callFillComplete (in) The new matrix is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_matrix has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.

      \return a reference counted pointer to the new redistributed matrix
  */
  Teuchos::RCP<Epetra_CrsMatrix>
     redistribute(const Epetra_CrsMatrix& input_matrix, bool callFillComplete= true);

  /** @ingroup partitioning_ptr_grp
      Method to accept a Epetra_CrsMatrix object, and
      return a redistributed Epetra_CrsMatrix object.

      \param inputMatrix (in) the matrix for which we want a new matrix that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param outputMatrix (out) pointer to the new redistributed matrix

      \param callFillComplete (in) The new matrix is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_matrix has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.
  */

  void redistribute(const Epetra_CrsMatrix& inputMatrix, Epetra_CrsMatrix * &outputMatrix, bool callFillComplete= true);

  /** @ingroup partitioning_rcp_grp
      Method to accept a Epetra_RowMatrix object, and
      return a redistributed Epetra_CrsMatrix object.

      \param input_matrix (in) the row matrix for which we want a new matrix that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param callFillComplete (in) The new matrix is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_matrix has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.

      \return a reference counted pointer to the new redistributed matrix
  */
  Teuchos::RCP<Epetra_CrsMatrix>
     redistribute(const Epetra_RowMatrix& input_matrix, bool callFillComplete= true);


  /** @ingroup partitioning_ptr_grp
      Method to accept a Epetra_RowMatrix object, and
      return a redistributed Epetra_CrsMatrix object.

      \param inputMatrix (in) the row matrix for which we want a new matrix that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param outputMatrix (out) pointer to the new redistributed matrix

      \param callFillComplete (in) The new matrix is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_matrix has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.
  */
     void 
     redistribute(const Epetra_RowMatrix& inputMatrix, Epetra_CrsMatrix * &outputMatrix, bool callFillComplete= true);

  /** @ingroup partitioning_rcp_grp
      Method to accept a Epetra_Vector object, and
      return a redistributed Epetra_Vector object.

      \param input_vector (in) the vector for which we want a new vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \return a reference counted pointer to the new redistributed vector
  */
  Teuchos::RCP<Epetra_Vector>
     redistribute(const Epetra_Vector& input_vector);

  /** @ingroup partitioning_ptr_grp
      Method to accept a Epetra_Vector object, and
      return a redistributed Epetra_Vector object.

      \param inputVector (in) the vector for which we want a new vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param outputVector (out) pointer to the new redistributed vector
  */
  void
  redistribute(const Epetra_Vector& inputVector, Epetra_Vector * &outputVector);

  /** @ingroup partitioning_rcp_grp 
      Method to accept a Epetra_MultiVector object, and
      return a redistributed Epetra_MultiVector object.

      \param input_vector (in) the multi vector for which we want a new multi vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \return a reference counted pointer to the new redistributed multi vector

  */
  Teuchos::RCP<Epetra_MultiVector>  
     redistribute(const Epetra_MultiVector& input_vector);


  /** @ingroup partitioning_ptr_grp 
      Method to accept a Epetra_MultiVector object, and
      return a redistributed Epetra_MultiVector object.

      \param inputVector (in) the multi vector for which we want a new multi vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param outputVector (out) a reference counted pointer to the new redistributed multi vector

  */
  void 
  redistribute(const Epetra_MultiVector& inputVector, Epetra_MultiVector * &outputVector);

  /** @ingroup partitioning_grp 
      Reverse redistribute an Epetra_Vector.

      \param input_vector (in) a vector that is distributed according to the partitioner that was used to create this Redistributor

      \param output_vector (out) a copy of the @c input_vector which has been redistributed according
                    to the reverse of the partitioner that was used to create this Redistributor

  */
  void
     redistribute_reverse(const Epetra_Vector& input_vector, Epetra_Vector& output_vector);

  /** @ingroup partitioning_grp 
      Reverse redistribute an Epetra_MultiVector.

      \param input_vector (in) a multi vector that is distributed according to the partitioner that was used to create this Redistributor

      \param output_vector (out) a copy of the @c input_vector which has been redistributed according
                    to the reverse of the partitioner that was used to create this Redistributor
  */
  void
     redistribute_reverse(const Epetra_MultiVector& input_vector, Epetra_MultiVector& output_vector);

  Epetra_Import &get_importer() { return *importer_;}

private:
  /** @ingroup partitioning_grp
      Create an importer object to be used in the redistribution
      \param src_map (in) the map describing the pattern of the import operation
   */
  void create_importer(const Epetra_BlockMap& src_map);

  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner_;
  Teuchos::RCP<Epetra_Import> importer_;
  Teuchos::RCP<Epetra_Map> target_map_;

}; //class Redistributor

}//namespace Epetra

}//namespace Isorropia

#endif //HAVE_EPETRA

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

