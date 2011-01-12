// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// Questions? Contact Lee Ann Riesen (lriesen@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_INPUTADAPTER_HPP_
#define _ZOLTAN2_INPUTADAPTER_HPP_

/*! \file Zoltan2_InputAdapter.hpp
    \brief The InputAdapter base class and derived classes.

*/

#include <Zoltan2_IdentifierMap.hpp>

namespace Z2
{

/*! Z2::InputAdapter
    \brief The adapter class for specific mesh, graph, or matrix sources.

    The InputAdapter presents a uniform interface to the Zoltan methods for 
    the distributed object to be partitioned, ordered, or colored.  The actual
    data may come from sources as diverse as a Tpetra_CrsMatrix, an Exodus II file,
    a view into the caller's arrays, or from legacy Zoltan1 query functions.

    The adapter maintains information about which one process owns each object in the
    input and can move the objects to new owners.

    The adapter has an IdentifierMap that maps each object to a process which owns
    it and to the value computed by Zoltan (part number, color number, order).
*/

template<typename Scalar, typename LNO , typename GNO, typename AppGID>
  class InputAdapter {

private:

  IdentifierMap<LNO, GNO, AppGID> ids;       /*!< The global IDs and maps */

public:

  /*! \name Constructors, Destructor, Copy and Assignment */

  /*! Constructor */
  InputAdapter();

  /*! Destructor */
  ~InputAdapter();

  /*! Copy Constructor */
  InputAdapter(const InputAdapter &os){
  }

  /*! Assignment operator */
  InputAdapter &operator=(const InputAdapter &os){
  }

  /*! \name The map */

  void set_id_map(IdentifierMap<LNO, GNO, AppGID> &map);

  IdentifierMap &get_id_map() { return ids;}

  /*! \name Traits of user input */

  /*! True if InputAdapter represents a matrix 
   **/
  virtual bool has_matrix();

  /*! True if InputAdapter is stored in dense matrix format.

      The function using the InputAdapter may be able to use
      it more efficiently if it know how it is stored.
   */
  virtual bool has_dense_matrix();

  /*! True if InputAdapter is stored in compressed sparse row format.

      The function using the InputAdapter may be able to use
      it more efficiently if it know how it is stored.
   */
  virtual bool has_csr_matrix();

  /*! True if InputAdapter is stored in compressed sparse column format.

      The function using the InputAdapter may be able to use
      it more efficiently if it know how it is stored.
   */
  virtual bool has_csc_matrix();

  /*! True if InputAdapter is stored as vertices and edges.
   */
  virtual bool has_graph();

  /*! True if we have vertices with coordinates.
   */
  virtual bool has_vertex_coordinates();

  /*! Returns the dimension of coordinates, or 0 if there are no coordinates.
   */
  virtual int coordinate_dimension();

  /*! Returns true if the matrix is symmetric (TODO - other important matrix properties?)
   */
  virtual bool symmetric_matrix();

  /*! Returns true if the graph is bipartite (TODO - other important graph properties?)
   */
  virtual bool bipartite_graph();

  /*! \name Access to user input for Zoltan methods */

  /*! Returns the count of objects.

     \param global_count will contain the global number of objects on return
     \param local_count will contain the local number of objects on return
   */
  virtual void get_number_of_objects(GNO &global_count, LNO &local_count);

  /*! Returns the global numbers of all the local objects.
   */
  virtual void get_global_numbers(std::vector<GNO> &gno);

  /*! Returns the coordinates of the requested local vertices.
   */
  virtual get_vertex_coordinates(std::vector<GNO> &gno, std::vector<Scalar> &x);
  virtual get_vertex_coordinates(std::vector<GNO> &gno, std::vector<Scalar> &x std::vector<Scalar> &y);
  virtual get_vertex_coordinates(std::vector<GNO> &gno, std::vector<Scalar> &x std::vector<Scalar> &y, std::vector<Scalar &z);

  /*! Returns the number of edges for each vertex 
   */
  virtual get_number_of_edges(std::vector<GNO> &gno, std::vector<int> &num_edges);

  /*! Returns the edges for each vertex 
   */
  virtual get_edges(std::vector<GNO> &gno, std::vector< std::vector<GNO> > &edgeGNO);

  /*! Returns the local and global number of hyperedges

      Each hyperedge is owned by only one process.
   */
  virtual get_number_of_hyperedges(GNO &global_count, LNO &local_count);

  /*! Returns the global numbers for the local hyperedges.
   */
  virtual void get_hyperedge_global_numbers(std::vector<GNO> &he_gno);

  /*! Returns the local hyperedges as vertex global ID lists.
   */
  virtual void get_hyperedges(std::vector<GNO> &he_gno, std::vector< std::vector<GNO> >&vtx_gno_list);

  /*! Returns the local and global number of matrix rows.

      On return, local_count can be invalid to indicate that rows are not 
      associated with processes.
   */
  virtual void get_number_of_matrix_rows(GNO &global_count, LNO &local_count);

  /*! Returns the local and global number of matrix columns.

      On return, local_count can be invalid to indicate that columns are not 
      associated with processes.
   */
  virtual void get_number_of_matrix_columns(GNO &global_count, LNO &local_count);

  /*! Returns the global numbers for the local rows.
   */
  virtual void get_matrix_row_global_numbers(std::vector<GNO> &row_gno);

  /*! Returns the global numbers for the local columns.
      associated with processes.
   */
  virtual void get_matrix_column_global_numbers(std::vector<GNO> &col_gno);

  /*! Returns the global numbers for the local columns.
      associated with processes.
   */
  virtual void get_matrix_row_non_zeros(std::vector<GNO> &row_gno, std::vector< std::vector<GNO> > &col_gno);

  virtual void get_matrix_column_non_zeros(std::vector<GNO> &col_gno, std::vector< std::vector<GNO> > &row_gno);

  /*! \name Methods to change ownership of objects and move them to the new owner */

  /*! Locally update the process owning the object.

      object_gno must be one belonging to local process.
   */
  virtual void update_process_owner(std::vector<GNO> &object_gno, std::vector<int> *procID);

  virtual void update_hyperedge_process_owner(std::vector<GNO> &object_gno, std::vector<int> *procID);

  /*! Move data to new owners.  This is a global call.
   */
  virtual void move_data(Teuchos::Comm &comm);
};

/*! Z2::Epetra_MultiVectorAdapter
    \brief The input is geometric coordinates in an Epetra_MultiVector
*/

template<typename Scalar, typename LNO , typename GNO, typename AppGID>
  class Epetra_MultiVectorAdapter : public InputAdapter<Scalar, LNO, GNO, AppGID>{

private:

  Epetra_MultiVector *mv;

public:

  /*! \name Constructors, Destructor, Copy and Assignment */

  /*! Constructor 

      Create the IdentifierMap mapping vector elements to processes and global numbers.
      This is easy for an Epetra_MultiVector, but much more work for something like a hypergraph
      provided in CSR format.
   */
  Epetra_MultiVectorAdapter(Epetra_MultiVector &emv) : mv=&emv {
    build_identifier_map(emv);
  }

  /*! Destructor */
  ~Epetra_MultiVectorAdapter();

  /*! Copy Constructor */
  Epetra_MultiVectorAdapter(const Epetra_MultiVectorAdapter &os){
  }

  /*! Assignment operator */
  Epetra_MultiVectorAdapter &operator=(const Epetra_MultiVectorAdapter &os){
  }

  /*! \name Traits of user input */

  bool has_matrix() {return false};

  /*! True if Epetra_MultiVectorAdapter is stored in dense matrix format.

      The function using the Epetra_MultiVectorAdapter may be able to use
      it more efficiently if it know how it is stored.
   */
  bool has_dense_matrix() {return false};

  /*! True if Epetra_MultiVectorAdapter is stored in compressed sparse row format.

      The function using the Epetra_MultiVectorAdapter may be able to use
      it more efficiently if it know how it is stored.
   */
  bool has_csr_matrix() {return false};

  /*! True if Epetra_MultiVectorAdapter is stored in compressed sparse column format.

      The function using the Epetra_MultiVectorAdapter may be able to use
      it more efficiently if it know how it is stored.
   */
  bool has_csc_matrix() {return false};

  /*! True if Epetra_MultiVectorAdapter is stored as vertices and edges.
   */
  bool has_graph() {return false};

  /*! True if vertices have coordinates.
   */
  bool has_vertex_coordinates() {return true;}

  /*! Returns the dimension of coordinates, or 0 if there are no coordinates.
   */
  int coordinate_dimension() { return mv->NumVectors(); }

  /*! Returns true if the matrix is symmetric (TODO - other important matrix properties?)
   */
  bool symmetric_matrix() {return false};

  /*! Returns true if the graph is bipartite (TODO - other important graph properties?)
   */
  bool bipartite_graph() {return false};

  /*! \name Access to user input for Zoltan methods */

  /*! Returns the count of objects.

     \param global_count will contain the global number of objects on return
     \param local_count will contain the local number of objects on return
   */
  void get_number_of_objects(GNO &global_count, LNO &local_count){
    global_count = mv->GlobalLength();
    local_count = mv->MyLength();
  }

  /*! Returns the global numbers of all the local objects.
   */
  void get_global_numbers(std::vector<GNO> &gno);

  /*! Returns the coordinates of the requested local vertices.
   */
  get_vertex_coordinates(std::vector<GNO> &gno, std::vector<Scalar> &x);
  get_vertex_coordinates(std::vector<GNO> &gno, std::vector<Scalar> &x std::vector<Scalar> &y);
  get_vertex_coordinates(std::vector<GNO> &gno, std::vector<Scalar> &x std::vector<Scalar> &y, std::vector<Scalar &z);

  /*! Returns the number of edges for each vertex 
   */
  get_number_of_edges(std::vector<GNO> &gno, std::vector<int> &num_edges) 
    { throw std::runtime_error("Undefined"); }

  /*! Returns the edges for each vertex 
   */
  get_edges(std::vector<GNO> &gno, std::vector< std::vector<GNO> > &edgeGNO)
    { throw std::runtime_error("Undefined"); }

  /*! Returns the local and global number of hyperedges

      Each hyperedge is owned by only one process.
   */
  get_number_of_hyperedges(GNO &global_count, LNO &local_count)
    { throw std::runtime_error("Undefined"); }

  /*! Returns the global numbers for the local hyperedges.
   */
  void get_hyperedge_global_numbers(std::vector<GNO> &he_gno)
    { throw std::runtime_error("Undefined"); }

  /*! Returns the local hyperedges as vertex global ID lists.
   */
  void get_hyperedges(std::vector<GNO> &he_gno, std::vector< std::vector<GNO> >&vtx_gno_list)
    { throw std::runtime_error("Undefined"); }

  /*! Returns the local and global number of matrix rows.

      On return, local_count can be invalid to indicate that rows are not 
      associated with processes.
   */
  void get_number_of_matrix_rows(GNO &global_count, LNO &local_count)
    { throw std::runtime_error("Undefined"); }

  /*! Returns the local and global number of matrix columns.

      On return, local_count can be invalid to indicate that columns are not 
      associated with processes.
   */
  void get_number_of_matrix_columns(GNO &global_count, LNO &local_count)
    { throw std::runtime_error("Undefined"); }

  /*! Returns the global numbers for the local rows.
   */
  void get_matrix_row_global_numbers(std::vector<GNO> &row_gno)
    { throw std::runtime_error("Undefined"); }

  /*! Returns the global numbers for the local columns.
      associated with processes.
   */
  void get_matrix_column_global_numbers(std::vector<GNO> &col_gno)
    { throw std::runtime_error("Undefined"); }

  /*! Returns the global numbers for the local columns.
      associated with processes.
   */
  void get_matrix_row_non_zeros(std::vector<GNO> &row_gno, std::vector< std::vector<GNO> > &col_gno)
    { throw std::runtime_error("Undefined"); }

  void get_matrix_column_non_zeros(std::vector<GNO> &col_gno, std::vector< std::vector<GNO> > &row_gno)
    { throw std::runtime_error("Undefined"); }

  /*! \name Methods to change ownership of objects and move them to the new owner */

  /*! Locally update the process owning the object.

      object_gno must be one belonging to local process.
   */
  void update_process_owner(std::vector<GNO> &object_gno, std::vector<int> *procID);

  void update_hyperedge_process_owner(std::vector<GNO> &object_gno, std::vector<int> *procID)
    { throw std::runtime_error("Undefined"); }

  /*! Move data to new owners.  This is a global call.
   */
  void move_data(Teuchos::Comm &comm);
};

} // namespace Z2

#endif /* _ZOLTAN2_INPUTADAPTER_HPP_ */
