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

#ifndef _ZOLTAN2_OBJECTSOURCE_HPP_
#define _ZOLTAN2_OBJECTSOURCE_HPP_

/*! \file Zoltan2_ObjectSource.hpp
    \brief The ObjectSource base class and derived classes.

  Detailed description of ObjectSource.
*/

/*
    TODO: allow request for the source to be randomized, or for 2D distribution
           of hypergraph.

    TODO: allow app to specify that not all info for source is needed.  So if graph
        vertices have coordinates, but we won't be using them, the source won't
        read them in.

    TODO: use Teuchos Array RCPs for persistent objects - views of our data that we
           supply to the caller

    TODO: This class does not support Zoltan1's global ID as is an arbitrary number
            of integers unless someone writes a class to represent that and
            uses it as the AppGID

    Weights are not inherent in the source.  They are part of the objective function.
    Same goes for number of parts and part sizes.  Fixed vertices are part of the
    constraints.  The ObjectSource is just the source, not the problem to be solved.

 ObjectSource traits:
    has_matrix()
    has_dense_matrix()
    has_csr_matrix()
    has_csc_matrix()
    has_graph()
    has_vertices()
    has_vertex_coordinates()
    coordinate_dimension()                 , etc.

 ObjectSource methods:
    get_number_of_vertices(),etc.  All methods to get local and global counts
    get_vertices()                 All methods to get local objects
    gid_to_gno()        All methods to convert back and forth between app GIDs 
    gid_to_lid()         and internal consecutive GNOs and LIDs
*/

class Tpetra::MultiVector;

namespace Zoltan2
{

/*! Zoltan2::ObjectSource
    \brief The adapter class for specific mesh, graph, hypergraph or matrix sources.

    The ObjectSource presents a uniform interface to the Zoltan methods for 
    the distributed object to be partitioned, ordered, or colored.  The actual
    data may come from sources as diverse as a Tpetra_CrsMatrix, an Exodus II file,
    a view into the caller's arrays, or from legacy Zoltan1 query functions.
*/

template<class Scalar, class LNO , class GNO, class AppGID>
  class ObjectSource {

private:

public:

  /*! Constructor */
  ObjectSource(){}

  /*! Destructor */
  ~ObjectSource(){}

  /*! Copy Constructor */
  ObjectSource(const ObjectSource &os){
  }

  /*! Assignment operator */
  ObjectSource &operator=(const ObjectSource &os){
  }

  /*! True if ObjectSource represents a matrix 
        TODO - declare const functions
   **/
  virtual bool has_matrix(){}

  /*! True if ObjectSource is stored in dense matrix format.

      The function using the ObjectSource may be able to use
      it more efficiently if it know how it is stored.
   */
  virtual bool has_dense_matrix(){}

  /*! True if ObjectSource is stored in compressed sparse row format.

      The function using the ObjectSource may be able to use
      it more efficiently if it know how it is stored.
   */
  virtual bool has_csr_matrix(){}

  /*! True if ObjectSource is stored in compressed sparse column format.

      The function using the ObjectSource may be able to use
      it more efficiently if it know how it is stored.
   */
  virtual bool has_csc_matrix(){}

  /*! True if ObjectSource is stored as vertices and edges.
   */
  virtual bool has_graph(){}

  /*! True if ObjectSource has a vertex list.
   */
  virtual bool has_vertices(){}

  /*! True if vertices have coordinates.
   */
  virtual bool has_vertex_coordinates(){}

  /*! Returns the dimension of coordinates, or 0 if there are no coordinates.
   */
  virtual int coordinate_dimension(){}

  /*! Returns true if the matrix is symmetric (TODO - other important matrix properties?)
   */
  virtual bool symmetric_matrix(){}

  /*! Returns true if the graph is bipartite (TODO - other important graph properties?)
   */
  virtual bool bipartite_graph(){}

  /*! Returns the global and local number of vertices.

      Vertices may also be abstract objects known only by an ID.  We call them
      "vertices" because the normal use case is geometric or a graph.

     \param global_count will contain the global number of vertices on return
     \param global_count will contain the local number of vertices on return
   */
  virtual void number_of_vertices(GNO &global_count, LNO &local_count){}

  /*! Returns the applications's object global IDs.

       TODO: we should be able to do partitioning without having global
             IDs from the application when there's no connectivity.
       The Zoltan methods shouldn't call this - they will work with GNOs.
       But when creating or reading the Result object, they need the AppGIDs.

     \param gids will contain a pointer to the application's global IDs on return 
   */
  virtual int caller_gid_list(AppGID * &gids){}

  /*! Returns the library's global numbering for the objects.

       Info returned by other methods is in the order corresponding
       to the order of the gnos in this list.

     \param gnos will contain a pointer to the ObjectSource's global numbers on return
   */
  virtual int global_number_list(GNO *gnos){}

  /*! Returns the coordinates of the vertices.
   */
  virtual int vertex_coordinates(Tpetra_MultiVector<Scalar, LNO, GNO> &coords){}

  // TODO: add methods to get matrix and graph info.

  /*! Return the local number corresponding to the global number.
    
      Local numbers are consecutive numbers beginning at 0 for each process.
      Global numbers are consecutive across the application and begin at 0.
    
      TODO - local/global/application ID management is likely to be the same across
        all sources.  Maybe this should be done in the base class.
    
      TODO - for efficiency, list versions of these calls.
   */
  virtual LNO gno_to_lno(GNO gno);

  /*! Return the global number corresponding to the local number. */

  virtual GNO lno_to_gno(LNO lno);

  /*! Return the application's GID corresponding to the global number. */
   
  virtual AppGID &gno_to_appGid(GNO gno);

  /*! Return the application's GID corresponding to the local number. */
   
  virtual AppGID &lno_to_appGid(GNO gno);
};

/*! Zoltan2::ExodusFileSource
    \brief This class represents a mesh read in from an Exodus file.
*/

template<class Scalar, class LNO , class GNO, class AppGID>
  class ExodusFileSource : public ObjectSource{

private:

  std::string file_name;
  int dimension;

public:

  /*! Constructor */
  ExodusFileSource():file_name("unset"){}

  /*! Constructor 
  
       TODO: A constructor and private variables
             containing the parallel file info in Zoltan's PARIO_INFO
   */
  ExodusFileSource(std::string &fileName){file_name = string(fileName);}

  /*! Destructor */
  ~ExodusFileSource(){}

  /*! Copy Constructor */
  ExodusFileSource(const ExodusFileSource &efs){
    this->file_name(efs.get_file_name());
  }

  /*! Assignment operator */
  ExodusFileSource &operator=(const ExodusFileSource &efs){
    if (this == &efs) return *this;

    ObjectSource *me = this;    // call base class operator
    ObjectSource *you = &efs;
    *me = *you

    file_name(efs.get_file_name());
    return *this;
  }

  bool has_matrix(){return false;}
  bool has_dense_matrix(){return false;}
  bool has_csr_matrix(){return false;}
  bool has_csc_matrix(){return false;}
  bool has_graph(){return false;}

  bool has_vertices(){return true;}
  bool has_vertex_coordinates(){return true;}
   
  int coordinate_dimension(){ return dimension;}

  bool symmetric_matrix(){return false;}
  bool bipartite_graph(){return false;}

  /*! Returns the global and local number of vertices.

   *  Vertices may also be abstract objects known only by an ID.  We call them
   *  "vertices" because the normal use case is geometric or a graph.
   */
  void number_of_vertices(GNO &global_count, LNO &local_count){}

  /*! Returns the applications's object global IDs.

   *   TODO: we should be able to do partitioning without having global
   *         IDs from the application when there's no connectivity.
   *   The Zoltan methods shouldn't call this - they will work with GNOs.
   *   But when creating or reading the Result object, they need the AppGIDs.
   */
  int caller_gid_list(AppGID *gids){}

  /*! Returns the library's global numbering for the objects.

   *   to the order of the gnos in this list.
   */
  int global_number_list(GNO *gnos){}

  /*! Returns the coordinates of the vertices. */

  int vertex_coordinates(Tpetra_MultiVector<Scalar, LNO, GNO> &coords){}

  // TODO: add methods to get matrix and graph info.

  /*! Return the local number corresponding to the global number.
    
      Local numbers are consecutive numbers beginning at 0 for each process.
      Global numbers are consecutive across the application and begin at 0.
    
      TODO - local/global/application ID management is likely to be the same across
        all sources.  Maybe this should be done in the base class.
    
      TODO - for efficiency, there should be list/vector versions of these calls.
   */
  LNO gno_to_lno(GNO gno);

  /*! Return the global number corresponding to the local number. */

  GNO lno_to_gno(LNO lno);

  /*! Return the application's GID corresponding to the global number. */

  AppGID &gno_to_appGid(GNO gno);

  /*! Return the application's GID corresponding to the local number. */

  AppGID &lno_to_appGid(GNO gno);
};

/*! Zoltan2::ExodusFileDualSource
    \brief This class represents a graph created from the mesh in an Exodus file.

   The mesh elements are the vertices and the edges represent shared faces.
*/

template<class Scalar, class LNO , class GNO, class AppGID>
  class ExodusFileDualSource : public ObjectSource{
};

/* TODO
    QueryCallerSource  (for legacy Zoltan code)
    ArrayViewSource   (pointers to caller's arrays)
    MMFileSource        matrix market file
    ChacoFileSource        chaco file
    TpetraCrsMatrixSource, etc.       Tpetra and Epetra objects
*/

} // namespace Zoltan2

#endif /* _ZOLTAN2_OBJECTSOURCE_HPP_ */
