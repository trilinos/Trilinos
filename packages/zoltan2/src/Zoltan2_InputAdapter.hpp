// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_InputAdapter.hpp

    \brief The abstract interface for an input adapter.
*/


#ifndef _ZOLTAN2_INPUTADAPTER_HPP_
#define _ZOLTAN2_INPUTADAPTER_HPP_

#include <Zoltan2_Standards.hpp>

namespace Zoltan2 {

/*! InputAdapter defines methods required by all InputAdapters
 *
 *  About local IDs:
 *    Applications are required to supply unique global IDs for
 *    their objects, such as vertices, coordinates, matrix rows,
 *    and so on.
 *
 *    Local IDs are optional.  Local IDs are symbols which are
 *    meaningful to the application and that reference the object.
 *    A local ID may be an array index or a pointer for example.
 *
 *    The impact on input adapter, models and solutions is this:
 *       1. If local IDs are supplied, they must appear in the
 *          solution if the application requests a solution.
 *       2. The "set" methods in the input adapter must accept
 *          local IDs and global IDs, and the application may
 *          use either or both.
 *       3. If local IDs are available, the model should query
 *          the input adapter using local IDs, since normally this 
 *          is more efficient.
 */

enum InputAdapterType {
  InvalidAdapterType = 0,
  MatrixAdapterType,
  MeshAdapterType,
  GraphAdapterType,
  CoordAdapterType,
  IdAdapterType,
  XpetraCrsMatrixAdapterType  // Special case for performance with Epetra/Tpetra
};

template <typename User>
struct InputTraits {
  // Input Adapter implementations must provide the following typedefs
  // for use in Zoltan2:
  //   scalar_t :  weights and coordinates
  //   lno_t    :  ordinal (e.g., int, long, int64_t) that can represent
  //               the number of local data items.
  //   gno_t    :  ordinal (e.g., int, long, int64_t) that can represent
  //               the number of global data items.
  //   lid_t    :  user type that represents a locally unique identifier 
  //               for data items.
  //   gid_t    :  user type that represents a globally unique identifier 
  //               for data items.
  //   node_t   :  Kokkos node.
  //
  // Default typedefs are included here. If a specialization of User is
  // not provided, these types will be used.
  typedef float scalar_t;
  typedef int   lno_t;
  typedef long  gno_t;
  typedef int   lid_t;
  typedef long  gid_t;
  typedef Kokkos::DefaultNode::DefaultNodeType node_t;
  static inline std::string name() {return "InputAdapter";}
};

template <typename User>
class InputAdapter {
private:
public:

  typedef User user_t;
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::lid_t    lid_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;

  // Return enumerated InputAdapterType for the input adapter.
  // This function is implemented in the MatrixAdapter, GraphAdapter,
  // MeshAdapter, CoordAdapter and IdAdapter subclasses.
  // Users do not have to implement this function for their adapters
  // as long as they inherit from one of the subclasses (which they must).
  virtual enum InputAdapterType inputAdapterType() = 0;

  /*! Pure virtual destructor
   */
  virtual ~InputAdapter() {};

  /*! Return a name that identifies the concrete adapter.
   *  Useful for debugging.
   */
  virtual std::string inputAdapterName() const = 0;

  /*! Returns true if input adapter uses local Ids for objects.
   */
  virtual bool haveLocalIds() const = 0;

  /*! Return true if local Ids are consecutive integral
   *   values and supply the base.  Providing this information
   *   can save memory, making local ID lists unneccesary.
   */
  virtual bool haveConsecutiveLocalIds(size_t &base) const = 0;

protected:
};
  
  
}  //namespace Zoltan2
  
#endif
