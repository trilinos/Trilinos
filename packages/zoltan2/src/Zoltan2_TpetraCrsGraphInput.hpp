// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_EpetraCrsGraphInput.hpp

    \brief An input adapter for a Epetra::CrsGraph.
*/

#ifndef _ZOLTAN2_TPETRACRSGRAPHINPUT_HPP_
#define _ZOLTAN2_TPETRACRSGRAPHINPUT_HPP_

#include <Tpetra_CrsGraph.hpp>
#include <Zoltan2_XpetraCrsGraphInput.hpp>
#include <Zoltan2_Standards.hpp>

namespace Zoltan2 {

/*! Zoltan2::TpetraCrsGraphInput
    \brief Provides access for Zoltan2 to Tpetra::CrsGraph data plus weights.

    TpetraCrsGraphInput provides access for Zoltan2 to a
    Tpetra::CrsGraph using an
    interface that is the same for all graph input.  Most of the methods used by
    Zoltan2 are in the XpetraGraphInput parent class.  The methods used
    by the Zoltan2 caller are in this class.

    The template parameter is the user's Tpetra::CrsGraph<lno, gno>.
*/

template<typename User>
class TpetraCrsGraphInput : public XpetraCrsGraphInput<User>
{
public:

  typedef typename InputAdapter<User>::scalar_t scalar_t;
  typedef typename InputAdapter<User>::lno_t    lno_t;
  typedef typename InputAdapter<User>::gno_t    gno_t;
  typedef typename InputAdapter<User>::lid_t    lid_t;
  typedef typename InputAdapter<User>::gid_t    gid_t;
  typedef typename InputAdapter<User>::node_t   node_t;

  /*! Name of input adapter type.
   */
  std::string inputAdapterName()const {return std::string("TpetraCrsGraph");}

  /*! Destructor
   */
  ~TpetraCrsGraphInput() { }

  /*! Constructor
   */
  TpetraCrsGraphInput(const RCP<const User> &graph): XpetraCrsGraphInput<User>(
    rcp(new Xpetra::TpetraCrsGraph<lno_t, gno_t, node_t>(
      rcp_const_cast<User>(graph)))) 
  {
    ingraph_ = graph;
  }


  /*! Access to the graph that instantiated the adapter
   */
  RCP<const User> getGraph()
  { 
    return ingraph_;
  }

private:
  RCP<const User > ingraph_;

};

} // namespace

#endif
