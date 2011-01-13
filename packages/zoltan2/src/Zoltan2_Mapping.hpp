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

#ifndef _ZOLTAN2_MAPPING_HPP_
#define _ZOLTAN2_MAPPING_HPP_

/*!
    \file Zoltan2_Mapping.hpp
    \brief The Mapping base class and derived classes.
*/

#include <Zoltan2_InputAdapter.hpp>
#include <Teuchos_ParameterList.hpp>

/*! Z2_Interface3
    \brief These classes are only used in interface 3.
*/

namespace Z2_Interface3
{

/*! Z2_Interface3::Mapping
    \brief A Mapping is an assignment of integers to the user's objects.

  The Mapping is the central object in Zoltan2.  It represents a mapping from
  the user's input (graph vertices, matrix rows, geometric coordinates, etc)
  to parts, colors, or an ordering.  The Mapping object has the user's input
  and objective and can compute and return the mapping.
*/

template<typename Scalar, typename LNO , typename GNO, typename AppGID>
  class Mapping {

private:

  Teuchos::ParameterList params;

public:

  /*! \name Constructors, Destructor, Copy and Assignment */

  /*! Constructors */
  Mapping(){}

  /*! Destructor */
  ~Mapping(){}

  /*! Copy Constructor 
  
    \param m is the initializer of this object.
  */
  Mapping(const Mapping &m){
  }

  /*! Assignment operator 

    \param m is the right hand side of the copy operator.
  */
  Mapping &operator=(const Mapping &m){
  }

  /*! \name Setting parameters for mapping. */

  void set_parameters(Teuchos::ParameterList &params);

  /*! \name the Mapping interface
   */
  
  /*! Compute the mapping.
   */

  virtual void solve();
  
  /*! Return the mapping.
   */

  virtual void get_mapping(std::vector<AppGID> &gid, std::vector<int> &value);

};

/*! Zoltan2::PartMapping
    \brief A PartMapping maps objects to part numbers.

  Because there are so many partitioning methods in Zoltan, it may
   make sense for PartMapping to be subclassed.
*/

template<class Scalar, class LNO , class GNO, class AppGID>
  class PartMapping : public Mapping<Scalar, LNO, GNO, AppGID> {

private:

  Z2::PartitioningObjective<Scalar, LNO, GNO, AppGID> objective;
  std::map<GNO, int> part_number;

public:

  /*! \name Constructors, Destructor, Copy and Assignment */

  PartMapping();

  PartMapping(Z2::PartitioningObjective &po);

  ~PartMapping();

  /*! \name Methods to set the objective.
   */

  void set_number_of_parts();

  void set_part_sizes();

  void set_this();

  void set_that();

  /*! \name Methods to get results.
   */

  void get_rcb_tree();

  void get_rcb_box();

  /*! \name Mapping interface
   */

  void solve();

  void get_mapping(std::vector<AppGID> &gid, std::vector<int> &value);
};

/*! Zoltan2::ColorMapping
    \brief A ColorMapping maps graph vertices to colors;
*/

template<class Scalar, class LNO , class GNO, class AppGID>
  class ColorMapping : public Mapping<Scalar, LNO, GNO, AppGID> {

private:

  Z2::ColoringObjective<Scalar, LNO, GNO, AppGID> objective;
  std::map<GNO, int> color;

public:

  /*! \name Constructors, Destructor, Copy and Assignment */

  ColorMapping();

  ColorMapping(Z2::ColoringObjective &co);

  ~ColorMapping();

  /*! \name Methods to set the objective.
   */

  void set_this();

  void set_that();

  /*! \name Methods to get results.
   */

  void get_number_of_colors();


  /*! \name Mapping interface
   */

  void solve();

  void get_mapping(std::vector<AppGID> &gid, std::vector<int> &value);
};

/*! Zoltan2::OrderMapping
    \brief A OrderMapping contains a matrix order problem.
*/

template<class Scalar, class LNO , class GNO, class AppGID>
  class OrderMapping : public Mapping<Scalar, LNO, GNO, AppGID> {

private:

  Z2::OrderingObjective<Scalar, LNO, GNO, AppGID> objective;
  std::map<GNO, int> order;

public:

  /*! \name Constructors, Destructor, Copy and Assignment */

  OrderMapping();

  OrderMapping(Z2::OrderingObjective &oo);

  ~OrderMapping();

  /*! \name Methods to set the objective.
   */

  void set_this();

  void set_that();

  /*! \name Mapping interface
   */

  void solve();

  void get_mapping(std::vector<AppGID> &gid, std::vector<int> &value);

};
};

} // namespace Z2_Interface3

#endif /* _ZOLTAN2_MAPPING_HPP_ */
