#ifndef XPETRA_REGION_UTILS_DECL_HPP_
#define XPETRA_REGION_UTILS_DECL_HPP_

// headers


// forward declarations
namespace Teuchos
{
  template<>
  class Comm<std::size_t>;
}

namespace Xpetra
{

/*! \brief Overloading the comparison operator for sorting of regions
 *
 * This compare class is used to run the sorting algorithm on the list of nodes
 * with associated regions they belong to.
 *
 * First, nodes are sorted in ascending order for region labels. Then, the
 * sorting shuffles the nodes in ascending node index for each given region.
 */
template<class GO>
bool compareRegions(const std::tuple<GO,GO>& lhs, ///< quantity on the left of the comparison operator
    const std::tuple<GO,GO>& rhs ///< quantity on the right of the comparison operator
    );

/*! \brief Overloading the comparison operator for sorting of nodes
 *
 * This compare is used to run the sorting algorithm where the nodes are ordered
 * in ascending order for the node indices, regardless of the associated region
 * index.
 */
template<class GO>
bool compareNodes(const std::tuple<GO,GO>& lhs, ///< quantity on the left of the comparison operator
    const std::tuple<GO,GO>& rhs ///< quantity on the right of the comparison operator
    );

/*! \class ChekcerNode
 *
 * Definition of the predicate for the node_ structure.
 * Given a tuple made of node index and a specific region it belongs to,
 * this predicate returns true if the node belongs to the region specified in input to the predicate.
 */
template<class GO>
class CheckerNode
{

public:

  //! @name Construction/Destruction
  //@{

  //! Constructor
  CheckerNode(GO regionIndex ///< ID of current region
      )
  {
    regionIndex_ = regionIndex;
  };

  //@}

  //! Unary Operator
  bool operator()(const std::tuple<GO,GO> &node)
  {
    return (std::get<1>(node) == regionIndex_);
  }

private:

  //! ID of region to be checked
  GO regionIndex_;

};

/*! \class CheckerNodesToRegion
 *
 *  Definition of the predicate for the nodesToRegion_ structure
 *  Given a tuple made of node index and a vector with labels of regions it belongs to,
 *  this predicate returns true if the node coincides with the node specified in input to the predicate.
 *
 *  \ author mayr.mt \date 09/2017
 */
template<class GO>
class CheckerNodesToRegion {

public:

  //! @name Construction/Destruction
  //@{

  //! Constructor
  CheckerNodesToRegion(GO node_index)
  {
    node_index_ = node_index;
  };


  //@}

  //! Unary Operator
  bool operator()(const std::tuple<GO,Teuchos::Array<GO> >& node)
  {
    return (std::get<0>(node) == node_index_);
  }

protected:

private:

  //! ID of node to be checked
  GO node_index_;

};

//Definition of the predicate for the regionToAll structure.
//Given a tuple made of node index and a specific region it belongs to,
//this predicate returns true if the node has composite index which coincides with the index specified in input.
template<class GO>
class CheckerRegionToAll {

public:

  //! @name Construction/Destruction
  //@{

  //! Constructor
  CheckerRegionToAll(GO node_index)
  {
    node_index_ = node_index;
  }

  //@}

  //Unary Operator
  bool operator()(const std::tuple<GO,GO> &node)
  {
    return (std::get<0>(node) == node_index_);
  }

private:

  GO node_index_;

};


//Definition of the predicate for the node_ structure.
//Given a tuple made of node index and a specific region it belongs to,
//this predicate returns true if the node has composite index which coincides with the index specified in input.
//It does the same thing as checkerRegionToAll but it works on a different data structure
template<class GO>
class CheckerAllToRegion {

public:

  //! @name Construction/Destruction
  //@{

  //! Constructor
  CheckerAllToRegion(GO node_index)
  {
    node_index_ = node_index;
  }
  ;

  //@}

  //Unary Operator
  bool operator()(const std::tuple<GO,GO> &node)
  {
    return (std::get < 1 > (node) == node_index_);
  }

private:

  GO node_index_;

};

//Definition of the predicate for the node_ structure.
//Given a tuple made of node index and a specific region it belongs to,
//this predicate returns true if the node has composite index which coincides with the index specified in input.
//This checker is specifically used only for nodes lying on the interface
template<class GO>
class CheckerInterfaceNodes {

public:

  //! @name Construction/Destruction
  //@{

  //! Constructor
  CheckerInterfaceNodes(GO node_index)
  {
    node_index_ = node_index;
  }

  //@}

  // Unary Operator
  bool operator()(const std::tuple<int,Teuchos::Array<GO> > &node)
  {
    return (std::get<0>(node) == node_index_);
  }

private:

  GO node_index_;

};

} /* namespace Xpetra */

#endif /* #ifndef XPETRA_REGION_UTILS_DECL_HPP_ */
