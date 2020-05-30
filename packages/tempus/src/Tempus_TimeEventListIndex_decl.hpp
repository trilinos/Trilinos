// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeEventListIndex_decl_hpp
#define Tempus_TimeEventListIndex_decl_hpp

#include <vector>

// Teuchos
#include "Teuchos_Time.hpp"

// Tempus
#include "Tempus_TimeEventBase.hpp"


namespace Tempus {


/** \brief TimeEventListIndex specifies a list of index events.
 *
 *
 */
template<class Scalar>
class TimeEventListIndex : virtual public TimeEventBase<Scalar>
{
public:

  /// Default constructor.
  TimeEventListIndex();

  /// Construct with full argument list of data members.
  TimeEventListIndex(std::string name, std::vector<int> indexList);

  /// Destructor
  virtual ~TimeEventListIndex() {}

  /// \name Basic methods
  //@{
    /// Test if index is near a index event (within tolerance).
    virtual bool isIndex(int index) const;

    /// How many indices until the next event. Negative indicating the last event is in the past.
    virtual int indexToNextEvent(int index) const;

    /// Index of the next event. Negative indicating the last event is in the past.
    virtual int indexOfNextEvent(int index) const;

    /// Test if an event occurs within the index range.
    virtual bool eventInRangeIndex(int index1, int index2) const;

    /// Describe member data.
    virtual void describe() const;
  //@}

  /// \name Accessor methods
  //@{
    virtual std::vector<int> getIndexList() const { return indexList_; }
    virtual void setIndexList(std::vector<int> indexList, bool sort = true);
    virtual void addIndex(int index);
    virtual void clearIndexList() { indexList_.clear(); }
  //@}


protected:

  std::vector<int> indexList_; // Sorted and unique list of index events.

};


} // namespace Tempus

#endif // Tempus_TimeEventListIndex_decl_hpp
