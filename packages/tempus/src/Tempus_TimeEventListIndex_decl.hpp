//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeEventListIndex_decl_hpp
#define Tempus_TimeEventListIndex_decl_hpp

#include <vector>

#include "Teuchos_Time.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Tempus_config.hpp"
#include "Tempus_TimeEventBase.hpp"

namespace Tempus {

/** \brief TimeEventListIndex specifies a list of index events.
 *
 *
 */
template <class Scalar>
class TimeEventListIndex : virtual public TimeEventBase<Scalar> {
 public:
  /// Default constructor.
  TimeEventListIndex();

  /// Construct with full argument list of data members.
  TimeEventListIndex(std::vector<int> indexList,
                     std::string name = "TimeEventListIndex");

  /// Destructor
  virtual ~TimeEventListIndex() {}

  /// \name Basic methods
  //@{
  /** \brief Test if index is a time event.
   *
   *  Return true if an event is the input index.
   *
   *  \param index [in] The input index.
   *  \return True if index is an event.
   */
  virtual bool isIndex(int index) const;

  /** \brief How many indices until the next event.
   *
   *  \param index [in] The input index.
   *  \return The number of steps (indices) to the next event.
   */
  virtual int indexToNextEvent(int index) const;

  /** \brief Return the index of the next event following the input index.
   *
   *  Returns the index of the next event that follows the input index.
   *  If the input index is before all events, the index of the first
   *  event is returned.  If the input index is after all events, the
   *  default index (an index in the distant future) is returned.  If the
   *  input index is an event index, the index of the next event is returned.
   *
   *  \param index     [in] Input index.
   *  \return Index of the next event.
   */
  virtual int indexOfNextEvent(int index) const;

  /** \brief Test if an event occurs within the index range.
   *
   *  Find if an event is within the input range,
   *  ( index1 < event <= index2 ).
   *
   *  \param index1 [in] Input index of one end of the range.
   *  \param index2 [in] Input index of the other end of the range.
   *  \return True if an index event is within the range.
   */
  virtual bool eventInRangeIndex(int index1, int index2) const;

  /// Describe member data.
  virtual void describe(Teuchos::FancyOStream &out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /// \name Accessor methods
  //@{
  /// \brief Return a vector of event indices.
  virtual std::vector<int> getIndexList() const { return indexList_; }

  /** \brief Set the vector of event indices.
   *
   *  This will completely replace the vector of event indices.
   *
   *  \param indexList [in] Vector of event indices.
   *  \param sort      [in] Sort vector into ascending order, if true.
   */
  virtual void setIndexList(std::vector<int> indexList, bool sort = true);

  /** \brief Add the index to event vector.
   *
   *  The input index will be inserted into the vector of
   *  events in ascending order.  If the index is already
   *  present, it is not added to keep the vector unique.
   *
   *  \param index [in] Index to insert to vector of events.
   */
  virtual void addIndex(int index);

  /// Clear the vector of all events.
  virtual void clearIndexList() { indexList_.clear(); }
  //@}

  /** \brief Return a valid ParameterList with current settings.
   *
   *  The returned ParameterList will contain the current parameters
   *  and can be used to reconstruct the exact same object using
   *  createTimeEventListIndex(...).
   *
   * \return Teuchos::ParameterList of TimeEventListIndex.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

 protected:
  std::vector<int> indexList_;  // Sorted and unique list of index events.
};

// Nonmember Contructors
// ------------------------------------------------------------------------

/** \brief Nonmember Constructor via ParameterList.
 *
 *  If the input ParameterList is Teuchos::null, return a default
 *  TimeEventListIndex.  A valid ParameterList can be obtained
 *  from getValidParameters().
 *
 *  \param pList [in] The input ParameterList to construct from.
 *  \return Constructed TimeEventListIndex.
 */
template <class Scalar>
Teuchos::RCP<TimeEventListIndex<Scalar> > createTimeEventListIndex(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_TimeEventListIndex_decl_hpp
