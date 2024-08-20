//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeEventIndexRange_decl_hpp
#define Tempus_TimeEventIndexRange_decl_hpp

#include <tuple>

#include "Teuchos_Time.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Tempus_config.hpp"
#include "Tempus_TimeEventBase.hpp"

namespace Tempus {

/** \brief TimeEventRangeIndex specifies a start, stop and stride index.
 *
 *
 */
template <class Scalar>
class TimeEventRangeIndex : virtual public TimeEventBase<Scalar> {
 public:
  /// Default constructor.
  TimeEventRangeIndex();

  /// Construct with full argument list of data members.
  TimeEventRangeIndex(int start, int stop, int stride, std::string name = "");

  /// Destructor
  virtual ~TimeEventRangeIndex() {}

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
   *  Find if an event is within the input range, inclusively
   *  ( index1 <= event <= index2 ).
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
  /** \brief Set the range of event indices.
   *
   *  This will completely replace the range
   *
   *  \param start  [in] The start of the index range.
   *  \param stop   [in] The stop of the index range.
   *  \param stride [in] The stride of the index range.
   */
  virtual void setIndexRange(int start, int stop, int stride)
  {
    setIndexStart(start);
    setIndexStop(stop);
    setIndexStride(stride);
  }

  /// Return the start of the index range.
  virtual int getIndexStart() const { return start_; }
  /// Set the start of the index range.
  virtual void setIndexStart(int start);

  /// Return the stop of the index range.
  virtual int getIndexStop() const { return stop_; }
  /// Set the stop of the index range.
  virtual void setIndexStop(int stop);

  /// Return the stride of the index range.
  virtual int getIndexStride() const { return stride_; }
  /// Set the stride of the index range.
  virtual void setIndexStride(int stride);

  /// Return the number of events.
  virtual int getNumEvents() const { return numEvents_; }

  /// Set the number of events from start_, stop_ and stride_.
  virtual void setNumEvents();
  //@}

  /** \brief Return a valid ParameterList with current settings.
   *
   *  The returned ParameterList will contain the current parameters
   *  and can be used to reconstruct the exact same object using
   *  createTimeEventRangeIndex(...).
   *
   * \return Teuchos::ParameterList of TimeEventRangeIndex.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

 protected:
  int start_;           ///< Start of index range.
  int stop_;            ///< Stop of index range.
  int stride_;          ///< Stride of index range.
  unsigned numEvents_;  ///< Number of events in index range.
};

// Nonmember Contructors
// ------------------------------------------------------------------------

/** \brief Nonmember Constructor via ParameterList.
 *
 *  If the input ParameterList is Teuchos::null, return a default
 *  TimeEventRangeIndex.  A valid ParameterList can be obtained
 *  from getValidParameters().
 *
 *  \param pList [in] The input ParameterList to construct from.
 *  \return Constructed TimeEventRangeIndex.
 */
template <class Scalar>
Teuchos::RCP<TimeEventRangeIndex<Scalar> > createTimeEventRangeIndex(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_TimeEventIndexRange_decl_hpp
