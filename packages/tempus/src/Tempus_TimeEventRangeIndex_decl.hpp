// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeEventIndexRange_decl_hpp
#define Tempus_TimeEventIndexRange_decl_hpp

#include <tuple>

// Teuchos
#include "Teuchos_Time.hpp"

// Tempus
#include "Tempus_TimeEventBase.hpp"


namespace Tempus {


/** \brief TimeEventRangeIndex specifies a start, stop and stride index.
 *
 *
 */
template<class Scalar>
class TimeEventRangeIndex : virtual public TimeEventBase<Scalar>
{
public:

  /// Default constructor.
  TimeEventRangeIndex();

  /// Construct with full argument list of data members.
  TimeEventRangeIndex(std::string name, int start, int stop, int stride);

  /// Destructor
  virtual ~TimeEventRangeIndex() {}

  /// \name Basic methods
  //@{
    /// Test if index is a time event.
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
    virtual void setIndexRange(int start, int stop, int stride)
    { setIndexStart(start); setIndexStop(stop); setIndexStride(stride); }

    virtual int getIndexStart() const { return start_; }
    virtual void setIndexStart(int start);

    virtual int getIndexStop() const { return stop_; }
    virtual void setIndexStop(int stop);

    virtual int getIndexStride() const { return stride_; }
    virtual void setIndexStride(int stride);

    virtual int getNumEvents() const { return numEvents_; }
    virtual void setNumEvents();
  //@}


protected:

  int start_;
  int stop_;
  int stride_;
  unsigned numEvents_;

};


} // namespace Tempus

#endif // Tempus_TimeEventIndexRange_decl_hpp
