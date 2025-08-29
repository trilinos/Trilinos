// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_Observable_hpp
#define percept_Observable_hpp

#include <iostream>
#include <vector>
#include <algorithm>

namespace percept
{

template<class DATA_TYPE>
class Observable;

template<class DATA_TYPE>
class Observer
{
public:

  Observer( Observable<DATA_TYPE>& observable)  : m_observable(observable)    
    {
      observable.addObserver(*this);
    }

  virtual ~Observer()
    {
      m_observable.removeObserver(*this);
    }

  virtual void notify(DATA_TYPE *data) = 0;

private:
   Observable<DATA_TYPE>& m_observable;
};

template<class DATA_TYPE>
class Observable
{
public:
  typedef std::vector<Observer<DATA_TYPE> *> Observers;

  Observable(){}

  void addObserver(Observer<DATA_TYPE>& observer)
    {
      m_observers.push_back(&observer);
    }

  void removeObserver(Observer<DATA_TYPE>& /*observer*/)
    {
#if 0
      std::vector<Observer*>::iterator found;
      found = std::find(m_observers.begin(), m_observers().end(), &observer);
      if (found != m_observers().end())
      {
        m_observers.erase(found);
      }
#endif
    }

  Observers& getObservers()
    {
      return m_observers;
    }

  void notifyObservers(DATA_TYPE *data) 
    {
      for (unsigned i = 0; i < m_observers.size(); i++)
      {
        m_observers[i]->notify(data);
      }
    }

private:
  Observers m_observers;
};




}

#endif
