#include <stk_percept/Observable.hpp>

namespace stk
{
namespace percept
{


#if 0

Observer::Observer(Observable& observable) : m_observable(observable)
{
}

Observer::~Observer()
{
  m_observable.removeObserver(*this);
}


void Observable::addObserver(Observer& observer)
{
  m_observers.push_back(&observer);
}

void Observable::removeObserver(Observer& observer)
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

Observable::Observers& Observable::getObservers()
{
  return m_observers;
}

#endif


}
}
