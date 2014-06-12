/**   ------------------------------------------------------------
 *    Copyright 2002-2007, 2009 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <stdexcept>
#include <vector>
#include <iostream>

#include <stk_util/diag/SignalHandler.hpp>

#include <signal.h>
#include <time.h>

extern "C" {
  static void signal_handler(int signal, siginfo_t *sip, void *ucp)
  {
    // This routine is called for all signals...
    // Just a C-callable wrapper around a call to the true signal handler...
    sierra::SignalHandler::instance().handle_signal(signal);
  }
}

namespace sierra {

namespace {

int
convert_name_to_signal(
  const String &	signal)
{
  if (signal == "SIGABRT" || signal == "SIGKILL") {
    return -1;
  }

#if defined(SIGILL)
  if (signal == "SIGILL") return SIGILL;
#endif
#if defined(SIGSEGV)
  if (signal == "SIGSEGV") return SIGSEGV;
#endif
#if defined(SIGALRM)
  if (signal == "SIGALRM") return SIGALRM;
#endif
#if defined(SIGFPE)
  if (signal == "SIGFPE")  return SIGFPE;
#endif
#if defined(SIGHUP)
  if (signal == "SIGHUP")  return SIGHUP;
#endif
#if defined(SIGINT)
  if (signal == "SIGINT")  return SIGINT;
#endif
#if defined(SIGPIPE)
  if (signal == "SIGPIPE") return SIGPIPE;
#endif
#if defined(SIGQUIT)
  if (signal == "SIGQUIT") return SIGQUIT;
#endif
#if defined(SIGTERM)
  if (signal == "SIGTERM") return SIGTERM;
#endif
#if defined(SIGUSR1)
  if (signal == "SIGUSR1") return SIGUSR1;
#endif
#if defined(SIGUSR2)
  if (signal == "SIGUSR2") return SIGUSR2;
#endif
  return -2;
}

} // namespace <unnamed>


SignalHandler &
SignalHandler::instance() {
  static SignalHandler signal_handler;

  return signal_handler;
}


void
SignalHandler::handle_signal(
  int			signal)
{
  typedef std::vector<const HandlerMap::value_type *> HandlerList;

  time_t now = ::time(NULL);

  std::cerr << "Sierra received signal " << signal << " at " << ::ctime(&now) << std::endl;

  HandlerList   handlers;

  std::pair<HandlerMap::const_iterator, HandlerMap::const_iterator> range = m_handlerMap.equal_range(signal);

  for (HandlerMap::const_iterator pos = range.first; pos != range.second; ++pos)
    handlers.push_back(&*pos);

  for (HandlerList::const_iterator it = handlers.begin(); it != handlers.end(); ++it) {
    CallbackBase &obj = *(*it)->second;
    obj();
  }
}


bool
SignalHandler::check_signal_name(
  const String &	signal)
{
  int isignal = convert_name_to_signal(signal);
  return (isignal >= 0);
}


void
SignalHandler::add_handler(
  const String &	signal,
  CallbackBase &	callback)
{
  int isignal = convert_name_to_signal(signal);
  if (isignal >= 0) {
    add_handler(isignal, callback);
  }
  else if (isignal == -1)
    throw std::runtime_error("signal cannot be handled");
  else if (isignal == -2)
    throw std::runtime_error("signal name invalid");
  else
    throw std::logic_error("invalid value from convert_node_to_signal()");
}


void
SignalHandler::add_handler(
  int			signal,
  CallbackBase &	callback)
{
  // See if already handling this signal...
  if (m_handlerMap.find(signal) == m_handlerMap.end()) {
    // Tell OS that we want to handle this signal...
    struct sigaction action;
    struct sigaction *old_action = new struct sigaction;

    action.sa_sigaction = signal_handler;
    sigemptyset(&action.sa_mask);
    action.sa_flags = SA_SIGINFO;
    ::sigaction(signal, &action, old_action);
    m_oldActionMap.insert(OldActionMap::value_type(signal, old_action));
  }
  m_handlerMap.insert(HandlerMap::value_type(signal, &callback));
}


void
SignalHandler::remove_handler(
  int			signal,
  CallbackBase &	callback)
{
  typedef std::pair<HandlerMap::iterator, HandlerMap::iterator> HandlerRange;

  HandlerRange handler_range = m_handlerMap.equal_range(signal);
  for (HandlerMap::iterator it = handler_range.first; it != handler_range.second; )
    if ((*it).second == &callback) {
      HandlerMap::iterator erase_it = it++;
      m_handlerMap.erase(erase_it);
    }
    else
      ++it;

  if (m_handlerMap.find(signal) == m_handlerMap.end()) {
    OldActionMap::iterator it = m_oldActionMap.find(signal);
    if (it != m_oldActionMap.end()) {
      ::sigaction(signal, (*it).second, NULL);
      delete (*it).second;
      m_oldActionMap.erase(it);
    }
  }
}


void
SignalHandler::remove_handler(
  const String &	signal,
  CallbackBase &	callback)
{
  int isignal = convert_name_to_signal(signal);
  if (isignal >= 0) {
    remove_handler(isignal, callback);
  }
  else if (isignal == -1)
    throw std::runtime_error("signal cannot be handled");
  else if (isignal == -2)
    throw std::runtime_error("signal name invalid");
  else
    throw std::logic_error("invalid value from convert_node_to_signal()");
}


void
SignalHandler::remove_all_handlers()
{
  m_handlerMap.clear();

  for (OldActionMap::iterator it = m_oldActionMap.begin(); it != m_oldActionMap.end(); ++it) {
    ::sigaction((*it).first, (*it).second, NULL);
    delete (*it).second;
  }
  m_oldActionMap.clear();
}


} // namespace sierra
