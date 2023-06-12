/*
// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
 */

#include "stk_util/util/Signal.hpp"
#include "stk_util/util/Callback.hpp"       // for Callback
#include "stk_util/util/FeatureTest.hpp"    // for SIERRA_MPI_ABORT_SIGNAL, SIERRA_USER_SHUTDOWN...
#include "stk_util/util/ReportHandler.hpp"  // for SignalHandler
#include "stk_util/util/SignalHandler.hpp"  // for SignalHandler
#include "stk_util/stk_config.h"
#include <csignal>                          // for SIGBUS, SIGILL, SIGSEGV, raise, SIGTERM
#include <sstream>

namespace sierra {
namespace Env {

namespace {

/**
 * @brief Class <b>EnvSignal</b> ...
 *
 */
class EnvSignal
{
public:
  static EnvSignal &instance() {
    static EnvSignal env_signal;

    return env_signal;
  }

private:
  EnvSignal()
    : m_enabled(false),
      m_hupReceived(false),
      hupCallback(*this, &EnvSignal::hupHandler),
      segvCallback(*this, &EnvSignal::segvHandler),
      illCallback(*this, &EnvSignal::illHandler),
      busCallback(*this, &EnvSignal::busHandler),
      termCallback(*this, &EnvSignal::termHandler)
  {}

public:
  /**
   * @brief Member function <b>getSigJmpBuf</b> establishes the signal handling.
   *
   * This function returns false under normal execution.  From that point forward, if any
   * signal occurs, the execution returns to this statement and returns value of true.  At
   * that point, there is no hope of recovery, but the traceback information is still
   * valid and, hopefully, the output streams are still functioning to receive the output.
   *
   * if (run_sierra.signalException())
   *   throw run_sierra.faultMessage();
   *
   * @return		a <b>bool</b> of <b>false</b> during normal
   *			execution.  Revisited and returns <b>true</b> when a signal
   *			occurs.
   */
  inline sigjmp_buf *getSigJmpBuf() {
    m_enabled = true;

    return &m_sigJmpBuf;
  }

  /**
   * @brief Member function <b>faultMessage</b> returns a
   * <b>std::runtime_error</b> with the message associated with the most recent
   * signal.
   *
   * @return		a <b>std::runtime_error</b> suitable for throwing.
   */
  const std::string &message() const {
    return m_message;
  }

  /**
   * @brief Member function <b>getHUPReceived</b> ...
   *
   */
  bool getHUPReceived() {
    return m_hupReceived;
  }

  /**
   * @brief Member function <b>activateSignals</b> ...
   *
   */
  void activateSignals();

  /**
   * @brief Member function <b>deactivateSignals</b> ...
   *
   */
  void deactivateSignals();

private:
  void doSignal(const char *message, int signal);
  void hupHandler();
  void segvHandler();
  void illHandler();
  void busHandler();
  void termHandler();

private:
  sigjmp_buf		m_sigJmpBuf;			///< setjmp/longjmp buffer
  bool			m_enabled;			///< Signal handling enabled state
  bool			m_hupReceived;			///< HUP has been received via handler
  std::string		m_message;			///< Message generated from signal

  Callback<EnvSignal>	hupCallback;
  Callback<EnvSignal>	segvCallback;
  Callback<EnvSignal>	illCallback;
  Callback<EnvSignal>	busCallback;
  Callback<EnvSignal>	termCallback;
};


void
EnvSignal::activateSignals()
{
#ifndef STK_ASAN_IS_ON
  SignalHandler::instance().add_handler(SIGSEGV, EnvSignal::segvCallback);
  SignalHandler::instance().add_handler(SIGILL, EnvSignal::illCallback);
  SignalHandler::instance().add_handler(SIGBUS, EnvSignal::busCallback);
//   SignalHandler::instance().add_handler(SIGINT, EnvSignal::intCallback);

#if defined(SIERRA_USER_SHUTDOWN_SIGNAL)
  SignalHandler::instance().add_handler(SIERRA_USER_SHUTDOWN_SIGNAL, EnvSignal::hupCallback);
#endif
#if defined(SIERRA_SHUTDOWN_SIGNAL)
  SignalHandler::instance().add_handler(SIERRA_SHUTDOWN_SIGNAL, EnvSignal::hupCallback);
#endif
#if defined(SIERRA_MPI_ABORT_SIGNAL)
  SignalHandler::instance().add_handler(SIERRA_MPI_ABORT_SIGNAL, EnvSignal::termCallback);
#endif
#endif
}


void
EnvSignal::deactivateSignals()
{
#ifndef STK_ASAN_IS_ON
  SignalHandler::instance().remove_handler(SIGSEGV, EnvSignal::segvCallback);
  SignalHandler::instance().remove_handler(SIGILL, EnvSignal::illCallback);
  SignalHandler::instance().remove_handler(SIGBUS, EnvSignal::busCallback);
//   SignalHandler::instance().add_handler(SIGINT, EnvSignal::intCallback);

#if defined(SIERRA_USER_SHUTDOWN_SIGNAL)
  SignalHandler::instance().remove_handler(SIERRA_USER_SHUTDOWN_SIGNAL, EnvSignal::hupCallback);
#endif
#if defined(SIERRA_SHUTDOWN_SIGNAL)
  SignalHandler::instance().remove_handler(SIERRA_SHUTDOWN_SIGNAL, EnvSignal::hupCallback);
#endif
#if defined(SIERRA_MPI_ABORT_SIGNAL)
  SignalHandler::instance().remove_handler(SIERRA_MPI_ABORT_SIGNAL, EnvSignal::termCallback);
#endif
#endif
}


void
EnvSignal::doSignal(
  const char *	message,
  int		signal)
{
  if (!m_enabled) {
//    std::cout << message << std::endl << "Signal exception handling not enabled" << std::endl;
    ::raise(signal);
    return;
  }
  else {
    m_enabled = false;
    std::ostringstream os;
    os << message << "\noccurred at:\n";
    stk::output_stacktrace(os) << "\n";
    m_message = os.str();
    ::siglongjmp(m_sigJmpBuf, signal);
  }
}


void
EnvSignal::hupHandler()
{
  m_hupReceived = true;
}


void
EnvSignal::segvHandler()
{
  SignalHandler::instance().remove_handler(SIGSEGV, EnvSignal::segvCallback);
  doSignal("Segmentation violation error", SIGSEGV);
}

void
EnvSignal::busHandler()
{
  SignalHandler::instance().remove_handler(SIGBUS, EnvSignal::busCallback);
  doSignal("Bus error", SIGBUS);
}


void
EnvSignal::illHandler()
{
  SignalHandler::instance().remove_handler(SIGILL, EnvSignal::illCallback);
  doSignal("Illegal instruction error", SIGILL);
}


void
EnvSignal::termHandler()
{
  SignalHandler::instance().remove_handler(SIGTERM, EnvSignal::termCallback);
  doSignal("Terminate signal received, likely due to an abort on another processor\n"
	   "Refer to standard output log for more information", SIGTERM);
}

} // namespace <unnamed>

void
activate_signals()
{
  EnvSignal::instance().activateSignals();
}

void
deactivate_signals()
{
  EnvSignal::instance().deactivateSignals();
}

sigjmp_buf *
get_sigjmpbuf()
{
  return EnvSignal::instance().getSigJmpBuf();
}

const std::string &
get_signal_message()
{
  return EnvSignal::instance().message();
}

bool
HUP_received()
{
  return EnvSignal::instance().getHUPReceived();
}

} // namespace Env
} // namespace sierra
