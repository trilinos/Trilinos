/**   ------------------------------------------------------------
 *    Copyright 2006-2007, 2009 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <iostream>

#include <signal.h>

#include <stk_util/diag/Env.hpp>
#include <stk_util/diag/Signal.hpp>
#include <stk_util/util/Callback.hpp>
#include <stk_util/util/FeatureTest.hpp>
#include <stk_util/diag/SignalHandler.hpp>

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
   * @brief Member function <b>disableSigLongJmp</b> disables the returning of
   * <b>true</b> to signalException when a signal is detected.
   *
   */
  void disableSigLongJmp() {
    m_enabled = false;
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
  void intHandler();

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
//  Callback<EnvSignal>	intCallback;

};


  void
EnvSignal::activateSignals()
{
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
}


void
EnvSignal::deactivateSignals()
{
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
    m_message = message;
    ::siglongjmp(m_sigJmpBuf, signal);
  }
}


void
EnvSignal::hupHandler()
{
  m_hupReceived = true;
}


void
EnvSignal::intHandler()
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


void
disable_siglongjmp()
{
  return EnvSignal::instance().disableSigLongJmp();
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
