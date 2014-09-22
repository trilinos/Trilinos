/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_ENVIRONMENT_LOGCONTROL_HPP
#define STK_UTIL_ENVIRONMENT_LOGCONTROL_HPP

#include <map>                          // for _Rb_tree_iterator, etc
#include <sstream>                      // for ostream, ostringstream, etc
#include <stk_util/util/string_case_compare.hpp>  // for LessCase
#include <string>                       // for string
#include <utility>                      // for pair


namespace stk {

/**
 * @brief Interface <code>LogControlRule</code> describes the interface to a log control rule.
 *
 */
struct LogControlRule
{
  /**
   * Destroys a <code>LogControlRule</code> instance.
   *
   */
  virtual ~LogControlRule()
  {}

  /**
   * @brief Member function <code>clone</code> creates a clone of the rule.
   *
   * @return			a <code>LogControlRule</code> pointer to newly created duplicate.
   */
  virtual LogControlRule *clone() const = 0;

  /**
   * @brief Member function <code>next</code> returns true if the log stream should write to the log
   * file, and false if the log stream should write to the cache.
   *
   * @return			a <code>bool</code> value of true if the log stream should write to
   *                            the log file, and false if the log stream should write to the cache.
   */
  virtual bool next() = 0;
};


/**
 * @brief Class <code>LogControlRuleAlways</code> is a log control rule that always wants to write to the log
 * file.
 *
 */
struct LogControlRuleAlways : public LogControlRule
{
  /**
   * Destroys a <code>LogControlRuleAlways</code> instance.
   *
   */
  virtual ~LogControlRuleAlways()
  {}

  /**
   * Creates a new <code>LogControlRuleAlways</code> instance.
   *
   */
  LogControlRuleAlways()
  {}

  /**
   * @brief Member function <code>clone</code> creates a duplicate LogControlRuleAlways object.
   *
   * @return			a <code>LogControlRule</code> pointer to the new duplicated always object.
   */
  virtual LogControlRule *clone() const {
    return new LogControlRuleAlways(*this);
  }

  /**
   * @brief Member function <code>next</code> returns true to indicate that the log stream should
   * write to the log file.
   *
   * @return			a <code>bool</code> returns true to indicate that the log stream
   *                            should write to the log file.
   */
  virtual bool next() {
    return true;
  }
};


struct LogControlRuleInterval : public LogControlRule
{
  /**
   * Creates a new <code>LogControlRuleInterval</code> instance.
   *
   * @param interval		an <code>int</code> interval to enable log output.
   *
   */
  LogControlRuleInterval(int interval);

  /**
   * Destroys a <code>LogControlRuleInterval</code> instance.
   *
   */
  virtual ~LogControlRuleInterval()
  {}

  /**
   * @brief Member function <code>clone</code> creates a duplicate LogControlRuleAlways object.
   *
   * @return			a <code>LogControlRule</code> pointer to the new duplicated always object.
   */
  virtual LogControlRule *clone() const {
    return new LogControlRuleInterval(*this);
  }

  /**
   * @brief Member function <code>next</code> returns true when the current count modulo the interval is zero.
   * whichs indicate that the log stream should write to the log file.
   *
   * @return			a <code>bool</code> returns true when the current count modulo the
   *                            interval is zero.  whichs indicate that the log stream should write
   *                            to the log file.
   */
  virtual bool next();

private:
  int           m_interval;
  int           m_count;
};


class RuleMap
{
public:
  typedef std::map<std::string, LogControlRule *, LessCase> Map;

  RuleMap()
    : m_ruleMap()
  {}

  ~RuleMap() {
    for (Map::iterator it = m_ruleMap.begin(); it != m_ruleMap.end(); ++it)
      delete (*it).second;
  }

  void addLogControlRule(const std::string &rule_name, const LogControlRule &rule) {
    Map::iterator it = m_ruleMap.find(rule_name);
    if (it != m_ruleMap.end())
      m_ruleMap.erase(it);

    m_ruleMap[rule_name] = rule.clone();
  }

  LogControlRule *getLogControlRule(const std::string &rule_name) {
    Map::iterator it = m_ruleMap.find(rule_name);

    if (it != m_ruleMap.end())
      return (*it).second;

    else {
      std::pair<Map::iterator, bool> result = m_ruleMap.insert(Map::value_type(rule_name, new LogControlRuleAlways));
      return (*result.first).second;
    }
  }

private:
  Map           m_ruleMap;
};


/**
 * @brief Enumeration <code>State</code> describes the current state of the caching for this
 * controller.
 *
 */
enum State {
  ON,           ///< Output is to be written to the log stream
  CACHE         ///< Output is to be written to the cache stream
};

/**
 * @brief Class <code>LogControl</code> provides a mechanism for reducing excessive output.  The
 * output is redirected to a cache where it can be written to the log stream where and error
 * condition arises.
 *
 * The controlling of the log stream is handled by creating a sentry which controls the stream
 * buffer of the specified stream using the specified rule.  The next() function executes the rule
 * and redirects the output to the log stream when the rule is true and to the cache when the rule
 * is false.
 *
 * LogControl sentries can be nested.  When nested, the current rule is if the parent is caching,
 * then child is forced to cache.  This behavior could change by passing parent state to next().
 *
 * It's important to note that LogControl sentries nearly always shared the same output stream.  So
 * the parent's original output stream buffer
 */
class LogControl
{
public:
  /**
   * Creates a new <code>LogControl</code> instance.
   *
   * @param log_stream		a <code>std::ostream</code> reference to the log stream to control.
   *
   * @param rule		a <code>LogControlRule</code> reference to the rule used to control the log
   *                            stream.
   *
   */
  LogControl(std::ostream &log_stream, const LogControlRule &rule);

  /**
   * Creates a new <code>LogControl</code> instance.
   *
   * @param log_stream		a <code>std::ostream</code> reference to the log stream to control.
   *
   * @param rule_name		a <code>std::string</code> constant reference to rule name used to
   *                            control the log stream.
   */
  LogControl(std::ostream &log_stream,const std::string &rule_name);

  /**
   * Destroys a <code>LogControl</code> instance.
   *
   */
  ~LogControl();

  /**
   * @brief Member function <code>next</code> executes the rule and sets the log stream to write to
   * the log file if true and to the cache if false.
   *
   */
  void next();

  /**
   * @brief Member function <code>fail</code> writes the cached output to the log stream due to an
   * error.
   *
   */
  void fail();

private:
  LogControl *          m_parent;               ///< Parent stream
  LogControlRule *      m_rule;                 ///< Rule to evaluate log destination

  State                 m_state;                ///< Current caching state

  std::ostream &        m_logStream;            ///< Log stream under control
  std::streambuf *      m_logStreambuf;         ///< Log stream original stream buffer
  std::ostringstream    m_cacheStream;          ///< Cache stream

  // Do not implement...
  LogControl(const LogControl&);
  LogControl & operator = (const LogControl&);
};

} // namespace stk

#endif //  STK_UTIL_ENVIRONMENT_LOGCONTROL_HPP
