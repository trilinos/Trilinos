#include <map>

#include <stk_util/environment/LogControl.hpp>

#include <stk_util/util/string_case_compare.hpp>

namespace stk {

namespace {

typedef std::map<std::ostream *, LogControl *> OStreamLogControlMap;

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


RuleMap &
get_rule_map() 
{
  static RuleMap s_ruleMap;

  return s_ruleMap;
}


OStreamLogControlMap &
get_ostream_log_control_map() 
{
  static OStreamLogControlMap s_ostreamLogControlMap;

  return s_ostreamLogControlMap;
}

} // namespace <unnamed>


LogControlRuleInterval::LogControlRuleInterval(
  int           interval) 
  : m_interval(interval),
    m_count(-1)
{}


bool
LogControlRuleInterval::next()
{
  ++m_count;
    
  if (m_count < 0) {
    return false;
  }
    
  else if (m_count == 0) {
    return true;
  }
    
  else {
    return m_count%m_interval == 0 ? true : false;
  }
}


void
addLogControlRule(
  const std::string &           rule_name,
  const LogControlRule &        rule)
{
  get_rule_map().addLogControlRule(rule_name, rule);
}


LogControl::LogControl(
  std::ostream &                log_ostream,
  const LogControlRule &        rule)
  : m_parent(0),
    m_rule(rule.clone()),
    m_state(ON),
    m_logStream(log_ostream),
    m_logStreambuf(log_ostream.rdbuf()),
    m_cacheStream()
{
  OStreamLogControlMap &ostream_log_control_map = get_ostream_log_control_map();

  // Append this as tail of linked list of LogControl's sharing this ostream.
  m_parent = ostream_log_control_map[&m_logStream];
  ostream_log_control_map[&m_logStream] = this;

  // Make sure log stream buffer is that of the root's.
  for (LogControl *parent = m_parent; parent != 0; parent = parent->m_parent)
    m_logStreambuf = parent->m_logStream.rdbuf();
}


LogControl::LogControl(
  std::ostream &        log_ostream,
  const std::string &   rule_name)
  : m_parent(0),
    m_rule(get_rule_map().getLogControlRule(rule_name)->clone()),
    m_state(ON),
    m_logStream(log_ostream),
    m_logStreambuf(log_ostream.rdbuf()),
    m_cacheStream()
{
  OStreamLogControlMap &ostream_log_control_map = get_ostream_log_control_map();

  // Append this as tail of linked list of LogControl's sharing this ostream.
  m_parent = ostream_log_control_map[&m_logStream];
  ostream_log_control_map[&m_logStream] = this;
}


LogControl::~LogControl()
{
  OStreamLogControlMap &ostream_log_control_map = get_ostream_log_control_map();

  // Reset tail pointer to this's parent.
  ostream_log_control_map[&m_logStream] = m_parent;

  // Reset log stream to either the parent's cache or the original log stream buffer. And,
  // concatenate cached text to the parent's cache or log stream.
  if (!m_parent || m_parent->m_state == ON) {           // Parent is writing
    m_logStream.rdbuf(m_logStreambuf);                  //   Set output stream back to real buffer
    if (m_state == CACHE)                               //   This is caching
      m_logStream << m_cacheStream.str();               //     Last cache is always written
  }
  else {                                                // Parent is caching
    m_logStream.rdbuf(m_parent->m_cacheStream.rdbuf()); //   Set output to parent's cache
    m_parent->m_cacheStream << m_cacheStream.str();     //   Append our cache to parent's
  }

  delete m_rule;
}


void
LogControl::fail()
{
  m_logStream.rdbuf(m_logStreambuf);
  m_logStream << m_cacheStream.str();
  m_cacheStream.str("");
  m_state = ON;
}


void
LogControl::next()
{
  m_cacheStream.str("");
  
  if (m_parent && m_parent->m_state == CACHE)
    m_state = CACHE;
  else 
    m_state = !m_rule || m_rule->next() ? ON : CACHE;
  
  if (m_state != CACHE) {
    if (m_logStream.rdbuf() != m_logStreambuf) {
      m_logStream.rdbuf(m_logStreambuf);
    }
  }
  else {
    if (m_logStream.rdbuf() != m_cacheStream.rdbuf())
      m_logStream.rdbuf(m_cacheStream.rdbuf());
  }
}


} // namespace stk
