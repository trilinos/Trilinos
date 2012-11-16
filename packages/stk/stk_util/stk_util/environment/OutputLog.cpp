/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/util/TeeStreambuf.hpp>
#include <stk_util/util/IndentStreambuf.hpp>


#include <map>
#include <list>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <cctype>

namespace stk {

namespace {

struct LogStream
{
  LogStream(const std::string &path, std::ostream *output_stream, std::ofstream *file_stream)
    : m_path(path),
      m_ostream(output_stream),
      m_ofstream(file_stream)
  {}

  ~LogStream();

  std::string           m_path;
  std::ostream *        m_ostream;
  std::ofstream *       m_ofstream;

  private:
  LogStream(const LogStream &);
  void operator = (const LogStream &);
};

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 444)
#endif
struct LogStreamMap : public std::map<std::string, LogStream *>
{
  LogStreamMap()
  {}

  ~LogStreamMap() {
    while (!empty()) {
      LogStream *log_stream = (*begin()).second;
      erase(begin());
      delete log_stream;
    }
  }
};
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

struct OStreamTeeStreambuf
{
  OStreamTeeStreambuf(std::ostream &output_stream)
    : m_ostream(&output_stream),
      m_origRdbuf(output_stream.rdbuf()),
      m_teeStreambuf(new tee_streambuf(&output_stream))
  {
    m_ostream->rdbuf(m_teeStreambuf);
  }

  ~OStreamTeeStreambuf();

  std::ostream *        m_ostream;
  std::streambuf *      m_origRdbuf;
  tee_streambuf *       m_teeStreambuf;

  private:
  OStreamTeeStreambuf(const OStreamTeeStreambuf &);
  void operator = (const OStreamTeeStreambuf &);
};

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 444)
#endif
struct OStreamTeeStreambufMap : public std::map<std::string, OStreamTeeStreambuf *>
{
  OStreamTeeStreambufMap()
  {}

  ~OStreamTeeStreambufMap() {
    while (!empty()) {
      OStreamTeeStreambuf *tee_streambuf = (*begin()).second;
      erase(begin());
      delete tee_streambuf;
    }
  }
};
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

LogStreamMap &
get_file_stream_map()
{
  static LogStreamMap s_logFileStreamMap;

  return s_logFileStreamMap;
}


OStreamTeeStreambufMap &
get_ostream_tee_streambuf_map()
{
  static OStreamTeeStreambufMap s_ostreamTeeStreambufMap;

  return s_ostreamTeeStreambufMap;
}


LogStream::~LogStream()
{
  m_ostream->flush();

  // If the output stream was created internally (via bind_output_stream), be sure to remove it from
  // all OStreamTeeStreamBuf's
  if (m_ofstream) {
    OStreamTeeStreambufMap &ostream_tee_streambuf_map = get_ostream_tee_streambuf_map();

    for (OStreamTeeStreambufMap::iterator it = ostream_tee_streambuf_map.begin(); it != ostream_tee_streambuf_map.end(); ++it)
      (*it).second->m_teeStreambuf->remove(m_ofstream);

    delete m_ofstream;
  }
}


OStreamTeeStreambuf::~OStreamTeeStreambuf()
{
  if (m_ostream) {
    m_ostream->flush();
    m_ostream->rdbuf(m_origRdbuf);
  }

  // Be sure to remove this from all OStreamTeeStreamBuf's
  OStreamTeeStreambufMap &ostream_tee_streambuf_map = get_ostream_tee_streambuf_map();

  for (OStreamTeeStreambufMap::iterator it = ostream_tee_streambuf_map.begin(); it != ostream_tee_streambuf_map.end(); ++it)
    (*it).second->m_teeStreambuf->remove(m_ostream);

  delete m_teeStreambuf;
}

} // namespace <empty>


void
create_log_file(
  const std::string &   name,
  const std::string &   path)
{
  LogStreamMap &file_stream_map = get_file_stream_map();

  close_log_file(name);

  std::ofstream *file_stream = new std::ofstream(path.c_str());

  if(!file_stream->good()) {

    std::ostringstream s;
    s << "Cannot open output log file '" << path << "' directory does not exist or is write protected.";

    throw std::runtime_error(s.str());

  }


  file_stream_map[name] = new LogStream(path, file_stream, file_stream);
}


void
close_log_file(
  const std::string &   name)
{
  LogStreamMap &file_stream_map = get_file_stream_map();

  LogStreamMap::iterator it = file_stream_map.find(name);

  if (it != file_stream_map.end()) {
    delete (*it).second;
    file_stream_map.erase(it);
  }
}


void
register_log_ostream(
  std::ostream &        os,
  const std::string &   name)
{
  LogStreamMap &file_stream_map = get_file_stream_map();

  LogStreamMap::iterator it = file_stream_map.find(name);

  if (it != file_stream_map.end()) {
    std::ostringstream s;
    s << "Log ostream " << name << " has already been registered";

    //Do we really want to throw if a stream is registered multiple times?
    //I don't think so... commenting this out.
    //throw std::runtime_error(s.str());
  }
  else {
    file_stream_map[name] = new LogStream(name, &os, 0);
  }
}


void
unregister_log_ostream(
  std::ostream &        os)
{
  LogStreamMap &file_stream_map = get_file_stream_map();

  for (LogStreamMap::iterator it = file_stream_map.begin(); it != file_stream_map.end(); ++it) {
    if ((*it).second->m_ostream == &os) {
      delete (*it).second;
      file_stream_map.erase(it);
      break;
    }
  }

  OStreamTeeStreambufMap &ostream_tee_streambuf_map = get_ostream_tee_streambuf_map();

  for (OStreamTeeStreambufMap::iterator it = ostream_tee_streambuf_map.begin(); it != ostream_tee_streambuf_map.end(); ++it)
    (*it).second->m_teeStreambuf->remove(&os);
}


const std::string &
get_log_path(
  const std::string &   name)
{
  static std::string not_found = "";

  LogStreamMap &file_stream_map = get_file_stream_map();

  LogStreamMap::iterator it = file_stream_map.find(name);

  return it == file_stream_map.end() ? not_found : (*it).second->m_path;
}


std::ostream *
get_log_ostream(
  const std::string &   name)
{
  LogStreamMap &file_stream_map = get_file_stream_map();

  LogStreamMap::iterator it = file_stream_map.find(name);

  return it == file_stream_map.end() ? 0 : (*it).second->m_ostream;
}


void
register_ostream(
  std::ostream &        os,
  const std::string &   name)
{
  OStreamTeeStreambufMap &ostream_tee_streambuf_map = get_ostream_tee_streambuf_map();

  unregister_ostream(os);

  OStreamTeeStreambufMap::iterator it = ostream_tee_streambuf_map.find(name);

  if (it != ostream_tee_streambuf_map.end()) {
//     delete (*it).second;
//     ostream_tee_streambuf_map.erase(it);
//   }
    std::ostringstream s;
    s << "Output stream " << name << " has already been registered";

    throw std::runtime_error(s.str());
  }

  ostream_tee_streambuf_map[name] = new OStreamTeeStreambuf(os);
}


void
unregister_ostream(
  std::ostream &        os)
{
  OStreamTeeStreambufMap &ostream_tee_streambuf_map = get_ostream_tee_streambuf_map();

  for (OStreamTeeStreambufMap::iterator it = ostream_tee_streambuf_map.begin(); it != ostream_tee_streambuf_map.end(); ++it)
    (*it).second->m_teeStreambuf->remove(&os);

  for (OStreamTeeStreambufMap::iterator it = ostream_tee_streambuf_map.begin(); it != ostream_tee_streambuf_map.end(); ++it) {
    if ((*it).second->m_ostream == &os) {
      delete (*it).second;
      ostream_tee_streambuf_map.erase(it);
      break;
    }
  }

}


std::ostream *
get_ostream_ostream(
  const std::string &   name)
{
  OStreamTeeStreambufMap &ostream_tee_streambuf_map = get_ostream_tee_streambuf_map();

  OStreamTeeStreambufMap::iterator it = ostream_tee_streambuf_map.find(name);

  return it == ostream_tee_streambuf_map.end() ? 0 : (*it).second->m_ostream;
}


tee_streambuf *
get_ostream_tee_streambuf(
  const std::string &   name)
{
  OStreamTeeStreambufMap &ostream_tee_streambuf_map = get_ostream_tee_streambuf_map();

  OStreamTeeStreambufMap::iterator it = ostream_tee_streambuf_map.find(name);

  return it == ostream_tee_streambuf_map.end() ? 0 : (*it).second->m_teeStreambuf;
}


std::ostream *
get_ostream_tee_ostream(
  const std::string &   name)
{
  OStreamTeeStreambufMap &ostream_tee_streambuf_map = get_ostream_tee_streambuf_map();

  OStreamTeeStreambufMap::iterator it = ostream_tee_streambuf_map.find(name);

  return it == ostream_tee_streambuf_map.end() ? 0 : (*it).second->m_ostream;
}


bool
is_registered_ostream(
  const std::string &   name)
{
  return get_ostream_ostream(name) != 0;
}


namespace {

struct Command
{
  virtual ~Command()
  {}

  virtual void execute() = 0;
};

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 444)
#endif
struct CommandList : public std::list<Command *>
{
  CommandList()
    : std::list<Command *>()
  {}

  ~CommandList()
  {
    for (std::list<Command *>::iterator it = begin(); it != end(); ++it)
      delete (*it);
  }
};
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

namespace {

tee_streambuf &
parse_tee_streambuf(
  const std::string &   tee_ostream_name)
{
  tee_streambuf *osb = get_ostream_tee_streambuf(tee_ostream_name);

  if (!osb) {
    std::ostringstream s;

    s << "Output stream " << tee_ostream_name << " has not been registered for output logging";
    throw std::runtime_error(s.str());
  }

  return *osb;
}


std::ostream *
parse_ostream(
  const std::string &   ostream_name)
{
  std::ostream *os = get_log_ostream(ostream_name);

  if (!os)
    os = get_ostream_tee_ostream(ostream_name);

  if (!os) {
    std::ostringstream s;

    s << "Log file '" << ostream_name << "' has not been registered";
    throw std::runtime_error(s.str());
  }

  return os;
}


struct OpenLog : public Command
{
  OpenLog(
    const std::string &name,
    const std::string &path)
    : m_name(name),
      m_path(path)
  {}

  virtual ~OpenLog()
  {}

  virtual void execute() {
    create_log_file(m_name, m_path);
  }

  std::string        m_name;
  std::string        m_path;
};


struct CloseLog : public Command
{
  CloseLog(
    const std::string &name)
    : m_name(name)
  {}

  virtual ~CloseLog()
  {}

  virtual void execute() {
    close_log_file(m_name);
  }

  std::string        m_name;
};


struct ClearTeeOStream : public Command
{
  ClearTeeOStream(
    const std::string & tee_ostream_name)
    : m_teeOStreamName(tee_ostream_name)
  {}

  virtual ~ClearTeeOStream()
  {}

  virtual void execute() {
    parse_tee_streambuf(m_teeOStreamName).clear();
  }

  std::string        m_teeOStreamName;
};


struct AddTeeOStream : public Command
{
  AddTeeOStream(
    const std::string & tee_ostream_name,
    const std::string & ostream_name)
    : m_teeOStreamName(tee_ostream_name),
      m_ostreamName(ostream_name)
  {}

  virtual ~AddTeeOStream()
  {}

  virtual void execute() {
    if (m_ostreamName != "null")
      parse_tee_streambuf(m_teeOStreamName).add(parse_ostream(m_ostreamName));
  }

  std::string        m_teeOStreamName;
  std::string        m_ostreamName;
};


struct RemoveTeeOStream : public Command
{
  RemoveTeeOStream(
    const std::string & tee_ostream_name,
    const std::string & ostream_name)
    : m_teeOStreamName(tee_ostream_name),
      m_ostreamName(ostream_name)
  {}

  virtual ~RemoveTeeOStream()
  {}

  virtual void execute() {
    parse_tee_streambuf(m_teeOStreamName).remove(parse_ostream(m_ostreamName));
  }

  std::string        m_teeOStreamName;
  std::string        m_ostreamName;
};

} // namespace <empty>


/*
 * Startup:     out > cout pout > cout dout > cout
 * Normal:      out > log-path+pout pout > null dout > out
 * Diagnostic:  out > out-path+pout pout > pout-path dout > out
 *
 * Modify:      out > +pout
 *              out > -pout
 */
void
parse_output_description(
  const std::string &   output_description,
  CommandList &         command_list)
{
  typedef std::pair<const char *, const char *>  Token;
  typedef std::list<Token> TokenList;

  command_list.clear();

  TokenList tokens;

  for (const char *c = output_description.c_str(); *c; ) {
    if (std::isspace(*c))
      ++c;

    else if (*c == '>' || *c == '+' || *c == '-' || *c == '=') {
      tokens.push_back(Token(c, c + 1));
      ++c;
    }

    else if (*c == '\"') {
      const char *d = c + 1;
      while (*d && *d != '\"')
        ++d;
      tokens.push_back(Token(c + 1, d));
      c = d + 1;
    }

    else {
      const char *d = c;
      while (std::isgraph(*d) && *d != '+' && *d != '-' &&*d != '=' && *d != '>')
        ++d;
      tokens.push_back(Token(c, d));
      c = d;
    }
  }

  for (TokenList::iterator it = tokens.begin(); it != tokens.end(); ) {
    std::string name((*it).first, (*it).second);

    ++it; if (it == tokens.end()) break;
    std::string operation((*it).first, (*it).second);

    if (operation == "=") {
      ++it;  if (it == tokens.end()) break;
      std::string path((*it).first, (*it).second);
      if (!path.empty())
        command_list.push_back(new OpenLog(name, path));
      else
        command_list.push_back(new CloseLog(name));
      ++it;  if (it == tokens.end()) break;
    }

    else if (operation == ">") {
      parse_tee_streambuf(name);

      ++it;  if (it == tokens.end()) break;
      std::string token(std::string((*it).first, (*it).second));
      if (token != "+" && token != "-") {
        std::string ostream_name(std::string((*it).first, (*it).second));

        command_list.push_back(new ClearTeeOStream(name));
        command_list.push_back(new AddTeeOStream(name, ostream_name));
        ++it;  if (it == tokens.end()) break;
      }

      while (it != tokens.end()) {
        token = std::string((*it).first, (*it).second);
        if (token == "+") {
          ++it;  if (it == tokens.end()) break;
          std::string ostream_name(std::string((*it).first, (*it).second));

          command_list.push_back(new AddTeeOStream(name, ostream_name));
          ++it;  if (it == tokens.end()) break;
        }

        else if (token == "-") {
          ++it;  if (it == tokens.end()) break;
          std::string ostream_name(std::string((*it).first, (*it).second));

          command_list.push_back(new RemoveTeeOStream(name, ostream_name));
          ++it;  if (it == tokens.end()) break;
        }
        else
          break;
      }
    }
  }
}

void
execute(
  const CommandList &   command_list)
{
  for (CommandList::const_iterator it = command_list.begin(); it != command_list.end(); ++it)
    (*it)->execute();
}

} // namespace <empty>

void
bind_output_streams(
  const std::string &   output_description)
{
  stk::CommandList command_list;

  parse_output_description(output_description, command_list);
  execute(command_list);
}

} // namespace stk

namespace sierra {

std::ostream &
out() {
  static std::ostream s_out(std::cout.rdbuf());

  return s_out;
}


std::ostream &
pout() {
  static std::ostream s_pout(std::cout.rdbuf());

  return s_pout;
}


std::ostream &
dout() {
  static std::ostream s_dout(std::cout.rdbuf());

  return s_dout;
}


std::ostream &
tout() {
  static std::ostream s_tout(std::cout.rdbuf());

  return s_tout;
}


std::ostream &
dwout() {
  static stk::indent_streambuf s_dwoutStreambuf(std::cout.rdbuf());
  static std::ostream s_dwout(&s_dwoutStreambuf);
  
  return s_dwout;
}

} // namespace sierra



