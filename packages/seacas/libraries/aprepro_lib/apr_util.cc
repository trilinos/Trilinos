#include <cstdio>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <errno.h>

#include "aprepro.h"
#include "aprepro_parser.h"

namespace {
  std::vector<char*> allocations;
}

namespace SEAMS {
extern Aprepro *aprepro;

void set_type(const SEAMS::Aprepro &apr, SEAMS::symrec* var, int type)
{
  if (var->name[0] == '_' || !apr.state_is_immutable()) {
    var->type = type;
  } else {
    if (type == Parser::token::VAR)
      var->type = Parser::token::IMMVAR;
    else if (type == Parser::token::SVAR)
      var->type = Parser::token::IMMSVAR;
    else
    var->type = type;
  }
}

void new_string(const char *from, char **to)
{
  int len=strlen(from);
  *to = new char[len+1];
  std::memcpy(*to, from, len+1);
  allocations.push_back(*to);
}

void concat_string(const char *from1, const char *from2, char **to)
{
  int len=strlen(from1) + strlen(from2);
  *to = new char[len+1];
  std::memcpy(*to, from1, len+1);
  std::strcat(*to, from2);
  allocations.push_back(*to);
}

  /* This function returns a pointer to a static character array.
   * If called multiple times, the name will change on each call,
   * so don't save the pointer; copy the data it points to if you
   * need to keep it for awhile.
   */
  char *get_temp_filename()
  {
    static char tmp_name[] = "./aprepro_temp_XXXXXX";
    int fd;

    std::strcpy(tmp_name, "./aprepro_temp_XXXXXX");
#if defined(__CYGWIN__) && defined(__NO_CYGWIN_OPTION__) 
    fd = mkstemps(tmp_name, 0);
#else
    fd = mkstemp(tmp_name);
#endif
    close(fd);
    return tmp_name;
  }  

void yyerror (const SEAMS::Aprepro &apr, const std::string &s)
{
  std::cerr << "Aprepro: ERROR:  '" << s << "' ("
	    << apr.ap_file_list.top().name << ", line "
	    << apr.ap_file_list.top().lineno + 1 << ")\n";
}

void immutable_modify(const SEAMS::Aprepro &apr, const SEAMS::symrec *var)
{
  std::cerr << "Aprepro: (IMMUTABLE) Variable " << var->name
	    << " is immutable and cannot be modified ("
	    << apr.ap_file_list.top().name << ", line "
	    << apr.ap_file_list.top().lineno + 1 << ")\n";
}

void undefined_warning (const SEAMS::Aprepro &apr, const std::string &var)
{
  if (apr.ap_options.warning_msg)
    std::cerr << "Aprepro: WARN: Undefined variable '"
	      << var << "' (" << apr.ap_file_list.top().name << ", line "
	      << apr.ap_file_list.top().lineno + 1 <<")\n";
}

void redefined_warning (const SEAMS::Aprepro &apr, const std::string &var)
{
  if (var[0] != '_' && apr.ap_options.warning_msg)
    std::cerr << "Aprepro: WARN: Variable '"
	      << var << "' redefined (" << apr.ap_file_list.top().name << ", line "
	      << apr.ap_file_list.top().lineno + 1 <<")\n";
}

void warning (const SEAMS::Aprepro &apr, const std::string &s)
{
  if (apr.ap_options.warning_msg)
    std::cerr << "Aprepro: WARN:  '" << s << "' ("
	      << apr.ap_file_list.top().name << ", line "
	      << apr.ap_file_list.top().lineno + 1 << ")\n";
}

void math_error(const SEAMS::Aprepro &apr, const char *function)
{
  if (errno != 0) {
    yyerror(apr, function);
  }
  if (errno == EDOM)
    perror("	DOMAIN error");
  else if (errno == ERANGE)
    perror("	RANGE error");
  else if (errno != 0)
    perror("	Unknown error");
}

void math_error(const char *function)
{
  if (errno != 0) {
    yyerror(*SEAMS::aprepro, function);
  }
  if (errno == EDOM)
    perror("	DOMAIN error");
  else if (errno == ERANGE)
    perror("	RANGE error");
  else if (errno != 0)
    perror("	Unknown error");
}

/* Convert string to all lower-case and replace all spaces with '_' */
void conv_string (char *string)
{
  char *p = string;
  while (*p != '\0')
    {
      if (*p == ' ')
	*p = '_';
      else if (isupper ((int)*p))
	*p = tolower ((int)*p);
      p++;
    }
}

void cleanup_memory()
{
  for (size_t i=0; i < allocations.size(); i++) {
    delete [] allocations[i];
  }
}
}
