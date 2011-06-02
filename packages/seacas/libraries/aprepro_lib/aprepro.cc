#include <fstream>
#include <cstring>
#include <vector>
#include <sstream>
#include <climits>
#include <iostream>
#include <ctime>

#include "aprepro.h"
#include "aprepro_parser.h"
#include "apr_scanner.h"
#include "apr_stats.h"
#include "apr_util.h"


namespace {
  const unsigned int HASHSIZE = 5939;
  const char* version_string = "3.08 (2011/05/05)";
  
  unsigned hash_symbol (const char *symbol)
  {
    unsigned hashval;
    for (hashval = 0; *symbol != '\0'; symbol++)
      hashval = *symbol + 65599 * hashval;
    return (hashval % HASHSIZE);
  }
}

namespace SEAMS {
  Aprepro *aprepro;  // A global for use in the library.  Clean this up...
  
  Aprepro::Aprepro()
    : sym_table(HASHSIZE)
  {
    ap_file_list.push(file_rec());
    init_table('#');
    ap_options.debugging = false;
    ap_options.trace_parsing = false;
    aprepro = this;

    // See the random number generator...
    time_t time_val = std::time ((time_t*)NULL);
    srand((unsigned)time_val);
  }

  Aprepro::~Aprepro()
  {
    for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
      for (symrec *ptr = sym_table[hashval]; ptr != NULL; ) {
	symrec *save = ptr;
	ptr = ptr->next;
	delete save;
      }
    }
    aprepro = NULL;
  }

  std::string Aprepro::version() const {return version_string;}
  
  bool Aprepro::parse_stream(std::istream& in, const std::string& in_name)
  {
    ap_file_list.top().name = in_name;

    Scanner scanner(*this, &in, &parsingResults);
    this->lexer = &scanner;

    Parser parser(*this);
    parser.set_debug_level(ap_options.trace_parsing);
    return (parser.parse() == 0);
  }

  bool Aprepro::parse_file(const std::string &filename)
  {
    std::ifstream in(filename.c_str());
    if (!in.good()) return false;
    return parse_stream(in, filename);
  }

  bool Aprepro::parse_string(const std::string &input, const std::string& sname)
  {
    std::istringstream iss(input);
    return parse_stream(iss, sname);
  }

  void Aprepro::error(const std::string& m) const
  {
    std::cerr << "Aprepro: ERROR: " << m << " ("
	      << ap_file_list.top().name<< ", line "
	      << ap_file_list.top().lineno + 1 << ")\n";
  }

  /* Two methods for opening files. In OPEN_FILE, the file must exist
     or else the code will exit. If CHECK_OPEN_FILE, it is OK if the
     file does not exist. A Null 'pointer' will be returned.
  */
  std::fstream* Aprepro::open_file(const std::string &file, const char *mode)
  {
    std::ios::openmode smode = std::ios::in;
    if (mode[0] == 'w')
      smode = std::ios::out;
  
    /* See if file exists in current directory (or as specified) */
    std::fstream *pointer = new std::fstream(file.c_str(), smode);
    if ((pointer == NULL || pointer->bad() || !pointer->good()) && !ap_options.include_path.empty()) {
      /* If there is an include path specified, try opening file there */
      std::string file_path(ap_options.include_path);
      file_path += "/";
      file_path += file;
      pointer = new std::fstream(file_path.c_str(), smode);
    }

    /* If pointer still null, print error message */
    if (pointer == NULL || pointer->bad() || !pointer->good()) {
      char tmpstr[128];
      sprintf(tmpstr, "Aprepro: ERROR:  Can't open '%s'",file.c_str()); 
      perror(tmpstr);
      exit(EXIT_FAILURE);
    }

    return pointer;
  }

  std::fstream *Aprepro::check_open_file(const std::string &file, const char *mode)
  {
    std::ios::openmode smode = std::ios::in;
    if (mode[0] == 'w')
      smode = std::ios::out;

    std::fstream *pointer = new std::fstream(file.c_str(), smode);

    if ((pointer == NULL || pointer->bad() || !pointer->good()) && !ap_options.include_path.empty()) {
      /* If there is an include path specified, try opening file there */
      std::string file_path(ap_options.include_path);
      file_path += "/";
      file_path += file;
      pointer = new std::fstream(file_path.c_str(), smode);
    }
    return pointer;
  }

  symrec *Aprepro::putsym (const char *sym_name, SYMBOL_TYPE sym_type, bool is_internal)
  {
    int parser_type = 0;
    switch (sym_type)
      {
      case VARIABLE:
	parser_type = Parser::token::VAR;
	break;
      case STRING_VARIABLE:
	parser_type = Parser::token::SVAR;
	break;
      case UNDEFINED_VARIABLE:
	parser_type = Parser::token::UNDVAR;
	break;
      case FUNCTION:
	parser_type = Parser::token::FNCT;
	break;
      case STRING_FUNCTION:
	parser_type = Parser::token::SFNCT;
	break;
      }
    symrec *ptr = new symrec(sym_name, parser_type, is_internal);
    if (ptr == NULL)
      return NULL;
  
    unsigned hashval = hash_symbol (ptr->name.c_str());
    ptr->next = sym_table[hashval];
    sym_table[hashval] = ptr;
    return ptr;
  }

  symrec *Aprepro::getsym (const char *sym_name) const
  {
    symrec *ptr = NULL;
    for (ptr = sym_table[hash_symbol (sym_name)]; ptr != NULL; ptr = ptr->next)
      if (strcmp (ptr->name.c_str(), sym_name) == 0)
	return ptr;
    return NULL;
  }

  void Aprepro::dumpsym (int type, bool doInternal) const
  {
    char comment = getsym("_C_")->value.svar[0];
  
    if (type == Parser::token::VAR || type == Parser::token::SVAR) {
      printf ("\n%c   Variable    = Value\n", comment);
      for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
	for (symrec *ptr = sym_table[hashval]; ptr != NULL; ptr = ptr->next) {
	  if ((doInternal && ptr->isInternal) || (!doInternal && !ptr->isInternal)) {
	    if (ptr->type == Parser::token::VAR)
	      printf ("%c  {%-10s\t= %.10g}\n", comment, ptr->name.c_str(), ptr->value.var);
	    else if (ptr->type == Parser::token::SVAR)
	      printf ("%c  {%-10s\t= \"%s\"}\n", comment, ptr->name.c_str(), ptr->value.svar);
	  }
	}
      }
    }
    else if (type == Parser::token::FNCT || type == Parser::token::SFNCT) {
      printf ("\nFunctions returning double:\n");
      for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
	for (symrec *ptr = sym_table[hashval]; ptr != NULL; ptr = ptr->next) {
	  if (ptr->type == Parser::token::FNCT) {
	    printf ("%-20s:  %s\n", ptr->syntax.c_str(), ptr->info.c_str());
	  }
	}
      }
      printf ("\nFunctions returning string:\n");
      for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
	for (symrec *ptr = sym_table[hashval]; ptr != NULL; ptr = ptr->next) {
	  if (ptr->type == Parser::token::SFNCT) {
	    printf ("%-20s:  %s\n", ptr->syntax.c_str(), ptr->info.c_str());
	  }
	}
      }
    }
  }

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

#define MAXLEN 16
  void Aprepro::statistics(std::ostream *out) const
  {
    std::ostream *output = out;
    if (output == NULL)
      output = &std::cout;
    
    symrec *ptr;
    unsigned entries = 0;
    int maxlen = 0;
    int minlen = INT_MAX;
    int lengths[MAXLEN];
    int longer = 0;
    double hash_ratio = 0.0;

    Stats stats;
    
    memset ((void *) lengths, 0, sizeof (lengths));

    for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
      int chain_len = 0;
      for (ptr = sym_table[hashval]; ptr != NULL; ptr = ptr->next)
	chain_len++;

      hash_ratio += chain_len * (chain_len + 1.0);
      entries += chain_len;
      if (chain_len >= MAXLEN)
	++longer;
      else
	++lengths[chain_len];

      minlen = min (minlen, chain_len);
      maxlen = max (maxlen, chain_len);

      if (chain_len > 0)
	stats.newsample (chain_len);
    }

    hash_ratio = hash_ratio / ((float) entries / HASHSIZE *
			       (float) (entries + 2.0 * HASHSIZE - 1.0));
    (*output) << entries << " entries in " << HASHSIZE << " element hash table, "
	      << lengths[0] << " (" << ((double) lengths[0] / HASHSIZE) * 100.0 << "%) empty.\n"
	      << "Hash ratio = " << hash_ratio << "\n"
	      << "Mean (nonempty) chain length = " << stats.mean()
	      << ", max = " << maxlen
	      << ", min = " << minlen
	      << ", deviation = " << stats.deviation() << "\n";

    for (int i = 0; i < MAXLEN; i++)
      if (lengths[i])
	(*output) << lengths[i] << " chain(s) of length " << i << "\n";

    if (longer)
      (*output) << longer << " chain(s) of length " << MAXLEN << " or longer\n";
  }

  void Aprepro::copyright(std::ostream *out) const
  {
    std::ostream *output = out;
    if (output == NULL)
      output = &std::cout;
    
    (*output) << "COPYRIGHT NOTICE\n";
  }
}
