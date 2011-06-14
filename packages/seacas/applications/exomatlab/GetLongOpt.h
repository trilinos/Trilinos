/* S Manoharan. Advanced Computer Research Institute. Lyon. France */

#ifndef _GetLongOpt_h_
#define _GetLongOpt_h_

#include <CodeTypes.h>
#include <iostream>

class GetLongOpt {
 public:
  enum OptType { 
    NoValue, OptionalValue, MandatoryValue
  };
 private:
  struct Cell {
    const char *option;	// option name
    OptType type;		// option type
    const char *description;	// a description of option
    const char *value;	// value of option (string)
    Cell *next;		// pointer to the next cell

    Cell() { option = description = value = 0; next = 0; }
  };
 private:
  Cell *table;		// option table
  const char *ustring;	// usage message
  char *pname;		// program basename
  Cell *last;			// last entry in option table 
  int enroll_done;		// finished enrolling
  char optmarker;		// option marker

 private:
  int setcell(Cell *c, char *valtoken, char *nexttoken, const char *p);
 public:
  explicit GetLongOpt(const char optmark = '-');
  ~GetLongOpt();

  static char *basename(char * const p);

  int parse(int argc, char * const *argv);
  int parse(char * const str, char * const p);

  int enroll(const char * const opt, const OptType t,
	     const char * const desc, const char * const val);
  const char *retrieve(const char * const opt) const;

  void usage(std::ostream &outfile = std::cout) const;
  void usage(const char *str)		{ ustring = str; }
};
#endif /* _GetLongOpt_h_ */
