// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef IOSS_Ioss_IOUtils_h
#define IOSS_Ioss_IOUtils_h

#include <Ioss_CodeTypes.h>

#include <string>

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cstdlib>

#define TOPTR(x) (x.empty() ? NULL : &x[0])
#define IOSS_ERROR(errmsg) throw std::runtime_error(errmsg.str())
#define IOSS_WARNING std::cerr

namespace Ioss {
  class GroupingEntity;
  class SideBlock;
  class Region;
  class Field;
  
  typedef std::vector<int> IntVector;

  class Utils {
  public:

    Utils();
    ~Utils() {};
    
    // Assignment operator
    // Copy constructor

    /*!
     * Fill time_string and date_string with current time and date
     * formatted as "HH:MM:SS" for time and "yy/mm/dd" or "yyyy/mm/dd"
     * for date
     */
    static void time_and_date(char* time_string, char* date_string, size_t length);

    static std::string decode_filename(const std::string &filename, int processor, int num_processors);
    static int    decode_entity_name(const std::string &entity_name);
    static std::string encode_entity_name(const std::string &entity_type, int id);

    /*!
     * Convert 'name' to lowercase and convert spaces to '_'
     */
    static void fixup_name(char *name);
    static void fixup_name(std::string &name);
    
    /*!
     * Process the base element type 'base' which has
     * 'nodes_per_element' nodes and a spatial dimension of 'spatial'
     * into a form that the IO system can (hopefully) recognize.
     * Lowercases the name; converts spaces to '_', adds
     * nodes_per_element at end of name (if not already there), and
     * does some other transformations to remove some exodusII ambiguity.
     */
    static std::string fixup_type(const std::string &base, int nodes_per_element, int spatial);

    static std::string uppercase(const std::string &name);
    static std::string lowercase(const std::string &name);

    static int case_strcmp(const std::string &s1, const std::string &s2);
    /*!
     * Return a string containing information about the current :
     * computing platform. This is used as information data in the
     * created results file to help in tracking when/where/... the file
     * was created.
     */
    static std::string platform_information();

    /*!
     * The following functions are wrappers around the sierra::Env functions
     * or similar substitutes to reduce coupling of the Sierra framewk to
     * the rest of the IO Subsystem. When compiled without the Framework code,
     * Only these functions need to be reimplemented.
     */
    static void abort();

    /*!
     * Return a filename relative to the specified working directory (if any)
     * of the current execution. Working_directory must end with '/' or be empty.
     */
    static std::string local_filename(const std::string &relative_filename,
				      const std::string &type,
				      const std::string &working_directory);

    static int field_warning(const Ioss::GroupingEntity *ge,
			     const Ioss::Field &field, const std::string& inout);

    static void calculate_sideblock_membership(IntVector &face_is_member,
					       const SideBlock *ef_blk,
					       const int *element, const int *sides,
					       int number_sides,
					       const Region *region);

    static unsigned int hash (const std::string& name);

    /*!
     * Return a vector of strings containing the lines of the input file.
     * Should only be called by a single processor or each processor will
     * be accessing the file at the same time...
     */
    static void input_file(const std::string &file_name,
			   std::vector<std::string> *lines, size_t max_line_length = 0);

    template <class T> static std::string to_string(const T & t)
      {
	std::ostringstream os;
	os << t;
	return os.str();
      }
    
    /*!
     * Many databases have a maximum length for variable names which can
     * cause a problem with variable name length.
     *
     * This routine tries to shorten long variable names to an acceptable
     * length ('max_var_len' characters max).  If the name is already less than this
     * length, it is returned unchanged...
     *
     * Since there is a (good) chance that two shortened names will match,
     * a 2-letter 'hash' code is appended to the end of the variable name.
     *
     * So, we shorten the name to a maximum of 'max_var_len'-3 characters and append a
     * 2 character hash+separator.
     *
     * It also converts name to lowercase and converts spaces to '_'
     */
    static std::string variable_name_kluge(const std::string &name,
					   size_t component_count, size_t copies,
					   size_t max_var_len);
      };
}
#endif
