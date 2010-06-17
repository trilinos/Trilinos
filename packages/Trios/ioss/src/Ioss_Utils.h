/*--------------------------------------------------------------------*/
/*    Copyright 2000-2008 Sandia Corporation.                         */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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

#define IOSS_ERROR(errmsg) throw std::runtime_error(errmsg.str())
#define IOSS_WARNING std::cerr

namespace Ioss {
  class GroupingEntity;
  class EntityBlock;
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

    static void calculate_faceblock_membership(IntVector &face_is_member,
					       const EntityBlock *ef_blk,
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
  };
}
#endif
