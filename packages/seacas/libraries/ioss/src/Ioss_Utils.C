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

#include <Ioss_Utils.h>
#include <assert.h>
#include <stddef.h>
#include <sys/select.h>
#include <time.h>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <stdint.h>

#ifndef _WIN32
#include <sys/utsname.h>
#endif

#include <Ioss_ElementBlock.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_Field.h>
#include <Ioss_GroupingEntity.h>
#include <Ioss_Region.h>
#include <Ioss_SideBlock.h>
#include <fstream>

#include "Ioss_Property.h"

namespace {
  inline int to_lower(int c) { return std::tolower(c); }
  inline int to_upper(int c) { return std::toupper(c); }
}
Ioss::Utils::Utils() {}
  
void Ioss::Utils::time_and_date(char* time_string, char* date_string,
				size_t length)
{
  time_t calendar_time = time(NULL);
  struct tm *local_time = localtime(&calendar_time);

  strftime(time_string, length, "%H:%M:%S", local_time);
  if (length == 8) {
    const char *fmt = "%y/%m/%d";
    strftime(date_string, length, fmt, local_time);
    date_string[8] = '\0';

  } else if (length >= 10) {
    strftime(date_string, length, "%Y/%m/%d", local_time);
    date_string[10] = '\0';
  }
  time_string[8] = '\0';
}

std::string Ioss::Utils::decode_filename(const std::string &filename, int processor, int num_processors)
{
  std::string decoded_filename(filename);
  // Current format for per-processor file names is:
  // PREFIX/basename.num_proc.cur_proc
  // the 'cur_proc' field is padded to be the same width as
  // the 'num_proc' field
  // Examples: basename.8.1, basename.64.03, basename.128.001

  // Create a std::string containing the total number of processors
  std::string num_proc = to_string(num_processors);
  size_t proc_width = num_proc.length();

  // Create a std::string containing the current processor number
  std::string cur_proc = to_string(processor);
  size_t cur_width = cur_proc.length();

  // Build the filename
  decoded_filename += ".";
  decoded_filename += num_proc;
  decoded_filename += ".";

  // Now, pad with zeros so that 'cur_proc' portion is same
  // width as 'num_proc' portion.
  while (cur_width++ < proc_width) {
    decoded_filename += "0";
  }

  decoded_filename += cur_proc;
  return decoded_filename;
}

int Ioss::Utils::decode_entity_name(const std::string &entity_name)
{
  // ExodusII stores block, nodeset, and sideset ids as integers
  // Sierra   stores these as std::strings. The string is created by
  // concatenating the type, the character '_' and the id.
  // This function reverses the process and returns the original id.

  const char *name = entity_name.c_str();
  while (*name != '_')
    name++;
  // Increment one more time to get past '_'
  assert(*name == '_');
  name++;

  // name now points to beginning of id portion of string
  assert(*name >= '0' && *name <= '9');
  int id = std::atoi(name);

  return id;
}

std::string Ioss::Utils::encode_entity_name(const std::string &entity_type, int64_t id)
{
  // ExodusII stores block, nodeset, and sideset ids as integers
  // Sierra   stores these as std::strings. The string is created by
  // concatenating the type, the character '_' and the id.

  std::string id_string = to_string(id);
  std::string entity_name = entity_type;
  entity_name += "_";
  entity_name += id_string;
  return entity_name;
}

std::string Ioss::Utils::fixup_type(const std::string &base, int nodes_per_element, int spatial)
{
  std::string type = base;
  Ioss::Utils::fixup_name(type); // Convert to lowercase; replace spaces with '_'

  // Fixup an exodusII kluge/ambiguity.
  // The element block type does not fully define the element. For
  // example, a block of type 'triangle' may have either 3 or 6
  // nodes.  To fix this, check the block type name and see if it
  // ends with a number.  If it does, assume it is OK; if not, append
  // the 'nodes_per_element'.
  if (!isdigit(*(type.rbegin()))) {
    if (nodes_per_element > 1) {
      type += Ioss::Utils::to_string(nodes_per_element);
    }
  }

  // Fixup an exodusII kludge.  For triangular elements, the same
  // name is used for 2D elements and 3D shell elements.  Convert
  // to unambiguous names for the IO Subsystem.  The 2D name
  // stays the same, the 3D name becomes 'trishell#'
  if (spatial == 3) {
    if      (type == "triangle3") type = "trishell3";
    else if (type == "triangle6") type = "trishell6";
    else if (type == "tri3")      type = "trishell3";
    else if (type == "tri6")      type = "trishell6";
  }

  if (spatial == 2) {
    if (type == "shell2")
      type = "shellline2d2";
    else if (type == "rod2" || type == "bar2" || type == "truss2")
      type = "rod2d2";
    else if (type == "shell3")
      type = "shellline2d3";
    else if (type == "bar3"  || type == "rod3"  || type == "truss3")
      type = "rod2d3";
  }

  if (std::strncmp(type.c_str(), "super", 5) == 0) {
    // A super element can have a varying number of nodes.  Create
    // an IO element type for this super element just so the IO
    // system can read a mesh containing super elements.  This
    // allows the "omit volume" command to be used in the Sierra
    // applications to skip creating a corresponding element block
    // in the application.
    type = "super" + Ioss::Utils::to_string(nodes_per_element);
  }
  return type;
}

void Ioss::Utils::abort()
{
  std::ostringstream errmsg("I/O abort");
  IOSS_ERROR(errmsg);
}

std::string Ioss::Utils::local_filename(const std::string& relative_filename,
					const std::string& type,
					const std::string& working_directory)
{
  if (relative_filename[0] == '/' || type == "generated" || working_directory.empty()) {
    return relative_filename;
  } else {
    std::string filename = working_directory;
    filename += relative_filename;
    return filename;
  }
}

int Ioss::Utils::field_warning(const Ioss::GroupingEntity *ge,
			       const Ioss::Field &field, const std::string& inout)
{
  IOSS_WARNING << ge->type_string() << " '" << ge->name()
	       << "'. Unknown " << inout << " field '"
	       << field.get_name() << "'\n";
  return -4;
}

std::string Ioss::Utils::platform_information()
{
  // Return a string containing the 'uname' output.
  // This is used as information data in the created results file
  // to help in tracking when/where/... the file was created
#ifndef _WIN32
  struct utsname sys_info;
  uname(&sys_info);

  std::string info = "Node: ";
  info += sys_info.nodename;
  info += ", OS: ";
  info += sys_info.sysname;
  info += " ";
  info += sys_info.release;
  info += ", ";
  info += sys_info.version;
  info += ", Machine: ";
  info += sys_info.machine;
#else
  std::string info = "Node: Unknown, OS: Unknown, Machine: Unknown";
#endif
  return info;
}

namespace {
  bool block_is_omitted(Ioss::ElementBlock *block) {
    bool omitted = false;
    if (block->property_exists("omitted"))
      omitted = (block->get_property("omitted").get_int() == 1);
    return omitted;
  }
}
void Ioss::Utils::calculate_sideblock_membership(IntVector &face_is_member,
						 const Ioss::SideBlock *ef_blk,
						 size_t int_byte_size,
						 const void *element, const void *sides,
						 int64_t number_sides,
						 const Ioss::Region *region)
{
  assert(ef_blk != NULL);
  
  face_is_member.reserve(number_sides);
  
  // Topology of faces in this face block...
  const ElementTopology *ftopo = ef_blk->topology();

  // Topology of parent element for faces in this face block
  const ElementTopology *parent_topo = ef_blk->parent_element_topology();

  // If split by element block then parent_block will be non-NULL
  const ElementBlock *parent_block = ef_blk->parent_element_block();

  // The element block containing the face we are working on...
  Ioss::ElementBlock *block = NULL;

  // Topology of face/edge in current element block
  const ElementTopology *common_ftopo = NULL;

  // Topology of elements in the element block containing this element
  const ElementTopology *block_topo = NULL;

  // Topology of the face we are currently working with...
  const ElementTopology *topo = NULL;

  // The element side that the current face is on the element...
  int64_t current_side = -1;

  if (number_sides > 0 && (element == NULL || sides == NULL)) {
    std::ostringstream errmsg;
    errmsg << "INTERNAL ERROR: null element or sides pointer passed to "
	   << "Ioss::Utils::calculate_sideblock_membership function.";
    IOSS_ERROR(errmsg);
  }

  for (int64_t iel = 0; iel < number_sides; iel++) {
    int64_t elem_id = 0;
    int64_t side_id = 0;
    if (int_byte_size == 4) {
      elem_id = ((int*)element)[iel];
      side_id = ((int*)sides)[iel];
    } else {
      elem_id = ((int64_t*)element)[iel];
      side_id = ((int64_t*)sides)[iel];
    }

    // Get the element block containing this face...
    if (block == NULL || !block->contains(elem_id)) {
      block = region->get_element_block(elem_id);
      block_topo = block->topology();
      // NULL if hetero face/edge on element
      common_ftopo = block->topology()->boundary_type(0);
      if (common_ftopo != NULL)
	topo = common_ftopo;
      current_side = -1;
    }

    // If the element topology of the element block containing this
    // face has heterogeneous topology (eg. wedge), then determine the
    // topology corresponding to the current side..
    if (common_ftopo == NULL && side_id != current_side) {
      current_side = side_id;
      topo = block->topology()->boundary_type(side_id);
    }

    // See if the face topology and the parent element topology for
    // the current face match the topology associated with this face block.
    if (topo == ftopo && block_topo == parent_topo &&
	(parent_block == NULL || parent_block == block )
	&& !block_is_omitted(block)) {
      // This face/edge  belongs in the face/edge block
      face_is_member.push_back(1);
    } else {
      face_is_member.push_back(0);
    }
  }
}

unsigned int Ioss::Utils::hash (const std::string& name)
{
  // Hash function from Aho, Sethi, Ullman "Compilers: Principles,
  // Techniques, and Tools.  Page 436

  const char* symbol = name.c_str();
  unsigned int hashval;
  unsigned int g;
  for (hashval = 0; *symbol != '\0'; symbol++) {
    hashval = (hashval << 4) + *symbol;
    g = hashval&0xf0000000;
    if (g != 0) {
      hashval = hashval ^ (g >> 24);
      hashval = hashval ^ g;
    }
  }
  return hashval;
}

void Ioss::Utils::input_file(const std::string &file_name,
			     std::vector<std::string> *lines,
			     size_t max_line_length)
{
  // Create an ifstream for the input file. This does almost the same
  // function as sierra::Env::input() except this is for a single
  // processor and the sierra::Env::input() is for parallel...

  if (file_name.length() != 0) {
    // Open the file and read into the vector...
    std::string input_line;
    std::ifstream infile(file_name.c_str());
    lines->push_back(file_name.substr(0,max_line_length));
    while (!std::getline(infile, input_line).fail()) {
      if (max_line_length == 0 || input_line.length() <= max_line_length) {
	lines->push_back(input_line);
      } else {
	// Split the line into pieces of length "max_line_length-1"
	// and append a "\" to all but the last. Don't worry about
	// splitting at whitespace...
	size_t ibeg = 0;
	do {
	  std::string sub = input_line.substr(ibeg, max_line_length-1);
	  if (ibeg+max_line_length-1 < input_line.length()) {
	    sub += "\\";
	  }
	  lines->push_back(sub);
	  ibeg += max_line_length-1;
	} while (ibeg < input_line.length());
      }
    }
  }
}

int Ioss::Utils::case_strcmp(const std::string &s1, const std::string &s2)
{
  const char *c1 = s1.c_str();
  const char *c2 = s2.c_str();
  for ( ; ; c1++, c2++) {
    if (std::tolower(*c1) != std::tolower(*c2))
      return (std::tolower(*c1) - std::tolower(*c2));
    if (*c1 == '\0')
      return 0;
  }
}

std::string Ioss::Utils::uppercase(const std::string &name)
{
  std::string s(name);
  std::transform(s.begin(), s.end(), s.begin(), to_upper);
  return s;
}

std::string Ioss::Utils::lowercase(const std::string &name)
{
  std::string s(name);
  std::transform(s.begin(), s.end(), s.begin(), to_lower);
  return s;
}

void Ioss::Utils::fixup_name(char *name)
{
  // Convert 'name' to lowercase and convert spaces to '_'
  assert(name != NULL);

  size_t len = std::strlen(name);
  for (size_t i=0; i < len; i++) {
    name[i] = static_cast<char>(tolower(name[i]));  // guaranteed(?) to be ascii...
    if (name[i] == ' ')
      name[i] = '_';
  }
}

void Ioss::Utils::fixup_name(std::string &name)
{
  // Convert 'name' to lowercase and convert spaces to '_'
  name = Ioss::Utils::lowercase(name);
  
  size_t len = name.length();
  for (size_t i=0; i < len; i++) {
    if (name[i] == ' ')
      name[i] = '_';
  }
}

namespace {
  std::string two_letter_hash (const char *symbol)
  {
    // Hash function from Aho, Sethi, Ullman "Compilers: Principles,
    // Techniques, and Tools.  Page 436
    const int HASHSIZE=673; // Largest prime less than 676 (26*26)
    char word[3];
    unsigned int hashval;
    unsigned int g;
    for (hashval = 0; *symbol != '\0'; symbol++) {
      hashval = (hashval << 4) + *symbol;
      g = hashval&0xf0000000;
      if (g != 0) {
	hashval = hashval ^ (g >> 24);
	hashval = hashval ^ g;
      }
    }

    // Convert to base-26 'number'
    hashval %= HASHSIZE;
    word[0] = char(hashval / 26) + 'a';
    word[1] = char(hashval % 26) + 'a';
    word[2] = '\0';
    return (std::string(word));
  }
}

std::string Ioss::Utils::variable_name_kluge(const std::string &name,
					     size_t component_count, size_t copies,
					     size_t max_var_len)
{
  // This routine tries to shorten long variable names to an acceptable
  // length ('max_var_len' characters max).  If the name is already less than this
  // length, it is returned unchanged except for the appending of the hash...
  //
  // Since there is a (good) chance that two shortened names will match,
  // a 2-letter 'hash' code is appended to the end of the variable
  // name. This can be treated as a 2-digit base 26 number
  //
  // So, we shorten the name to a maximum of 'max_var_len-3' characters and append a
  // dot ('.') and 2 character hash.
  //
  // But, we also have to deal with the suffices that Ioex_DatabaseIO
  // appends on non-scalar values.  For the 'standard' types, the
  // maximum suffix is 4 characters (underscore + 1, 2, or 3 characters).
  // So...shorten name to maximum of 'max_var_len-3-{3|4|n}' characters depending on the
  // number of components.
  //
  // This function alo converts name to lowercase and converts spaces
  // to '_'

  // Width = 'max_var_len'.
  // Reserve space for suffix '_00...'
  // Reserve 3 for hash   '.xx'
  int hash_len = 3;
  int comp_len = 3; // _00
  int copy_len = 0;

  if (copies > 1) {
    assert(component_count % copies == 0);
    component_count /= copies;
  }

  if (component_count     <=      1)
    comp_len = 0;
  else if (component_count <    100)
    comp_len = 3;
  else if (component_count <   1000)
    comp_len = 4; //   _000
  else if (component_count <  10000)
    comp_len = 5; //  _0000
  else if (component_count < 100000)
    comp_len = 6; // _00000
  else {
    std::ostringstream errmsg;
    errmsg << "Variable '" << name << "' has " << component_count
	   << " components which is larger than the current maximum"
	   << " of 100,000. Please contact developer.";
    IOSS_ERROR(errmsg);
  }

  if (copies     <=      1)
    copy_len = 0;
  else if (copies <     10)
    copy_len = 2;
  else if (copies <    100)
    copy_len = 3; 
  else if (copies <   1000)
    copy_len = 4; //   _000
  else if (copies <  10000)
    copy_len = 5; //  _0000
  else if (copies < 100000)
    copy_len = 6; // _00000
  else {
    std::ostringstream errmsg;
    errmsg << "Variable '" << name << "' has " << copies
	   << " copies which is larger than the current maximum"
	   << " of 100,000. Please contact developer.";
    IOSS_ERROR(errmsg);
  }

  size_t maxlen = max_var_len - comp_len - copy_len;

  std::string new_str = name;
  if (name.length() <= maxlen) {
    // If name fits without kluging, then just use name as it is
    // without adding on the hash...
    std::transform(new_str.begin(), new_str.end(), new_str.begin(), to_lower);
    return new_str;
  } else {
    // Know that the name is too long, try to shorten. Need room for
    // hash now.
    maxlen -= hash_len;
    int len = name.length();

    // Take last 'maxlen' characters.  Motivation for this is that the
    // beginning of the composed (or generated) variable name is the
    // names of the mechanics and mechanics instances in which this
    // variable is nested, so they will be similar for all variables at
    // the same scope and the differences will occur at the variable
    // name level...
    //
    // However, there will likely be variables at the
    // same level but in different scope that have the same name which
    // would cause a clash, so we *hope* that the hash will make those
    // scope names unique...
    std::string s = std::string(name.c_str()).substr(len-maxlen,len);
    assert(s.length() <= maxlen);
    new_str = s;
  }

  // NOTE: The hash is not added if the name is not shortened. 
  std::string hash_string = two_letter_hash(name.c_str());
  new_str += std::string(".");
  new_str += hash_string;
  std::transform(new_str.begin(), new_str.end(), new_str.begin(), to_lower);
  return new_str;
}
