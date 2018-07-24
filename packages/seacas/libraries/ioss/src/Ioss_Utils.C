// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//
//     * Neither the name of NTESS nor the names of its
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

#include <Ioss_CodeTypes.h>
#include <Ioss_MeshCopyOptions.h>
#include <Ioss_Utils.h>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>
#include <sys/select.h>
#include <tokenize.h>
#include <vector>

#ifndef _WIN32
#include <sys/utsname.h>
#endif

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/kern_return.h> // for kern_return_t
#include <mach/mach.h>
#include <mach/message.h> // for mach_msg_type_number_t
#include <mach/task_info.h>
#endif

#if defined(BGQ_LWK) && defined(__linux__)
#include <spi/include/kernel/location.h>
#include <spi/include/kernel/memory.h>
#endif

#include <Ioss_SubSystem.h>

#include <fstream>

// For copy_database...
namespace {
  size_t MAX(size_t a, size_t b) { return b ^ ((a ^ b) & -static_cast<int>(a > b)); }
  size_t max_field_size = 0;

  struct DataPool
  {
    // Data space shared by most field input/output routines...
    std::vector<char>    data;
    std::vector<int>     data_int;
    std::vector<int64_t> data_int64;
    std::vector<double>  data_double;
    std::vector<Complex> data_complex;
#ifdef SEACAS_HAVE_KOKKOS
    Kokkos::View<char *>    data_view_char;
    Kokkos::View<int *>     data_view_int;
    Kokkos::View<int64_t *> data_view_int64;
    Kokkos::View<double *>  data_view_double;
    // Kokkos::View<Kokkos_Complex *> data_view_complex cannot be a global variable,
    // Since Kokkos::initialize() has not yet been called. Also, a Kokkos:View cannot
    // have type std::complex entities.
    Kokkos::View<char **>    data_view_2D_char;
    Kokkos::View<int **>     data_view_2D_int;
    Kokkos::View<int64_t **> data_view_2D_int64;
    Kokkos::View<double **>  data_view_2D_double;
    // Kokkos::View<Kokkos_Complex **> data_view_2D_complex cannot be a global variable,
    // Since Kokkos::initialize() has not yet been called. Also, a Kokkos:View cannot
    // have type std::complex entities.
    Kokkos::View<char **, Kokkos::LayoutRight, Kokkos::HostSpace> data_view_2D_char_layout_space;
    Kokkos::View<int **, Kokkos::LayoutRight, Kokkos::HostSpace>  data_view_2D_int_layout_space;
    Kokkos::View<int64_t **, Kokkos::LayoutRight, Kokkos::HostSpace>
        data_view_2D_int64_layout_space;
    Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::HostSpace>
        data_view_2D_double_layout_space;
    // Kokkos::View<Kokkos_Complex **, Kokkos::LayoutRight, Kokkos::HostSpace>
    // data_view_2D_complex_layout_space cannot be a global variable,
    // Since Kokkos::initialize() has not yet been called. Also, a Kokkos:View cannot
    // have type std::complex entities.
#endif
  };

  void show_step(int istep, double time, bool verbose, int rank);

  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, DataPool &pool,
                          bool debug, bool verbose, int rank);
  void transfer_structuredblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                                 bool verbose, int rank);
  void transfer_elementblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                              bool verbose, int rank);
  void transfer_edgeblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                           bool verbose, int rank);
  void transfer_faceblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                           bool verbose, int rank);
  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                         bool verbose, int rank);
  void transfer_edgesets(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                         bool verbose, int rank);
  void transfer_facesets(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                         bool verbose, int rank);
  void transfer_elemsets(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                         bool verbose, int rank);
  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                         bool verbose, int rank);
  void transfer_commsets(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                         bool verbose, int rank);
  void transfer_coordinate_frames(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                                  bool verbose, int rank);

  template <typename T>
  void transfer_fields(const std::vector<T *> &entities, Ioss::Region &output_region,
                       Ioss::Field::RoleType role, const Ioss::MeshCopyOptions &options, int rank);

  void transfer_fields(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                       Ioss::Field::RoleType role, const std::string &prefix = "");

  template <typename T>
  void transfer_field_data(const std::vector<T *> &entities, Ioss::Region &output_region,
                           DataPool &pool, Ioss::Field::RoleType role,
                           const Ioss::MeshCopyOptions &options);

  void transfer_field_data(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge, DataPool &pool,
                           Ioss::Field::RoleType role, const Ioss::MeshCopyOptions &options,
                           const std::string &prefix = "");

  void transfer_properties(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge);

  void transfer_qa_info(Ioss::Region &in, Ioss::Region &out);

  void transfer_field_data_internal(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                                    DataPool &pool, const std::string &field_name,
                                    const Ioss::MeshCopyOptions &options);

  template <typename INT>
  void set_owned_node_count(Ioss::Region &region, int my_processor, INT dummy);

  ////////////////////////////////////////////////////////////////////////
  inline int to_lower(int c) { return std::tolower(c); }
  inline int to_upper(int c) { return std::toupper(c); }

  bool is_separator(const char separator, const char value) { return separator == value; }

  size_t match(const char *name1, const char *name2)
  {
    size_t l1  = std::strlen(name1);
    size_t l2  = std::strlen(name2);
    size_t len = l1 < l2 ? l1 : l2;
    for (size_t i = 0; i < len; i++) {
      if (name1[i] != name2[i]) {
        while (i > 0 && (isdigit(name1[i - 1]) != 0) && (isdigit(name2[i - 1]) != 0)) {
          i--;
          // Back up to first non-digit so to handle "evar0000, evar0001, ..., evar 1123"
        }
        return i;
      }
    }
    return len;
  }

  // Split 'str' into 'tokens' based on the 'separator' character.
  // If 'str' starts with 1 or more 'separator', they are part of the
  // first token and not used for splitting.  If there are multiple
  // 'separator' characters in a row, then the first is used to split
  // and the subsequent 'separator' characters are put as leading
  // characters of the next token.
  // __this___is_a_string__for_tokens will split to 6 tokens:
  // __this __is a string _for tokens
  void field_tokenize(const std::string &str, const char separator,
                      std::vector<std::string> &tokens)
  {
    std::string curr_token;
    // Skip leading separators...
    size_t i = 0;
    while (i < str.length() && is_separator(separator, str[i])) {
      curr_token += str[i++];
    }
    for (; i < str.length(); ++i) {
      char curr_char = str[i];

      // determine if current character is a separator
      bool is_sep = is_separator(separator, curr_char);
      if (is_sep && curr_token != "") {
        // we just completed a token
        tokens.push_back(curr_token);
        curr_token.clear();
        while (i++ < str.length() && is_separator(separator, str[i])) {
          curr_token += str[i];
        }
        i--;
      }
      else if (!is_sep) {
        curr_token += curr_char;
      }
    }
    if (curr_token != "") {
      tokens.push_back(curr_token);
    }
  }

} // namespace

/** \brief Get formatted time and date strings.
 *
 *  Fill time_string and date_string with current time and date
 *  formatted as "HH:MM:SS" for time and "yy/mm/dd" or "yyyy/mm/dd"
 *  for date.
 *
 *  \param[out] time_string The formatted time string.
 *  \param[out] date_string The formatted date string.
 *  \param[in] length Use 8 for short-year date format, or 10 for long-year date format.
 */
void Ioss::Utils::time_and_date(char *time_string, char *date_string, size_t length)
{
  time_t     calendar_time = time(nullptr);
  struct tm *local_time    = localtime(&calendar_time);

  strftime(time_string, length, "%H:%M:%S", local_time);
  if (length == 8) {
    const char *fmt = "%y/%m/%d";
    strftime(date_string, length, fmt, local_time);
    date_string[8] = '\0';
  }
  else if (length >= 10) {
    strftime(date_string, length, "%Y/%m/%d", local_time);
    date_string[10] = '\0';
  }
  time_string[8] = '\0';
}

void Ioss::Utils::check_non_null(void *ptr, const char *type, const std::string &name,
                                 const std::string &func)
{
  if (ptr == nullptr) {
    std::ostringstream errmsg;
    errmsg << "INTERNAL ERROR: Could not find " << type << " '" << name << "'."
           << " Something is wrong in " << func << ". Please report.\n";
    IOSS_ERROR(errmsg);
  }
}

std::string Ioss::Utils::decode_filename(const std::string &filename, int processor,
                                         int num_processors)
{
  std::string decoded_filename(filename);
  // Current format for per-processor file names is:
  // PREFIX/basename.num_proc.cur_proc
  // the 'cur_proc' field is padded to be the same width as
  // the 'num_proc' field
  // Examples: basename.8.1, basename.64.03, basename.128.001

  // Create a std::string containing the total number of processors
  std::string num_proc   = std::to_string(num_processors);
  size_t      proc_width = num_proc.length();

  // Create a std::string containing the current processor number
  std::string cur_proc  = std::to_string(processor);
  size_t      cur_width = cur_proc.length();

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

size_t Ioss::Utils::get_number(const std::string &suffix)
{
  int  N       = 0;
  bool all_dig = suffix.find_first_not_of("0123456789") == std::string::npos;
  if (all_dig) {
    N = std::stoi(suffix);
  }
  return N;
}

int64_t Ioss::Utils::extract_id(const std::string &name_id)
{
  int64_t id = 0;

  std::vector<std::string> tokens = Ioss::tokenize(name_id, "_");
  if (tokens.size() > 1) {
    // Check whether last token is an integer...
    std::string str_id = tokens.back();
    id                 = get_number(str_id);
  }
  return id;
}

std::string Ioss::Utils::encode_entity_name(const std::string &entity_type, int64_t id)
{
  // ExodusII stores block, nodeset, and sideset ids as integers
  // Sierra   stores these as std::strings. The string is created by
  // concatenating the type, the character '_' and the id.

  std::string id_string   = std::to_string(id);
  std::string entity_name = entity_type;
  entity_name += "_";
  entity_name += id_string;
  return entity_name;
}

char **Ioss::Utils::get_name_array(size_t count, int size)
{
  auto names = new char *[count];
  for (size_t i = 0; i < count; i++) {
    names[i] = new char[size + 1];
    std::memset(names[i], '\0', size + 1);
  }
  return names;
}

void Ioss::Utils::delete_name_array(char **names, int count)
{
  for (int i = 0; i < count; i++) {
    delete[] names[i];
  }
  delete[] names;
}

/** \brief Process the base element type 'base' which has
 *         'nodes_per_element' nodes and a spatial dimension of 'spatial'
 *         into a form that the IO system can (hopefully) recognize.
 *
 *  Lowercases the name; converts spaces to '_', adds
 *  nodes_per_element at end of name (if not already there), and
 *  does some other transformations to remove some exodusII ambiguity.
 *
 *  \param[in] base The element base name.
 *  \param[in] nodes_per_element The number of nodes per element.
 *  \param[in] spatial The spatial dimension of the element.
 *  \returns The Ioss-formatted element name.
 */
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
  if (isdigit(*(type.rbegin())) == 0) {
    if (nodes_per_element > 1) {
      type += std::to_string(nodes_per_element);
    }
  }

  // Fixup an exodusII kludge.  For triangular elements, the same
  // name is used for 2D elements and 3D shell elements.  Convert
  // to unambiguous names for the IO Subsystem.  The 2D name
  // stays the same, the 3D name becomes 'trishell#'
  if (spatial == 3) {
    if (type == "triangle3") {
      type = "trishell3";
    }
    else if (type == "triangle4") {
      type = "trishell4";
    }
    else if (type == "triangle6") {
      type = "trishell6";
    }
    else if (type == "tri3") {
      type = "trishell3";
    }
    else if (type == "tri4") {
      type = "trishell4";
    }
    else if (type == "tri6") {
      type = "trishell6";
    }
  }

  if (spatial == 2) {
    if (type == "shell2") {
      type = "shellline2d2";
    }
    else if (type == "rod2" || type == "bar2" || type == "truss2") {
      type = "rod2d2";
    }
    else if (type == "shell3") {
      type = "shellline2d3";
    }
    else if (type == "bar3" || type == "rod3" || type == "truss3") {
      type = "rod2d3";
    }
  }

  if (std::strncmp(type.c_str(), "super", 5) == 0) {
    // A super element can have a varying number of nodes.  Create
    // an IO element type for this super element just so the IO
    // system can read a mesh containing super elements.  This
    // allows the "omit volume" command to be used in the Sierra
    // applications to skip creating a corresponding element block
    // in the application.
    type = "super" + std::to_string(nodes_per_element);
  }
  return type;
}

/* \brief Throw a runtime exception with message "I/O abort".
 */
void Ioss::Utils::abort()
{
  std::ostringstream errmsg("I/O abort");
  IOSS_ERROR(errmsg);
}

/** \brief Get a filename relative to the specified working directory (if any)
 *         of the current execution.
 *
 *  Working_directory must end with '/' or be empty.
 *
 *  \param[in] relative_filename The file path to be appended to the working directory path.
 *  \param[in] type The file type. "generated" file types are treated differently.
 *  \param[in] working_directory the path to which the relative_filename path is appended.
 *  \returns The full path (working_directory + relative_filename)
 */
std::string Ioss::Utils::local_filename(const std::string &relative_filename,
                                        const std::string &type,
                                        const std::string &working_directory)
{
  if (relative_filename[0] == '/' || type == "generated" || working_directory.empty()) {
    return relative_filename;
  }
  std::string filename = working_directory;
  filename += relative_filename;
  return filename;
}

int Ioss::Utils::field_warning(const Ioss::GroupingEntity *ge, const Ioss::Field &field,
                               const std::string &inout)
{
  IOSS_WARNING << ge->type_string() << " '" << ge->name() << "'. Unknown " << inout << " field '"
               << field.get_name() << "'\n";
  return -4;
}

namespace {
  const Ioss::VariableType *match_composite_field(char **names, Ioss::IntVector &which_names,
                                                  const char suffix_separator)
  {
    // ASSUME: Fields are in order...
    // The field we are trying to match will be a composite field of
    // type base_x_1, base_y_1, base_z_1, ...., base_y_3, base_z_3.
    // The composite field type currently always has a numeric Real[N]
    // type field as the last suffix and the other field as the first
    // suffix.
    // If we take the last suffix of the last name, it should give us
    // the 'N' in the Real[N] field.  Dividing 'which_names.size()' by
    // 'N' will give the number of components in the inner field.

    char suffix[2];
    suffix[0] = suffix_separator;
    suffix[1] = 0;

    std::vector<std::string> tokens =
        Ioss::tokenize(names[which_names[which_names.size() - 1]], suffix);

    if (tokens.size() <= 2) {
      return nullptr;
    }

    assert(tokens.size() > 2);

    // Check that suffix is a number -- all digits
    size_t N = Ioss::Utils::get_number(tokens[tokens.size() - 1]);

    if (N == 0) {
      return nullptr;
    }

    if (which_names.size() % N != 0) {
      return nullptr;
    }

    size_t inner_token = tokens.size() - 2;
    size_t inner_comp  = which_names.size() / N;

    // Gather the first 'inner_ccomp' inner field suffices...
    std::vector<Ioss::Suffix> suffices;
    for (size_t i = 0; i < inner_comp; i++) {
      std::vector<std::string> ltokens = Ioss::tokenize(names[which_names[i]], suffix);
      // The second-last token is the suffix for this component...
      Ioss::Suffix tmp(ltokens[inner_token]);
      suffices.push_back(tmp);
    }

    // check that the suffices on the next copies of the inner field
    // match the first copy...
    size_t j = inner_comp;
    for (size_t copy = 1; copy < N; copy++) {
      for (size_t i = 0; i < inner_comp; i++) {
        std::vector<std::string> ltokens = Ioss::tokenize(names[which_names[j++]], suffix);
        // The second-last token is the suffix for this component...
        if (suffices[i] != ltokens[inner_token]) {
          return nullptr;
        }
      }
    }

    // All 'N' copies of the inner field match, now see the
    // suffices actually defines a field...
    const Ioss::VariableType *type = Ioss::VariableType::factory(suffices);
    if (type != nullptr) {
      type = Ioss::VariableType::factory(type->name(), N);
    }
    return type;
  }

  const Ioss::VariableType *match_single_field(char **names, Ioss::IntVector &which_names,
                                               const char suffix_separator)
  {
    // Strip off the suffix from each name indexed in 'which_names'
    // and see if it defines a valid type...
    std::vector<Ioss::Suffix> suffices;

    char suffix[2];
    suffix[0] = suffix_separator;
    suffix[1] = 0;

    for (int which_name : which_names) {
      std::vector<std::string> tokens     = Ioss::tokenize(names[which_name], suffix);
      size_t                   num_tokens = tokens.size();

      // The last token is the suffix for this component...
      Ioss::Suffix tmp(tokens[num_tokens - 1]);
      suffices.push_back(tmp);
    }
    const Ioss::VariableType *type = Ioss::VariableType::factory(suffices);
    return type;
  }

  Ioss::Field get_next_field(char **names, int num_names, size_t count,
                             Ioss::Field::RoleType fld_role, const char suffix_separator,
                             int *truth_table)
  {
    // NOTE: 'names' are all lowercase at this point.

    // Assumption 1: To convert to a non-SCALAR type, the variable name
    // must have an field_suffix_sep in the name separating the suffixes from
    // the main name.

    // Find first unused name (used names have '\0' as first character...
    int  index       = 0;
    bool found_valid = false;
    for (index = 0; index < num_names; index++) {
      assert(truth_table == nullptr || truth_table[index] == 1 || truth_table[index] == 0);
      if ((truth_table == nullptr || truth_table[index] == 1) && names[index][0] != '\0') {
        found_valid = true;
        break;
      }
    }

    if (!found_valid) {
      // Return an invalid field...
      return Ioss::Field("", Ioss::Field::INVALID, IOSS_SCALAR(), fld_role, 1);
    }

    // At this point, name[index] should be a valid potential field
    // name and all names[i] with i < index are either already used or
    // not valid for this grouping entity (truth_table entry == 0).
    assert(index < num_names && names[index][0] != '\0' &&
           (truth_table == nullptr || truth_table[index] == 1));
    char *name = names[index];

    // Split the name up into tokens separated by the
    // 'suffix_separator'.  Note that the basename itself could
    // contain a suffix_separator (back_stress_xx or
    // back_stress_xx_01). Need to ignore embedded separators
    // (back_stress) and also recognize composite variable types
    // (back_stress_xx_01). At the current time, a composite variable
    // type can only contain two non-composite variable types, so we
    // only need to look to be concerned with the last 1 or 2 tokens...
    std::vector<std::string> tokens;
    field_tokenize(name, suffix_separator, tokens);
    size_t num_tokens = tokens.size();

    // Check that tokenizer did not return empty tokens...
    bool invalid = tokens[0].empty() || tokens[num_tokens - 1].empty();
    if (num_tokens == 1 || invalid) {
      // It is not a (Sierra-generated) name for a non-SCALAR variable
      // Return a SCALAR field
      Ioss::Field field(name, Ioss::Field::REAL, IOSS_SCALAR(), fld_role, count);
      field.set_index(index);
      names[index][0] = '\0';
      return field;
    }

    // KNOW: num_tokens > 1 at this point.  Possible that we still
    // just have a scalar with one or more embedded separator characters...
    int suffix_size = 1;
    if (num_tokens > 2) {
      suffix_size = 2;
    }

    // If num_tokens > 2, then we can potentially have a composite
    // variable type which would have a double suffix (_xx_01).

    // Gather all names which match in the first
    // (num_tokens-suffix_size) tokens and see if their suffices form
    // a valid variable type...
    while (suffix_size > 0) {
      Ioss::IntVector which_names; // Contains index of names that
      // potentially match as components
      // of a higher-order type.

      std::string base_name = tokens[0];
      for (size_t i = 1; i < num_tokens - suffix_size; i++) {
        base_name += suffix_separator;
        base_name += tokens[i];
      }
      base_name += suffix_separator;
      size_t bn_len = base_name.length(); // Length of basename portion only
      size_t length = std::strlen(name);  // Length of total name (with suffix)

      // Add the current name...
      which_names.push_back(index);

      // Gather all other names that are valid for this entity, and
      // have the same overall length and match in the first 'bn_len'
      // characters.
      //
      // Check that they have the same number of tokens,
      // It is possible that the first name(s) that match with two
      // suffices have a basename that match other names with only a
      // single suffix lc_cam_x, lc_cam_y, lc_sfarea.
      for (int i = index + 1; i < num_names; i++) {
        char *                   tst_name = names[i];
        std::vector<std::string> subtokens;
        field_tokenize(tst_name, suffix_separator, subtokens);
        if ((truth_table == nullptr || truth_table[i] == 1) && // Defined on this entity
            std::strlen(tst_name) == length &&                 // names must be same length
            std::strncmp(name, tst_name, bn_len) == 0 &&       // base portion must match
            subtokens.size() == num_tokens) {
          which_names.push_back(i);
        }
      }

      const Ioss::VariableType *type = nullptr;
      if (suffix_size == 2) {
        if (which_names.size() > 1) {
          type = match_composite_field(names, which_names, suffix_separator);
        }
      }
      else {
        assert(suffix_size == 1);
        type = match_single_field(names, which_names, suffix_separator);
      }

      if (type != nullptr) {
        // A valid variable type was recognized.
        // Mark the names which were used so they aren't used for another field on this entity.
        // Create a field of that variable type.
        assert(type->component_count() == static_cast<int>(which_names.size()));
        Ioss::Field field(base_name.substr(0, bn_len - 1), Ioss::Field::REAL, type, fld_role,
                          count);
        field.set_index(index);
        for (const auto &which_name : which_names) {
          names[which_name][0] = '\0';
        }
        return field;
      }
      if (suffix_size == 1) {
        Ioss::Field field(name, Ioss::Field::REAL, IOSS_SCALAR(), fld_role, count);
        field.set_index(index);
        names[index][0] = '\0';
        return field;
      }

      suffix_size--;
    }
    return Ioss::Field("", Ioss::Field::INVALID, IOSS_SCALAR(), fld_role, 1);
  }

  // common
  bool define_field(size_t nmatch, size_t match_length, char **names,
                    std::vector<Ioss::Suffix> &suffices, size_t entity_count,
                    Ioss::Field::RoleType fld_role, std::vector<Ioss::Field> &fields)
  {
    // Try to define a field of size 'nmatch' with the suffices in 'suffices'.
    // If this doesn't define a known field, then assume it is a scalar instead
    // and return false.
    if (nmatch > 1) {
      const Ioss::VariableType *type = Ioss::VariableType::factory(suffices);
      if (type == nullptr) {
        nmatch = 1;
      }
      else {
        char *name         = names[0];
        name[match_length] = '\0';
        Ioss::Field field(name, Ioss::Field::REAL, type, fld_role, entity_count);
        if (field.is_valid()) {
          fields.push_back(field);
        }
        for (size_t j = 0; j < nmatch; j++) {
          names[j][0] = '\0';
        }
        return true;
      }
    }

    // NOTE: nmatch could be reset inside previous if block.
    // This is not an 'else' block, it is a new if block.
    if (nmatch == 1) {
      Ioss::Field field(names[0], Ioss::Field::REAL, IOSS_SCALAR(), fld_role, entity_count);
      if (field.is_valid()) {
        fields.push_back(field);
      }
      names[0][0] = '\0';
      return false;
    }
    return false; // Can't get here...  Quiet the compiler
  }
} // namespace
// Read scalar fields off an input database and determine whether
// they are components of a higher order type (vector, tensor, ...).
// This routine is used if there is no field component separator.  E.g.,
// fieldx, fieldy, fieldz instead of field_x field_y field_z

void Ioss::Utils::get_fields(int64_t entity_count, // The number of objects in this entity.
                             char ** names,        // Raw list of field names from exodus
                             size_t  num_names,    // Number of names in list
                             Ioss::Field::RoleType fld_role, // Role of field
                             bool enable_field_recognition, const char suffix_separator,
                             int *local_truth, // Truth table for this entity;
                             // null if not applicable.
                             std::vector<Ioss::Field> &fields) // The fields that were found.
{
  if (!enable_field_recognition) {
    // Create a separate field for each name.
    for (size_t i = 0; i < num_names; i++) {
      if (local_truth == nullptr || local_truth[i] == 1) {
        Ioss::Field field(names[i], Ioss::Field::REAL, IOSS_SCALAR(), fld_role, entity_count);
        fields.push_back(field);
        names[i][0] = '\0';
      }
    }
  }
  else if (suffix_separator != 0) {
    while (true) {
      // NOTE: 'get_next_field' determines storage type (vector, tensor,...)
      Ioss::Field field =
          get_next_field(names, num_names, entity_count, fld_role, suffix_separator, local_truth);
      if (field.is_valid()) {
        fields.push_back(field);
      }
      else {
        break;
      }
    }
  }
  else {
    size_t                    nmatch = 1;
    size_t                    ibeg   = 0;
    size_t                    pmat   = 0;
    std::vector<Ioss::Suffix> suffices;
  top:

    while (ibeg + nmatch < num_names) {
      if (local_truth != nullptr) {
        while (ibeg < num_names && local_truth[ibeg] == 0) {
          ibeg++;
        }
      }
      for (size_t i = ibeg + 1; i < num_names; i++) {
        size_t mat = match(names[ibeg], names[i]);
        if (local_truth != nullptr && local_truth[i] == 0) {
          mat = 0;
        }

        // For all fields, the total length of the name is the same
        // for all components of that field.  The 'basename' of the
        // field will also be the same for all cases.
        //
        // It is possible that the length of the match won't be the
        // same for all components of a field since the match may
        // include a portion of the suffix; (sigxx, sigxy, sigyy
        // should match only 3 characters of the basename (sig), but
        // sigxx and sigxy will match 4 characters) so consider a
        // valid match if the match length is >= previous match length.
        if ((std::strlen(names[ibeg]) == std::strlen(names[i])) && mat > 0 &&
            (pmat == 0 || mat >= pmat)) {
          nmatch++;
          if (nmatch == 2) {
            // Get suffix for first field in the match
            pmat = mat;
            Ioss::Suffix tmp(&names[ibeg][pmat]);
            suffices.push_back(tmp);
          }
          // Get suffix for next fields in the match
          Ioss::Suffix tmp(&names[i][pmat]);
          suffices.push_back(tmp);
        }
        else {

          bool multi_component =
              define_field(nmatch, pmat, &names[ibeg], suffices, entity_count, fld_role, fields);
          if (!multi_component) {
            // Although we matched multiple suffices, it wasn't a
            // higher-order field, so we only used 1 name instead of
            // the 'nmatch' we thought we might use.
            i = ibeg + 1;
          }

          // Cleanout the suffices vector.
          clear(suffices);

          // Reset for the next time through the while loop...
          nmatch = 1;
          pmat   = 0;
          ibeg   = i;
          break;
        }
      }
    }
    // We've gone through the entire list of names; see if what we
    // have forms a multi-component field; if not, then define a
    // scalar field and jump up to the loop again to handle the others
    // that had been gathered.
    if (ibeg < num_names) {
      if (local_truth == nullptr || local_truth[ibeg] == 1) {
        bool multi_component =
            define_field(nmatch, pmat, &names[ibeg], suffices, entity_count, fld_role, fields);
        clear(suffices);
        if (nmatch > 1 && !multi_component) {
          ibeg++;
          goto top;
        }
      }
      else {
        ibeg++;
        goto top;
      }
    }
  }
}

/** \brief Get a string containing 'uname' output.
 *
 *  This output contains information about the current computing platform.
 *  This is used as information data in the created results file to help
 *  in tracking when/where/... the file was created.
 *
 *  \returns The platform information string.
 */
std::string Ioss::Utils::platform_information()
{
#ifndef _WIN32
  struct utsname sys_info
  {
  };
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

/** \brief Return amount of memory being used on this processor */
size_t Ioss::Utils::get_memory_info()
{
  size_t memory_usage = 0;
#if defined(__APPLE__) && defined(__MACH__)
  static size_t               original = 0;
  kern_return_t               error;
  mach_msg_type_number_t      outCount;
  mach_task_basic_info_data_t taskinfo{};

  taskinfo.virtual_size = 0;
  outCount              = MACH_TASK_BASIC_INFO_COUNT;
  error                 = task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                    reinterpret_cast<task_info_t>(&taskinfo), &outCount);
  if (error == KERN_SUCCESS) {
    // type is mach_vm_size_t
    if (original == 0) {
      original = taskinfo.virtual_size;
    }
    memory_usage = taskinfo.virtual_size - original;
  }
#elif __linux__
#if defined(BGQ_LWK)
  uint64_t    heap;
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  memory_usage      = heap;
#else
  std::string line(128, '\0');

  /* Read memory size data from /proc/self/status
   * run "man proc" to get info on the contents of /proc/self/status
   */
  std::ifstream proc_status("/proc/self/status");
  if (!proc_status) {
    return memory_usage;
  }

  while (1) {
    if (!std::getline(proc_status, line)) {
      return memory_usage;
    }

    if (line.substr(0, 6) == "VmRSS:") {
      std::string        vmrss = line.substr(7);
      std::istringstream iss(vmrss);
      iss >> memory_usage;
      memory_usage *= 1024;
      break;
    }
  }
  proc_status.close();
#endif
#endif
  return memory_usage;
}

size_t Ioss::Utils::get_hwm_memory_info()
{
  size_t memory_usage = 0;
#if defined(__linux__)
#if defined(BGQ_LWK)

#else
  std::string line(128, '\0');

  /* Read memory size data from /proc/self/status
   * run "man proc" to get info on the contents of /proc/self/status
   */
  std::ifstream proc_status("/proc/self/status");
  if (!proc_status)
    return memory_usage;

  while (1) {

    if (!std::getline(proc_status, line))
      return memory_usage;
    if (line.substr(0, 6) == "VmHWM:") {
      std::string        vmrss = line.substr(7);
      std::istringstream iss(vmrss);
      iss >> memory_usage;
      memory_usage *= 1024;
      break;
    }
  }
  proc_status.close();
#endif
#endif
  return memory_usage;
}

/** \brief Determine whether an entity has the property "omitted."
 *
 *  \param[in] block The entity.
 *  \returns True if the entity has the property "omitted."
 */
bool Ioss::Utils::block_is_omitted(Ioss::GroupingEntity *block)
{
  bool omitted = false;
  if (block->property_exists("omitted")) {
    omitted = (block->get_property("omitted").get_int() == 1);
  }
  return omitted;
}

void Ioss::Utils::calculate_sideblock_membership(IntVector &            face_is_member,
                                                 const Ioss::SideBlock *ef_blk,
                                                 size_t int_byte_size, const void *element,
                                                 const void *sides, int64_t number_sides,
                                                 const Ioss::Region *region)
{
  assert(ef_blk != nullptr);

  face_is_member.reserve(number_sides);

  const ElementTopology *unknown = Ioss::ElementTopology::factory("unknown");

  // Topology of faces in this face block...
  const ElementTopology *ftopo = ef_blk->topology();

  // Topology of parent element for faces in this face block
  const ElementTopology *parent_topo = ef_blk->parent_element_topology();

  // If split by element block then parent_block will be non-nullptr
  const ElementBlock *parent_block = ef_blk->parent_element_block();

  // The element block containing the face we are working on...
  Ioss::ElementBlock *block = nullptr;

  // Topology of face/edge in current element block
  const ElementTopology *common_ftopo = nullptr;

  // Topology of elements in the element block containing this element
  const ElementTopology *block_topo = nullptr;

  // Topology of the face we are currently working with...
  const ElementTopology *topo = nullptr;

  // The element side that the current face is on the element...
  int64_t current_side = -1;

  if (number_sides > 0 && (element == nullptr || sides == nullptr)) {
    std::ostringstream errmsg;
    errmsg << "INTERNAL ERROR: null element or sides pointer passed to "
           << "Ioss::Utils::calculate_sideblock_membership function.";
    IOSS_ERROR(errmsg);
  }

  for (int64_t iel = 0; iel < number_sides; iel++) {
    int64_t elem_id = 0;
    int64_t side_id = 0;
    if (int_byte_size == 4) {
      elem_id = ((int *)element)[iel];
      side_id = ((int *)sides)[iel];
    }
    else {
      elem_id = ((int64_t *)element)[iel];
      side_id = ((int64_t *)sides)[iel];
    }

    // Get the element block containing this face...
    if (block == nullptr || !block->contains(elem_id)) {
      block      = region->get_element_block(elem_id);
      block_topo = block->topology();
      // nullptr if hetero face/edge on element
      common_ftopo = block->topology()->boundary_type(0);
      if (common_ftopo != nullptr) {
        topo = common_ftopo;
      }
      current_side = -1;
    }

    // If the element topology of the element block containing this
    // face has heterogeneous topology (eg. wedge), then determine the
    // topology corresponding to the current side..
    if (common_ftopo == nullptr && side_id != current_side) {
      current_side = side_id;
      topo         = block->topology()->boundary_type(side_id);
    }

    bool face_topo_match  = ftopo == unknown || topo == ftopo;
    bool block_topo_match = parent_topo == unknown || block_topo == parent_topo;
    // See if the face topology and the parent element topology for
    // the current face match the topology associated with this face block.
    if (face_topo_match && block_topo_match && (parent_block == nullptr || parent_block == block) &&
        !block_is_omitted(block)) {
      // This face/edge  belongs in the face/edge block
      face_is_member.push_back(1);
    }
    else {
      face_is_member.push_back(0);
    }
  }
}

/** \brief Get the appropriate index offset for the sides of elements in a SideBlock.
 *
 *  And yet another idiosyncrasy of sidesets...
 *  The side of an element (especially shells) can be
 *  either a face or an edge in the same sideset.  The
 *  ordinal of an edge is (local_edge_number+numfaces) on the
 *  database, but needs to be (local_edge_number) for Sierra...
 *
 *  If the sideblock has a "parent_element_topology" and a
 *  "topology", then we can determine whether to offset the
 *  side ordinals...
 *
 *  \param[in] sb Compute the offset for element sides in this SideBlock
 *  \returns The offset.
 */
int64_t Ioss::Utils::get_side_offset(const Ioss::SideBlock *sb)
{

  const Ioss::ElementTopology *side_topo   = sb->topology();
  const Ioss::ElementTopology *parent_topo = sb->parent_element_topology();
  int64_t                      side_offset = 0;
  if ((side_topo != nullptr) && (parent_topo != nullptr)) {
    int side_topo_dim = side_topo->parametric_dimension();
    int elem_topo_dim = parent_topo->parametric_dimension();
    int elem_spat_dim = parent_topo->spatial_dimension();

    if (side_topo_dim + 1 < elem_spat_dim && side_topo_dim < elem_topo_dim) {
      side_offset = parent_topo->number_faces();
    }
  }
  return side_offset;
}

unsigned int Ioss::Utils::hash(const std::string &name)
{
  // Hash function from Aho, Sethi, Ullman "Compilers: Principles,
  // Techniques, and Tools.  Page 436

  const char * symbol = name.c_str();
  unsigned int hashval;
  unsigned int g;
  for (hashval = 0; *symbol != '\0'; symbol++) {
    hashval = (hashval << 4) + *symbol;
    g       = hashval & 0xf0000000;
    if (g != 0) {
      hashval = hashval ^ (g >> 24);
      hashval = hashval ^ g;
    }
  }
  return hashval;
}

double Ioss::Utils::timer()
{
#ifdef SEACAS_HAVE_MPI
  return MPI_Wtime();
#else
  static auto begin = std::chrono::high_resolution_clock::now();

  auto now = std::chrono::high_resolution_clock::now();
  return std::chrono::duration<double>(now - begin).count();
#endif
}

/** \brief Convert an input file to a vector of strings containing one string for each line of the
 * file.
 *
 *  Should only be called by a single processor or each processor will be accessing the file
 *  at the same time...
 *
 *  \param[in] file_name The name of the file.
 *  \param[out] lines The vector of strings containing the lines of the file
 *  \param[in] max_line_length The maximum number of characters in any line of the file.
 */
void Ioss::Utils::input_file(const std::string &file_name, std::vector<std::string> *lines,
                             size_t max_line_length)
{
  // Create an ifstream for the input file. This does almost the same
  // function as sierra::Env::input() except this is for a single
  // processor and the sierra::Env::input() is for parallel...

  if (file_name.length() != 0) {
    // Open the file and read into the vector...
    std::string   input_line;
    std::ifstream infile(file_name);
    lines->push_back(file_name.substr(0, max_line_length));
    while (!std::getline(infile, input_line).fail()) {
      if (max_line_length == 0 || input_line.length() <= max_line_length) {
        lines->push_back(input_line);
      }
      else {
        // Split the line into pieces of length "max_line_length-1"
        // and append a "\" to all but the last. Don't worry about
        // splitting at whitespace...
        size_t ibeg = 0;
        do {
          std::string sub = input_line.substr(ibeg, max_line_length - 1);
          if (ibeg + max_line_length - 1 < input_line.length()) {
            sub += "\\";
          }
          lines->push_back(sub);
          ibeg += max_line_length - 1;
        } while (ibeg < input_line.length());
      }
    }
  }
}

/** \brief Case-insensitive string comparison.
 *
 *  \param[in] s1 First string
 *  \param[in] s2 Second string
 *  \returns 0 if strings are equal, nonzero otherwise.
 */
int Ioss::Utils::case_strcmp(const std::string &s1, const std::string &s2)
{
  const char *c1 = s1.c_str();
  const char *c2 = s2.c_str();
  for (;; c1++, c2++) {
    if (std::tolower(*c1) != std::tolower(*c2)) {
      return (std::tolower(*c1) - std::tolower(*c2));
    }
    if (*c1 == '\0') {
      return 0;
    }
  }
}

/** \brief Convert a string to upper case.
 *
 *  \param[in] name The string to convert.
 *  \returns The converted string.
 */
std::string Ioss::Utils::uppercase(std::string name)
{
  std::transform(name.begin(), name.end(), name.begin(), to_upper);
  return name;
}

/** \brief Convert a string to lower case.
 *
 *  \param[in] name The string to convert.
 *  \returns The converted string.
 */
std::string Ioss::Utils::lowercase(std::string name)
{
  std::transform(name.begin(), name.end(), name.begin(), to_lower);
  return name;
}

/** \brief Check whether property 'prop_name' exists and if so, set 'prop_value'
 *
 * based on the property value.  Either "TRUE", "YES", "ON", or nonzero for true;
 * or "FALSE", "NO", "OFF", or 0 for false.
 * \param[in] properties the Ioss::PropertyManager containing the properties to be checked.
 * \param[in] prop_name the name of the property to check whether it exists and if so, set its
 * value.
 * \param[out] prop_value if prop_name exists and has a valid value, set prop_value accordingly.
 * \returns true/false depending on whether property found and value set.
 */

bool Ioss::Utils::check_set_bool_property(const Ioss::PropertyManager &properties,
                                          const std::string &prop_name, bool &prop_value)
{
  bool found_property = false;
  if (properties.exists(prop_name)) {
    found_property = true;
    if (properties.get(prop_name).get_type() == Ioss::Property::INTEGER) {
      prop_value = properties.get(prop_name).get_int() != 0;
    }
    else {
      std::string yesno = Ioss::Utils::uppercase(properties.get(prop_name).get_string());
      if (yesno == "TRUE" || yesno == "YES" || yesno == "ON") {
        prop_value = true;
      }
      else if (yesno == "FALSE" || yesno == "NO" || yesno == "OFF") {
        prop_value = false;
      }
      else {
        std::ostringstream errmsg;
        errmsg << "ERROR: Unrecognized value found for " << prop_name << ". Found '" << yesno
               << "' which is not one of TRUE|FALSE|YES|NO|ON|OFF";
        IOSS_ERROR(errmsg);
      }
    }
  }
  return found_property;
}

/** \brief Convert a string to lower case, and convert spaces to '_'.
 *
 *  The conversion is performed in place.
 *
 *  \param[in,out] name On input, the string to convert. On output, the converted string.
 *
 */
void Ioss::Utils::fixup_name(char *name)
{
  assert(name != nullptr);

  size_t len = std::strlen(name);
  for (size_t i = 0; i < len; i++) {
    name[i] = static_cast<char>(tolower(name[i])); // guaranteed(?) to be ascii...
    if (name[i] == ' ') {
      name[i] = '_';
    }
  }
}

/** \brief Convert a string to lower case, and convert spaces to '_'.
 *
 *  The conversion is performed in place.
 *
 *  \param[in,out] name On input, the string to convert. On output, the converted string.
 *
 */
void Ioss::Utils::fixup_name(std::string &name)
{
  name = Ioss::Utils::lowercase(name);

  size_t len = name.length();
  for (size_t i = 0; i < len; i++) {
    if (name[i] == ' ') {
      name[i] = '_';
    }
  }
}

namespace {

  /** \brief Hash function from Aho, Sethi, Ullman "Compilers: Principles,
   *         Techniques, and Tools.  Page 436
   */
  std::string two_letter_hash(const char *symbol)
  {
    const int    HASHSIZE = 673; // Largest prime less than 676 (26*26)
    char         word[3];
    unsigned int hashval;
    unsigned int g;
    for (hashval = 0; *symbol != '\0'; symbol++) {
      hashval = (hashval << 4) + *symbol;
      g       = hashval & 0xf0000000;
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
} // namespace

/** \brief Tries to shorten long variable names to an acceptable length, and converts to
 *         lowercase and spaces to '_'
 *
 *   Many databases have a maximum length for variable names which can
 *   cause a problem with variable name length.
 *
 *   This routine tries to shorten long variable names to an acceptable
 *   length ('max_var_len' characters max).  If the name is already less than
 *   this length, it is returned unchanged except for the appending of the hash...
 *
 *   Since there is a (good) chance that two shortened names will match,
 *   a 2-letter 'hash' code is appended to the end of the variable
 *   name. This can be treated as a 2-digit base 26 number
 *
 *   So, we shorten the name to a maximum of 'max_var_len-3' characters and
 *   append a dot ('.') and 2 character hash.
 *
 *   But, we also have to deal with the suffices that Ioex_DatabaseIO
 *   appends on non-scalar values.  For the 'standard' types, the
 *   maximum suffix is 4 characters (underscore + 1, 2, or 3 characters).
 *   So...shorten name to maximum of 'max_var_len-3-{3|4|n}' characters
 *   depending on the number of components.
 *
 *   This function also converts name to lowercase and converts spaces
 *   to '_'.
 */
std::string Ioss::Utils::variable_name_kluge(const std::string &name, size_t component_count,
                                             size_t copies, size_t max_var_len)
{

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

  if (component_count <= 1) {
    comp_len = 0;
  }
  else if (component_count < 100) {
    comp_len = 3;
  }
  else if (component_count < 1000) {
    comp_len = 4; //   _000
  }
  else if (component_count < 10000) {
    comp_len = 5; //  _0000
  }
  else if (component_count < 100000) {
    comp_len = 6; // _00000
  }
  else {
    std::ostringstream errmsg;
    errmsg << "Variable '" << name << "' has " << component_count
           << " components which is larger than the current maximum"
           << " of 100,000. Please contact developer.";
    IOSS_ERROR(errmsg);
  }

  if (copies <= 1) {
    copy_len = 0;
  }
  else if (copies < 10) {
    copy_len = 2;
  }
  else if (copies < 100) {
    copy_len = 3;
  }
  else if (copies < 1000) {
    copy_len = 4; //   _000
  }
  else if (copies < 10000) {
    copy_len = 5; //  _0000
  }
  else if (copies < 100000) {
    copy_len = 6; // _00000
  }
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
  }
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
  std::string s = std::string(name).substr(len - maxlen, len);
  assert(s.length() <= maxlen);
  new_str = s;

  // NOTE: The hash is not added if the name is not shortened.
  std::string hash_string = two_letter_hash(name.c_str());
  new_str += std::string(".");
  new_str += hash_string;
  std::transform(new_str.begin(), new_str.end(), new_str.begin(), to_lower);
  return new_str;
}

/** \brief Create a nominal mesh for use in history databases.
 *
 *  The model for a history file is a single sphere element (1 node, 1 element).
 *  This is needed for some applications that read this file that require a
 *  "mesh" even though a history file is just a collection of global variables
 *  with no real mesh. This routine will add the mesh portion to a history file.
 *
 *  \param[in,out] region The region on which the nominal mesh is to be defined.
 */
void Ioss::Utils::generate_history_mesh(Ioss::Region *region)
{
  Ioss::DatabaseIO *db = region->get_database();
  if (db->parallel_rank() == 0) {
    region->begin_mode(Ioss::STATE_DEFINE_MODEL);

    // Node Block
    Ioss::NodeBlock *nb = new Ioss::NodeBlock(db, "nodeblock_1", 1, 3);
    region->add(nb);

    // Element Block
    Ioss::ElementBlock *eb = new Ioss::ElementBlock(db, "e1", "sphere", 1);
    eb->property_add(Ioss::Property("id", 1));
    eb->property_add(Ioss::Property("guid", 1));
    region->add(eb);
    region->end_mode(Ioss::STATE_DEFINE_MODEL);

    region->begin_mode(Ioss::STATE_MODEL);
    static double coord[3] = {1.1, 2.2, 3.3};
    static int    ids[1]   = {1};
    nb->put_field_data("ids", ids, sizeof(int));
    nb->put_field_data("mesh_model_coordinates", coord, 3 * sizeof(double));

    static int connect[1] = {1};
    eb->put_field_data("ids", ids, sizeof(int));
    eb->put_field_data("connectivity", connect, 1 * sizeof(int));

    region->end_mode(Ioss::STATE_MODEL);
  }
}

namespace {
  const int tab64[64] = {63, 0,  58, 1,  59, 47, 53, 2,  60, 39, 48, 27, 54, 33, 42, 3,
                         61, 51, 37, 40, 49, 18, 28, 20, 55, 30, 34, 11, 43, 14, 22, 4,
                         62, 57, 46, 52, 38, 26, 32, 41, 50, 36, 17, 19, 29, 10, 13, 21,
                         56, 45, 25, 31, 35, 16, 9,  12, 44, 24, 15, 8,  23, 7,  6,  5};
}

int Ioss::Utils::log_power_2(uint64_t value)
{
  assert(value > 0);
  value = (value << 1) - 1;
  value |= value >> 1;
  value |= value >> 2;
  value |= value >> 4;
  value |= value >> 8;
  value |= value >> 16;
  value |= value >> 32;
  return tab64[((uint64_t)((value - (value >> 1)) * 0x07EDD5E59A4E28C2)) >> 58];
}

void Ioss::Utils::copy_database(Ioss::Region &region, Ioss::Region &output_region,
                                Ioss::MeshCopyOptions &options)
{
  DataPool data_pool;

  bool              memory_stats = options.memory_statistics;
  Ioss::DatabaseIO *dbi          = region.get_database();

  int rank = dbi->util().parallel_rank();

  bool appending = output_region.get_database()->open_create_behavior() == Ioss::DB_APPEND;

  if (!appending) {

    if (options.debug && rank == 0) {
      std::cerr << "DEFINING MODEL ... \n";
    }
    if (memory_stats) {
      dbi->util().progress("DEFINING MODEL");
    }
    if (!output_region.begin_mode(Ioss::STATE_DEFINE_MODEL)) {
      if (options.verbose && rank == 0) {
        std::cerr << "ERROR: Could not put output region into define model state\n";
      }
      std::exit(EXIT_FAILURE);
    }

    // Get all properties of input database...
    transfer_properties(&region, &output_region);
    transfer_qa_info(region, output_region);

    transfer_nodeblock(region, output_region, data_pool, options.debug, options.verbose, rank);

#ifdef SEACAS_HAVE_MPI
    // This also assumes that the node order and count is the same for input
    // and output regions... (This is checked during nodeset output)
    if (output_region.get_database()->needs_shared_node_information()) {
      if (options.ints_64_bit)
        set_owned_node_count(region, rank, (int64_t)0);
      else
        set_owned_node_count(region, rank, (int)0);
    }
#endif

    transfer_edgeblocks(region, output_region, options.debug, options.verbose, rank);
    transfer_faceblocks(region, output_region, options.debug, options.verbose, rank);
    transfer_elementblocks(region, output_region, options.debug, options.verbose, rank);
    transfer_structuredblocks(region, output_region, options.debug, options.verbose, rank);

    transfer_nodesets(region, output_region, options.debug, options.verbose, rank);
    transfer_edgesets(region, output_region, options.debug, options.verbose, rank);
    transfer_facesets(region, output_region, options.debug, options.verbose, rank);
    transfer_elemsets(region, output_region, options.debug, options.verbose, rank);

    transfer_sidesets(region, output_region, options.debug, options.verbose, rank);
    transfer_commsets(region, output_region, options.debug, options.verbose, rank);

    transfer_coordinate_frames(region, output_region, options.debug, options.verbose, rank);

    if (options.debug && rank == 0) {
      std::cerr << "END STATE_DEFINE_MODEL... " << '\n';
    }
    if (memory_stats) {
      dbi->util().progress("END STATE_DEFINE_MODEL");
    }

    output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

    if (options.verbose && rank == 0) {
      std::cerr << "Maximum Field size = " << max_field_size << " bytes.\n";
    }
    data_pool.data.resize(max_field_size);
    if (options.verbose && rank == 0) {
      std::cerr << "Resize finished...\n";
    }
    if (options.debug && rank == 0) {
      std::cerr << "TRANSFERRING MESH FIELD DATA ... " << '\n';
    }
    if (memory_stats) {
      dbi->util().progress("TRANSFERRING MESH FIELD DATA ... ");
    }

    // Model defined, now fill in the model data...
    output_region.begin_mode(Ioss::STATE_MODEL);

    // Transfer MESH field_data from input to output...
    // Transfer MESH field_data from input to output...
    bool node_major = output_region.node_major();

    if (!node_major) {
      transfer_field_data(region.get_element_blocks(), output_region, data_pool, Ioss::Field::MESH,
                          options);
      transfer_field_data(region.get_element_blocks(), output_region, data_pool,
                          Ioss::Field::ATTRIBUTE, options);
    }

    if (region.mesh_type() != Ioss::MeshType::STRUCTURED) {
      transfer_field_data(region.get_node_blocks(), output_region, data_pool, Ioss::Field::MESH,
                          options);
      transfer_field_data(region.get_node_blocks(), output_region, data_pool,
                          Ioss::Field::ATTRIBUTE, options);
    }

    if (node_major) {
      transfer_field_data(region.get_element_blocks(), output_region, data_pool, Ioss::Field::MESH,
                          options);
      transfer_field_data(region.get_element_blocks(), output_region, data_pool,
                          Ioss::Field::ATTRIBUTE, options);
    }

    transfer_field_data(region.get_structured_blocks(), output_region, data_pool, Ioss::Field::MESH,
                        options);
    transfer_field_data(region.get_structured_blocks(), output_region, data_pool,
                        Ioss::Field::ATTRIBUTE, options);

    transfer_field_data(region.get_edge_blocks(), output_region, data_pool, Ioss::Field::MESH,
                        options);
    transfer_field_data(region.get_edge_blocks(), output_region, data_pool, Ioss::Field::ATTRIBUTE,
                        options);

    transfer_field_data(region.get_face_blocks(), output_region, data_pool, Ioss::Field::MESH,
                        options);
    transfer_field_data(region.get_face_blocks(), output_region, data_pool, Ioss::Field::ATTRIBUTE,
                        options);

    transfer_field_data(region.get_nodesets(), output_region, data_pool, Ioss::Field::MESH,
                        options);
    transfer_field_data(region.get_nodesets(), output_region, data_pool, Ioss::Field::ATTRIBUTE,
                        options);

    transfer_field_data(region.get_edgesets(), output_region, data_pool, Ioss::Field::MESH,
                        options);
    transfer_field_data(region.get_edgesets(), output_region, data_pool, Ioss::Field::ATTRIBUTE,
                        options);

    transfer_field_data(region.get_facesets(), output_region, data_pool, Ioss::Field::MESH,
                        options);
    transfer_field_data(region.get_facesets(), output_region, data_pool, Ioss::Field::ATTRIBUTE,
                        options);

    transfer_field_data(region.get_elementsets(), output_region, data_pool, Ioss::Field::MESH,
                        options);
    transfer_field_data(region.get_elementsets(), output_region, data_pool, Ioss::Field::ATTRIBUTE,
                        options);

    transfer_field_data(region.get_commsets(), output_region, data_pool, Ioss::Field::MESH,
                        options);
    transfer_field_data(region.get_commsets(), output_region, data_pool, Ioss::Field::ATTRIBUTE,
                        options);
    transfer_field_data(region.get_commsets(), output_region, data_pool, Ioss::Field::COMMUNICATION,
                        options);

    // Side Sets
    if (region.mesh_type() == Ioss::MeshType::UNSTRUCTURED) {
      const auto &fss = region.get_sidesets();
      for (const auto &ifs : fss) {
        std::string name = ifs->name();
        if (options.debug && rank == 0) {
          std::cerr << name << ", ";
        }
        // Find matching output sideset
        Ioss::SideSet *ofs = output_region.get_sideset(name);

        if (ofs != nullptr) {
          transfer_field_data(ifs, ofs, data_pool, Ioss::Field::MESH, options);
          transfer_field_data(ifs, ofs, data_pool, Ioss::Field::ATTRIBUTE, options);

          const auto &fbs = ifs->get_side_blocks();
          for (const auto &ifb : fbs) {

            // Find matching output sideblock
            std::string fbname = ifb->name();
            if (options.debug && rank == 0) {
              std::cerr << fbname << ", ";
            }
            Ioss::SideBlock *ofb = ofs->get_side_block(fbname);

            if (ofb != nullptr) {
              transfer_field_data(ifb, ofb, data_pool, Ioss::Field::MESH, options);
              transfer_field_data(ifb, ofb, data_pool, Ioss::Field::ATTRIBUTE, options);
            }
          }
        }
      }
      if (options.debug && rank == 0) {
        std::cerr << '\n';
      }
    }
    if (options.debug && rank == 0) {
      std::cerr << "END STATE_MODEL... " << '\n';
    }
    if (memory_stats) {
      dbi->util().progress("END STATE_MODEL... ");
    }
    output_region.end_mode(Ioss::STATE_MODEL);

    if (options.delete_timesteps) {
      Ioss::Utils::clear(data_pool.data);
      return;
    }
  } // appending

  if (options.debug && rank == 0) {
    std::cerr << "DEFINING TRANSIENT FIELDS ... " << '\n';
  }
  if (memory_stats) {
    dbi->util().progress("DEFINING TRANSIENT FIELDS ... ");
  }

  if (region.property_exists("state_count") && region.get_property("state_count").get_int() > 0) {
    if (options.verbose && rank == 0) {
      std::cerr << "\n Number of time steps on database     =" << std::setw(12)
                << region.get_property("state_count").get_int() << "\n\n";
    }

    output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

    // For each 'TRANSIENT' field in the node blocks and element
    // blocks, transfer to the output node and element blocks.
    transfer_fields(&region, &output_region, Ioss::Field::TRANSIENT);

    transfer_fields(region.get_node_blocks(), output_region, Ioss::Field::TRANSIENT, options, rank);
    transfer_fields(region.get_edge_blocks(), output_region, Ioss::Field::TRANSIENT, options, rank);
    transfer_fields(region.get_face_blocks(), output_region, Ioss::Field::TRANSIENT, options, rank);
    transfer_fields(region.get_element_blocks(), output_region, Ioss::Field::TRANSIENT, options,
                    rank);
    transfer_fields(region.get_structured_blocks(), output_region, Ioss::Field::TRANSIENT, options,
                    rank);

    transfer_fields(region.get_nodesets(), output_region, Ioss::Field::TRANSIENT, options, rank);
    transfer_fields(region.get_edgesets(), output_region, Ioss::Field::TRANSIENT, options, rank);
    transfer_fields(region.get_facesets(), output_region, Ioss::Field::TRANSIENT, options, rank);
    transfer_fields(region.get_elementsets(), output_region, Ioss::Field::TRANSIENT, options, rank);

    // Side Sets
    {
      const auto &fss = region.get_sidesets();
      for (const auto &ifs : fss) {
        std::string name = ifs->name();
        if (options.debug && rank == 0) {
          std::cerr << name << ", ";
        }

        // Find matching output sideset
        Ioss::SideSet *ofs = output_region.get_sideset(name);
        if (ofs != nullptr) {
          transfer_fields(ifs, ofs, Ioss::Field::TRANSIENT);

          const auto &fbs = ifs->get_side_blocks();
          for (const auto &ifb : fbs) {

            // Find matching output sideblock
            std::string fbname = ifb->name();
            if (options.debug && rank == 0) {
              std::cerr << fbname << ", ";
            }

            Ioss::SideBlock *ofb = ofs->get_side_block(fbname);
            if (ofb != nullptr) {
              transfer_fields(ifb, ofb, Ioss::Field::TRANSIENT);
            }
          }
        }
      }
      if (options.debug && rank == 0) {
        std::cerr << '\n';
      }
    }
    if (options.debug && rank == 0) {
      std::cerr << "END STATE_DEFINE_TRANSIENT... " << '\n';
    }
    if (memory_stats) {
      dbi->util().progress("END STATE_DEFINE_TRANSIENT... ");
    }
    output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
  }

  if (options.debug && rank == 0) {
    std::cerr << "TRANSFERRING TRANSIENT FIELDS ... " << '\n';
  }
  if (memory_stats) {
    dbi->util().progress("TRANSFERRING TRANSIENT FIELDS... ");
  }

  output_region.begin_mode(Ioss::STATE_TRANSIENT);
  // Get the timesteps from the input database.  Step through them
  // and transfer fields to output database...

  int step_count = region.get_property("state_count").get_int();

  for (int istep = 1; istep <= step_count; istep++) {
    double time = region.get_state_time(istep);
    if (time < options.minimum_time) {
      continue;
    }
    if (options.maximum_time != 0.0 && time > options.maximum_time) {
      break;
    }

    int ostep = output_region.add_state(time);
    show_step(istep, time, options.verbose, rank);

    output_region.begin_state(ostep);
    region.begin_state(istep);

    transfer_field_data(&region, &output_region, data_pool, Ioss::Field::TRANSIENT, options);

    if (region.mesh_type() != Ioss::MeshType::STRUCTURED) {
      transfer_field_data(region.get_node_blocks(), output_region, data_pool,
                          Ioss::Field::TRANSIENT, options);
    }
    transfer_field_data(region.get_edge_blocks(), output_region, data_pool, Ioss::Field::TRANSIENT,
                        options);
    transfer_field_data(region.get_face_blocks(), output_region, data_pool, Ioss::Field::TRANSIENT,
                        options);
    transfer_field_data(region.get_element_blocks(), output_region, data_pool,
                        Ioss::Field::TRANSIENT, options);
    transfer_field_data(region.get_structured_blocks(), output_region, data_pool,
                        Ioss::Field::TRANSIENT, options);

    transfer_field_data(region.get_nodesets(), output_region, data_pool, Ioss::Field::TRANSIENT,
                        options);
    transfer_field_data(region.get_edgesets(), output_region, data_pool, Ioss::Field::TRANSIENT,
                        options);
    transfer_field_data(region.get_facesets(), output_region, data_pool, Ioss::Field::TRANSIENT,
                        options);
    transfer_field_data(region.get_elementsets(), output_region, data_pool, Ioss::Field::TRANSIENT,
                        options);

    // Side Sets
    {
      const auto &fss = region.get_sidesets();
      for (const auto &ifs : fss) {
        std::string name = ifs->name();
        if (options.debug && rank == 0) {
          std::cerr << name << ", ";
        }

        // Find matching output sideset
        Ioss::SideSet *ofs = output_region.get_sideset(name);
        if (ofs != nullptr) {
          transfer_field_data(ifs, ofs, data_pool, Ioss::Field::TRANSIENT, options);

          const auto &fbs = ifs->get_side_blocks();
          for (const auto &ifb : fbs) {

            // Find matching output sideblock
            std::string fbname = ifb->name();
            if (options.debug && rank == 0) {
              std::cerr << fbname << ", ";
            }

            Ioss::SideBlock *ofb = ofs->get_side_block(fbname);
            if (ofb != nullptr) {
              transfer_field_data(ifb, ofb, data_pool, Ioss::Field::TRANSIENT, options);
            }
          }
        }
      }
    }
    region.end_state(istep);
    output_region.end_state(ostep);
    if (options.delay > 0.0) {
      struct timespec delay;
      delay.tv_sec  = (int)options.delay;
      delay.tv_nsec = (options.delay - delay.tv_sec) * 1000000000L;
      nanosleep(&delay, nullptr);
    }
  }
  if (options.debug && rank == 0) {
    std::cerr << "END STATE_TRANSIENT... " << '\n';
  }
  if (memory_stats) {
    dbi->util().progress("END STATE_TRANSIENT ... ");
  }
  output_region.end_mode(Ioss::STATE_TRANSIENT);
  Ioss::Utils::clear(data_pool.data);
}

namespace {
  void transfer_nodeblock(Ioss::Region &region, Ioss::Region &output_region, DataPool &pool,
                          bool debug, bool verbose, int rank)
  {
    const auto &nbs = region.get_node_blocks();
    size_t      id  = 1;
    for (const auto &inb : nbs) {
      std::string name = inb->name();
      if (debug && rank == 0) {
        std::cerr << name << ", ";
      }
      size_t num_nodes = inb->entity_count();
      size_t degree    = inb->get_property("component_degree").get_int();
      if (verbose && rank == 0) {
        std::cerr << " Number of  Coordinates per Node        =" << std::setw(12) << degree << "\n";
      }
      if (verbose && rank == 0) {
        std::cerr << " Number of             Nodes            =" << std::setw(12) << num_nodes
                  << "\n";
      }
      auto nb = new Ioss::NodeBlock(output_region.get_database(), name, num_nodes, degree);
      output_region.add(nb);

      transfer_properties(inb, nb);

      if (output_region.get_database()->needs_shared_node_information()) {
        // If the "owning_processor" field exists on the input
        // nodeblock, transfer it and the "ids" field to the output
        // nodeblock at this time since it is used to determine
        // per-processor sizes of nodeblocks and nodesets.
        if (inb->field_exists("owning_processor")) {
          size_t isize = inb->get_field("ids").get_size();
          pool.data.resize(isize);
          inb->get_field_data("ids", pool.data.data(), isize);
          nb->put_field_data("ids", pool.data.data(), isize);

          isize = inb->get_field("owning_processor").get_size();
          pool.data.resize(isize);
          inb->get_field_data("owning_processor", pool.data.data(), isize);
          nb->put_field_data("owning_processor", pool.data.data(), isize);
        }
      }

      transfer_fields(inb, nb, Ioss::Field::MESH);
      transfer_fields(inb, nb, Ioss::Field::ATTRIBUTE);
      ++id;
    }
    if (debug && rank == 0) {
      std::cerr << '\n';
    }
  }

  template <typename T>
  void transfer_fields(const std::vector<T *> &entities, Ioss::Region &output_region,
                       Ioss::Field::RoleType role, const Ioss::MeshCopyOptions &options, int rank)
  {
    for (const auto &entity : entities) {
      std::string name = entity->name();
      if (options.debug && rank == 0) {
        std::cerr << name << ", ";
      }

      // Find the corresponding output node_block...
      Ioss::GroupingEntity *oeb = output_region.get_entity(name, entity->type());
      if (oeb != nullptr) {
        transfer_fields(entity, oeb, role);
      }
    }
    if (options.debug && rank == 0) {
      std::cerr << '\n';
    }
  }

  template <typename T>
  void transfer_field_data(const std::vector<T *> &entities, Ioss::Region &output_region,
                           DataPool &pool, Ioss::Field::RoleType role,
                           const Ioss::MeshCopyOptions &options)
  {
    for (const auto &entity : entities) {
      std::string name = entity->name();

      // Find the corresponding output block...
      Ioss::GroupingEntity *output = output_region.get_entity(name, entity->type());
      if (output != nullptr) {
        transfer_field_data(entity, output, pool, role, options);
      }
    }
  }

  template <typename T>
  void transfer_blocks(const std::vector<T *> &blocks, Ioss::Region &output_region, bool debug,
                       bool verbose, int rank)
  {
    if (!blocks.empty()) {
      size_t total_entities = 0;
      for (const auto &iblock : blocks) {
        std::string name = iblock->name();
        if (debug && rank == 0) {
          std::cerr << name << ", ";
        }
        std::string type  = iblock->get_property("topology_type").get_string();
        size_t      count = iblock->entity_count();
        total_entities += count;

        auto block = new T(output_region.get_database(), name, type, count);
        if (iblock->property_exists("original_block_order")) {
          block->property_add(iblock->get_property("original_block_order"));
        }
        output_region.add(block);
        transfer_properties(iblock, block);
        transfer_fields(iblock, block, Ioss::Field::MESH);
        transfer_fields(iblock, block, Ioss::Field::ATTRIBUTE);
      }
      if (verbose && rank == 0) {
        std::cerr << " Number of " << std::setw(16) << (*blocks.begin())->type_string()
                  << "s            =" << std::setw(12) << blocks.size() << "\t"
                  << "Length of entity list   =" << std::setw(12) << total_entities << "\n";
      }
      if (debug && rank == 0) {
        std::cerr << '\n';
      }
    }
  }

  void transfer_structuredblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                                 bool verbose, int rank)
  {
    auto blocks = region.get_structured_blocks();
    if (!blocks.empty()) {
      size_t total_entities = 0;
      for (const auto &iblock : blocks) {
        std::string name = iblock->name();
        if (debug && rank == 0) {
          std::cerr << name << ", ";
        }
        size_t count = iblock->entity_count();
        total_entities += count;

        auto block = iblock->clone(output_region.get_database());
        output_region.add(block);
        transfer_properties(iblock, block);
        transfer_fields(iblock, block, Ioss::Field::MESH);
        transfer_fields(iblock, block, Ioss::Field::ATTRIBUTE);
      }
      if (verbose && rank == 0) {
        std::cerr << " Number of " << std::setw(16) << (*blocks.begin())->type_string()
                  << "s            =" << std::setw(12) << blocks.size() << "\t"
                  << "Length of entity list   =" << std::setw(12) << total_entities << "\n";
      }
      if (debug && rank == 0) {
        std::cerr << '\n';
      }
    }
  }

  void transfer_elementblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                              bool verbose, int rank)
  {
    const auto &ebs = region.get_element_blocks();
    transfer_blocks(ebs, output_region, debug, verbose, rank);
  }

  void transfer_edgeblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                           bool verbose, int rank)
  {
    const auto &ebs = region.get_edge_blocks();
    transfer_blocks(ebs, output_region, debug, verbose, rank);
  }

  void transfer_faceblocks(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                           bool verbose, int rank)
  {
    const auto &ebs = region.get_face_blocks();
    transfer_blocks(ebs, output_region, debug, verbose, rank);
  }

  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                         bool verbose, int rank)
  {
    const auto &fss         = region.get_sidesets();
    size_t      total_sides = 0;
    for (const auto &ss : fss) {
      std::string name = ss->name();
      if (debug && rank == 0) {
        std::cerr << name << ", ";
      }

      auto        surf = new Ioss::SideSet(output_region.get_database(), name);
      const auto &fbs  = ss->get_side_blocks();
      for (const auto &fb : fbs) {
        std::string fbname = fb->name();
        if (debug && rank == 0) {
          std::cerr << fbname << ", ";
        }
        std::string fbtype   = fb->get_property("topology_type").get_string();
        std::string partype  = fb->get_property("parent_topology_type").get_string();
        size_t      num_side = fb->entity_count();
        total_sides += num_side;

        auto block =
            new Ioss::SideBlock(output_region.get_database(), fbname, fbtype, partype, num_side);
        surf->add(block);
        transfer_properties(fb, block);
        transfer_fields(fb, block, Ioss::Field::MESH);
        transfer_fields(fb, block, Ioss::Field::ATTRIBUTE);
        if (fb->parent_block() != nullptr) {
          auto               fb_name = fb->parent_block()->name();
          Ioss::EntityBlock *parent =
              dynamic_cast<Ioss::EntityBlock *>(output_region.get_entity(fb_name));
          block->set_parent_block(parent);
        }
      }
      transfer_properties(ss, surf);
      transfer_fields(ss, surf, Ioss::Field::MESH);
      transfer_fields(ss, surf, Ioss::Field::ATTRIBUTE);
      output_region.add(surf);
    }
    if (verbose && rank == 0) {
      std::cerr << " Number of          SideSets            =" << std::setw(12) << fss.size()
                << "\t"
                << "Number of element sides =" << std::setw(12) << total_sides << "\n";
    }
    if (debug && rank == 0) {
      std::cerr << '\n';
    }
  }

  template <typename T>
  void transfer_sets(const std::vector<T *> &sets, Ioss::Region &output_region, bool debug,
                     bool verbose, int rank)
  {
    if (!sets.empty()) {
      size_t total_entities = 0;
      for (const auto &set : sets) {
        std::string name = set->name();
        if (debug && rank == 0) {
          std::cerr << name << ", ";
        }
        size_t count = set->entity_count();
        total_entities += count;
        auto o_set = new T(output_region.get_database(), name, count);
        output_region.add(o_set);
        transfer_properties(set, o_set);
        transfer_fields(set, o_set, Ioss::Field::MESH);
        transfer_fields(set, o_set, Ioss::Field::ATTRIBUTE);
      }

      if (verbose && rank == 0) {
        std::cerr << " Number of " << std::setw(16) << (*sets.begin())->type_string()
                  << "s            =" << std::setw(12) << sets.size() << "\t"
                  << "Length of entity list   =" << std::setw(12) << total_entities << "\n";
      }
      if (debug && rank == 0) {
        std::cerr << '\n';
      }
    }
  }

  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                         bool verbose, int rank)
  {
    const auto &nss = region.get_nodesets();
    transfer_sets(nss, output_region, debug, verbose, rank);
  }

  void transfer_edgesets(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                         bool verbose, int rank)
  {
    const auto &nss = region.get_edgesets();
    transfer_sets(nss, output_region, debug, verbose, rank);
  }

  void transfer_facesets(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                         bool verbose, int rank)
  {
    const auto &nss = region.get_facesets();
    transfer_sets(nss, output_region, debug, verbose, rank);
  }

  void transfer_elemsets(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                         bool verbose, int rank)
  {
    const auto &nss = region.get_elementsets();
    transfer_sets(nss, output_region, debug, verbose, rank);
  }

  void transfer_commsets(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                         bool verbose, int rank)
  {
    const auto &css = region.get_commsets();
    for (const auto &ics : css) {
      std::string name = ics->name();
      if (debug && rank == 0) {
        std::cerr << name << ", ";
      }
      std::string type  = ics->get_property("entity_type").get_string();
      size_t      count = ics->entity_count();
      auto        cs    = new Ioss::CommSet(output_region.get_database(), name, type, count);
      output_region.add(cs);
      transfer_properties(ics, cs);
      transfer_fields(ics, cs, Ioss::Field::MESH);
      transfer_fields(ics, cs, Ioss::Field::ATTRIBUTE);
      transfer_fields(ics, cs, Ioss::Field::COMMUNICATION);
    }
    if (debug && rank == 0) {
      std::cerr << '\n';
    }
  }

  void transfer_coordinate_frames(Ioss::Region &region, Ioss::Region &output_region, bool debug,
                                  bool verbose, int rank)
  {
    Ioss::CoordinateFrameContainer cf = region.get_coordinate_frames();
    for (const auto &frame : cf) {
      output_region.add(frame);
    }
  }

  void transfer_fields(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                       Ioss::Field::RoleType role, const std::string &prefix)
  {
    // Check for transient fields...
    Ioss::NameList fields;
    ige->field_describe(role, &fields);

    // Iterate through results fields and transfer to output
    // database...  If a prefix is specified, only transfer fields
    // whose names begin with the prefix
    for (const auto &field_name : fields) {
      Ioss::Field field = ige->get_field(field_name);
      max_field_size    = MAX(max_field_size, field.get_size());
      if (field_name != "ids" && !oge->field_exists(field_name) &&
          (prefix.length() == 0 ||
           std::strncmp(prefix.c_str(), field_name.c_str(), prefix.length()) == 0)) {
        // If the field does not already exist, add it to the output node block
        oge->field_add(field);
      }
    }
  }

  void transfer_field_data(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge, DataPool &pool,
                           Ioss::Field::RoleType role, const Ioss::MeshCopyOptions &options,
                           const std::string &prefix)
  {
    // Iterate through the TRANSIENT-role fields of the input
    // database and transfer to output database.
    Ioss::NameList state_fields;
    ige->field_describe(role, &state_fields);

    // Complication here is that if the 'role' is 'Ioss::Field::MESH',
    // then the 'ids' field must be transferred first...
    if (role == Ioss::Field::MESH && ige->field_exists("ids")) {
      assert(oge->field_exists("ids"));
      transfer_field_data_internal(ige, oge, pool, "ids", options);
    }

    for (const auto &field_name : state_fields) {
      // All of the 'Ioss::EntityBlock' derived classes have a
      // 'connectivity' field, but it is only interesting on the
      // Ioss::ElementBlock class. On the other classes, it just
      // generates overhead...
      if (field_name == "connectivity" && ige->type() != Ioss::ELEMENTBLOCK) {
        continue;
      }
      if (field_name == "ids") {
        continue;
      }
      if ((prefix.length() == 0 ||
           std::strncmp(prefix.c_str(), field_name.c_str(), prefix.length()) == 0)) {
        assert(oge->field_exists(field_name));
        transfer_field_data_internal(ige, oge, pool, field_name, options);
      }
    }
  }

  void transfer_field_data_internal(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                                    DataPool &pool, const std::string &field_name,
                                    const Ioss::MeshCopyOptions &options)
  {

    size_t isize = ige->get_field(field_name).get_size();
    assert(isize == oge->get_field(field_name).get_size());

    int basic_type = ige->get_field(field_name).get_type();

    if (field_name == "mesh_model_coordinates_x") {
      return;
    }
    if (field_name == "mesh_model_coordinates_y") {
      return;
    }
    if (field_name == "mesh_model_coordinates_z") {
      return;
    }
    if (field_name == "connectivity_raw") {
      return;
    }
    if (field_name == "element_side_raw") {
      return;
    }
    if (field_name == "ids_raw") {
      return;
    }
    if (field_name == "implicit_ids") {
      return;
    }
    if (field_name == "node_connectivity_status") {
      return;
    }
    if (field_name == "owning_processor") {
      return;
    }
    if (field_name == "entity_processor_raw") {
      return;
    }
    if (field_name == "ids" && ige->type() == Ioss::SIDEBLOCK) {
      return;
    }
    if (field_name == "ids" && ige->type() == Ioss::STRUCTUREDBLOCK) {
      return;
    }
    if (field_name == "cell_ids" && ige->type() == Ioss::STRUCTUREDBLOCK) {
      return;
    }
    if (field_name == "cell_node_ids" && ige->type() == Ioss::STRUCTUREDBLOCK) {
      return;
    }

    if (options.data_storage_type == 1 || options.data_storage_type == 2) {
      if (pool.data.size() < isize) {
        pool.data.resize(isize);
      }
    }
    else {
    }

    assert(pool.data.size() >= isize);

    switch (options.data_storage_type) {
    case 1: ige->get_field_data(field_name, pool.data.data(), isize); break;
    case 2:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        ige->get_field_data(field_name, pool.data);
      }
      else if (basic_type == Ioss::Field::INT32) {
        ige->get_field_data(field_name, pool.data_int);
      }
      else if (basic_type == Ioss::Field::INT64) {
        ige->get_field_data(field_name, pool.data_int64);
      }
      else if (basic_type == Ioss::Field::REAL) {
        ige->get_field_data(field_name, pool.data_double);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        ige->get_field_data(field_name, pool.data_complex);
      }
      else {
      }
      break;
#ifdef SEACAS_HAVE_KOKKOS
    case 3:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        ige->get_field_data<char>(field_name, pool.data_view_char);
      }
      else if (basic_type == Ioss::Field::INT32) {
        ige->get_field_data<int>(field_name, pool.data_view_int);
      }
      else if (basic_type == Ioss::Field::INT64) {
        ige->get_field_data<int64_t>(field_name, pool.data_view_int64);
      }
      else if (basic_type == Ioss::Field::REAL) {
        ige->get_field_data<double>(field_name, pool.data_view_double);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        // Since data_view_complex cannot be a global variable.
        ige->get_field_data(field_name, pool.data.data(), isize);
      }
      else {
      }
      break;
    case 4:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        ige->get_field_data<char>(field_name, pool.data_view_2D_char);
      }
      else if (basic_type == Ioss::Field::INT32) {
        ige->get_field_data<int>(field_name, pool.data_view_2D_int);
      }
      else if (basic_type == Ioss::Field::INT64) {
        ige->get_field_data<int64_t>(field_name, pool.data_view_2D_int64);
      }
      else if (basic_type == Ioss::Field::REAL) {
        ige->get_field_data<double>(field_name, pool.data_view_2D_double);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        // Since data_view_complex cannot be a global variable.
        ige->get_field_data(field_name, pool.data.data(), isize);
      }
      else {
      }
      break;
    case 5:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        ige->get_field_data<char, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_char_layout_space);
      }
      else if (basic_type == Ioss::Field::INT32) {
        ige->get_field_data<int, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_int_layout_space);
      }
      else if (basic_type == Ioss::Field::INT64) {
        ige->get_field_data<int64_t, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_int64_layout_space);
      }
      else if (basic_type == Ioss::Field::REAL) {
        ige->get_field_data<double, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_double_layout_space);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        // Since data_view_complex cannot be a global variable.
        ige->get_field_data(field_name, pool.data.data(), isize);
      }
      else {
      }
      break;
#endif
    default:
      if (field_name == "mesh_model_coordinates") {
        std::cerr << "data_storage option not recognized.";
      }
      return;
    }

    switch (options.data_storage_type) {
    case 1: oge->put_field_data(field_name, pool.data.data(), isize); break;
    case 2:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        oge->put_field_data(field_name, pool.data);
      }
      else if (basic_type == Ioss::Field::INT32) {
        oge->put_field_data(field_name, pool.data_int);
      }
      else if (basic_type == Ioss::Field::INT64) {
        oge->put_field_data(field_name, pool.data_int64);
      }
      else if (basic_type == Ioss::Field::REAL) {
        oge->put_field_data(field_name, pool.data_double);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        oge->put_field_data(field_name, pool.data_complex);
      }
      else {
      }
      break;
#ifdef SEACAS_HAVE_KOKKOS
    case 3:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        oge->put_field_data<char>(field_name, pool.data_view_char);
      }
      else if (basic_type == Ioss::Field::INT32) {
        oge->put_field_data<int>(field_name, pool.data_view_int);
      }
      else if (basic_type == Ioss::Field::INT64) {
        oge->put_field_data<int64_t>(field_name, pool.data_view_int64);
      }
      else if (basic_type == Ioss::Field::REAL) {
        oge->put_field_data<double>(field_name, pool.data_view_double);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        // Since data_view_complex cannot be a global variable.
        oge->put_field_data(field_name, pool.data.data(), isize);
      }
      else {
      }
      break;
    case 4:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        oge->put_field_data<char>(field_name, pool.data_view_2D_char);
      }
      else if (basic_type == Ioss::Field::INT32) {
        oge->put_field_data<int>(field_name, pool.data_view_2D_int);
      }
      else if (basic_type == Ioss::Field::INT64) {
        oge->put_field_data<int64_t>(field_name, pool.data_view_2D_int64);
      }
      else if (basic_type == Ioss::Field::REAL) {
        oge->put_field_data<double>(field_name, pool.data_view_2D_double);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        // Since data_view_complex cannot be a global variable.
        oge->put_field_data(field_name, pool.data.data(), isize);
      }
      else {
      }
      break;
    case 5:
      if ((basic_type == Ioss::Field::CHARACTER) || (basic_type == Ioss::Field::STRING)) {
        oge->put_field_data<char, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_char_layout_space);
      }
      else if (basic_type == Ioss::Field::INT32) {
        oge->put_field_data<int, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_int_layout_space);
      }
      else if (basic_type == Ioss::Field::INT64) {
        oge->put_field_data<int64_t, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_int64_layout_space);
      }
      else if (basic_type == Ioss::Field::REAL) {
        oge->put_field_data<double, Kokkos::LayoutRight, Kokkos::HostSpace>(
            field_name, pool.data_view_2D_double_layout_space);
      }
      else if (basic_type == Ioss::Field::COMPLEX) {
        // Since data_view_complex cannot be a global variable.
        oge->put_field_data(field_name, pool.data.data(), isize);
      }
      else {
      }
      break;
#endif
    default: return;
    }
    return;
  }

  void transfer_qa_info(Ioss::Region &in, Ioss::Region &out)
  {
    out.add_information_records(in.get_information_records());

    const std::vector<std::string> &qa = in.get_qa_records();
    for (size_t i = 0; i < qa.size(); i += 4) {
      out.add_qa_record(qa[i + 0], qa[i + 1], qa[i + 2], qa[i + 3]);
    }
  }

  void transfer_properties(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge)
  {
    Ioss::NameList properties;
    ige->property_describe(&properties);

    // Iterate through properties and transfer to output database...
    for (const auto &property : properties) {
      if (!oge->property_exists(property)) {
        oge->property_add(ige->get_property(property));
      }
    }
  }

  void show_step(int istep, double time, bool verbose, int rank)
  {
    if (verbose && rank == 0) {
      std::cerr.setf(std::ios::scientific);
    }
    if (verbose && rank == 0) {
      std::cerr.setf(std::ios::showpoint);
    }
    if (verbose && rank == 0) {
      std::cerr << "     Time step " << std::setw(5) << istep << " at time " << std::setprecision(5)
                << time << '\n';
    }
  }

  template <typename INT>
  void set_owned_node_count(Ioss::Region &region, int my_processor, INT /*dummy*/)
  {
    Ioss::NodeBlock *nb = region.get_node_block("nodeblock_1");
    if (nb->field_exists("owning_processor")) {
      std::vector<int> my_data;
      nb->get_field_data("owning_processor", my_data);

      INT owned = std::count(my_data.begin(), my_data.end(), my_processor);
      nb->property_add(Ioss::Property("locally_owned_count", owned));

      // Set locally_owned_count property on all nodesets...
      Ioss::NodeSetContainer nss = region.get_nodesets();
      for (auto ns : nss) {

        std::vector<INT> ids;
        ns->get_field_data("ids_raw", ids);
        owned = 0;
        for (size_t n = 0; n < ids.size(); n++) {
          INT id = ids[n];
          if (my_data[id - 1] == my_processor) {
            owned++;
          }
        }
        ns->property_add(Ioss::Property("locally_owned_count", owned));
      }
    }
  }
} // namespace
