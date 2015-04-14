#include <exodus/Ioex_Utils.h>
#include <Ioss_Utils.h>
#include <Ioss_Region.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_VariableType.h>
#include <tokenize.h>
#include <cstring>
#include <algorithm>


namespace {
  const size_t max_line_length   = MAX_LINE_LENGTH;

  const std::string SEP() {return std::string("@");} // Separator for attribute offset storage
  const std::string SCALAR()     {return std::string("scalar");}
  const std::string VECTOR3D()   {return std::string("vector_3d");}
  const std::string SYM_TENSOR() {return std::string("sym_tensor_33");}

  template <typename INT>
  void internal_write_coordinate_frames(int exoid, const Ioss::CoordinateFrameContainer &frames, INT /*dummy*/)
  {
    // Query number of coordinate frames...
    int nframes = (int)frames.size();
    if (nframes > 0) {
      std::vector<char> tags(nframes);
      std::vector<double> coordinates(nframes*9);
      std::vector<INT>  ids(nframes);

      for (size_t i=0; i < frames.size(); i++) {
	ids[i]  = frames[i].id();
	tags[i] = frames[i].tag();
	const double *coord = frames[i].coordinates();
	for (size_t j=0; j < 9; j++) {
	  coordinates[9*i+j] = coord[j];
	}
      }
      ex_put_coordinate_frames(exoid, nframes, TOPTR(ids), TOPTR(coordinates), TOPTR(tags));
    }
  }

  void update_last_time_attribute(int exodusFilePtr, double value)
  {
    char errmsg[MAX_ERR_LENGTH];
    const char *routine = "Ioex::Utils::update_last_time_attribute()";
    
    double tmp = 0.0;
    int rootid = (unsigned)exodusFilePtr & EX_FILE_ID_MASK;
    int status = nc_get_att_double(rootid, NC_GLOBAL, "last_written_time", &tmp);
    if (status == NC_NOERR && value > tmp) {
      status=nc_put_att_double(rootid, NC_GLOBAL, "last_written_time",
			       NC_DOUBLE, 1, &value);
      if (status != NC_NOERR) {
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to define 'last_written_time' attribute to file id %d",
		exodusFilePtr);
	ex_err(routine,errmsg,status);
      }
    }
  }

  bool read_last_time_attribute(int exodusFilePtr, double *value)
  {
    // Check whether the "last_written_time" attribute exists.  If it does,
    // return the value of the attribute in 'value' and return 'true'.
    // If not, don't change 'value' and return 'false'.
    bool found = false;

    int rootid = (unsigned)exodusFilePtr & EX_FILE_ID_MASK;
    nc_type att_type = NC_NAT;
    size_t att_len = 0;
    int status = nc_inq_att(rootid, NC_GLOBAL, "last_written_time", &att_type, &att_len);
    if (status == NC_NOERR && att_type == NC_DOUBLE) {
      // Attribute exists on this database, read it...
      double tmp = 0.0;
      status = nc_get_att_double(rootid, NC_GLOBAL, "last_written_time", &tmp);
      if (status == NC_NOERR) {
	*value = tmp;
	found = true;
      } else {
	char errmsg[MAX_ERR_LENGTH];
	const char *routine = "Ioex::Utils::read_last_time_attribute()";
	ex_opts(EX_VERBOSE);
	sprintf(errmsg,
		"Error: failed to read last_written_time attribute from file id %d", exodusFilePtr);
	ex_err(routine,errmsg,status);
	found = false;
      }
    }
    return found;
  }

  template <typename INT>
  void internal_add_coordinate_frames(int exoid, Ioss::Region *region, INT /*dummy*/)
  {
    // Query number of coordinate frames...
    int nframes = 0;
    ex_get_coordinate_frames(exoid, &nframes, NULL, NULL, NULL);

    if (nframes > 0) {
      std::vector<char> tags(nframes);
      std::vector<double> coord(nframes*9);
      std::vector<INT>  ids(nframes);
      ex_get_coordinate_frames(exoid, &nframes, TOPTR(ids), TOPTR(coord), TOPTR(tags));

      for (int i=0; i<nframes; i++) {
	Ioss::CoordinateFrame cf(ids[i], tags[i], &coord[9*i]);
	region->add(cf);
      }
    }
  }

  size_t match(const char *name1, const char *name2)
  {
    size_t l1 = std::strlen(name1);
    size_t l2 = std::strlen(name2);
    size_t len = l1 < l2 ? l1 : l2;
    for (size_t i=0; i < len; i++) {
      if (name1[i] != name2[i]) {
	while (i > 0 && isdigit(name1[i-1]) && isdigit(name2[i-1])) {
	  i--;
	  // Back up to first non-digit so to handle "evar0000, evar0001, ..., evar 1123"
	}
	return i;
      }
    }
    return len;
  }

  size_t get_number(const std::string &suffix)
  {
    int N = 0;
    bool all_dig = suffix.find_first_not_of("0123456789") == std::string::npos;
    if (all_dig) {
      N = std::strtol(suffix.c_str(), NULL, 10);
    }
    return N;
  }

}

namespace Ioex {
  const char *Version() {return "Ioex_DatabaseIO.C 2015/04/13";}

  bool type_match(const std::string& type, const char *substring)
  {
    // Returns true if 'substring' is a sub-string of 'type'.
    // The comparisons are case-insensitive
    // 'substring' is required to be in all lowercase.
    const char *s = substring;
    const char *t = type.c_str();

    assert(s != NULL && t != NULL);
    while (*s != '\0' && *t != '\0') {
      if (*s++ != tolower(*t++)) {
	return false;
      }
    }
    return true;
  }

  void decode_surface_name(Ioex::SideSetMap &fs_map, Ioex::SideSetSet &fs_set, const std::string &name)
  {
    std::vector<std::string> tokens;
    Ioss::tokenize(name, "_", tokens);
    if (tokens.size() >= 4) {
      // Name of form: "name_eltopo_sidetopo_id" or
      // "name_block_id_sidetopo_id" "name" is typically "surface".
      // The sideset containing this should then be called "name_id"

      // Check whether the second-last token is a side topology and
      // the third-last token is an element topology.
      const Ioss::ElementTopology *side_topo = Ioss::ElementTopology::factory(tokens[tokens.size()-2], true);
      if (side_topo != NULL) {
	const Ioss::ElementTopology *element_topo = Ioss::ElementTopology::factory(tokens[tokens.size()-3], true);
	if (element_topo != NULL || tokens[tokens.size()-4] == "block") {
	  // The remainder of the tokens will be used to create
	  // a side set name and then this sideset will be
	  // a side block in that set.
	  std::string fs_name;
	  size_t last_token = tokens.size()-3;
	  if (element_topo == NULL)
	    last_token--;
	  for (size_t tok=0; tok < last_token; tok++) {
	    fs_name += tokens[tok];
	  }
	  fs_name += "_";
	  fs_name += tokens[tokens.size()-1]; // Add on the id.

	  fs_set.insert(fs_name);
	  fs_map.insert(Ioex::SideSetMap::value_type(name,fs_name));
	}
      }
    }
  }

  bool set_id(const Ioss::GroupingEntity *entity, ex_entity_type type, Ioex::EntityIdSet *idset)
  {
    // See description of 'get_id' function.  This function just primes
    // the idset with existing ids so that when we start generating ids,
    // we don't overwrite an existing one.

    // Avoid a few string constructors/destructors
    static std::string prop_name("name");
    static std::string id_prop("id");

    bool succeed = false;
    if (entity->property_exists(id_prop)) {
      int64_t id = entity->get_property(id_prop).get_int();

      // See whether it already exists...
      succeed = idset->insert(std::make_pair((int)type,id)).second;
      if (!succeed) {
	// Need to remove the property so it doesn't cause problems
	// later...
	Ioss::GroupingEntity *new_entity = const_cast<Ioss::GroupingEntity*>(entity);
	new_entity->property_erase(id_prop);
	assert(!entity->property_exists(id_prop));
      }
    }
    return succeed;
  }

  // Potentially extract the id from a name possibly of the form name_id.
  // If not of this form, return 0;
  int64_t extract_id(const std::string &name_id)
  {
    std::vector<std::string> tokens;
    Ioss::tokenize(name_id,"_",tokens);

    if (tokens.size() == 1)
      return 0;

    // Check whether last token is an integer...
    std::string str_id = tokens[tokens.size()-1];
    size_t len = str_id.length();
    bool is_int = true;
    for (size_t i=0; i < len; i++) {
      if (str_id[i] < '0' || str_id[i] > '9') {
	is_int = false;
	break;
      }
    }      
    if (is_int)
      return std::atoi(str_id.c_str());

    return 0;
  }

  int64_t get_id(const Ioss::GroupingEntity *entity, ex_entity_type type, Ioex::EntityIdSet *idset)
  {
    // Sierra uses names to refer to grouping entities; however,
    // exodusII requires integer ids.  When reading an exodusII file,
    // the DatabaseIO creates a name by concatenating the entity
    // type (e.g., 'block') and the id separated by an underscore.  For
    // example, an exodusII element block with an id of 100 would be
    // encoded into "block_100"

    // This routine tries to determine the id of the entity using 3
    // approaches:
    //
    // 1. If the entity contains a property named 'id', this is used.
    // The DatabaseIO actually stores the id in the "id" property;
    // however, other grouping entity creators are not required to do
    // this so the property is not guaranteed to exist.
    //
    // 2.If property does not exist, it tries to decode the entity name
    // based on the above encoding.  Again, it is not required that the
    // name follow this convention so success is not guaranteed.
    //
    // 3. If all other schemes fail, the routine picks an id for the entity
    // and returns it.  It also stores this id in the "id" property so an
    // entity will always return the same id for multiple calls.
    // Note that this violates the 'const'ness of the entity so we use
    // a const-cast.

    // Avoid a few string constructors/destructors
    static std::string prop_name("name");
    static std::string id_prop("id");

    int64_t id = 1;

    if (entity->property_exists(id_prop)) {
      id = entity->get_property(id_prop).get_int();
      return id;

    } else {
      // Try to decode an id from the name.
      std::string name_string = entity->get_property(prop_name).get_string();
      id = extract_id(name_string);
      if (id <= 0) id = 1;
    }

    // At this point, we either have an id equal to '1' or we have an id
    // extracted from the entities name. Increment it until it is
    // unique...
    while (idset->find(std::make_pair(int(type), id)) != idset->end()) {
      ++id;
    }

    // 'id' is a unique id for this entity type...
    idset->insert(std::make_pair((int)type,id));
    Ioss::GroupingEntity *new_entity = const_cast<Ioss::GroupingEntity*>(entity);
    new_entity->property_add(Ioss::Property(id_prop, id));
    return id;
  }

  bool find_displacement_field(Ioss::NameList &fields,
			       const Ioss::GroupingEntity *block,
			       int ndim,
			       std::string *disp_name)
  {
    // This is a kluge to work with many of the SEACAS codes.  The
    // convention used (in Blot and others) is that the first 'ndim'
    // nodal variables are assumed to be displacements *if* the first
    // character of the names is 'D' and the last characters match the
    // coordinate labels (typically 'X', 'Y', and 'Z').  This routine
    // looks for the field that has the longest match with the string
    // "displacement" and is of the correct storage type (VECTOR_2D or
    // VECTOR_3D).  If found, it returns the name.
    //

    static char displace[] = "displacement";

    Ioss::NameList::const_iterator IF;
    Ioss::NameList::const_iterator IFend = fields.end();
    size_t max_span = 0;

    for (IF = fields.begin(); IF != IFend; ++IF) {
      const char *name = (*IF).c_str();
      std::string lc_name(name);

      Ioss::Utils::fixup_name(lc_name);
      size_t span = match(lc_name.c_str(), displace);
      if (span > max_span) {
	const Ioss::VariableType *var_type =
	  block->get_field((*IF)).transformed_storage();
	int comp_count = var_type->component_count();
	if (comp_count == ndim) {
	  max_span  = span;
	  *disp_name = *IF;
	}
      }
    }
    return max_span > 0 ? true : false;
  }

  void fix_bad_name(char* name)
  {
    assert(name != NULL);

    size_t len = std::strlen(name);
    for (size_t i=0; i < len; i++) {
      if (name[i] < 32 || name[i] > 126) {
	// Zero out entire name if a bad character found anywhere in the name.
	for (size_t j=0; j < len; j++) {
	  name[j] = '\0';
	}
	return;
      }
    }
  }

  std::string get_entity_name(int exoid, ex_entity_type type, int64_t id,
                              const std::string &basename, int length,
                              bool &db_has_name)
  {
    std::vector<char> buffer(length+1);
    buffer[0] = '\0';
    int error = ex_get_name(exoid, type, id, TOPTR(buffer));
    if (error < 0)
      exodus_error(exoid, __LINE__, -1);
    if (buffer[0] != '\0') {
      Ioss::Utils::fixup_name(TOPTR(buffer));
      // Filter out names of the form "basename_id" if the name
      // id doesn't match the id in the name...
      size_t base_size = basename.size();
      if (std::strncmp(basename.c_str(), &buffer[0], base_size) == 0) {
        int64_t name_id = extract_id(TOPTR(buffer));
        if (name_id > 0 && name_id != id) {
          // See if name is truly of form "basename_name_id"
          std::string tmp_name =  Ioss::Utils::encode_entity_name(basename, name_id);
          if (tmp_name == TOPTR(buffer)) {
            std::string new_name =  Ioss::Utils::encode_entity_name(basename, id);
            IOSS_WARNING << "WARNING: The entity named '" << TOPTR(buffer) << "' has the id " << id
			 << " which does not match the embedded id " << name_id
			 << ".\n         This can cause issues later on; the entity will be renamed to '"
			 << new_name << "' (IOSS)\n\n";
            db_has_name = false;
            return new_name;
          }
        }
      }
      db_has_name = true;
      return (std::string(TOPTR(buffer)));
    } else {
      db_has_name = false;
      return Ioss::Utils::encode_entity_name(basename, id);
    }
  }

  char ** get_exodus_names(size_t count, int size)
  {
    char **names = new char* [count];
    for (size_t i=0; i < count; i++) {
      names[i] = new char [size+1];
      std::memset(names[i], '\0', size+1);
    }
    return names;
  }

  void delete_exodus_names(char **names, int count)
  {
    for (int i=0; i < count; i++) {delete [] names[i];}
    delete [] names;
  }

  void exodus_error(int exoid, int lineno, int /* processor */) {
    std::ostringstream errmsg;
    // Create errmsg here so that the exerrval doesn't get cleared by
    // the ex_close call.
    errmsg << "Exodus error (" << exerrval << ")" << nc_strerror(exerrval) << " at line " << lineno
	   << " in file '" << Version()
	   << "' Please report to gdsjaar@sandia.gov if you need help.";

    ex_err(NULL, NULL, EX_PRTLASTMSG);
    if (exoid > 0)
      ex_close(exoid);
    IOSS_ERROR(errmsg);
  }

  // common
  const Ioss::VariableType *match_composite_field(char** names, Ioss::IntVector &which_names,
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

    std::vector<std::string> tokens;
    Ioss::tokenize(names[which_names[which_names.size()-1]] ,suffix, tokens);

    if (tokens.size() <= 2)
      return NULL;

    assert(tokens.size() > 2);

    // Check that suffix is a number -- all digits
    size_t N = get_number(tokens[tokens.size()-1]);

    if (N == 0)
      return NULL;

    if (which_names.size() % N != 0) {
      return NULL;
    }

    size_t inner_token = tokens.size() - 2;
    size_t inner_comp = which_names.size() / N;

    // Gather the first 'inner_ccomp' inner field suffices...
    std::vector<Ioss::Suffix> suffices;
    for (size_t i=0; i < inner_comp; i++) {
      std::vector<std::string> ltokens;
      Ioss::tokenize(names[which_names[i]], suffix, ltokens);
      // The second-last token is the suffix for this component...
      Ioss::Suffix tmp(ltokens[inner_token]);
      suffices.push_back(tmp);
    }

    // check that the suffices on the next copies of the inner field
    // match the first copy...
    size_t j = inner_comp;
    for (size_t copy = 1; copy < N; copy++) {
      for (size_t i=0; i < inner_comp; i++) {
	std::vector<std::string> ltokens;
	Ioss::tokenize(names[which_names[j++]], suffix, ltokens);
	// The second-last token is the suffix for this component...
	if (suffices[i] != ltokens[inner_token]) {
	  return NULL;
	}
      }
    }

    // All 'N' copies of the inner field match, now see the
    // suffices actually defines a field...
    const Ioss::VariableType *type = Ioss::VariableType::factory(suffices);
    if (type != NULL) {
      type = Ioss::VariableType::factory(type->name(), N);
    }
    return type;
  }

  const Ioss::VariableType *match_single_field(char** names, Ioss::IntVector &which_names,
					       const char suffix_separator)
  {
    // Strip off the suffix from each name indexed in 'which_names'
    // and see if it defines a valid type...
    std::vector<Ioss::Suffix> suffices;

    char suffix[2];
    suffix[0] = suffix_separator;
    suffix[1] = 0;

    for (size_t i=0; i < which_names.size(); i++) {
      std::vector<std::string> tokens;
      Ioss::tokenize(names[which_names[i]], suffix, tokens);
      size_t num_tokens = tokens.size();
      // The last token is the suffix for this component...
      Ioss::Suffix tmp(tokens[num_tokens-1]);
      suffices.push_back(tmp);
    }
    const Ioss::VariableType *type = Ioss::VariableType::factory(suffices);
    return type;
  }

  Ioss::Field get_next_field(char** names, int num_names, size_t count,
			     Ioss::Field::RoleType fld_role,
			     const char suffix_separator, int *truth_table)
  {
    // NOTE: 'names' are all lowercase at this point.

    // Assumption 1: To convert to a non-SCALAR type, the variable name
    // must have an field_suffix_sep in the name separating the suffixes from
    // the main name.

    // Find first unused name (used names have '\0' as first character...
    int index = 0;
    bool found_valid = false;
    for (index = 0; index < num_names; index++) {
      assert(truth_table == NULL || truth_table[index] == 1 || truth_table[index] == 0);
      if ((truth_table == NULL || truth_table[index] == 1) && names[index][0] != '\0') {
	found_valid = true;
	break;
      }
    }

    if (!found_valid) {
      // Return an invalid field...
      return Ioss::Field("", Ioss::Field::INVALID, SCALAR(), fld_role, 1);
    }

    // At this point, name[index] should be a valid potential field
    // name and all names[i] with i < index are either already used or
    // not valid for this grouping entity (truth_table entry == 0).
    assert (index < num_names && names[index][0] != '\0' && (truth_table == NULL || truth_table[index] == 1));
    char *name = names[index];

    // Split the name up into tokens separated by the
    // 'suffix_separator'.  Note that the basename itself could
    // contain a suffix_separator (back_stress_xx or
    // back_stress_xx_01). Need to ignore embedded separators
    // (back_stress) and also recognize composite variable types
    // (back_stress_xx_01). At the current time, a composite variable
    // type can only contain two non-composite variable types, so we
    // only need to look to be concerned with the last 1 or 2 tokens...
    char suffix[2];
    suffix[0] = suffix_separator;
    suffix[1] = 0;
    std::vector<std::string> tokens;
    Ioss::tokenize(name, suffix, tokens);
    size_t num_tokens = tokens.size();

    // Check that separator is not first or last character of the name...
    bool invalid = tokens[0].empty() || tokens[num_tokens-1].empty();
    if (num_tokens == 1 || invalid) {
      // It is not a (Sierra-generated) name for a non-SCALAR variable
      // Return a SCALAR field
      Ioss::Field field(name, Ioss::Field::REAL, SCALAR(), fld_role, count);
      names[index][0] = '\0';
      return field;
    }

    // KNOW: The field_suffix_sep is not in first or last position.
    // KNOW: num_tokens > 1 at this point.  Possible that we still
    // just have a scalar with an embedded separator character...
    int suffix_size = 1;
    if (num_tokens > 2)
      suffix_size = 2;

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
      for (size_t i=1; i < num_tokens-suffix_size; i++) {
	base_name += suffix_separator;
	base_name += tokens[i];
      }
      base_name += suffix_separator;
      size_t bn_len = base_name.length(); // Length of basename portion only
      size_t length = std::strlen(name); // Length of total name (with suffix)

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
      for (int i = index+1; i < num_names; i++) {
	char *tst_name = names[i];
	std::vector<std::string> subtokens;
	Ioss::tokenize(tst_name,suffix,subtokens);
	if ((truth_table == NULL || truth_table[i] == 1) &&  // Defined on this entity
	    std::strlen(tst_name) == length &&              // names must be same length
	    std::strncmp(name, tst_name, bn_len) == 0 &&   // base portion must match
	    subtokens.size() == num_tokens) {
	  which_names.push_back(i);
	}
      }

      const Ioss::VariableType *type = NULL;
      if (suffix_size == 2) {
	if (which_names.size() > 1)
	  type = match_composite_field(names, which_names, suffix_separator);
      } else {
	assert(suffix_size == 1);
	type = match_single_field(names, which_names, suffix_separator);
      }

      if (type != NULL) {
	// A valid variable type was recognized.
	// Mark the names which were used so they aren't used for another field on this entity.
	// Create a field of that variable type.
	assert(type->component_count() == static_cast<int>(which_names.size()));
	Ioss::Field field(base_name.substr(0,bn_len-1), Ioss::Field::REAL, type, fld_role, count);
	for (size_t i=0; i < which_names.size(); i++) {
	  names[which_names[i]][0] = '\0';
	}
	return field;
      } else {
	if (suffix_size == 1) {
	  Ioss::Field field(name, Ioss::Field::REAL, SCALAR(), fld_role, count);
	  names[index][0] = '\0';
	  return field;
	}
      }
      suffix_size--;
    }
    return Ioss::Field("", Ioss::Field::INVALID, SCALAR(), fld_role, 1);
  }

  // common
  bool define_field(size_t nmatch, size_t match_length,
		    char **names, std::vector<Ioss::Suffix> &suffices,
		    size_t entity_count, Ioss::Field::RoleType fld_role,
		    std::vector<Ioss::Field> &fields)
  {
    // Try to define a field of size 'nmatch' with the suffices in 'suffices'.
    // If this doesn't define a known field, then assume it is a scalar instead
    // and return false.
    if (nmatch > 1) {
      const Ioss::VariableType *type = Ioss::VariableType::factory(suffices);
      if (type == NULL) {
	nmatch = 1;
      } else {
	char *name = names[0];
	name[match_length] = '\0';
	Ioss::Field field(name, Ioss::Field::REAL, type, fld_role, entity_count);
	if (field.is_valid()) {
	  fields.push_back(field);
	}
	for (size_t j = 0; j < nmatch; j++)
	  names[j][0] = '\0';
	return true;
      }
    }

    // NOTE: nmatch could be reset inside previous if block.
    // This is not an 'else' block, it is a new if block.
    if (nmatch == 1) {
      Ioss::Field field(names[0], Ioss::Field::REAL, SCALAR(), fld_role, entity_count);
      if (field.is_valid()) {
	fields.push_back(field);
      }
      names[0][0] = '\0';
      return false;
    }
    return false; // Can't get here...  Quiet the compiler
  }

  // Read scalar fields off an input database and determine whether
  // they are components of a higher order type (vector, tensor, ...).
  // This routine is used if there is no field component separator.  E.g.,
  // fieldx, fieldy, fieldz instead of field_x field_y field_z

  void get_fields(int64_t entity_count, // The number of objects in this entity.
		  char** names,     // Raw list of field names from exodus
		  size_t num_names, // Number of names in list
		  Ioss::Field::RoleType fld_role, // Role of field
		  const char suffix_separator,
		  int *local_truth, // Truth table for this entity;
		  // null if not applicable.
		  std::vector<Ioss::Field> &fields) // The fields that were found.
  {
    if (suffix_separator != 0) {
      while (1) {
	// NOTE: 'get_next_field' determines storage type (vector, tensor,...)
	Ioss::Field field = get_next_field(names, num_names, entity_count, fld_role,
					   suffix_separator, local_truth);
	if (field.is_valid()) {
	  fields.push_back(field);
	} else {
	  break;
	}
      }
    } else {
      size_t nmatch = 1;
      size_t ibeg   = 0;
      size_t pmat   = 0;
      std::vector<Ioss::Suffix> suffices;
    top:

      while (ibeg+nmatch < num_names) {
	if (local_truth != NULL) {
	  while (ibeg < num_names && local_truth[ibeg] == 0)
	    ibeg++;
	}
	for (size_t i=ibeg+1; i < num_names; i++) {
	  size_t mat = match(names[ibeg], names[i]);
	  if (local_truth != NULL && local_truth[i] == 0)
	    mat = 0;

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
	  if ((std::strlen(names[ibeg]) == std::strlen(names[i])) &&
	      mat > 0 && (pmat == 0 || mat >= pmat)) {
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
	  } else {

	    bool multi_component = define_field(nmatch, pmat, &names[ibeg], suffices,
						entity_count, fld_role, fields);
	    if (!multi_component) {
	      // Although we matched multiple suffices, it wasn't a
	      // higher-order field, so we only used 1 name instead of
	      // the 'nmatch' we thought we might use.
	      i = ibeg + 1;
	    }

	    // Cleanout the suffices vector.
	    std::vector<Ioss::Suffix>().swap(suffices);

	    // Reset for the next time through the while loop...
	    nmatch=1;
	    pmat = 0;
	    ibeg=i;
	    break;
	  }
	}
      }
      // We've gone through the entire list of names; see if what we
      // have forms a multi-component field; if not, then define a
      // scalar field and jump up to the loop again to handle the others
      // that had been gathered.
      if (ibeg < num_names) {
	if (local_truth == NULL || local_truth[ibeg] == 1) {
	  bool multi_component = define_field(nmatch, pmat, &names[ibeg], suffices,
					      entity_count, fld_role, fields);
	  std::vector<Ioss::Suffix>().swap(suffices);
	  if (nmatch > 1 && !multi_component) {
	    ibeg++;
	    goto top;
	  }
	} else {
	  ibeg++;
	  goto top;
	}
      }
    }
  }

  void check_non_null(void *ptr, const char *type, const std::string &name)
  {
    if (ptr == NULL) {
      std::ostringstream errmsg;
      errmsg << "INTERNAL ERROR: Could not find " << type << " '" << name << "'."
	     << " Something is wrong in the Ioex::DatabaseIO class. Please report.\n";
      IOSS_ERROR(errmsg);
    }
  }

  void add_map_fields(int exoid, Ioss::ElementBlock *block, int64_t my_element_count,
		      size_t name_length)
  {
    // Check for optional element maps...
    int map_count = ex_inquire_int(exoid, EX_INQ_ELEM_MAP);
    if (map_count <= 0)
      return;

    // Get the names of the maps...
    char **names = Ioex::get_exodus_names(map_count, name_length);
    int ierr = ex_get_names(exoid, EX_ELEM_MAP, names);
    if (ierr < 0)
      Ioex::exodus_error(exoid, __LINE__, -1);

    // Convert to lowercase.
    for (int i=0; i < map_count; i++) {
      Ioss::Utils::fixup_name(names[i]);
    }

    if (map_count == 2 && std::strncmp(names[0], "skin:", 5) == 0 && std::strncmp(names[1], "skin:", 5) == 0) {
      // Currently, only support the "skin" map -- It will be a 2
      // component field consisting of "parent_element":"local_side"
      // pairs.  The parent_element is an element in the original mesh,
      // not this mesh.
      block->field_add(Ioss::Field("skin", block->field_int_type(), "Real[2]",
				   Ioss::Field::MESH, my_element_count));
    }
    Ioex::delete_exodus_names(names, map_count);
  }

  void write_coordinate_frames(int exoid, const Ioss::CoordinateFrameContainer &frames) {
    if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
      internal_write_coordinate_frames(exoid, frames, (int64_t)0);
    }
    else {
      internal_write_coordinate_frames(exoid, frames, (int)0);
    }
  }

  void add_coordinate_frames(int exoid, Ioss::Region *region)
  {
    if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
      internal_add_coordinate_frames(exoid, region, (int64_t)0);
    }
    else {
      internal_add_coordinate_frames(exoid, region, (int)0);
    }
  }

  bool filter_node_list(Ioss::Int64Vector &nodes,
			const std::vector<unsigned char> &node_connectivity_status)
  {
    // Iterate through 'nodes' and determine which of the nodes are
    // not connected to any non-omitted blocks. The index of these
    // nodes is then put in the 'nodes' list.
    // Assumes that there is at least one omitted element block.  The
    // 'nodes' list on entry contains 1-based local node ids, not global.
    // On return, the nodes list contains indices.  To filter a nodeset list:
    // for (size_t i = 0; i < nodes.size(); i++) {
    //    active_values[i] = some_nset_values[nodes[i]];
    // }

    size_t orig_size = nodes.size();
    size_t active = 0;
    for (size_t i=0; i < orig_size; i++) {
      if (node_connectivity_status[nodes[i]-1] >= 2) {
	// Node is connected to at least 1 active element...
	nodes[active++] = i;
      }
    }
    nodes.resize(active);
    std::vector<int64_t>(nodes).swap(nodes); // shrink to fit
    return (active != orig_size);
  }

  void filter_element_list(Ioss::Region *region, 
			   Ioss::Int64Vector &elements,
			   Ioss::Int64Vector &sides,
			   bool remove_omitted_elements)
  {
    // Iterate through 'elements' and remove the elements which are in an omitted block.
    // Precondition is that there is at least one omitted element block.
    // The 'elements' list contains local element ids, not global.
    // Since there are typically a small number of omitted blocks, do
    // the following:
    // For each omitted block, determine the min and max element id in
    // that block.  Iterate 'elements' vector and set the id to zero if
    // min <= id <= max.  Once all omitted blocks have been processed,
    // then iterate the vector and compress out all zeros.  Keep 'sides'
    // array consistent.

    // Get all element blocks in region...
    bool omitted = false;
    Ioss::ElementBlockContainer element_blocks = region->get_element_blocks();
    for (size_t blk=0; blk < element_blocks.size(); blk++) {
      Ioss::ElementBlock *block = element_blocks[blk];
      if (Ioss::Utils::block_is_omitted(block)) {
	ssize_t min_id = block->get_offset() + 1;
	ssize_t max_id = min_id + block->get_property("entity_count").get_int() - 1;
	for (size_t i=0; i < elements.size(); i++) {
	  if (min_id <= elements[i] && elements[i] <= max_id) {
	    omitted = true;
	    elements[i] = 0;
	    sides[i]    = 0;
	  }
	}
      }
    }
    if (remove_omitted_elements && omitted) {
      elements.erase(std::remove(elements.begin(), elements.end(), 0), elements.end());
      sides.erase(   std::remove(sides.begin(),    sides.end(),    0), sides.end());
    }
  }

  void separate_surface_element_sides(Ioss::Int64Vector &element,
				      Ioss::Int64Vector &sides,
				      Ioss::Region *region,
				      Ioex::TopologyMap &topo_map,
				      Ioex::TopologyMap &side_map,
				      Ioss::SurfaceSplitType split_type)
  {
    if (!element.empty()) {
      Ioss::ElementBlock *block = NULL;
      // Topology of sides in current element block
      const Ioss::ElementTopology *common_ftopo = NULL;
      const Ioss::ElementTopology *topo = NULL; // Topology of current side
      int64_t current_side = -1;

      for (size_t iel = 0; iel < element.size(); iel++) {
	int64_t elem_id = element[iel];
	if (block == NULL || !block->contains(elem_id)) {
	  block = region->get_element_block(elem_id);
	  assert(block != NULL);
	  assert(!Ioss::Utils::block_is_omitted(block)); // Filtered out above.

	  // NULL if hetero sides on element
	  common_ftopo = block->topology()->boundary_type(0);
	  if (common_ftopo != NULL)
	    topo = common_ftopo;
	  current_side = -1;
	}

	if (common_ftopo == NULL && sides[iel] != current_side) {
	  current_side = sides[iel];
	  assert(current_side > 0 && current_side <= block->topology()->number_boundaries());
	  topo = block->topology()->boundary_type(sides[iel]);
	  assert(topo != NULL);
	}
	std::pair<std::string, const Ioss::ElementTopology*> name_topo;
	if (split_type == Ioss::SPLIT_BY_TOPOLOGIES) {
	  name_topo = std::make_pair(block->topology()->name(), topo);
	} else if (split_type == Ioss::SPLIT_BY_ELEMENT_BLOCK) {
	  name_topo = std::make_pair(block->name(), topo);
	}
	topo_map[name_topo]++;
	if (side_map[name_topo] == 0)
	  side_map[name_topo] = sides[iel];
	else if (side_map[name_topo] != sides[iel]) {
	  // Not a consistent side for all sides in this
	  // sideset. Set to large number. Note that maximum
	  // sides/element is 6, so don't have to worry about
	  // a valid element having 999 sides (unless go to
	  // arbitrary polyhedra some time...) Using a large
	  // number instead of -1 makes it easier to check the
	  // parallel consistency...
	  side_map[name_topo] = 999;
	}
      }
    }
  }

}
