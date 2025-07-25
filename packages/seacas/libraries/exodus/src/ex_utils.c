/*
 * Copyright(C) 1999-2025 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * exutils - utility routines
 *
 *****************************************************************************/

#define _GNU_SOURCE
#include <stdlib.h>

#include <errno.h>
#include <stdbool.h>

#include "exodusII.h"
#include "exodusII_int.h"

struct exi_obj_stats *exoII_eb  = NULL;
struct exi_obj_stats *exoII_ed  = NULL;
struct exi_obj_stats *exoII_fa  = NULL;
struct exi_obj_stats *exoII_ns  = NULL;
struct exi_obj_stats *exoII_es  = NULL;
struct exi_obj_stats *exoII_fs  = NULL;
struct exi_obj_stats *exoII_ss  = NULL;
struct exi_obj_stats *exoII_els = NULL;
struct exi_obj_stats *exoII_em  = NULL;
struct exi_obj_stats *exoII_edm = NULL;
struct exi_obj_stats *exoII_fam = NULL;
struct exi_obj_stats *exoII_nm  = NULL;
struct exi_obj_stats *exoII_ass = NULL;

/*****************************************************************************
 *
 * utility routines for string conversions
 * exi_catstr  - concatenate  string/number (where number is converted to ASCII)
 * exi_catstr2 - concatenate  string1/number1/string2/number2   "
 *
 * NOTE: these routines reuse the same storage over and over to build
 *        concatenated strings, because the strings are just passed to netCDF
 *        routines as names used to look up variables.  if the strings returned
 *        by these routines are needed for any other purpose, they should
 *        immediately be copied into other storage.
 *****************************************************************************/

static char  ret_string[10 * (MAX_VAR_NAME_LENGTH + 1)];
static char *cur_string = &ret_string[0];

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
#define WINDOWS__ 1
#endif

#ifndef WINDOWS__
#if (NC_VERSION_MAJOR > 4 || (NC_VERSION_MAJOR == 4 && NC_VERSION_MINOR > 4)) && NC_HAS_HDF5
extern int NC4_hdf5get_libversion(unsigned *, unsigned *, unsigned *);
#endif
#endif

#if 0
#if NC_HAS_PNETCDF
extern char *ncmpi_inq_libvers();
#endif
#endif

/*!
  \ingroup Utilities
  \undoc
*/
const char *ex_config(void)
{
  static char buffer[2048];
  int         buffer_size = sizeof(buffer) / sizeof(buffer[0]);

  int j = snprintf(buffer, buffer_size, "\tExodus Version %s, Released %s\n", EXODUS_VERSION,
                   EXODUS_RELEASE_DATE);
#if defined(PARALLEL_AWARE_EXODUS)
  j += snprintf(buffer + j, buffer_size - j, "\t\tParallel enabled\n");
#else
  j += snprintf(buffer + j, buffer_size - j, "\t\tParallel NOT enabled\n");
#endif
#if defined(EXODUS_THREADSAFE)
  j += snprintf(buffer + j, buffer_size - j, "\t\tThread Safe enabled\n");
#else
  j += snprintf(buffer + j, buffer_size - j, "\t\tThread Safe NOT enabled\n");
#endif
#if defined(SEACAS_HIDE_DEPRECATED_CODE)
  j += snprintf(buffer + j, buffer_size - j, "\t\tDeprecated Functions NOT built\n\n");
#else
  j += snprintf(buffer + j, buffer_size - j, "\t\tDeprecated Functions available\n\n");
#endif
#if defined(NC_VERSION)
  j += snprintf(buffer + j, buffer_size - j, "\tNetCDF Version %s\n", NC_VERSION);
#else
  j += snprintf(buffer + j, buffer_size - j, "\tNetCDF Version < 4.3.3\n");
#endif
#if NC_HAS_CDF5
  j += snprintf(buffer + j, buffer_size - j, "\t\tCDF5 enabled\n");
#endif
#ifndef WINDOWS__
#if NC_HAS_HDF5
  {
#if NC_VERSION_MAJOR > 4 || (NC_VERSION_MAJOR == 4 && NC_VERSION_MINOR > 4)
    unsigned major, minor, release;
    NC4_hdf5get_libversion(&major, &minor, &release);
    j += snprintf(buffer + j, buffer_size - j, "\t\tHDF5 enabled (%u.%u.%u)\n", major, minor,
                  release);
#endif
  }
  j += snprintf(buffer + j, buffer_size - j, "\t\tZlib Compression (read/write) enabled\n");
#if NC_HAS_SZIP_WRITE == 1
  j += snprintf(buffer + j, buffer_size - j, "\t\tSZip Compression (read/write) enabled\n");
#else
  j += snprintf(buffer + j, buffer_size - j, "\t\tSZip Compression (read/write) NOT enabled\n");
#endif
#if NC_HAS_ZSTD == 1
  j += snprintf(buffer + j, buffer_size - j, "\t\tZstd Compression enabled\n");
#else
  j += snprintf(buffer + j, buffer_size - j, "\t\tZstd Compression NOT enabled\n");
#endif
#if NC_HAS_QUANTIZE == 1
  j += snprintf(buffer + j, buffer_size - j, "\t\tQuantization support enabled\n");
#else
  j += snprintf(buffer + j, buffer_size - j, "\t\tQuantization support NOT enabled\n");
#endif
#endif
#endif
#if defined(PARALLEL_AWARE_EXODUS)
#if NC_HAS_PARALLEL
  j += snprintf(buffer + j, buffer_size - j, "\t\tParallel IO enabled via HDF5 and/or PnetCDF\n");
#else
  j +=
      snprintf(buffer + j, buffer_size - j,
               "\t\tParallel IO *NOT* enabled via HDF5 and/or PnetCDF (PROBABLY A BUILD ERROR!)\n");
#endif
#if NC_HAS_PARALLEL4
  j += snprintf(buffer + j, buffer_size - j, "\t\tParallel IO enabled via HDF5\n");
#else
  j += snprintf(buffer + j, buffer_size - j, "\t\tParallel IO *NOT* enabled via HDF5\n");
#endif
#if NC_HAS_PAR_FILTERS
  j += snprintf(buffer + j, buffer_size - j, "\t\tParallel IO supports filters\n");
#endif
#if NC_HAS_PNETCDF
  {
#if 0
    char *libver = ncmpi_inq_libvers();
    j += snprintf(buffer + j, buffer_size - j, "\t\tParallel IO enabled via PnetCDF (%s)\n", libver);
#else
    j += snprintf(buffer + j, buffer_size - j, "\t\tParallel IO enabled via PnetCDF.\n");
#endif
  }
#else
  j += snprintf(buffer + j, buffer_size - j, "\t\tParallel IO *NOT* enabled via PnetCDF\n");
#endif
#endif /* PARALLEL_AWARE_EXODUS */

#if NC_HAS_ERANGE_FILL
  j += snprintf(buffer + j, buffer_size - j, "\t\tERANGE_FILL support\n");
#endif
#if NC_RELAX_COORD_BOUND
  j += snprintf(buffer + j, buffer_size - j, "\t\tRELAX_COORD_BOUND defined\n");
#endif
#if defined(NC_COMPACT)
  j += snprintf(buffer + j, buffer_size - j, "\t\tNC_COMPACT defined\n");
#endif
#if defined(NC_HAVE_META_H)
  j += snprintf(buffer + j, buffer_size - j, "\t\tNC_HAVE_META_H defined\n");
#endif
#if NC_HAS_NC2
  j += snprintf(buffer + j, buffer_size - j, "\t\tAPI Version 2 support enabled\n");
#else
  j += snprintf(buffer + j, buffer_size - j, "\t\tAPI Version 2 support NOT enabled\n");
#endif
  j += snprintf(buffer + j, buffer_size - j, "\n");

  assert(j < 2048);
  return buffer;
}
void ex_print_config(void) { fprintf(stderr, "%s\n", ex_config()); }

/*!
  \ingroup Utilities
  \undoc
*/
int exi_check_file_type(const char *path, int *type)
{
  /* Based on (stolen from?) NC_check_file_type from netcdf sources.

     Type is set to:
     1 if this is a netcdf classic file,
     2 if this is a netcdf 64-bit offset file,
     4 pnetcdf cdf5 file.
     5 if this is an hdf5 file
  */

#define MAGIC_NUMBER_LEN 4

  char magic[MAGIC_NUMBER_LEN + 1];
  EX_FUNC_ENTER();

  *type = 0;

  /* Get the 4-byte magic from the beginning of the file. */
  {
    FILE *fp;
    if (!(fp = fopen(path, "r"))) {
      char errmsg[MAX_ERR_LENGTH];
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: Could not open file '%s',\n\t\terror = %s.", path,
               strerror(errno));
      ex_err(__func__, errmsg, EX_WRONGFILETYPE);
      EX_FUNC_LEAVE(EX_FATAL);
    }
    int i                   = fread(magic, 1, MAGIC_NUMBER_LEN, fp);
    magic[MAGIC_NUMBER_LEN] = '\0';
    fclose(fp);
    if (i != MAGIC_NUMBER_LEN) {
      char errmsg[MAX_ERR_LENGTH];
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: Could not read magic data from file '%s',\n\t\terror = %s.", path,
               strerror(errno));
      ex_err(__func__, errmsg, EX_WRONGFILETYPE);
      EX_FUNC_LEAVE(EX_FATAL);
    }
  }

  /* Ignore the first byte for HDF */
  if (magic[1] == 'H' && magic[2] == 'D' && magic[3] == 'F') {
    *type = 5;
  }
  else if (magic[0] == 'C' && magic[1] == 'D' && magic[2] == 'F') {
    if (magic[3] == '\001') {
      *type = 1;
    }
    else if (magic[3] == '\002') {
      *type = 2;
    }
    else if (magic[3] == '\005') {
      *type = 4; /* cdf5 (including pnetcdf) file */
    }
  }
  if (*type == 0) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(
        errmsg, MAX_ERR_LENGTH,
        "ERROR: Could not recognize %s as a valid Exodus/NetCDF file variant.  Magic value is '%s'",
        path, magic);
    ex_err(__func__, errmsg, EX_WRONGFILETYPE);
    EX_FUNC_LEAVE(EX_FATAL);
  }
  EX_FUNC_LEAVE(EX_NOERR);
}

/*!
  \ingroup Utilities
  \undoc
*/
int ex_set_max_name_length(int exoid, int length)
{
  char errmsg[MAX_ERR_LENGTH];

  EX_FUNC_ENTER();
  if (exi_check_valid_file_id(exoid, __func__) == EX_FATAL) {
    EX_FUNC_LEAVE(EX_FATAL);
  }
  if (length <= 0) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: Max name length must be positive.");
    ex_err_fn(exoid, __func__, errmsg, NC_EMAXNAME);
    EX_FUNC_LEAVE(EX_FATAL);
  }
  if (length > EX_MAX_NAME) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: Max name length (%d) exceeds netcdf max name size (%d).", length, EX_MAX_NAME);
    ex_err_fn(exoid, __func__, errmsg, NC_EMAXNAME);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  ex_set_option(exoid, EX_OPT_MAX_NAME_LENGTH, length);

  EX_FUNC_LEAVE(EX_NOERR);
}

/*!
  \ingroup Utilities
  \undoc
*/
void exi_update_max_name_length(int exoid, int length)
{
  EX_FUNC_ENTER();
  if (exi_check_valid_file_id(exoid, __func__) == EX_FATAL) {
    EX_FUNC_VOID();
  }

  /* Get current value of the maximum_name_length attribute... */
  int status;
  int db_length = 0;
  int rootid    = exoid & EX_FILE_ID_MASK;
  if ((status = nc_get_att_int(rootid, NC_GLOBAL, ATT_MAX_NAME_LENGTH, &db_length)) != EX_NOERR) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to update 'max_name_length' attribute in file id %d", exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
  }

  if (length > db_length) {
    /* Update with new value... */
    ex_set_max_name_length(exoid, length);
    if ((status = nc_put_att_int(rootid, NC_GLOBAL, ATT_MAX_NAME_LENGTH, NC_INT, 1, &length)) !=
        EX_NOERR) {
      char errmsg[MAX_ERR_LENGTH];
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to update 'max_name_length' attribute with new value in file id %d",
               exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
    }
    nc_sync(rootid);
  }
  EX_FUNC_VOID();
}

/*!
  \internal
  \undoc
*/
int exi_put_names(int exoid, int varid, size_t num_names, char *const *names,
                  ex_entity_type obj_type, const char *subtype, const char *routine)
{
  EX_FUNC_ENTER();
  if (exi_check_valid_file_id(exoid, __func__) == EX_FATAL) {
    EX_FUNC_LEAVE(EX_FATAL);
  }
  /* inquire previously defined dimensions  */
  size_t name_length = ex_inquire_int(exoid, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH) + 1;

  char *int_names = NULL;
  if (!(int_names = calloc(num_names * name_length, 1))) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to allocate memory for internal int_names "
             "array in file id %d",
             exoid);
    ex_err_fn(exoid, __func__, errmsg, EX_MEMFAIL);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  size_t idx          = 0;
  int    max_name_len = 0;
  int    found_name   = 0;
  for (size_t i = 0; i < num_names; i++) {
    if (names != NULL && *names != NULL && *names[i] != '\0') {
      found_name = 1;
      ex_copy_string(&int_names[idx], names[i], name_length);
      size_t length = strlen(names[i]) + 1;
      if (length > (size_t)name_length) {
        fprintf(stderr,
                "Warning: The %s %s name '%s' is too long.\n\tIt will "
                "be truncated from %d to %d characters. [Called from %s]\n",
                ex_name_of_object(obj_type), subtype, names[i], (int)length - 1,
                (int)name_length - 1, routine);
        length = name_length;
      }

      if (length > (size_t)max_name_len) {
        max_name_len = length;
      }
    }
    idx += name_length;
  }

  int status;
  if ((status = nc_put_var_text(exoid, varid, int_names)) != EX_NOERR) {
    free(int_names);
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to store %s names in file id %d",
             ex_name_of_object(obj_type), exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  if (found_name) {

    /* Update the maximum_name_length attribute on the file. */
    exi_update_max_name_length(exoid, max_name_len - 1);
  }
  free(int_names);

  EX_FUNC_LEAVE(EX_NOERR);
}

/*!
  \internal
  \undoc
*/
int exi_put_name(int exoid, int varid, size_t index, const char *name, ex_entity_type obj_type,
                 const char *subtype, const char *routine)
{
  if (exi_check_valid_file_id(exoid, __func__) == EX_FATAL) {
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* inquire previously defined dimensions  */
  size_t name_length = ex_inquire_int(exoid, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH) + 1;

  if (name != NULL && *name != '\0') {
    size_t start[2], count[2];
    int    too_long = 0;
    start[0]        = index;
    start[1]        = 0;

    count[0] = 1;
    count[1] = strlen(name) + 1;

    if (count[1] > name_length) {
      fprintf(stderr,
              "Warning: The %s %s name '%s' is too long.\n\tIt will be "
              "truncated from %d to %d characters. [Called from %s]\n",
              ex_name_of_object(obj_type), subtype, name, (int)strlen(name), (int)name_length - 1,
              routine);
      count[1] = name_length;
      too_long = 1;
    }

    int status;
    if ((status = nc_put_vara_text(exoid, varid, start, count, name)) != EX_NOERR) {
      char errmsg[MAX_ERR_LENGTH];
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to store %s name in file id %d",
               ex_name_of_object(obj_type), exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      return EX_FATAL;
    }

    /* Add the trailing null if the variable name was too long */
    if (too_long) {
      start[1] = name_length - 1;
      nc_put_var1_text(exoid, varid, start, "\0");
    }

    /* Update the maximum_name_length attribute on the file. */
    exi_update_max_name_length(exoid, count[1] - 1);
  }
  return EX_NOERR;
}

/*!
  \internal
  \undoc
*/
int exi_get_names(int exoid, int varid, size_t num_names, char **names, ex_entity_type obj_type,
                  const char *routine)
{
  /* Query size of names on file
   * Use the smaller of the size on file or user-specified length
   */
  int db_name_size  = ex_inquire_int(exoid, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH);
  int api_name_size = ex_inquire_int(exoid, EX_INQ_MAX_READ_NAME_LENGTH);
  int name_size     = db_name_size < api_name_size ? db_name_size : api_name_size;

  for (size_t i = 0; i < num_names; i++) {
    int status = exi_get_name(exoid, varid, i, names[i], name_size, obj_type, routine);
    if (status != EX_NOERR) {
      return status;
    }
  }
  return EX_NOERR;
}

/*!
  \internal
  \undoc
*/
int exi_get_name(int exoid, int varid, size_t index, char *name, int name_size,
                 ex_entity_type obj_type, const char *routine)
{

  /* read the name */
  size_t start[2], count[2];
  start[0] = index;
  count[0] = 1;
  start[1] = 0;
  count[1] = name_size + 1;

  int status = nc_get_vara_text(exoid, varid, start, count, name);
  if (status != EX_NOERR) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to get %s name at index %d from file id %d [Called from %s]",
             ex_name_of_object(obj_type), (int)index, exoid, routine);
    ex_err_fn(exoid, __func__, errmsg, status);
    return EX_FATAL;
  }

  int api_name_size   = ex_inquire_int(exoid, EX_INQ_MAX_READ_NAME_LENGTH);
  name[api_name_size] = '\0';

  exi_trim(name);
  return EX_NOERR;
}

/*!
  \internal
  \undoc
*/
void exi_trim(char *name)
{
  /* Thread-safe, reentrant */
  /* Trim trailing spaces... */
  if (name == NULL) {
    return;
  }

  size_t size = strlen(name);
  if (size == 0) {
    return;
  }

  char *end = name + size - 1;
  while (end >= name && isspace(*end)) {
    end--;
  }

  *(end + 1) = '\0';
}

/** exi_catstr  - concatenate  string/number (where number is converted to ASCII)
 */
/*!
  \internal
  \undoc
*/
char *exi_catstr(const char *string, int num)
{
  /* Only called from an already locked function */
  char *tmp_string = cur_string;
  cur_string += snprintf(cur_string, MAX_VAR_NAME_LENGTH + 1, "%s%d", string, num) + 1;
  if (cur_string - ret_string > 9 * (MAX_VAR_NAME_LENGTH + 1)) {
    cur_string = ret_string;
  }
  return tmp_string;
}

/** exi_catstr2 - concatenate  string1num1string2num2   */
/*!
  \internal
  \undoc
*/
char *exi_catstr2(const char *string1, int num1, const char *string2, int num2)
{
  /* Only called from an already locked function */
  char *tmp_string = cur_string;
  cur_string +=
      snprintf(cur_string, MAX_VAR_NAME_LENGTH + 1, "%s%d%s%d", string1, num1, string2, num2) + 1;
  if (cur_string - ret_string > 9 * (MAX_VAR_NAME_LENGTH + 1)) {
    cur_string = ret_string;
  }
  return tmp_string;
}

/*!
  \internal
  \undoc
*/
char *ex_name_of_object(ex_entity_type obj_type)
{
  /* Thread-safe and reentrant */
  switch (obj_type) {
  case EX_ASSEMBLY: return "assembly";
  case EX_BLOB: return "blob";
  case EX_COORDINATE: /* kluge so some wrapper functions work */ return "coordinate";
  case EX_NODAL: return "nodal";
  case EX_EDGE_BLOCK: return "edge block";
  case EX_FACE_BLOCK: return "face block";
  case EX_ELEM_BLOCK: return "element block";
  case EX_NODE_SET: return "node set";
  case EX_EDGE_SET: return "edge set";
  case EX_FACE_SET: return "face set";
  case EX_SIDE_SET: return "side set";
  case EX_ELEM_SET: return "element set";
  case EX_ELEM_MAP: return "element map";
  case EX_NODE_MAP: return "node map";
  case EX_EDGE_MAP: return "edge map";
  case EX_FACE_MAP: return "face map";
  case EX_GLOBAL: return "global";
  default: return "invalid type";
  }
}

/*!
  \internal
  \undoc
*/
ex_entity_type ex_var_type_to_ex_entity_type(char var_type)
{
  /* Thread-safe and reentrant */
  char var_lower = tolower(var_type);
  if (var_lower == 'n') {
    return EX_NODAL;
  }
  if (var_lower == 'l') {
    return EX_EDGE_BLOCK;
  }
  if (var_lower == 'f') {
    return EX_FACE_BLOCK;
  }
  if (var_lower == 'e') {
    return EX_ELEM_BLOCK;
  }
  if (var_lower == 'm') {
    return EX_NODE_SET;
  }
  else if (var_lower == 'd') {
    return EX_EDGE_SET;
  }
  else if (var_lower == 'a') {
    return EX_FACE_SET;
  }
  else if (var_lower == 's') {
    return EX_SIDE_SET;
  }
  else if (var_lower == 't') {
    return EX_ELEM_SET;
  }
  else if (var_lower == 'g') {
    return EX_GLOBAL;
  }
  else {
    return EX_INVALID;
  }
}

/*!
  \internal
  \undoc
*/
char *exi_dim_num_objects(ex_entity_type obj_type)
{
  switch (obj_type) {
  case EX_NODAL: return DIM_NUM_NODES;
  case EX_ASSEMBLY: return DIM_NUM_ASSEMBLY;
  case EX_BLOB: return DIM_NUM_BLOB;
  case EX_ELEM_BLOCK: return DIM_NUM_EL_BLK;
  case EX_EDGE_BLOCK: return DIM_NUM_ED_BLK;
  case EX_FACE_BLOCK: return DIM_NUM_FA_BLK;
  case EX_NODE_SET: return DIM_NUM_NS;
  case EX_EDGE_SET: return DIM_NUM_ES;
  case EX_FACE_SET: return DIM_NUM_FS;
  case EX_ELEM_SET: return DIM_NUM_ELS;
  case EX_SIDE_SET: return DIM_NUM_SS;
  case EX_ELEM_MAP: return DIM_NUM_EM;
  case EX_FACE_MAP: return DIM_NUM_FAM;
  case EX_EDGE_MAP: return DIM_NUM_EDM;
  case EX_NODE_MAP: return DIM_NUM_NM;
  default: {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: object type %d not supported in call to %s", obj_type,
             __func__);
    ex_err(__func__, errmsg, EX_BADPARAM);
    return NULL;
  }
  }
}

/*!
  \internal
  \undoc
*/
char *exi_dim_num_entries_in_object(ex_entity_type obj_type, int idx)
{
  switch (obj_type) {
  case EX_NODAL: return DIM_NUM_NODES;
  case EX_ASSEMBLY: return DIM_NUM_ENTITY_ASSEMBLY(idx);
  case EX_BLOB: return DIM_NUM_VALUES_BLOB(idx);
  case EX_EDGE_BLOCK: return DIM_NUM_ED_IN_EBLK(idx);
  case EX_FACE_BLOCK: return DIM_NUM_FA_IN_FBLK(idx);
  case EX_ELEM_BLOCK: return DIM_NUM_EL_IN_BLK(idx);
  case EX_NODE_SET: return DIM_NUM_NOD_NS(idx);
  case EX_EDGE_SET: return DIM_NUM_EDGE_ES(idx);
  case EX_FACE_SET: return DIM_NUM_FACE_FS(idx);
  case EX_SIDE_SET: return DIM_NUM_SIDE_SS(idx);
  case EX_ELEM_SET: return DIM_NUM_ELE_ELS(idx);
  default: {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: object type %d not supported in call to %s", obj_type,
             __func__);
    ex_err(__func__, errmsg, EX_BADPARAM);
    return NULL;
  }
  }
}

/*!
  \internal
  \undoc
*/
char *exi_name_var_of_object(ex_entity_type obj_type, int i, int j)
{
  switch (obj_type) {
  case EX_ASSEMBLY: return VAR_ASSEMBLY_VAR(i, j);
  case EX_BLOB: return VAR_BLOB_VAR(i, j);
  case EX_EDGE_BLOCK: return VAR_EDGE_VAR(i, j);
  case EX_FACE_BLOCK: return VAR_FACE_VAR(i, j);
  case EX_ELEM_BLOCK: return VAR_ELEM_VAR(i, j);
  case EX_NODE_SET: return VAR_NS_VAR(i, j);
  case EX_EDGE_SET: return VAR_ES_VAR(i, j);
  case EX_FACE_SET: return VAR_FS_VAR(i, j);
  case EX_SIDE_SET: return VAR_SS_VAR(i, j);
  case EX_ELEM_SET: return VAR_ELS_VAR(i, j);
  default: {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: object type %d not supported in call to %s", obj_type,
             __func__);
    ex_err(__func__, errmsg, EX_BADPARAM);
    return NULL;
  }
  }
}

/*!
  \internal
  \undoc
*/
char *exi_name_red_var_of_object(ex_entity_type obj_type, int id)
{
  switch (obj_type) {
  case EX_ASSEMBLY: return VAR_ASSEMBLY_RED_VAR(id);
  case EX_BLOB: return VAR_BLOB_RED_VAR(id);
  case EX_EDGE_BLOCK: return VAR_EDGE_RED_VAR(id);
  case EX_FACE_BLOCK: return VAR_FACE_RED_VAR(id);
  case EX_ELEM_BLOCK: return VAR_ELEM_RED_VAR(id);
  case EX_NODE_SET: return VAR_NS_RED_VAR(id);
  case EX_EDGE_SET: return VAR_ES_RED_VAR(id);
  case EX_FACE_SET: return VAR_FS_RED_VAR(id);
  case EX_SIDE_SET: return VAR_SS_RED_VAR(id);
  case EX_ELEM_SET: return VAR_ELS_RED_VAR(id);
  default: {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: object type %d not supported in call to %s", obj_type,
             __func__);
    ex_err(__func__, errmsg, EX_BADPARAM);
    return NULL;
  }
  }
}

/*!
  \internal
  \undoc
*/
char *exi_name_of_map(ex_entity_type map_type, int map_index)
{
  switch (map_type) {
  case EX_NODE_MAP: return VAR_NODE_MAP(map_index);
  case EX_EDGE_MAP: return VAR_EDGE_MAP(map_index);
  case EX_FACE_MAP: return VAR_FACE_MAP(map_index);
  case EX_ELEM_MAP: return VAR_ELEM_MAP(map_index);
  default: return NULL;
  }
}

/*****************************************************************************
*
* exi_id_lkup - look up id
*
* entry conditions -
*   input parameters:
*       int            exoid             exodus file id
*       ex_entity_type id_type           id type name:
*                                         elem_ss
*                                         node_ns
2*                                         side_ss
*       int     num                     id value
*
* exit conditions -
*       int     return                  index into table (1-based)
*
*****************************************************************************/

/*!
  \internal
  \undoc
*/
int exi_id_lkup(int exoid, ex_entity_type id_type, ex_entity_id num)
{
  char    *id_table   = NULL;
  char    *id_dim     = NULL;
  char    *stat_table = NULL;
  int      varid, dimid;
  size_t   dim_len;
  int64_t  i;
  int64_t *id_vals   = NULL;
  int     *stat_vals = NULL;

  static bool           sequential = false;
  struct exi_obj_stats *tmp_stats;
  int                   status;
  char                  errmsg[MAX_ERR_LENGTH];

  switch (id_type) {
  case EX_NODAL: return 0;
  case EX_GLOBAL: return 0;
  case EX_ASSEMBLY: return num;
  case EX_BLOB: return num;
  case EX_ELEM_BLOCK:
    id_table   = VAR_ID_EL_BLK;   /* id array name */
    id_dim     = DIM_NUM_EL_BLK;  /* id array dimension name*/
    stat_table = VAR_STAT_EL_BLK; /* id status array name */
    tmp_stats  = exi_get_stat_ptr(exoid, &exoII_eb);
    break;
  case EX_NODE_SET:
    id_table   = VAR_NS_IDS;
    id_dim     = DIM_NUM_NS;
    stat_table = VAR_NS_STAT;
    tmp_stats  = exi_get_stat_ptr(exoid, &exoII_ns);
    break;
  case EX_SIDE_SET:
    id_table   = VAR_SS_IDS;
    id_dim     = DIM_NUM_SS;
    stat_table = VAR_SS_STAT;
    tmp_stats  = exi_get_stat_ptr(exoid, &exoII_ss);
    break;
  case EX_EDGE_BLOCK:
    id_table   = VAR_ID_ED_BLK;
    id_dim     = DIM_NUM_ED_BLK;
    stat_table = VAR_STAT_ED_BLK;
    tmp_stats  = exi_get_stat_ptr(exoid, &exoII_ed);
    break;
  case EX_FACE_BLOCK:
    id_table   = VAR_ID_FA_BLK;
    id_dim     = DIM_NUM_FA_BLK;
    stat_table = VAR_STAT_FA_BLK;
    tmp_stats  = exi_get_stat_ptr(exoid, &exoII_fa);
    break;
  case EX_EDGE_SET:
    id_table   = VAR_ES_IDS;
    id_dim     = DIM_NUM_ES;
    stat_table = VAR_ES_STAT;
    tmp_stats  = exi_get_stat_ptr(exoid, &exoII_es);
    break;
  case EX_FACE_SET:
    id_table   = VAR_FS_IDS;
    id_dim     = DIM_NUM_FS;
    stat_table = VAR_FS_STAT;
    tmp_stats  = exi_get_stat_ptr(exoid, &exoII_fs);
    break;
  case EX_ELEM_SET:
    id_table   = VAR_ELS_IDS;
    id_dim     = DIM_NUM_ELS;
    stat_table = VAR_ELS_STAT;
    tmp_stats  = exi_get_stat_ptr(exoid, &exoII_els);
    break;
  case EX_NODE_MAP:
    id_table   = VAR_NM_PROP(1);
    id_dim     = DIM_NUM_NM;
    stat_table = "";
    tmp_stats  = exi_get_stat_ptr(exoid, &exoII_nm);
    break;
  case EX_EDGE_MAP:
    id_table   = VAR_EDM_PROP(1);
    id_dim     = DIM_NUM_EDM;
    stat_table = "";
    tmp_stats  = exi_get_stat_ptr(exoid, &exoII_edm);
    break;
  case EX_FACE_MAP:
    id_table   = VAR_FAM_PROP(1);
    id_dim     = DIM_NUM_FAM;
    stat_table = "";
    tmp_stats  = exi_get_stat_ptr(exoid, &exoII_fam);
    break;
  case EX_ELEM_MAP:
    id_table   = VAR_EM_PROP(1);
    id_dim     = DIM_NUM_EM;
    stat_table = "";
    tmp_stats  = exi_get_stat_ptr(exoid, &exoII_em);
    break;
  default:
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: unsupported id array type %d for file id %d", id_type,
             exoid);
    ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
    return EX_FATAL;
  }

  if ((tmp_stats->id_vals == NULL) || (!(tmp_stats->valid_ids))) {

    /* first time through or id arrays haven't been completely filled yet */

    /* get size of id array */

    /* First get dimension id of id array */
    if ((status = nc_inq_dimid(exoid, id_dim, &dimid)) != EX_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to locate id array dimension in file id %d",
               exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      return EX_FATAL;
    }

    /* Next get value of dimension */
    if ((status = nc_inq_dimlen(exoid, dimid, &dim_len)) != EX_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to locate %s array length in file id %d",
               id_table, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      return EX_FATAL;
    }

    /* get variable id of id array */
    if ((status = nc_inq_varid(exoid, id_table, &varid)) != EX_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to locate %s array in file id %d", id_table,
               exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      return EX_FATAL;
    }

    /* allocate space for id array and initialize to zero to ensure
       that the higher bits don't contain garbage while copy from ints */
    if (!(id_vals = calloc(dim_len, sizeof(int64_t)))) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to allocate memory for %s array for file id %d", id_table, exoid);
      ex_err_fn(exoid, __func__, errmsg, EX_MEMFAIL);
      return EX_FATAL;
    }

    if (ex_int64_status(exoid) & EX_IDS_INT64_API) {
      status = nc_get_var_longlong(exoid, varid, (long long *)id_vals);
    }
    else {
      int *id_vals_int;
      if (!(id_vals_int = malloc(dim_len * sizeof(int)))) {
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: failed to allocate memory for temporary array "
                 "id_vals_int for file id %d",
                 exoid);
        ex_err_fn(exoid, __func__, errmsg, EX_MEMFAIL);
        free(id_vals);
        return EX_FATAL;
      }
      status = nc_get_var_int(exoid, varid, id_vals_int);
      if (status == EX_NOERR) {
        for (size_t iii = 0; iii < dim_len; iii++) {
          id_vals[iii] = (int64_t)id_vals_int[iii];
        }
      }
      free(id_vals_int);
    }

    if (status != EX_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to get %s array from file id %d", id_table,
               exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      free(id_vals);
      return EX_FATAL;
    }

    /* check if values in stored arrays are filled with non-zeroes */
    bool filled = true;
    sequential  = true;
    for (size_t iii = 0; iii < dim_len; iii++) {
      if (id_vals[iii] != (int64_t)iii + 1) {
        sequential = false;
      }
      if (id_vals[iii] == EX_INVALID_ID || id_vals[iii] == NC_FILL_INT) {
        filled     = false;
        sequential = false;
        break; /* id array hasn't been completely filled with valid ids yet */
      }
    }

    if (filled) {
      tmp_stats->valid_ids  = true;
      tmp_stats->sequential = sequential;
      tmp_stats->num        = dim_len;
      tmp_stats->id_vals    = id_vals;
    }
  }
  else {
    id_vals    = tmp_stats->id_vals;
    dim_len    = tmp_stats->num;
    sequential = tmp_stats->sequential;
  }

  if (sequential && (size_t)num < dim_len) {
    i = num - 1;
  }
  else {
    /* Do a linear search through the id array to find the array value
       corresponding to the passed index number */
    for (i = 0; i < (int64_t)dim_len; i++) {
      if (id_vals[i] == num) {
        break; /* found the id requested */
      }
    }
  }
  if (i >= (int64_t)dim_len) /* failed to find id number */
  {
    if (!(tmp_stats->valid_ids)) {
      free(id_vals);
    }
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to locate id %" PRId64 " for file id %d", num,
             exoid);
    ex_set_err(__func__, errmsg, EX_LOOKUPFAIL);
    return -EX_LOOKUPFAIL; /*if we got here, the id array value doesn't exist */
  }

  /* Now check status array to see if object is null */
  if ((tmp_stats->stat_vals == NULL) || (!(tmp_stats->valid_stat))) {

    /* allocate space for new status array */
    if (!(stat_vals = malloc(dim_len * sizeof(int)))) {
      free(id_vals);
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to allocate memory for %s array for file id %d", id_table, exoid);
      ex_err_fn(exoid, __func__, errmsg, EX_MEMFAIL);
      return EX_FATAL;
    }

    /* first time through or status arrays haven't been filled yet */
    if (nc_inq_varid(exoid, stat_table, &varid) == EX_NOERR) {
      /* get variable id of status array */
      /* if status array exists, use it, otherwise assume object exists
         to be backward compatible */

      if ((status = nc_get_var_int(exoid, varid, stat_vals)) != EX_NOERR) {
        free(id_vals);
        free(stat_vals);
        snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to get %s array from file id %d",
                 stat_table, exoid);
        ex_err_fn(exoid, __func__, errmsg, status);
        return EX_FATAL;
      }
    }
    else {
      for (size_t j = 0; j < dim_len; j++) {
        stat_vals[j] = 1;
      }
    }

    if (tmp_stats->valid_ids) {
      /* status array is valid only if ids are valid */
      tmp_stats->valid_stat = true;
      tmp_stats->stat_vals  = stat_vals;
    }
  }
  else {
    stat_vals = tmp_stats->stat_vals;
  }

  if (!(tmp_stats->valid_ids)) {
    free(id_vals);
  }

  if (stat_vals[i] == 0) /* is this object null? */ {
    ex_err_fn(exoid, __func__, "", EX_NULLENTITY);
    if (!(tmp_stats->valid_stat)) {
      free(stat_vals);
    }
    return -((int)i + 1); /* return index into id array (1-based) */
  }
  if (!(tmp_stats->valid_stat)) {
    free(stat_vals);
  }
  return i + 1; /* return index into id array (1-based) */
}

/******************************************************************************
 *
 * exi_get_stat_ptr - returns a pointer to a structure of object ids
 *
 *****************************************************************************/

/*! this routine returns a pointer to a structure containing the ids of
 * element blocks, node sets, or side sets according to exoid;  if there
 * is not a structure that matches the exoid, one is created.
 *
 * NOTE: If this file contains `groups` or `change sets`, then each
 *       group will have its own stat_ptr set of lists...
 * \internal
 */

struct exi_obj_stats *exi_get_stat_ptr(int exoid, struct exi_obj_stats **obj_ptr)
{
  struct exi_obj_stats *tmp_ptr = *obj_ptr;

  while (tmp_ptr) {
    if ((tmp_ptr)->exoid == exoid) {
      break;
    }
    tmp_ptr = (tmp_ptr)->next;
  }

  if (!tmp_ptr) { /* exoid not found */
    tmp_ptr             = (struct exi_obj_stats *)calloc(1, sizeof(struct exi_obj_stats));
    tmp_ptr->exoid      = exoid;
    tmp_ptr->next       = *obj_ptr;
    tmp_ptr->id_vals    = 0;
    tmp_ptr->stat_vals  = 0;
    tmp_ptr->num        = 0;
    tmp_ptr->valid_ids  = 0;
    tmp_ptr->valid_stat = 0;
    *obj_ptr            = tmp_ptr;
  }
  return tmp_ptr;
}

/******************************************************************************
 *
 * exi_rm_stat_ptr - removes a pointer to a structure of object ids
 *
 *****************************************************************************/

/*! this routine removes a pointer to a structure containing the ids of
 * element blocks, node sets, or side sets according to exoid;  this
 * is necessary to clean up because netCDF reuses file ids;  should be
 * called from ex_close
 *
 * NOTE: If this file contains `groups` or `change sets`, then each
 *       group will have its own stat_ptr set of lists.  However,
 *       this routine is called from ex_close which only closes the
 *       root id, so we need to iterate the entire list to see if there
 *       are any subgroups of the root group.
 * \internal
 */

void exi_rm_stat_ptr(int exoid, struct exi_obj_stats **obj_ptr)
{
  struct exi_obj_stats *tmp_ptr            = *obj_ptr;
  struct exi_obj_stats *last_head_list_ptr = *obj_ptr; /* save last head pointer */

  int root_id = exoid & EX_FILE_ID_MASK;
  while (tmp_ptr) {                                      /* Walk linked list of file ids/vals */
    if (root_id == (tmp_ptr->exoid & EX_FILE_ID_MASK)) { /* linear search for exodus file id */
      if (tmp_ptr == *obj_ptr) {                         /* Are we at the head of the list? */
        *obj_ptr = (*obj_ptr)->next;                     /*   yes, reset ptr to head of list */
      }
      else { /*   no, remove this record from chain*/
        last_head_list_ptr->next = tmp_ptr->next;
      }
      struct exi_obj_stats *tmp = tmp_ptr;
      tmp_ptr                   = tmp_ptr->next;
      free(tmp->id_vals); /* free up memory */
      free(tmp->stat_vals);
      free(tmp);
    }
    else {
      last_head_list_ptr = tmp_ptr;       /* save last head pointer */
      tmp_ptr            = tmp_ptr->next; /* Loop back if not */
    }
  }
}

/* structures to hold number of blocks of that type for each file id */
static struct exi_list_item *ed_ctr_list = NULL; /* edge blocks */
static struct exi_list_item *fa_ctr_list = NULL; /* face blocks */
static struct exi_list_item *eb_ctr_list = NULL; /* element blocks */

/* structures to hold number of sets of that type for each file id */
static struct exi_list_item *ns_ctr_list  = NULL; /* node sets */
static struct exi_list_item *es_ctr_list  = NULL; /* edge sets */
static struct exi_list_item *fs_ctr_list  = NULL; /* face sets */
static struct exi_list_item *ss_ctr_list  = NULL; /* side sets */
static struct exi_list_item *els_ctr_list = NULL; /* element sets */

/* structures to hold number of blobs/assemblies for each file id */
static struct exi_list_item *assm_ctr_list = NULL; /* assemblies */
static struct exi_list_item *blob_ctr_list = NULL; /* blobs */

/* structures to hold number of maps of that type for each file id */
static struct exi_list_item *nm_ctr_list  = NULL; /* node maps */
static struct exi_list_item *edm_ctr_list = NULL; /* edge maps */
static struct exi_list_item *fam_ctr_list = NULL; /* face maps */
static struct exi_list_item *em_ctr_list  = NULL; /* element maps */

/*!
  \internal
  \undoc
*/
struct exi_list_item **exi_get_counter_list(ex_entity_type obj_type)
{
  /* Thread-safe, but is dealing with globals */
  /* Only called from a routine which will be using locks */
  switch (obj_type) {
  case EX_ASSEMBLY: return &assm_ctr_list;
  case EX_BLOB: return &blob_ctr_list;
  case EX_ELEM_BLOCK: return &eb_ctr_list;
  case EX_NODE_SET: return &ns_ctr_list;
  case EX_SIDE_SET: return &ss_ctr_list;
  case EX_ELEM_MAP: return &em_ctr_list;
  case EX_NODE_MAP: return &nm_ctr_list;
  case EX_EDGE_BLOCK: return &ed_ctr_list;
  case EX_FACE_BLOCK: return &fa_ctr_list;
  case EX_EDGE_SET: return &es_ctr_list;
  case EX_FACE_SET: return &fs_ctr_list;
  case EX_ELEM_SET: return &els_ctr_list;
  case EX_EDGE_MAP: return &edm_ctr_list;
  case EX_FACE_MAP: return &fam_ctr_list;
  default: return NULL;
  }
}

/******************************************************************************
 *
 * exi_inc_file_item - increment file item
 *
 *****************************************************************************/

/*! this routine sets up a structure to track and increment a counter for
 * each open exodus file.  it is designed to be used by the routines
 * ex_put_elem_block() and ex_put_set_param(),
 * to keep track of the number of element blocks, and each type of set,
 * respectively, for each open exodus II file.
 *
 * The list structure is used as follows:
 *
 *   ptr -----------> list item structure
 *                    -------------------
 *                    exodus file id
 *                    item value (int)
 *                    ptr to next (NULL if last)
 *
 *
 * NOTE: since netCDF reuses its file ids, and a user may open and close any
 *       number of files in one application, items must be taken out of the
 *       linked lists in each of the above routines.  these should be called
 *       after ncclose().
 * \internal
 */

/* NOTE: If this is done for a file which contains `groups` or `change sets`, then
 * the exoid refers to the specific group id and not the `root_id`...
 */
int exi_inc_file_item(int                    exoid,    /* file id */
                      struct exi_list_item **list_ptr) /* ptr to ptr to list_item */
{
  struct exi_list_item *tlist_ptr = *list_ptr; /* use temp list ptr to walk linked list */
  while (tlist_ptr) {                          /* Walk linked list of file ids/vals */
    if (exoid == tlist_ptr->exo_id) {          /* linear search for exodus file id */
      break;                                   /* Quit if found */
    }
    tlist_ptr = tlist_ptr->next; /* Loop back if not */
  }

  if (!tlist_ptr) { /* ptr NULL? */
    /* allocate space for new structure record */
    tlist_ptr         = (struct exi_list_item *)calloc(1, sizeof(struct exi_list_item));
    tlist_ptr->exo_id = exoid;     /* insert file id */
    tlist_ptr->next   = *list_ptr; /* insert into head of list */
    *list_ptr         = tlist_ptr; /* fix up new head of list  */
  }
  return tlist_ptr->value++;
}

/*****************************************************************************
 *
 * exi_get_file_item - return count
 *
 *****************************************************************************/

/*! this routine accesses a structure to track and increment a counter for
 * each open exodus file.  it is designed to be used by the routines
 * ex_put_elem_block(), and ex_put_set_param(),
 * to get the number of element blocks, or a type of set,
 * respectively, for an open exodus II file.
 *
 * The list structure is used as follows:
 *
 *   ptr -----------> list item structure
 *                    -------------------
 *                    exodus file id
 *                    item value (int)
 *                    ptr to next (NULL if last)
 *
 *
 * NOTE: since netCDF reuses its file ids, and a user may open and close any
 *       number of files in one application, items must be taken out of the
 *       linked lists in each of the above routines.  these should be called
 *       after nc_close().
 * \internal
 */

int exi_get_file_item(int                    exoid,    /* file id */
                      struct exi_list_item **list_ptr) /* ptr to ptr to list_item */
{
  /* Not thread-safe: list_ptr passed in is a global
   * Would probably work ok with multiple threads since read-only,
   * but possible that list_ptr will be modified while being used
   */
  struct exi_list_item *tlist_ptr = *list_ptr; /* use temp list ptr to walk linked list */
  while (tlist_ptr) {                          /* Walk linked list of file ids/vals */
    if (exoid == tlist_ptr->exo_id) {          /* linear search for exodus file id */
      break;                                   /* Quit if found */
    }
    tlist_ptr = tlist_ptr->next; /* Loop back if not */
  }

  if (!tlist_ptr) { /* ptr NULL? */
    return -1;
  }

  return tlist_ptr->value;
}

/*****************************************************************************
 *
 * exi_rm_file_item - remove file item
 *
 *****************************************************************************/

/*! this routine removes a structure to track and increment a counter for
 * each open exodus file.
 *
 * The list structure is used as follows:
 *
 *   ptr -----------> list item structure
 *                    -------------------
 *                    exodus file id
 *                    item value (int)
 *                    ptr to next (NULL if last)
 *
 *
 * NOTE: since netCDF reuses its file ids, and a user may open and close any
 *       number of files in one application, items must be taken out of the
 *       linked lists in each of the above routines.  these should be called
 *       after ncclose().
 *
 * NOTE: this is called from `ex_close` which only closes the root file and
 *       not each group  or change set within that `root` file. If the file
 *       contains groups, then we need to check the entire list...
 * \internal
 */

void exi_rm_file_item(int                    exoid,    /* file id */
                      struct exi_list_item **list_ptr) /* ptr to ptr to list_item */

{
  struct exi_list_item *last_head_list_ptr = *list_ptr; /* save last head pointer */

  int root_id = exoid & EX_FILE_ID_MASK;

  struct exi_list_item *tlist_ptr = *list_ptr;
  while (tlist_ptr) {                                       /* Walk linked list of file ids/vals */
    if (root_id == (tlist_ptr->exo_id & EX_FILE_ID_MASK)) { /* linear search for exodus file id */
      if (tlist_ptr == *list_ptr) {                         /* Are we at the head of the list? */
        *list_ptr = (*list_ptr)->next;                      /*   yes, reset ptr to head of list */
      }
      else { /*   no, remove this record from chain*/
        last_head_list_ptr->next = tlist_ptr->next;
      }
      struct exi_list_item *temp = tlist_ptr;
      tlist_ptr                  = tlist_ptr->next;
      free(temp); /* free up memory */
    }
    else {
      last_head_list_ptr = tlist_ptr;       /* save last head pointer */
      tlist_ptr          = tlist_ptr->next; /* Loop back if not */
    }
  }
}

/*!
  \ingroup Utilities
  \undoc
*/
int ex_get_num_props(int exoid, ex_entity_type obj_type)
{
  int   cntr, varid;
  char *var_name;
  char  errmsg[MAX_ERR_LENGTH];

  EX_FUNC_ENTER();
  cntr = 0;

  /* loop until there is not a property variable defined; the name of */
  /* the variables begin with an increment of 1 ("xx_prop1") so use cntr+1 */
  while (true) {
    switch (obj_type) {
    case EX_ELEM_BLOCK: var_name = VAR_EB_PROP(cntr + 1); break;
    case EX_EDGE_BLOCK: var_name = VAR_ED_PROP(cntr + 1); break;
    case EX_FACE_BLOCK: var_name = VAR_FA_PROP(cntr + 1); break;
    case EX_NODE_SET: var_name = VAR_NS_PROP(cntr + 1); break;
    case EX_EDGE_SET: var_name = VAR_ES_PROP(cntr + 1); break;
    case EX_FACE_SET: var_name = VAR_FS_PROP(cntr + 1); break;
    case EX_SIDE_SET: var_name = VAR_SS_PROP(cntr + 1); break;
    case EX_ELEM_SET: var_name = VAR_ELS_PROP(cntr + 1); break;
    case EX_ASSEMBLY: var_name = VAR_ASSEMBLY_PROP(cntr + 1); break;
    case EX_BLOB: var_name = VAR_BLOB_PROP(cntr + 1); break;
    case EX_ELEM_MAP: var_name = VAR_EM_PROP(cntr + 1); break;
    case EX_FACE_MAP: var_name = VAR_FAM_PROP(cntr + 1); break;
    case EX_EDGE_MAP: var_name = VAR_EDM_PROP(cntr + 1); break;
    case EX_NODE_MAP: var_name = VAR_NM_PROP(cntr + 1); break;
    default:
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: object type %d not supported; file id %d", obj_type,
               exoid);
      ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
      EX_FUNC_LEAVE(EX_FATAL);
    }

    if (nc_inq_varid(exoid, var_name, &varid) != EX_NOERR) {
      /*   no variable with this name; return cntr which is now the number of */
      /*   properties for this type of entity */
      EX_FUNC_LEAVE(cntr);
    }
    cntr++;
  }
  EX_FUNC_LEAVE(EX_FATAL);
}

/*!
  \ingroup Utilities
  \undoc
*/
int exi_get_cpu_ws(void) { return sizeof(float); }

/* swap - interchange v[i] and v[j] */
/*!
  \internal
  \undoc
*/
static void ex_swap(int v[], int64_t i, int64_t j)
{
  /* Thread-safe, reentrant */
  int temp = v[i];
  v[i]     = v[j];
  v[j]     = temp;
}

/*!
  \internal
  \undoc
*/
static void ex_swap64(int64_t v[], int64_t i, int64_t j)
{
  /* Thread-safe, reentrant */
  int64_t temp = v[i];
  v[i]         = v[j];
  v[j]         = temp;
}

#define EX_QSORT_CUTOFF 12

/*!
  \internal
  \undoc
*/
static int ex_int_median3(int v[], int iv[], int64_t left, int64_t right)
{
  /* Thread-safe, reentrant */
  int64_t center = (left + right) / 2;

  if (v[iv[left]] > v[iv[center]]) {
    ex_swap(iv, left, center);
  }
  if (v[iv[left]] > v[iv[right]]) {
    ex_swap(iv, left, right);
  }
  if (v[iv[center]] > v[iv[right]]) {
    ex_swap(iv, center, right);
  }

  ex_swap(iv, center, right - 1);
  return iv[right - 1];
}

/*!
  \internal
  \undoc
*/
static int64_t ex_int_median3_64(int64_t v[], int64_t iv[], int64_t left, int64_t right)
{
  /* Thread-safe, reentrant */
  int64_t center = (left + right) / 2;

  if (v[iv[left]] > v[iv[center]]) {
    ex_swap64(iv, left, center);
  }
  if (v[iv[left]] > v[iv[right]]) {
    ex_swap64(iv, left, right);
  }
  if (v[iv[center]] > v[iv[right]]) {
    ex_swap64(iv, center, right);
  }

  ex_swap64(iv, center, right - 1);
  return iv[right - 1];
}

/*!
  \internal
  \undoc
*/
static void ex_int_iqsort(int v[], int iv[], int left, int right)
{
  /* Thread-safe, reentrant */
  if (left + EX_QSORT_CUTOFF <= right) {
    int pivot = ex_int_median3(v, iv, left, right);
    int i     = left;
    int j     = right - 1;

    for (;;) {
      while (v[iv[++i]] < v[pivot]) {
        ;
      }
      while (v[iv[--j]] > v[pivot]) {
        ;
      }
      if (i < j) {
        ex_swap(iv, i, j);
      }
      else {
        break;
      }
    }

    ex_swap(iv, i, right - 1);
    ex_int_iqsort(v, iv, left, i - 1);
    ex_int_iqsort(v, iv, i + 1, right);
  }
}

/*!
  \internal
  \undoc
*/
static void ex_int_iqsort64(int64_t v[], int64_t iv[], int64_t left, int64_t right)
{
  /* Thread-safe, reentrant */
  if (left + EX_QSORT_CUTOFF <= right) {
    int64_t pivot = ex_int_median3_64(v, iv, left, right);
    int64_t i     = left;
    int64_t j     = right - 1;

    for (;;) {
      while (v[iv[++i]] < v[pivot]) {
        ;
      }
      while (v[iv[--j]] > v[pivot]) {
        ;
      }
      if (i < j) {
        ex_swap64(iv, i, j);
      }
      else {
        break;
      }
    }

    ex_swap64(iv, i, right - 1);
    ex_int_iqsort64(v, iv, left, i - 1);
    ex_int_iqsort64(v, iv, i + 1, right);
  }
}

/*!
  \internal
  \undoc
*/
static void ex_int_iisort(int v[], int iv[], int N)
{
  /* Thread-safe, reentrant */
  int ndx   = 0;
  int small = v[iv[0]];
  for (int i = 1; i < N; i++) {
    if (v[iv[i]] < small) {
      small = v[iv[i]];
      ndx   = i;
    }
  }
  /* Put smallest value in slot 0 */
  ex_swap(iv, 0, ndx);

  for (int i = 1; i < N; i++) {
    int tmp = iv[i];
    int j;
    for (j = i; v[tmp] < v[iv[j - 1]]; j--) {
      iv[j] = iv[j - 1];
    }
    iv[j] = tmp;
  }
}

/*!
  \internal
  \undoc
*/
static void ex_int_iisort64(int64_t v[], int64_t iv[], int64_t N)
{
  /* Thread-safe, reentrant */
  int64_t ndx   = 0;
  int64_t small = v[iv[0]];
  for (int64_t i = 1; i < N; i++) {
    if (v[iv[i]] < small) {
      small = v[iv[i]];
      ndx   = i;
    }
  }
  /* Put smallest value in slot 0 */
  ex_swap64(iv, 0, ndx);

  for (int64_t i = 1; i < N; i++) {
    int64_t tmp = iv[i];
    int64_t j;
    for (j = i; v[tmp] < v[iv[j - 1]]; j--) {
      iv[j] = iv[j - 1];
    }
    iv[j] = tmp;
  }
}

/*!
 * \ingroup Utilities
 * \internal
 * The following 'indexed qsort' routine is modified from Sedgewicks
 * algorithm It selects the pivot based on the median of the left,
 * right, and center values to try to avoid degenerate cases occurring
 * when a single value is chosen.  It performs a quicksort on
 * intervals down to the #EX_QSORT_CUTOFF size and then performs a final
 * insertion sort on the almost sorted final array.  Based on data in
 * Sedgewick, the #EX_QSORT_CUTOFF value should be between 5 and 20.
 *
 * See Sedgewick for further details
 * Define `DEBUG_QSORT` at the top of this file and recompile to compile
 * in code that verifies that the array is sorted.
 *
 * NOTE: The 'int' implementation below assumes that *both* the items
 *       being sorted and the *number* of items being sorted are both
 *       representable as 'int'.
 * \internal
 */
void exi_iqsort(int v[], int iv[], size_t N)
{
  /* Thread-safe, reentrant */
  ex_int_iqsort(v, iv, 0, N - 1);
  ex_int_iisort(v, iv, N);

#if defined(DEBUG_QSORT)
  fprintf(stderr, "Checking sort of %zu values\n", N + 1);
  int i;
  for (i = 1; i < N; i++) {
    assert(v[iv[i - 1]] <= v[iv[i]]);
  }
#endif
}

/*! \sa exi_iqsort() */
void exi_iqsort64(int64_t v[], int64_t iv[], int64_t N)
{
  /* Thread-safe, reentrant */
  ex_int_iqsort64(v, iv, 0, N - 1);
  ex_int_iisort64(v, iv, N);

#if defined(DEBUG_QSORT)
  fprintf(stderr, "Checking sort of %" PRId64 " values\n", N + 1);
  int i;
  for (i = 1; i < N; i++) {
    assert(v[iv[i - 1]] <= v[iv[i]]);
  }
#endif
}

/*!
 * Determine whether the new large model storage is being used in this
 * file, or old method. Basically, the difference is whether the
 * coordinates and nodal variables are stored in a blob (xyz
 * components together) or as a variable per component per
 * nodal_variable.
 *
 * \ingroup Utilities
 */
int ex_large_model(int exoid)
{
  if (exoid < 0) {
    return EXODUS_DEFAULT_SIZE; /* Specified in exodusII_int.h */
  }

  /* See if the ATT_FILESIZE attribute is defined in the file */
  EX_FUNC_ENTER();
  int file_size = 0;
  int rootid    = exoid & EX_FILE_ID_MASK;
  if (nc_get_att_int(rootid, NC_GLOBAL, ATT_FILESIZE, &file_size) != EX_NOERR) {
    /* Variable not found; default is 0 */
    file_size = 0;
  }
  EX_FUNC_LEAVE(file_size);
}

/*!
  \internal
  \undoc
*/
int exi_get_dimension(int exoid, const char *DIMENSION, const char *label, size_t *count,
                      int *dimid, const char *routine)
{
  int status;

  *count = 0;
  *dimid = -1;

  if ((status = nc_inq_dimid(exoid, DIMENSION, dimid)) != EX_NOERR) {
    if (routine != NULL) {
      if (status == NC_EBADDIM) {
        char errmsg[MAX_ERR_LENGTH];
        snprintf(errmsg, MAX_ERR_LENGTH, "Warning: no dimension defining '%s' found in file id %d",
                 label, exoid);
        ex_err_fn(exoid, __func__, errmsg, status);
      }
      else {
        char errmsg[MAX_ERR_LENGTH];
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: failed to locate dimension defining number of '%s' in file id %d", label,
                 exoid);
        ex_err_fn(exoid, __func__, errmsg, status);
      }
    }
    return status;
  }

  if ((status = nc_inq_dimlen(exoid, *dimid, count)) != EX_NOERR) {
    if (routine != NULL) {
      char errmsg[MAX_ERR_LENGTH];
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to get length of dimension defining number of '%s' in file id %d",
               label, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
    }
    return status;
  }
  return status;
}

/*!
  \deprecated
*/
size_t ex_header_size(int exoid)
{
  EX_UNUSED(exoid);
  return 0;
}

void exi_set_compact_storage(int exoid, int varid)
{
  /* Capability was released in version 4.7.4
     Only applicable to netcdf-4 files, but will
     succeed on other files; just won't do anything
  */
#if defined(NC_COMPACT)
  nc_def_var_chunking(exoid, varid, NC_COMPACT, NULL);
#else
  EX_UNUSED(exoid);
  EX_UNUSED(varid);
#endif
}

/*
 * type = 1 for integer, 2 for real, 3 for character
 * If type < 0, then don't compress, but do set collective on parallel
 */

/*!
  \internal
  \undoc
*/

void exi_compress_variable(int exoid, int varid, int type)
{
#if NC_HAS_HDF5
  struct exi_file_item *file = exi_find_file_item(exoid);

  if (!file) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: unknown file id %d for exi_compress_variable().",
             exoid);
    ex_err_fn(exoid, __func__, errmsg, EX_BADFILEID);
  }
  else {
    /* Compression only supported on HDF5 (NetCDF-4) files; Do not try to compress character data */
    int status = EX_NOERR;
    if ((type == 1 || type == 2) && file->is_hdf5) {
      if (file->compression_algorithm == EX_COMPRESS_GZIP) {
        int deflate_level = file->compression_level;
        if (deflate_level > 0) {
          int compress = 1;
          int shuffle  = file->shuffle;
          status       = nc_def_var_deflate(exoid, varid, shuffle, compress, deflate_level);
        }
      }
      else if (file->compression_algorithm == EX_COMPRESS_SZIP) {
#if NC_HAS_SZIP_WRITE == 1
        /* See: https://support.hdfgroup.org/doc_resource/SZIP/ and
                https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-SetSzip
           for details on SZIP library and parameters.
        */

#if !defined(NC_SZIP_EC)
        const int NC_SZIP_EC = 4; /* Selects entropy coding method for szip. */
#endif
#if !defined(NC_SZIP_NN)
        const int NC_SZIP_NN = 32; /* Selects nearest neighbor coding method for szip. */
#endif
        /* Even and between 4 and 32; typical values are 8, 10, 16, 32 */
        const int SZIP_PIXELS_PER_BLOCK =
            file->compression_level == 0 ? 32 : file->compression_level;
        status = nc_def_var_szip(exoid, varid, NC_SZIP_NN, SZIP_PIXELS_PER_BLOCK);
#else
        char errmsg[MAX_ERR_LENGTH];
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: Compression algorithm SZIP is not supported yet (EXPERIMENTAL).");
        ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
#endif
      }
      else if (file->compression_algorithm == EX_COMPRESS_ZSTD) {
#if NC_HAS_ZSTD == 1
        status = nc_def_var_zstandard(exoid, varid, file->compression_level);
#else
        char errmsg[MAX_ERR_LENGTH];
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: Compression algorithm ZSTANDARD is not supported in this version of the "
                 "netCDF library.");
        ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
#endif
      }
      else if (file->compression_algorithm == EX_COMPRESS_BZ2) {
#if NC_HAS_BZ2 == 1
        status = nc_def_var_bzip2(exoid, varid, file->compression_level);
#else
        char errmsg[MAX_ERR_LENGTH];
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: Compression algorithm BZIP2 / BZ2 is not supported in this version of the "
                 "netCDF library.");
        ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
#endif
      }
      if (status != EX_NOERR) {
        char errmsg[MAX_ERR_LENGTH];
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: failed to set compression attribute on variable in file id %d", exoid);
        ex_err_fn(exoid, __func__, errmsg, status);
      }

      if (type == 2 && file->quantize_nsd > 0) {
#if NC_HAS_QUANTIZE == 1
        // Lossy compression using netCDF quantize methods.
        if ((status = nc_def_var_quantize(exoid, varid, NC_QUANTIZE_GRANULARBR,
                                          file->quantize_nsd)) != EX_NOERR) {
          char errmsg[MAX_ERR_LENGTH];
          snprintf(errmsg, MAX_ERR_LENGTH,
                   "ERROR: failed to set quanitzation method on variable in file id %d", exoid);
          ex_err_fn(exoid, __func__, errmsg, status);
        }
#else
        char errmsg[MAX_ERR_LENGTH];
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: Quanitzation is not supported in this version of netCDF library.");
        ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
#endif
      }
    }
#if defined(PARALLEL_AWARE_EXODUS)
    if (file->is_parallel) {
      nc_var_par_access(exoid, varid, NC_COLLECTIVE);
    }
#endif
  }
#endif
}

/*!
  \internal
  \undoc
*/
int exi_leavedef(int exoid, const char *call_func)
{
  int status;

  struct exi_file_item *file = exi_find_file_item(exoid);
  if (!file) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: unknown file id %d for %s called from %s.", exoid,
             __func__, call_func);
    ex_err_fn(exoid, __func__, errmsg, EX_BADFILEID);
    return EX_FATAL;
  }

  if (!file->persist_define_mode && file->in_define_mode) {
    if ((status = nc_enddef(exoid)) != EX_NOERR) {
      char errmsg[MAX_ERR_LENGTH];
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to complete definition for file id %d",
               exoid);
      ex_err_fn(exoid, call_func, errmsg, status);

      return EX_FATAL;
    }
    file->in_define_mode = 0;
  }
  return EX_NOERR;
}

int exi_redef(int exoid, const char *call_func)
{
  int status;

  struct exi_file_item *file = exi_find_file_item(exoid);

  if (!file) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: unknown file id %d for %s called from %s.", exoid,
             __func__, call_func);
    ex_err_fn(exoid, __func__, errmsg, EX_BADFILEID);
    return EX_FATAL;
  }

  if (!file->in_define_mode) {
    if ((status = nc_redef(exoid)) != EX_NOERR) {
      char errmsg[MAX_ERR_LENGTH];
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to put file %d into definition mode in exi_redef called from %s",
               exoid, call_func);
      ex_err_fn(exoid, __func__, errmsg, status);
      return EX_FATAL;
    }
    file->in_define_mode = 1;
  }
  return EX_NOERR;
}

int exi_persist_redef(int exoid, const char *call_func)
{
  int status;

  struct exi_file_item *file = exi_find_file_item(exoid);

  if (!file) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: unknown file id %d for %s called from %s.", exoid,
             __func__, call_func);
    ex_err_fn(exoid, __func__, errmsg, EX_BADFILEID);
    return EX_FATAL;
  }

  if ((++file->persist_define_mode == 1) && !file->in_define_mode) {
    if ((status = nc_redef(exoid)) != EX_NOERR) {
      char errmsg[MAX_ERR_LENGTH];
      snprintf(
          errmsg, MAX_ERR_LENGTH,
          "ERROR: failed to put file %d into definition mode in exi_persist_redef called from %s",
          exoid, call_func);
      ex_err_fn(exoid, __func__, errmsg, status);
      return EX_FATAL;
    }
    file->in_define_mode = 1;
  }
  return EX_NOERR;
}

int exi_persist_leavedef(int exoid, const char *call_func)
{
  int status;

  struct exi_file_item *file = exi_find_file_item(exoid);
  if (!file) {
    char errmsg[MAX_ERR_LENGTH];
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: unknown file id %d for %s called from %s.", exoid,
             __func__, call_func);
    ex_err_fn(exoid, __func__, errmsg, EX_BADFILEID);
    return EX_FATAL;
  }

  if ((file->persist_define_mode-- == 1) && file->in_define_mode) {
    if ((status = nc_enddef(exoid)) != EX_NOERR) {
      char errmsg[MAX_ERR_LENGTH];
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to complete definition for file id %d",
               exoid);
      ex_err_fn(exoid, call_func, errmsg, status);

      return EX_FATAL;
    }
    file->in_define_mode      = 0;
    file->persist_define_mode = 0;
  }
  return EX_NOERR;
}

static int warning_output = 0;

int exi_check_version(int run_version)
{
  if (run_version != EX_API_VERS_NODOT && warning_output == 0) {
    int run_version_major = run_version / 100;
    int run_version_minor = run_version % 100;
    int lib_version_major = EXODUS_VERSION_MAJOR;
    int lib_version_minor = EXODUS_VERSION_MINOR;
    fprintf(stderr,
            "EXODUS: Warning: This code was compiled with exodus "
            "version %d.%02d,\n          but was linked with exodus "
            "library version %d.%02d\n          This is probably an "
            "error in the build process of this code.\n",
            run_version_major, run_version_minor, lib_version_major, lib_version_minor);
    warning_output = 1;
  }
  return warning_output;
}

/*!
  \internal
  \undoc
*/
int exi_handle_mode(unsigned int my_mode, int is_parallel, int run_version)
{
  char       errmsg[MAX_ERR_LENGTH];
  int        nc_mode      = 0;
  static int netcdf4_mode = -1;
#if NC_HAS_CDF5
  static int netcdf5_mode = -1;
#endif

  int filesiz = 1;
  int int64_status;
  int pariomode = 0;

  /* Contains a 1 in all bits corresponding to file modes */
  /* Do not include EX_64BIT_DATA in this list */
  static unsigned int all_modes = EX_NORMAL_MODEL | EX_64BIT_OFFSET | EX_NETCDF4 | EX_PNETCDF;

  exi_check_version(run_version);

/*
 * See if specified mode is supported in the version of netcdf we
 * are using
 */
#if !NC_HAS_HDF5
  if (my_mode & EX_NETCDF4) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "EXODUS: ERROR: File format specified as netcdf-4, but the "
             "NetCDF library being used was not configured to enable "
             "this format\n");
    ex_err(__func__, errmsg, EX_BADPARAM);
    EX_FUNC_LEAVE(EX_FATAL);
  }
#endif

#if !NC_HAS_CDF5
  if (my_mode & EX_64BIT_DATA) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "EXODUS: ERROR: File format specified as 64bit_data, but "
             "the NetCDF library being used does not support this "
             "format\n");
    ex_err(__func__, errmsg, EX_BADPARAM);
    EX_FUNC_LEAVE(EX_FATAL);
  }
#endif

  /* EX_64_BIT_DATA is 64-bit integer version of EX_PNETCDF.  If
     EX_64_BIT_DATA and EX_PNETCDF is not set, then set EX_PNETCDF... */
  if (my_mode & EX_64BIT_DATA) {
    my_mode |= EX_PNETCDF;
  }

  /* Check that one and only one format mode is specified... */
  {
    unsigned int set_modes = all_modes & my_mode;

    if (set_modes == 0) {
      my_mode |= EX_64BIT_OFFSET; /* Default if nothing specified */
    }
    else {
      /* Checks that only a single bit is set */
      set_modes = !(set_modes & (set_modes - 1));
      if (!set_modes) {
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "EXODUS: ERROR: More than 1 file format "
                 "(EX_NORMAL_MODEL, EX_LARGE_MODEL, EX_64BIT_OFFSET, "
                 "or EX_NETCDF4)\nwas specified in the "
                 "mode argument of the ex_create call. Only a single "
                 "format can be specified.\n");
        ex_err(__func__, errmsg, EX_BADPARAM);
        EX_FUNC_LEAVE(EX_FATAL);
      }
    }
  }

  /*
   * See if any integer data is to be stored as int64 (long long). If
   * so, then need to set NC_NETCDF4 and unset NC_CLASSIC_MODEL (or
   * set EX_NOCLASSIC.  Output meaningful error message if the library
   * is not NetCDF-4 enabled...
   *
   * As of netcdf-4.4.0, can also use NC_64BIT_DATA (CDF5) mode for this...
   */
  int64_status = my_mode & (EX_ALL_INT64_DB | EX_ALL_INT64_API);

  if ((int64_status & EX_ALL_INT64_DB) != 0) {
#if NC_HAS_HDF5 || NC_HAS_CDF5
    /* Library DOES support netcdf4 and/or cdf5 ... See if user
     * specified either of these and use that one; if not, pick
     * netcdf4, non-classic as default.
     */
    if (my_mode & EX_NETCDF4) {
      my_mode |= EX_NOCLASSIC;
    }
#if NC_HAS_CDF5
    else if (my_mode & EX_64BIT_DATA) {
      ; /* Do nothing, already set */
    }
    else if (my_mode & EX_PNETCDF) {
      my_mode |= EX_64BIT_DATA;
    }
#endif
    else {
      /* Unset the current mode so we don't have multiples specified */
      /* ~all_modes sets to 1 all bits not associated with file format */
      my_mode &= ~all_modes;
#if NC_HAS_HDF5
      /* Pick netcdf4 as default mode for 64-bit integers */
      my_mode |= EX_NOCLASSIC;
      my_mode |= EX_NETCDF4;
#else
      /* Pick 64bit_data as default mode for 64-bit integers */
      my_mode |= EX_64BIT_DATA;
#endif
    }
#else
    /* Library does NOT support netcdf4 or cdf5 */
    snprintf(errmsg, MAX_ERR_LENGTH,
             "EXODUS: ERROR: 64-bit integer storage requested, but the "
             "netcdf library does not support the required netcdf-4 or "
             "64BIT_DATA extensions.\n");
    ex_err(__func__, errmsg, EX_BADPARAM);
    EX_FUNC_LEAVE(EX_FATAL);
#endif
  }

#if defined(PARALLEL_AWARE_EXODUS)
  /* Check parallel io mode.  Valid is NC_MPIPOSIX or NC_MPIIO or NC_PNETCDF
   * Exodus uses different flag values; map to netcdf values
   *
   * NOTE: In current versions of NetCDF, MPIPOSIX and MPIIO are ignored and the
   *       underlying format is either NC_PNETCDF or NC_NETCDF4 (hdf5-based)
   *       They map NC_MPIIO to NC_PNETCDF, but in the past, exodus mapped EX_MPIIO
   *       to EX_NETCDF4.
   */
  if (is_parallel) {
    int tmp_mode = 0;
    if (my_mode & EX_MPIPOSIX) {
      pariomode = NC_MPIIO;
      tmp_mode  = EX_NETCDF4;
#if !NC_HAS_HDF5
      snprintf(errmsg, MAX_ERR_LENGTH,
               "EXODUS: ERROR: EX_MPIPOSIX parallel output requested "
               "which requires NetCDF-4 support, but the library does "
               "not have that option enabled.\n");
      ex_err(__func__, errmsg, EX_BADPARAM);
      EX_FUNC_LEAVE(EX_FATAL);
#endif
    }
    else if (my_mode & EX_MPIIO) {
      pariomode = NC_MPIIO;
      tmp_mode  = EX_NETCDF4;
#if !NC_HAS_HDF5
      snprintf(errmsg, MAX_ERR_LENGTH,
               "EXODUS: ERROR: EX_MPIIO parallel output requested which "
               "requires NetCDF-4 support, but the library does not "
               "have that option enabled.\n");
      ex_err(__func__, errmsg, EX_BADPARAM);
      EX_FUNC_LEAVE(EX_FATAL);
#endif
    }
    else if (my_mode & EX_NETCDF4) {
      pariomode = NC_MPIIO;
      tmp_mode  = EX_NETCDF4;
#if !NC_HAS_HDF5
      snprintf(errmsg, MAX_ERR_LENGTH,
               "EXODUS: ERROR: EX_NETCDF4 parallel output requested which "
               "requires NetCDF-4 support, but the library does not "
               "have that option enabled.\n");
      ex_err(__func__, errmsg, EX_BADPARAM);
      EX_FUNC_LEAVE(EX_FATAL);
#endif
    }
    else if (my_mode & EX_PNETCDF) {
      pariomode = NC_PNETCDF;
      /* See if client specified 64-bit or not... */
      if ((my_mode & EX_64BIT_DATA) || (int64_status & EX_ALL_INT64_DB)) {
        tmp_mode = EX_64BIT_DATA;
      }
      else {
        if (my_mode & EX_64BIT_DATA) {
          tmp_mode = EX_64BIT_DATA;
        }
        else {
          tmp_mode = EX_64BIT_OFFSET;
        }
      }
#if !NC_HAS_PNETCDF
      snprintf(errmsg, MAX_ERR_LENGTH,
               "EXODUS: ERROR: EX_PNETCDF parallel output requested "
               "which requires PNetCDF support, but the library does "
               "not have that option enabled.\n");
      ex_err(__func__, errmsg, EX_BADPARAM);
      EX_FUNC_LEAVE(EX_FATAL);
#endif
    }

    /* If tmp_mode was set here, then need to clear any other mode that
       was potentially already set in my_mode... */
    my_mode &= ~all_modes;
    my_mode |= tmp_mode;
  }
#else
  EX_UNUSED(is_parallel);
#endif /* PARALLEL_AWARE_EXODUS */

  if (my_mode & EX_NETCDF4) {
    nc_mode |= NC_NETCDF4;
  }
  else {
    if (netcdf4_mode == -1) {
      char *option = getenv("EXODUS_NETCDF4");
      if (option != NULL) {
        netcdf4_mode = NC_NETCDF4;
        if (option[0] != 'q') {
          fprintf(stderr, "EXODUS: Using netcdf version 4 selected via "
                          "EXODUS_NETCDF4 environment variable\n");
        }
      }
      else {
        netcdf4_mode = 0;
      }
    }
    nc_mode |= netcdf4_mode;
  }

  if (!(my_mode & EX_NOCLASSIC)) {
    nc_mode |= NC_CLASSIC_MODEL;
  }

#if NC_HAS_CDF5
  if (my_mode & EX_64BIT_DATA) {
    nc_mode |= (NC_64BIT_DATA);
  }
  else {
    if (netcdf5_mode == -1) {
      char *option = getenv("EXODUS_NETCDF5");
      if (option != NULL) {
        netcdf5_mode = NC_64BIT_DATA;
        if (option[0] != 'q') {
          fprintf(stderr, "EXODUS: Using netcdf version 5 (CDF5) selected via "
                          "EXODUS_NETCDF5 environment variable\n");
        }
      }
      else {
        netcdf5_mode = 0;
      }
    }
    nc_mode |= netcdf5_mode;
  }
#endif

  /*
   * Hardwire filesiz to 1 for all created files. Reduce complexity in nodal output routines.
   * has been default for a decade or so, but still support it on read...
   */
  if (
#if NC_HAS_HDF5
      !(nc_mode & NC_NETCDF4) &&
#endif
#if NC_HAS_CDF5
      !(nc_mode & NC_64BIT_DATA) &&
#endif
      filesiz == 1) {
    nc_mode |= NC_64BIT_OFFSET;
  }

  if (my_mode & EX_SHARE) {
    nc_mode |= NC_SHARE;
  }

  /*
   * set error handling mode to no messages, non-fatal errors
   * unless specified differently via environment.
   */
  {
    if (exoptval != EX_VERBOSE) {
      /* Avoid getenv call if already in verbose mode */
      char *option = getenv("EXODUS_VERBOSE");
      if (option != NULL) {
        exoptval = EX_VERBOSE;
        if (option[0] != 'q') {
          fprintf(stderr, "EXODUS: Setting EX_VERBOSE mode since EXODUS_VERBOSE environment "
                          "variable is set.\n");
        }
      }
    }
    ex_opts(exoptval); /* call required to set ncopts first time through */
  }

  if (my_mode & EX_CLOBBER) {
    nc_mode |= NC_CLOBBER;
  }
  else {
    nc_mode |= NC_NOCLOBBER;
  }

#if NC_HAS_DISKLESS
  /* Use of diskless (in-memory) and parallel is not tested... */
  if (my_mode & EX_DISKLESS) {
    nc_mode |= NC_DISKLESS;
    nc_mode |= NC_WRITE;
#if defined NC_PERSIST
    nc_mode |= NC_PERSIST;
#endif
  }
#endif
  return nc_mode | pariomode;
}

/*!
  \internal
  \undoc
*/
int exi_populate_header(int exoid, const char *path, int my_mode, int is_parallel, int *comp_ws,
                        int *io_ws)
{
  int  status;
  int  old_fill;
  int  lio_ws;
  int  filesiz    = 1;
  bool is_hdf5    = false;
  bool is_pnetcdf = false;

  float version;
  char  errmsg[MAX_ERR_LENGTH];
  int   int64_status = my_mode & (EX_ALL_INT64_DB | EX_ALL_INT64_API);

  int format = 0;
  int mode;

  /* turn off automatic filling of netCDF variables */
  if ((status = nc_set_fill(exoid, NC_NOFILL, &old_fill)) != EX_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to set nofill mode in file id %d", exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    return EX_FATAL;
  }

  /* Verify that there is not an existing file_item struct for this
     exoid This could happen (and has) when application calls
     ex_open(), but then closes file using nc_close() and then reopens
     file.  NetCDF will possibly reuse the exoid which results in
     internal corruption in exodus data structures since exodus does
     not know that file was closed and possibly new file opened for
     this exoid
  */
  if (exi_find_file_item(exoid) != NULL) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: There is an existing file already using the file "
             "id %d which was also assigned to file %s.\n\tWas "
             "nc_close() called instead of ex_close() on an open Exodus "
             "file?\n",
             exoid, path);
    ex_err_fn(exoid, __func__, errmsg, EX_BADFILEID);
    nc_close(exoid);
    return EX_FATAL;
  }

  /* initialize floating point size conversion.  since creating new file,
   * i/o wordsize attribute from file is zero.
   */
  /* Determine format being used for underlying NetCDF file */
  nc_inq_format_extended(exoid, &format, &mode);

  if (format & NC_FORMAT_PNETCDF) {
    is_pnetcdf = true;
    ;
  }

  if (format & NC_FORMAT_NC_HDF5) {
    is_hdf5 = true;
  }

  if (exi_conv_init(exoid, comp_ws, io_ws, 0, int64_status, is_parallel, is_hdf5, is_pnetcdf,
                    my_mode & EX_WRITE) != EX_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to init conversion routines in file id %d",
             exoid);
    ex_err_fn(exoid, __func__, errmsg, EX_LASTERR);
    return EX_FATAL;
  }

  /* put the EXODUS version number, and i/o floating point word size as
   * netcdf global attributes
   */

  /* store Exodus API version # as an attribute */
  {
    float version_major = EXODUS_VERSION_MAJOR;
    float version_minor = EXODUS_VERSION_MINOR;
    version             = version_major + version_minor / 100.0;
    if ((status = nc_put_att_float(exoid, NC_GLOBAL, ATT_API_VERSION, NC_FLOAT, 1, &version)) !=
        EX_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to store Exodus II API version attribute in file id %d", exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      return EX_FATAL;
    }
  }

  /* store Exodus file version # as an attribute */
  {
    float version_major = EXODUS_VERSION_MAJOR;
    float version_minor = EXODUS_VERSION_MINOR;
    version             = version_major + version_minor / 100.0;
    if ((status = nc_put_att_float(exoid, NC_GLOBAL, ATT_VERSION, NC_FLOAT, 1, &version)) !=
        EX_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to store Exodus II file version attribute in file id %d", exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      return EX_FATAL;
    }
  }

  /* store Exodus file float word size  as an attribute */
  lio_ws = (*io_ws);
  if ((status = nc_put_att_int(exoid, NC_GLOBAL, ATT_FLT_WORDSIZE, NC_INT, 1, &lio_ws)) !=
      EX_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to store Exodus II file float word size "
             "attribute in file id %d",
             exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    return EX_FATAL;
  }

  /* store Exodus file size (1=large, 0=normal) as an attribute */
  if ((status = nc_put_att_int(exoid, NC_GLOBAL, ATT_FILESIZE, NC_INT, 1, &filesiz)) != EX_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to store Exodus II file size attribute in file id %d", exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    return EX_FATAL;
  }

  {
    int max_so_far = 32;
    if ((status = nc_put_att_int(exoid, NC_GLOBAL, ATT_MAX_NAME_LENGTH, NC_INT, 1, &max_so_far)) !=
        EX_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to add maximum_name_length attribute in file id %d", exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      return EX_FATAL;
    }
  }

  {
    int int64_db_status = int64_status & EX_ALL_INT64_DB;
    if ((status = nc_put_att_int(exoid, NC_GLOBAL, ATT_INT64_STATUS, NC_INT, 1,
                                 &int64_db_status)) != EX_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to add int64_status attribute in file id %d",
               exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      return EX_FATAL;
    }
  }

#if 0
  /* Testing to see if can eliminate some nc_enddef movement of vars/recs */
  if ((status = nc__enddef(exoid, 10000, 4, 10000, 4)) != EX_NOERR) {
#else
  if ((status = nc_enddef(exoid)) != EX_NOERR) {
#endif
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to complete definition for file id %d", exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    return EX_FATAL;
  }
  return EX_NOERR;
}

/*!
  \internal
  \undoc
  Safer than strncpy -- guarantees null termination
*/
char *ex_copy_string(char *dest, char const *source, size_t elements)
{
  char *d;
  for (d = dest; d + 1 < dest + elements && *source; d++, source++) {
    *d = *source;
  }
  *d = '\0';
  return d;
}

/* Used by the entity attribute code `ex_get_attribute(s)` and `ex_put_attribute(s)` */
int exi_get_varid(int exoid, ex_entity_type obj_type, ex_entity_id id)
{
  if (exi_check_valid_file_id(exoid, __func__) == EX_FATAL) {
    EX_FUNC_LEAVE(EX_FATAL);
  }

  if (obj_type == EX_GLOBAL) {
    return NC_GLOBAL;
  }

  int  status = 0;
  char errmsg[MAX_ERR_LENGTH];

  if (obj_type == EX_NODAL) {
    /* For the nodal entity attributes, we store it on the
    "coor_names" variable.  Not exactly logical, but it exists in any
    model that has nodes, so it at least gives us a place to put the
    attribute. Another possibility would be the nodal x-coordinate
    variable...
    */
    int varid = 0;
    if ((status = nc_inq_varid(exoid, VAR_NAME_COOR, &varid)) != EX_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to locate node block variable id in file id %d", exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      return EX_FATAL;
    }
    return varid;
  }

  /* Everything else ... */
  /* First, locate index of this objects id `obj_type` id array */
  /* First, locate index of this objects id `obj_type` id array */
  int id_ndx = exi_id_lkup(exoid, obj_type, id);
  if (id_ndx <= 0) {
    ex_get_err(NULL, NULL, &status);
    if (status != 0) {
      if (status == EX_NULLENTITY) { /* NULL object?    */
        return EX_NOERR;
      }
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to locate %s id  %" PRId64 " in id array in file id %d",
               ex_name_of_object(obj_type), id, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      return EX_FATAL;
    }
  }

  const char *entryptr = NULL;
  switch (obj_type) {
  case EX_ASSEMBLY: entryptr = VAR_ENTITY_ASSEMBLY(id_ndx); break;
  case EX_BLOB: entryptr = VAR_ENTITY_BLOB(id_ndx); break;
  case EX_NODE_SET: entryptr = VAR_NODE_NS(id_ndx); break;
  case EX_EDGE_SET: entryptr = VAR_EDGE_ES(id_ndx); break;
  case EX_FACE_SET: entryptr = VAR_FACE_FS(id_ndx); break;
  case EX_SIDE_SET: entryptr = VAR_ELEM_SS(id_ndx); break;
  case EX_ELEM_SET: entryptr = VAR_ELEM_ELS(id_ndx); break;
  case EX_EDGE_BLOCK: entryptr = VAR_EBCONN(id_ndx); break;
  case EX_FACE_BLOCK: entryptr = VAR_FBCONN(id_ndx); break;
  case EX_ELEM_BLOCK: entryptr = VAR_CONN(id_ndx); break;
  default:
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: object type %d not supported in call to %s", obj_type,
             __func__);
    ex_err(__func__, errmsg, EX_BADPARAM);
    return EX_FATAL;
  }

  int varid = 0;
  if ((status = nc_inq_varid(exoid, entryptr, &varid)) != EX_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to locate entity list array for %s %" PRId64 " in file id %d",
             ex_name_of_object(obj_type), id, exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    return EX_FATAL;
  }
  return varid;
}

/*
 * Code from:
 * https://stackoverflow.com/questions/11034002/how-to-get-absolute-path-of-file-or-directory-that-does-not-exist
 *
 * Return the input path in a canonical form. This is achieved by
 * expanding all symbolic links, resolving references to "." and "..",
 * and removing duplicate "/" characters.
 *
 * If the file exists, its path is canonicalized and returned. If the file,
 * or parts of the containing directory, do not exist, path components are
 * removed from the end until an existing path is found. The remainder of the
 * path is then appended to the canonical form of the existing path,
 * and returned. Consequently, the returned path may not exist. The portion
 * of the path which exists, however, is represented in canonical form.
 *
 * If successful, this function returns a C-string, which needs to be freed by
 * the caller using free().
 *
 * ARGUMENTS:
 *   file_path
 *   File path, whose canonical form to return.
 *
 * RETURNS:
 *   On success, returns the canonical path to the file, which needs to be freed
 *   by the caller.
 *
 *   On failure, returns NULL.
 */
char *exi_canonicalize_filename(char const *file_path)
{
#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
  return _fullpath(NULL, file_path, _MAX_PATH);
#else
  char        *canonical_file_path = NULL;
  unsigned int file_path_len       = strlen(file_path);

  if (file_path_len > 0) {
    canonical_file_path = realpath(file_path, NULL);
    if (canonical_file_path == NULL && errno == ENOENT) {
      // The file was not found. Back up to a segment which exists,
      // and append the remainder of the path to it.
      char *file_path_copy = NULL;
      if (file_path[0] == '/' || (strncmp(file_path, "./", 2) == 0) ||
          (strncmp(file_path, "../", 3) == 0)) {
        // Absolute path, or path starts with "./" or "../"
        file_path_copy = strdup(file_path);
      }
      else {
        // Relative path
        file_path_copy = (char *)malloc(strlen(file_path) + 3);
        strcpy(file_path_copy, "./");
        strcat(file_path_copy, file_path);
      }

      // Remove path components from the end, until an existing path is found
      for (int char_idx = strlen(file_path_copy) - 1; char_idx >= 0 && canonical_file_path == NULL;
           --char_idx) {
        if (file_path_copy[char_idx] == '/') {
          // Remove the slash character
          file_path_copy[char_idx] = '\0';

          canonical_file_path = realpath(file_path_copy, NULL);
          if (canonical_file_path != NULL) {
            // An existing path was found. Append the remainder of the path
            // to a canonical form of the existing path.
            char *combined_file_path = (char *)malloc(strlen(canonical_file_path) +
                                                      strlen(file_path_copy + char_idx + 1) + 2);
            strcpy(combined_file_path, canonical_file_path);
            strcat(combined_file_path, "/");
            strcat(combined_file_path, file_path_copy + char_idx + 1);
            free(canonical_file_path);
            canonical_file_path = combined_file_path;
          }
          else {
            // The path segment does not exist. Replace the slash character
            // and keep trying by removing the previous path component.
            file_path_copy[char_idx] = '/';
          }
        }
      }

      free(file_path_copy);
    }
  }
  return canonical_file_path;
#endif
}
