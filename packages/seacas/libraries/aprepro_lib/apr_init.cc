// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

/***
   NAME
     init
   PURPOSE
     Initialize variables and functions Aprepro
***/
#include "apr_builtin.h"
#include "apr_tokenize.h"
#include "aprepro.h"      // for symrec, Aprepro, etc
#include "init_structs.h" // for array_a_init, array_c_init, etc
#include <array>
#include <string> // for string

namespace SEAMS {
  init arith_0_fncts[] = {
      {"seconds", do_time, "seconds()", "Seconds since epoch (useful for srand())."},
      {nullptr, nullptr, nullptr, nullptr}};

  init_d arith_fncts[] = {
      {"abs", do_fabs, "abs(x)", "Absolute value of x. |x|."},
      {"acos", do_acos, "acos(x)", "Inverse cosine of x, returns radians."},
      {"acosd", do_acosd, "acosd(x)", "Inverse cosine of x, returns degrees."},
      {"acosh", do_acosh, "acosh(x)", "Inverse hyperbolic cosine of x."},
      {"asin", do_asin, "asin(x)", "Inverse sine of x, returns radians."},
      {"asind", do_asind, "asind(x)", "Inverse sine of x, returns degrees."},
      {"asinh", do_asinh, "asinh(x)", "Inverse hyperbolic sine of x."},
      {"atan", do_atan, "atan(x)", "Inverse tangent of x, returns radians."},
      {"atand", do_atand, "atand(x)", "Inverse tangent of x, returns degrees."},
      {"atanh", do_atanh, "atanh(x)", "Inverse hyperbolic tangent of x."},
      {"cbrt", do_cbrt, "cbrt(x)", "Cube root of x. "},
      {"ceil", do_ceil, "ceil(x)", "Smallest integer not less than x."},
      {"cos", do_cos, "cos(x)", "Cosine of x, with x in radians"},
      {"cosd", do_cosd, "cosd(x)", "Cosine of x, with x in degrees"},
      {"cosh", do_cosh, "cosh(x)", "Hyperbolic cosine of x."},
      {"d2r", do_d2r, "d2r(x)", "Degrees to radians."},
      {"erf", do_erf, "erf(x)", "Error Function of x"},
      {"erfc", do_erfc, "erfc(x)", "Complementary Error Function of x"},
      {"exp", do_exp, "exp(x)", "Exponential: e^x"},
      {"expm1", do_expm1, "expm1(x)", "Exponential: Accurate version of e^x - 1.0 for small x"},
      {"floor", do_floor, "floor(x)", "Largest integer not greater than x."},
      {"int", do_int, "int(x), [x]", "Integer part of x truncated toward 0."},
      {"lgamma", do_lgamma, "lgamma(x)", "log(Gamma(x))."},
      {"tgamma", do_tgamma, "tgamma(x)", "Gamma(x)."},
      {"ln", do_log, "ln(x)", "Natural (base e) logarithm of x."},
      {"log", do_log, "log(x)", "Natural (base e) logarithm of x."},
      {"log10", do_log10, "log10(x)", "Base 10 logarithm of x. "},
      {"log1p", do_log1p, "log1p(x)", "log(1+x) "},
      {"nint", do_nint, "nint(x)", "Rounds x to nearest integer. <0.5 down; >=0.5 up."},
      {"r2d", do_r2d, "r2d(x)", "Radians to degrees. "},
      {"sin", do_sin, "sin(x)", "Sine of x, with x in radians. "},
      {"sind", do_sind, "sind(x)", "Sine of x, with x in degrees. "},
      {"sinh", do_sinh, "sinh(x)", "Hyperbolic sine of x "},
      {"srand", do_srand, "srand(seed)",
       "Seed the random number generator with the given integer value. "},
      {"sqrt", do_sqrt, "sqrt(x)", "Square root of x. "},
      {"tan", do_tan, "tan(x)", "Tangent of x, with x in radians. "},
      {"tand", do_tand, "tand(x)", "Tangent of x, with x in radians. "},
      {"tanh", do_tanh, "tanh(x)", "Hyperbolic tangent of x. "},
      {"FtoC", do_FtoC, "FtoC(x)",
       "Convert temperature x from degrees F to degrees C (212F -> 100C)"},
      {"CtoF", do_CtoF, "CtoF(x)",
       "Convert temperature x from degrees C to degrees F (100C -> 212F)"},
      {nullptr, nullptr, nullptr, nullptr}};

  init_a arith_a_fncts[] = {
      {"rows", do_rows, "rows(array)", "Returns the number of rows in the array. "},
      {"cols", do_cols, "cols(array)", "Returns the number of columns in the array. "},
      {nullptr, nullptr, nullptr, nullptr}};

  init_dd arith_dd_fncts[] = {
      {"atan2", do_atan2, "atan2(y,x)", "Inverse tangent of y/x, returns radians."},
      {"atan2d", do_atan2d, "atan2d(y,x)", "Inverse tangent of y/x, returns degrees."},
      {"dim", do_dim, "dim(x,y)", "x - min(x,y)"},
      {"fmod", do_fmod, "fmod(x,y)", "Floating-point remainder of x/y."},
      {"hypot", do_hypot, "hypot(x,y)", "sqrt(x^2+y^2)."},
      {"max", do_max, "max(x,y)", "Maximum of x and y. "},
      {"min", do_min, "min(x,y)", "Minimum of x and y. "},
      {"polarX", do_polarX, "polarX(r,a)", "r * cos(a), a is in degrees "},
      {"polarY", do_polarY, "polarY(r,a)", "r * sin(a), a is in degrees "},
      {"pow", do_pow, "pow(x,y)", "x^y "},
      {"rand", do_rand, "rand(xl,xh)", "Random value between xl and xh; uniformly distributed. "},
      {"rand_normal", do_rand_normal, "rand_normal(m,s)",
       "Random value normally distributed with mean m and stddev s."},
      {"rand_lognormal", do_rand_lognormal, "rand_lognormal(m,s)",
       "Random value with lognormal distribution with mean m and stddev s."},
      {"rand_weibull", do_rand_weibull, "rand_weibull(a, b)",
       "Random value with weibull distribution with alpha=a and beta=b. "},
      {"sign", do_sign, "sign(x,y)", "x * sgn(y)"},
      {nullptr, nullptr, nullptr, nullptr}};

  init_dddd arith_dddd_fncts[] = {
      {"Vangle", do_angle, "Vangle(x1,y1,x2,y2)",
       "Angle (radians) between vector x1_i+y1_j and x2_i+y2_j."},
      {"Vangled", do_angled, "Vangled(x1,y1,x2,y2)",
       "Angle (degrees) between vector x1_i+y1_j and x2_i+y2_j."},
      {"dist", do_dist, "dist(x1,y1, x2,y2)", "sqrt((x1-x2)^2 + (y1-y2)^2)"},
      {nullptr, nullptr, nullptr, nullptr}};

  init_cc arith_cc_fncts[] = {{"word_count", do_word_count, "word_count(svar,del)",
                               "Number of words in svar. Words are separated by one or more of the "
                               "characters\n\t\t\tin the "
                               "string variable 'del'."},
                              {nullptr, nullptr, nullptr, nullptr}};

  init_ccc arith_ccc_fncts[] = {
      {"find_word", do_find_word, "find_word(w,s,d)",
       "Find the 1-based index of word 'w' in variable 's'. Words are separated "
       "by one or more of the\n\t\t\tcharacters in the "
       "string variable 'd'. Returns 0 if not found."},
      {nullptr, nullptr, nullptr, nullptr}};

  init_c arith_c_fncts[] = {{"strtod", do_strtod, "strtod(svar)",
                             "Returns a double-precision floating-point number "
                             "equal to the value represented by the\n\t\t\tcharacter "
                             "string pointed to by svar."},
                            {nullptr, nullptr, nullptr, nullptr}};

  init_cd arith_cd_fncts[] = {{"option", do_option, "option(?,?)", "Internal"},
                              {nullptr, nullptr, nullptr, nullptr}};

  init_ddd arith_ddd_fncts[] = {
      {"julday", do_julday, "julday(mm, dd, yy)", "Julian day corresponding to mm/dd/yy. "},
      {nullptr, nullptr, nullptr, nullptr}};

  init_dddddd arith_dddddd_fncts[] = {{"juldayhms", do_juldayhms, "juldayhms(m,d,y,h,m,s)",
                                       "Julian day corresponding to m/d/y at h:m:s "},
                                      {nullptr, nullptr, nullptr, nullptr}};

  str_init string_fncts[] = {
      {"DUMP", do_dumpsym, "DUMP()",
       "Output a list of all user-defined variables and their value."},
      {"DUMP_JSON", do_dumpsym_json, "DUMP_JSON()",
       "Output a list of all user-defined variables and their value in JSON format."},
      {"DUMP_FUNC", do_dumpfunc, "DUMP_FUNC()",
       "Output a list of all double and string functions recognized by aprepro."},
      {"DUMP_PREVAR", do_dumpvar, "DUMP_PREVAR()",
       "Output a list of all predefined variables and their value."},
      {"get_date", do_get_date, "get_date()",
       "Returns a string representing the current date in the form YYYY/MM/DD."},
      {"get_iso_date", do_get_iso_date, "get_iso_date()",
       "Returns a string representing the current date in the form YYYYMMDD."},
      {"get_time", do_get_time, "get_time()",
       "Returns a string representing the current time in the form HH:MM:SS."},
      {"get_temp_filename", do_get_temp_filename, "get_temp_filename()",
       "Returns a string which can be used for a temporary filename without conflicting with any "
       "other filenames."},
      {"version", do_version, "version()",
       "Return the version string (See also _VERSION variable)"},
      {nullptr, nullptr, nullptr, nullptr}};

  str_c_init string_c_fncts[] = {
      {"DUMP", do_dumpsym1, "DUMP(str)",
       "Output a list of all defined variables and their value if name contains 'str'."},
      {"DUMP_FUNC", do_dumpfunc1, "DUMP_FUNC(str)",
       "Output a list of all double and string functions recognized by aprepro if name contains "
       "'str'."},
      {"DUMP_PREVAR", do_dumpvar1, "DUMP_PREVAR()",
       "Output a list of all predefined variables and their value if name contains 'str'."},
      {"tolower", do_tolower, "tolower(svar)",
       "Translates all uppercase characters in svar to "
       "lowercase. It modifies svar and returns the "
       "resulting string.  "},
      {"toupper", do_toupper, "toupper(svar)",
       "Translates all lowercase character in svar to "
       "uppercase. It modifies svar and returns the "
       "resulting string. "},
      {"to_lower", do_tolower, "to_lower(svar)",
       "Translates all uppercase characters in svar to "
       "lowercase. It modifies svar and returns the "
       "resulting string.  "},
      {"to_upper", do_toupper, "to_upper(svar)",
       "Translates all lowercase character in svar to "
       "uppercase. It modifies svar and returns the "
       "resulting string. "},
      {"import", do_import, "import(svar)",
       "Include/import the file pointed to by the string variable or expression 'svar'"},
      {"getenv", do_getenv, "getenv(svar)",
       "Returns a string containing the value of the environment variable svar. If the environment "
       "\n\t\t\tvariable is not defined, an empty string is returned. "},
      {"file_to_string", do_file_to_string, "file_to_string(fn)",
       "Opens the file specified by fn and returns the contents as a multi-line string."},
      {"error", do_error, "error(svar)",
       "Outputs the string svar to stderr and then terminates the code with an error exit status."},
      {"execute", do_execute, "execute(svar)",
       "svar is parsed and executed as if it were a line read from the input file."},
      {"output", do_output, "output(filename)",
       "Creates the file specified by filename and sends "
       "\n\t\t\tall subsequent output from aprepro to that file."},
      {"output_append", do_append, "output_append(fn)",
       "If file with name fn exists, append output to it; otherwise create the file and "
       "send\n\t\t\tall "
       "subsequent output from aprepro to that file."},
      {"rescan", do_rescan, "rescan(svar)",
       "The difference between execute(sv1) and rescan(sv2) "
       "is that sv1 must be a valid expression,\n\t\t\tbut sv2 can "
       "contain zero or more expressions. "},
      {"include_path", do_include_path, "include_path(path)",
       "Specify an optional path to be prepended to a filename when opening a file.\n\t\t\tCan "
       "also be "
       "specified via the -I command line option when executing aprepro."},
      {"Units", do_Units, "Units(svar)",
       "See manual. svar is one of the defined units "
       "systems:\n\t\t\t'si', 'cgs', 'cgs-ev', 'shock', 'swap', "
       "'ft-lbf-s', 'ft-lbm-s', 'in-lbf-s'"},
      {"delete", do_delete, "delete(var_name)", "Delete the variable with name 'var_name'."},
      {"if", do_str_if, "if(x)",
       "Handles the if statements. x can be any valid expression; nonzero is true"},
      {"If", do_str_if, "If(x)",
       "Handles the if statements. x can be any valid expression; nonzero is true"},
      {"elseif", do_str_elseif, "elseif(x)",
       "Handles the if statements. x can be any valid expression; nonzero is true"},
      {"Elseif", do_str_elseif, "Elseif(x)",
       "Handles the if statements. x can be any valid expression; nonzero is true"},
      {"_ifdef", do_str_if, "ifdef(x)",
       "Handles the if statements. x can be any valid expression; "
       "nonzero is true (deprecated, use if)"},
      {"_ifndef", do_str_notif, "ifndef(x)",
       "Handles the if statements. x can be any valid "
       "expression; nonzero is true (deprecated, use if)"},
#if defined(EXODUS_SUPPORT)
      {"exodus_meta", do_exodus_meta, "exodus_meta(filename)",
       "Creates several variables and arrays related to the exodus metadata in the specified "
       "file. "},
#endif
      {nullptr, nullptr, nullptr, nullptr}};

  str_d_init string_d_fncts[] = {
      {"IO", do_intout, "IO(x)", "Convert x to an integer and then to a string. "},
      {"to_string", do_tostring, "to_string(x)",
       "Returns a string representation of the numerical variable x. The variable x is unchanged."},
      {"tostring", do_tostring, "tostring(x)",
       "Returns a string representation of the numerical variable x. The variable x is unchanged."},
      {"if", do_if, "if(x)",
       "Handles the if statements. x can be any valid expression; nonzero is true"},
      {"If", do_if, "If(x)",
       "Handles the if statements. x can be any valid expression; nonzero is true"},
      {"elseif", do_elseif, "elseif(x)",
       "Handles the if statements. x can be any valid expression; nonzero is true"},
      {"Elseif", do_elseif, "Elseif(x)",
       "Handles the if statements. x can be any valid expression; nonzero is true"},
      {"_ifdef", do_if, "ifdef(x)",
       "Handles the if statements. x can be any valid expression; "
       "nonzero is true.\n\t\t\tEats leading whitespace. (deprecated, use "
       "if)"},
      {"_ifndef", do_notif, "ifndef(x)",
       "Handles the if statements. x can be any valid "
       "expression; nonzero is true.\n\t\t\tEats leading whitespace. "
       "(deprecated, use if)"},
      {"switch", do_switch, "switch(x)",
       "Switch statement. Select from the following case "
       "statements which matches 'x' and execute that one.\n\t\t\tEnd "
       "with endswitch"},
      {"Switch", do_switch, "Switch(x)",
       "Switch statement. Select from the following case "
       "statements which matches 'x' and execute that one.\n\t\t\tEnd "
       "with endswitch"},
      {"case", do_case, "case(x)",
       "Switch statement. A case used in a containing switch statement."},
      {"Case", do_case, "Case(x)",
       "Switch statement. A case used in a containing switch statement."},
      {nullptr, nullptr, nullptr, nullptr}};

  str_dcc_init string_dcc_fncts[] = {
      {"get_word", do_get_word, "get_word(n,svar,del)",
       "Returns a string containing the nth word of svar. The words are separated by one or more "
       "\n\t\t\tcharacters in the string variable del "},
      {nullptr, nullptr, nullptr, nullptr}};

  str_ccc_init string_ccc_fncts[] = {
      {"extract", do_extract, "extract(s, b, e)",
       "Return substring [b,e). 'b' is included; 'e' is not. If 'b' not found, return empty; If "
       "'e' not found,\n\t\t\treturn rest of string. If 'b' empty, start at beginning; if 'e' "
       "empty, "
       "return rest of string."},
#if defined(EXODUS_SUPPORT)
      {"exodus_info", do_exodus_info_range, "exodus_info(filename, beg, end)",
       "Parses the info records starting after 'beg' and ending before 'end'"},
#endif
      {nullptr, nullptr, nullptr, nullptr}};

  str_cc_init string_cc_fncts[] = {
#if defined(EXODUS_SUPPORT)
      {"exodus_info", do_exodus_info, "exodus_info(filename, prefix)",
       "Parses the info records that begin with 'prefix' extracted from the exodus file 'ex_fn'"},
#endif
      {nullptr, nullptr, nullptr, nullptr}};

  str_a_init string_a_fncts[] = {
      {"print_array", do_print_array, "print_array(array)", "Prints the data in the array."},
      {nullptr, nullptr, nullptr, nullptr}};

  array_c_init array_c_fncts[] = {{"csv_array", do_csv_array1, "csv_array(filename)",
                                   "Create a 2D array from the data in a csv file."},
                                  {nullptr, nullptr, nullptr, nullptr}};

  array_cd_init array_cd_fncts[] = {
      {"csv_array", do_csv_array, "csv_array(filename, [skip])",
       "Create a 2D array from the data in a csv file optionally skipping rows."
       " If skip is integer\n\t\t\tskip that many rows; if skip is a character, skip lines "
       "beginning with "
       "that character"},
      {nullptr, nullptr, nullptr, nullptr}};

  array_cc_init array_cc_fncts[] = {
      {"csv_array", do_csv_array2, "csv_array(filename, [skip])",
       "Create a 2D array from the data in a csv file optionally skipping rows."
       " If skip is integer skip that many rows; if skip is a character, skip lines beginning with "
       "that character"},
      {"array_from_string", do_array_from_string, "array_from_string(string, delim)",
       "Create a 1D array from the data in a delimited string."
       " The array double values are\n\t\t\tseparated by one or more of the characters in the "
       "string "
       "variable delim."},
      {"sym_tensor_from_string", do_sym_tensor_from_string, "sym_tensor_from_string(string, delim)",
       "Create a 3x3 symmetric array from the data in a delimited string."
       " The six array values are\n\t\t\tseparated by one or more of the characters in the "
       "string variable delim. Order is xx, yy, zz, xy, yz, xz."},
      {nullptr, nullptr, nullptr, nullptr}};

  array_ddd_init array_ddd_fncts[] = {
      {"linear_array", do_linear_array, "linear_array(init, final, count)",
       "Create a 1D array of 'count' rows. Values linearly spaced from 'init' to 'final'."},
      {"make_array", do_make_array_init, "make_array(rows, cols, init=0)",
       "Create a 2D array of size 'rows' by 'cols' initialized to 'init'. 0 if not specified."},
      {nullptr, nullptr, nullptr, nullptr}};

  array_dd_init array_dd_fncts[] = {
      {"make_array", do_make_array, "make_array(rows, cols)",
       "Create a 2D array of size 'rows' by 'cols' initialized to zero."},
      {nullptr, nullptr, nullptr, nullptr}};

  array_d_init array_d_fncts[] = {
      {"identity", do_identity, "identity(size)",
       "Create a 2D identity array with 'size' rows and columns. Diagonal = 1.0"},
      {nullptr, nullptr, nullptr, nullptr}};

  array_a_init array_a_fncts[] = {
      {"transpose", do_transpose, "transpose(array)", "Return the transpose of input array"},
      {"principal_stress", do_principal, "principal_stress(array)",
       "Calculate principal stresses of symmetric 3x3 stress tensor (array)."},
      {nullptr, nullptr, nullptr, nullptr}};

  // clang-format off
  var_init variables[] = {
      {"DEG",  57.29577951308232087680},  /* 180/pi, degrees per radian */
      {"RAD",   0.01745329251994329576},  /* pi/180, radians per degree */
      {"E",     2.71828182845904523536},  /* e, base of natural log     */
      {"GAMMA", 0.57721566490153286060},  /* euler-mascheroni constant  */
      {"PHI",   1.61803398874989484820},  /* golden ratio               */
      {"TAU",   6.28318530717958623200},  /* 2*PI see Tau Manifesto, http://tauday.com */
      {"PI",    3.14159265358979323846},  /* pi                         */
      {"PI_2",  1.57079632679489661923},  /* pi / 2                      */
      {"SQRT2", 1.41421356237309504880},  /* square root of 2            */
      {"TRUE",  1},
      {"FALSE", 0},
      {nullptr, 0}
  };
  // clang-format on

  svar_init svariables[] = {{"_FORMAT", "%.10g"}, /* Default output format */
                            {"_UNITS_SYSTEM", "none"},
                            {nullptr, nullptr}};
  /* NOTE: The current comment is stored in "_C_"
   *     Since it can be changed by user on command line, we
   *     initialize is differently than the other string variables.
   */

#define internal_init_table(functions, func_type, sym_type)                                        \
  do {                                                                                             \
    for (int i = 0; functions[i].fname != nullptr; i++) {                                          \
      symrec *ptr          = putsym(functions[i].fname, sym_type, true);                           \
      ptr->value.func_type = functions[i].fnct;                                                    \
      ptr->info            = functions[i].description;                                             \
      ptr->syntax          = functions[i].syntax;                                                  \
    }                                                                                              \
  } while (0)

  extern SEAMS::Aprepro *aprepro;

  void Aprepro::init_table(const char *comment)
  {
    // clang-format off
    internal_init_table(arith_0_fncts,      fnctptr,        SYMBOL_TYPE::FUNCTION);
    internal_init_table(arith_fncts,        fnctptr_d,      SYMBOL_TYPE::FUNCTION);
    internal_init_table(arith_dd_fncts,     fnctptr_dd,     SYMBOL_TYPE::FUNCTION);
    internal_init_table(arith_a_fncts,      fnctptr_a,      SYMBOL_TYPE::FUNCTION);
    internal_init_table(arith_dddd_fncts,   fnctptr_dddd,   SYMBOL_TYPE::FUNCTION);
    internal_init_table(arith_ccc_fncts,    fnctptr_ccc,    SYMBOL_TYPE::FUNCTION);
    internal_init_table(arith_cc_fncts,     fnctptr_cc,     SYMBOL_TYPE::FUNCTION);
    internal_init_table(arith_c_fncts,      fnctptr_c,      SYMBOL_TYPE::FUNCTION);
    internal_init_table(arith_cd_fncts,     fnctptr_cd,     SYMBOL_TYPE::FUNCTION);
    internal_init_table(arith_ddd_fncts,    fnctptr_ddd,    SYMBOL_TYPE::FUNCTION);
    internal_init_table(arith_dddddd_fncts, fnctptr_dddddd, SYMBOL_TYPE::FUNCTION);

    internal_init_table(string_fncts,       strfnct,        SYMBOL_TYPE::STRING_FUNCTION);
    internal_init_table(string_c_fncts,     strfnct_c,      SYMBOL_TYPE::STRING_FUNCTION);
    internal_init_table(string_d_fncts,     strfnct_d,      SYMBOL_TYPE::STRING_FUNCTION);
    internal_init_table(string_dcc_fncts,   strfnct_dcc,    SYMBOL_TYPE::STRING_FUNCTION);
    internal_init_table(string_ccc_fncts,   strfnct_ccc,    SYMBOL_TYPE::STRING_FUNCTION);
    internal_init_table(string_cc_fncts,    strfnct_cc,     SYMBOL_TYPE::STRING_FUNCTION);
    internal_init_table(string_a_fncts,     strfnct_a,      SYMBOL_TYPE::STRING_FUNCTION);

    internal_init_table(array_c_fncts,      arrfnct_c,      SYMBOL_TYPE::ARRAY_FUNCTION);
    internal_init_table(array_cc_fncts,     arrfnct_cc,     SYMBOL_TYPE::ARRAY_FUNCTION);
    internal_init_table(array_cd_fncts,     arrfnct_cd,     SYMBOL_TYPE::ARRAY_FUNCTION);
    internal_init_table(array_d_fncts,      arrfnct_d,      SYMBOL_TYPE::ARRAY_FUNCTION);
    internal_init_table(array_dd_fncts,     arrfnct_dd,     SYMBOL_TYPE::ARRAY_FUNCTION);
    internal_init_table(array_ddd_fncts,    arrfnct_ddd,    SYMBOL_TYPE::ARRAY_FUNCTION);
    internal_init_table(array_a_fncts,      arrfnct_a,      SYMBOL_TYPE::ARRAY_FUNCTION);
    // clang-format on

    for (int i = 0; variables[i].vname != nullptr; i++) {
      // These should be immutable, but possible user is using them for some other purpose...
      // For backward compatibility, keep as mutable.
      add_variable(variables[i].vname, variables[i].value, false, true);
    }

    for (int i = 0; svariables[i].vname != nullptr; i++) {
      add_variable(svariables[i].vname, svariables[i].value, false, true);
    }

    add_variable("_C_", comment, false, true);

    std::string ap_version = SEAMS::Aprepro::version();
    auto        tokens     = tokenize(ap_version, " ");
    double      ver        = std::stod(tokens[0]);
    add_variable("_VERSION_", ver, true, true);
  }
} // namespace SEAMS
