/* define if we want to include experimental code */
#define INCLUDE_ZOLTAN2_EXPERIMENTAL

/* define if we have MPI */
#define HAVE_ZOLTAN2_MPI

/* define if we want to use OpenMP */
/* #undef HAVE_ZOLTAN2_OMP */

/* define if we have the METIS library*/
/* #undef HAVE_ZOLTAN2_METIS */

/* define if we have the PaToH library*/
/* #undef HAVE_ZOLTAN2_PATOH */

/* define if we have the Scotch library*/
/* #undef HAVE_ZOLTAN2_SCOTCH */

/* define if we have the ParMETIS library*/
/* #undef HAVE_ZOLTAN2_PARMETIS */

/* define if we have the AMD library*/
/* #undef HAVE_ZOLTAN2_AMD */

/* define if we have the Ovis library*/
/* #undef HAVE_ZOLTAN2_OVIS */

/* define if we have the Ovis library*/
/* #undef HAVE_ZOLTAN2_ANASAZI */

/* define to turn on runtime display of
 * total memory committed at points throughout
 * the execution of zoltan2
 */
/* #undef SHOW_ZOLTAN2_LINUX_MEMORY */

/* define to turn on runtime display of
 * maximum memory used by scotch
 * (requires a hack to the scotch library)
 */
/* #undef SHOW_ZOLTAN2_SCOTCH_MEMORY */

/* Define if Teuchos was compiled to support "long long".
 */
/* #undef HAVE_ZOLTAN2_LONG_LONG */

/* Tpetra explicit instantiation:
 *    scalar type is float
 *    local id type is int
 *    global id type is long
 */
/* #undef HAVE_ZOLTAN2_INST_FLOAT_INT_LONG */

/* Tpetra explicit instantiation:
 *    scalar type is double
 *    local id type is int
 *    global id type is long
 */
/* #undef HAVE_ZOLTAN2_INST_DOUBLE_INT_LONG */

/* Tpetra explicit instantiation:
 *    scalar type is float
 *    local id type is int
 *    global id type is int
 */
/* #undef HAVE_ZOLTAN2_INST_FLOAT_INT_INT */

/* Tpetra explicit instantiation:
 *    scalar type is double
 *    local id type is int
 *    global id type is int
 */
/* #undef HAVE_ZOLTAN2_INST_DOUBLE_INT_INT */

/* Zoltan (1) is included in the Trilinos build.
 */
#define HAVE_ZOLTAN2_ZOLTAN
