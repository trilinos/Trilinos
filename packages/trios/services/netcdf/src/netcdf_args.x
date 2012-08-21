/* -------------------------------------------------------------------------- */
/**
 *   @file netcdf_args.x
 *
 *   @brief XDR argument structures for the netcdf proxy.
 *
 *   @author Ron Oldfield (raoldfi\@cs.sandia.gov).
 *   $Revision: 342 $.
 *   $Date: 2005-05-01 23:30:57 -0600 (Sun, 01 May 2005) $.
 *
 */

/* Extra stuff to put at the beginning of the header file */
#ifdef RPC_HDR
%#include "Trios_xdr.h"
#endif

/* Extra stuff to put at the beginning of the C file */
#ifdef RPC_XDR
%#include <pnetcdf.h>
%#include "Trios_xdr.h"
#endif

/**
 * Operation codes for the netcdf proxy.
 */
enum netcdf_opcode {
    NETCDF_NULL_OP = 10000,
    NETCDF_CREATE_OP,
    NETCDF_OPEN_OP,
    NETCDF_DEF_DIM_OP,
    NETCDF_DEF_VAR_OP,
    NETCDF_GET_ATT_OP,
    NETCDF_PUT_ATT_OP,
    NETCDF_REDEF_OP,
    NETCDF_ENDDEF_OP,
    NETCDF_PUT_VARS_OP,
    NETCDF_GET_VARS_OP,
    NETCDF_SYNC_OP,
    NETCDF_CLOSE_OP,
    NETCDF_BEGIN_INDEP_OP,
    NETCDF_END_INDEP_OP,
    NETCDF_SET_FILL_OP
};


typedef int64_t nc_size_t;

enum arg_type {
    NC_ARG_NAT,
    NC_ARG_VOID,
    NC_ARG_TEXT,
    NC_ARG_UCHAR,
    NC_ARG_SCHAR,
    NC_ARG_SHORT,
    NC_ARG_INT,
    NC_ARG_LONG,
    NC_ARG_FLOAT,
    NC_ARG_DOUBLE,
    NC_ARG_UBYTE,
    NC_ARG_USHORT,
    NC_ARG_UINT,
    NC_ARG_LONGLONG,
    NC_ARG_ULONGLONG
};

enum extra_errcodes {
    NC_ENOTSUPP = -99
};

enum write_type {
    WRITE_DIRECT,
    WRITE_AGGREGATE_INDEPENDENT,
    WRITE_AGGREGATE_COLLECTIVE,
    WRITE_CACHING_INDEPENDENT,
    WRITE_CACHING_COLLECTIVE
};

const NC_PATH_MAX = 256;

/* ********* ARGUMENTS FOR STUB FUNCTIONS ************* */

/**
 * Argument structure for nc_create.
 */
struct nc_create_args {
    string path<NC_PATH_MAX>;
    int32_t cmode;
    nc_size_t initialsz;
    nc_size_t chunksizehint;
    enum write_type write_type;
    int32_t num_participants;
};

/**
 * Result structure for nc_create.
 */
struct nc_create_res {
    int32_t ncid;
    nc_size_t chunksizehint;
};



/**
 * Argument structure for nc__open
 */
struct nc_open_args {
    string path<NC_PATH_MAX>;
    int32_t mode;
    nc_size_t chunksizehint;
    enum write_type write_type;
    int32_t num_participants;
};

struct nc_dim {
    int32_t dimid;
    string name<NC_MAX_NAME>;
    nc_size_t len;
};

struct nc_att {
    string name<NC_MAX_NAME>;
    int32_t xtype;
    nc_size_t len;
};

struct nc_var {
    int32_t varid;
    string name<NC_MAX_NAME>;
    int32_t xtype;
    int32_t dimids<NC_MAX_DIMS>;
    struct nc_att atts<NC_MAX_ATTRS>;
};

/**
 * A group is essentially a netcdf directory.
 */
struct nc_group {
    string name<NC_MAX_NAME>;
    struct nc_dim dims<NC_MAX_DIMS>;
    struct nc_var vars<NC_MAX_VARS>;
    struct nc_att atts<NC_MAX_ATTRS>;
    struct nc_group groups<NC_MAX_DIMS>;
    int32_t ncid;
    int32_t parent_ncid;
    int32_t unlimdimid;
};


/**
 * Structure for nc__open result.
 *
 * As an added optimization, we also return all variable
 * and attribute metadata along with the chunksizehint and ncid.
 */
struct nc_open_res {
    nc_size_t chunksizehint;
    int32_t unlimdimid;
    int32_t format;
    struct nc_att gatts<NC_MAX_ATTRS>;
    struct nc_group root_group;
};

/**
 * Argument structure for nc_def_dim
 */
struct nc_def_dim_args {
    int32_t ncid;
    string name<NC_MAX_NAME>;
    nc_size_t len;
};

/**
 * Marshalled argument structure for nc_def_var
 */
struct nc_def_var_args {
    int32_t ncid;
    int32_t xtype;
    string name<NC_MAX_NAME>;
    int32_t dimids<NC_MAX_DIMS>;
};


/**
 * Marshaled arguments for nc_put_att
 */
struct nc_put_att_args {
    int32_t ncid;
    int32_t varid;
    string name<NC_MAX_NAME>;
    int32_t xtype;
    arg_type atype;
    opaque data<>;
};

/**
 * Marshalled arguments for nc_get_att
 */
struct nc_get_att_args {
    int32_t ncid;
    int32_t varid;
    string name<NC_MAX_NAME>;
};


/**
 * Marshaled arguments for nc_put_var
 */
struct nc_put_vars_args {
    int32_t ncid;
    int32_t varid;
    nc_size_t start<NC_MAX_DIMS>;
    nc_size_t count<NC_MAX_DIMS>;
    nc_size_t stride<NC_MAX_DIMS>;
    arg_type atype;
    int32_t   buftype;
    nc_size_t element_count;
    nc_size_t len;
};

typedef struct nc_put_vars_args nc_get_vars_args;

/**
 * Marshaled arguments for nc_put_var (i.e., write)
 */
struct netcdf_put_args {
    int32_t unused;
};

/**
 * Marshaled arguments for nc_set_fill
 */
struct nc_set_fill_args {
    int32_t ncid;
    int32_t new_fill_mode;
};

/**
 * Result structure for nc_set_fill
 */
struct nc_set_fill_res {
    int32_t old_fill_mode;
};

/*
program NETCDF_PROG {
        version NETCDF_VERS {
                data_t xfer(data_array_t) = 1;
        } = 1;
} = 0x23451111;
 */
