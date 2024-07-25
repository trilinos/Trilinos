"""
exodus.py v 1.21.3 (seacas-py3) is a python wrapper of some of the exodus library
(Python 3 Version)

Exodus is a common database for multiple application codes (mesh
generators, analysis codes, visualization software, etc.) affording
flexibility and robustness for both
the application code developer and application code user. By using the
Exodus data model, a user inherits the flexibility of using a large
array of application codes (including vendor-supplied codes) which
access this common data file directly or via translators.

The uses of the Exodus data model include the following:

* problem definition -- mesh generation, specification
  of locations of boundary conditions and load application, specification
  of material types.

* simulation -- model input and results output.

* visualization -- model verification, results postprocessing,
  data interrogation, and analysis tracking.

The data in Exodus files can be divided into three primary categories:
* initialization data,
* model data, and
* results data.

* Initialization data includes sizing parameters (number of
  nodes, number of elements, etc.), optional quality assurance
  information (names of codes that have operated on the data),
  and optional informational text.

* The model is described by data which are static (do not change
  through time). This data includes nodal coordinates, element
  connectivity (node lists for each element), element attributes,
  and node sets and side sets (used to aid in applying loading
  conditions and boundary constraints).

* The results are optional and include five types of variables -- nodal,
  element, nodeset, sideset, and global -- each of which is stored
  through time.

  * Nodal results are output (at each time step) for all the
    nodes in the model. An example of a nodal variable is displacement in
    the X direction.
  * Element, nodeset, and sideset results are output (at
    each time step) for all entities (elements, nodes, sides) in one or
    more entity block. For example, stress may be an element
    variable.
  * Another use of element variables is to record element status
    (a binary flag indicating whether each element is "alive" or "dead")
    through time.
  * Global results are output (at each time step) for a
    single element or node, or for a single property. Linear momentum of a
    structure and the acceleration at a particular point are both examples
    of global variables.
  * Although these examples correspond to typical FE
    applications, the data format is flexible enough to accommodate a
    spectrum of uses.

Copyright(C) 1999-2023 National Technology & Engineering Solutions
of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
NTESS, the U.S. Government retains certain rights in this software.

See packages/seacas/LICENSE for details

"""

import sys
if sys.version_info[0] < 3:
    raise Exception("Python-3 version. If using python-2, try `import exodus2 as exodus`")

import ctypes
import os
import locale
from enum import Enum

EXODUS_PY_COPYRIGHT_AND_LICENSE = __doc__

EXODUS_PY_VERSION = "1.21.3 (seacas-py3)"

EXODUS_PY_COPYRIGHT = """
You are using exodus.py v 1.21.3 (seacas-py3), a python wrapper of some of the exodus library.

Copyright (c) 2013-2023 National Technology &
Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of
Contract DE-NA0003525 with NTESS, the U.S. Government retains certain
rights in this software.
"""

EXODUS_PY_CONTACTS = """
Authors:
  Greg Sjaardema   (gdsjaar@sandia.gov)
  Mario LoPrinzi   (mvlopri@sandia.gov)
  Timothy Shelton  (trshelt@sandia.gov)
  Michael Veilleux (mgveill@sandia.gov)
  David Littlewood (djlittl@sandia.gov)
"""

# show the banner on first use
SHOW_BANNER = True

sys.dont_write_bytecode = True

ONELINE = "Gather from or export to Exodus files using the Exodus library"


def basename(file_name):
    """
    Extract base name from file_name.
    `basename("test.e") -> "test"`
    """
    return os.path.splitext(file_name)[0]


def getExodusVersion():
    """
    Parse the exodusII.h header file and return the version number or 0 if not
    found.
    """

    return _parse_exodus_version('@EXODUS_VERSION@')


def _parse_exodus_version(version_string):
    if version_string:
        assert version_string.startswith("#define EXODUS_VERSION "), "Received a incorrectly formatted verstion string. Please check the CMakeLists.txt"
        return int(version_string.strip().split()[-1].strip('"').replace('.', ''))
    else:
        return 0


try:
    locale.setlocale(locale.LC_ALL, 'en_US.utf-8')
except locale.Error:
    locale.setlocale(locale.LC_ALL, 'C')


class ex_options(Enum):
    """
    `ex_opts()` function codes - codes are OR'ed into exopts

    Attributes
    ----------
    EX_DEFAULT
         Application responsible for calling `ex_err()` to get error and warning messages to output; library is quiet
    EX_VERBOSE
         Verbose mode -- output all error and warning messages
    EX_DEBUG
         Output Debug messages
    EX_ABORT
         If an error is detected, library will abort instead of letting application decide
    EX_NULLVERBOSE
         Output error and warning messages for NULL Entity errors and warnings
    """
    EX_DEFAULT = 0
    EX_VERBOSE = 1
    EX_DEBUG = 2
    EX_ABORT = 4
    EX_NULLVERBOSE = 8


if os.name == 'nt':
    so_prefix = ''
    so_suffix = 'dll'
else:
    if os.uname()[0] == 'Darwin':
        so_prefix = 'lib'
        so_suffix = 'dylib'
    else:
        so_prefix = 'lib'
        so_suffix = 'so'
pip_path = os.path.dirname(__file__)
pip_so_path = os.path.join(pip_path, f"{so_prefix}exodus.{so_suffix}")
try:
    EXODUS_LIB = ctypes.cdll.LoadLibrary(pip_so_path)
except Exception:
    ACCESS = os.getenv('ACCESS', '@ACCESSDIR@')
    EXODUS_SO = f"{ACCESS}/@SEACAS_LIBDIR@/{so_prefix}exodus.{so_suffix}"
    EXODUS_LIB = ctypes.cdll.LoadLibrary(EXODUS_SO)

MAX_STR_LENGTH = 32      # match exodus default
MAX_NAME_LENGTH = 256     # match exodus default
MAX_LINE_LENGTH = 80      # match exodus default

EX_API_VERSION_NODOT = getExodusVersion()
EX_READ = 0x0002 if EX_API_VERSION_NODOT >= 608 else 0x0000
EX_WRITE = 0x0001  # ex_open(): open existing file for appending.
EX_NOCLOBBER = 0x0004  # does not overwrite existing exodus file
EX_CLOBBER = 0x0008  # overwrites existing exodus file

EX_NORMAL_MODEL = 0x0010  # disable mods that permit storage of larger models
EX_64BIT_OFFSET = 0x0020  # enable mods that permit storage of larger models
# enable mods that permit storage of larger models
EX_LARGE_MODEL = EX_64BIT_OFFSET
EX_64BIT_DATA = 0x400000  # CDF-5 format: classic model but 64 bit dimensions and sizes
EX_NETCDF4 = 0x0040  # use the hdf5-based netcdf4 output
EX_NOSHARE = 0x0080  # Do not open netcdf file in "share" mode
EX_SHARE = 0x0100  # Do open netcdf file in "share" mode
EX_NOCLASSIC = 0x0200  # Do not force netcdf to classic mode in netcdf4 mode

EX_DISKLESS = 0x100000  # Experimental
EX_MMAP = 0x200000  # Experimental

# Need to distinguish between storage on database (DB in name) and
# passed through the API functions (API in name).
EX_MAPS_INT64_DB = 0x0400  # All maps (id, order, ...) store int64_t values
# All entity ids (sets, blocks, maps) are int64_t values
EX_IDS_INT64_DB = 0x0800
# All integer bulk data (local indices, counts, maps); not ids
EX_BULK_INT64_DB = 0x1000
# All of the above...
EX_ALL_INT64_DB = EX_MAPS_INT64_DB + EX_IDS_INT64_DB + EX_BULK_INT64_DB

EX_MAPS_INT64_API = 0x2000  # All maps (id, order, ...) store int64_t values
# All entity ids (sets, blocks, maps) are int64_t values
EX_IDS_INT64_API = 0x4000
# All integer bulk data (local indices, counts, maps); not ids
EX_BULK_INT64_API = 0x8000
EX_INQ_INT64_API = 0x10000  # Integers passed to/from ex_inquire are int64_t
# All of the above...
EX_ALL_INT64_API = EX_MAPS_INT64_API + EX_IDS_INT64_API + \
    EX_BULK_INT64_API + EX_INQ_INT64_API

# Parallel IO mode flags...
EX_MPIIO = 0x20000
EX_MPIPOSIX = 0x40000  # \deprecated As of libhdf5 1.8.13.
EX_PNETCDF = 0x80000

# set exodus error output option
exErrPrintMode = ctypes.c_int(ex_options.EX_VERBOSE.value)
EXODUS_LIB.ex_opts(exErrPrintMode)


class ex_inquiry(Enum):
    EX_INQ_FILE_TYPE = 1                    # inquire EXODUS file type
    EX_INQ_API_VERS = 2                     # inquire API version number
    EX_INQ_DB_VERS = 3                      # inquire database version number
    EX_INQ_TITLE = 4                        # inquire database title
    EX_INQ_DIM = 5                          # inquire number of dimensions
    EX_INQ_NODES = 6                        # inquire number of nodes
    EX_INQ_ELEM = 7                         # inquire number of elements
    EX_INQ_ELEM_BLK = 8                     # inquire number of element blocks
    EX_INQ_NODE_SETS = 9                    # inquire number of node sets
    EX_INQ_NS_NODE_LEN = 10                 # inquire length of node set node list
    EX_INQ_SIDE_SETS = 11                   # inquire number of side sets
    EX_INQ_SS_NODE_LEN = 12                 # inquire length of side set node list
    EX_INQ_SS_ELEM_LEN = 13                 # inquire length of side set element list
    EX_INQ_QA = 14                          # inquire number of QA records
    EX_INQ_INFO = 15                        # inquire number of info records
    EX_INQ_TIME = 16                        # inquire number of time steps in the database
    EX_INQ_EB_PROP = 17                     # inquire number of element block properties
    EX_INQ_NS_PROP = 18                     # inquire number of node set properties
    EX_INQ_SS_PROP = 19                     # inquire number of side set properties
    # inquire length of node set distribution factor list
    EX_INQ_NS_DF_LEN = 20
    # inquire length of side set distribution factor list
    EX_INQ_SS_DF_LEN = 21
    EX_INQ_LIB_VERS = 22                    # inquire API Lib vers number
    EX_INQ_EM_PROP = 23                     # inquire number of element map properties
    EX_INQ_NM_PROP = 24                     # inquire number of node map properties
    EX_INQ_ELEM_MAP = 25                    # inquire number of element maps
    EX_INQ_NODE_MAP = 26                    # inquire number of node maps
    EX_INQ_EDGE = 27                        # inquire number of edges
    EX_INQ_EDGE_BLK = 28                    # inquire number of edge blocks
    EX_INQ_EDGE_SETS = 29                   # inquire number of edge sets
    # inquire length of concat edge set edge list
    EX_INQ_ES_LEN = 30
    # inquire length of concat edge set dist factor list
    EX_INQ_ES_DF_LEN = 31
    # inquire number of properties stored per edge block
    EX_INQ_EDGE_PROP = 32
    # inquire number of properties stored per edge set
    EX_INQ_ES_PROP = 33
    EX_INQ_FACE = 34                        # inquire number of faces
    EX_INQ_FACE_BLK = 35                    # inquire number of face blocks
    EX_INQ_FACE_SETS = 36                   # inquire number of face sets
    # inquire length of concat face set face list
    EX_INQ_FS_LEN = 37
    # inquire length of concat face set dist factor list
    EX_INQ_FS_DF_LEN = 38
    # inquire number of properties stored per face block
    EX_INQ_FACE_PROP = 39
    # inquire number of properties stored per face set
    EX_INQ_FS_PROP = 40
    EX_INQ_ELEM_SETS = 41                   # inquire number of element sets
    # inquire length of concat element set element list
    EX_INQ_ELS_LEN = 42
    # inquire length of concat element set dist factor list
    EX_INQ_ELS_DF_LEN = 43
    # inquire number of properties stored per elem set
    EX_INQ_ELS_PROP = 44
    EX_INQ_EDGE_MAP = 45                    # inquire number of edge maps
    EX_INQ_FACE_MAP = 46                    # inquire number of face maps
    EX_INQ_COORD_FRAMES = 47                # inquire number of coordinate frames
    # inquire size of MAX_NAME_LENGTH dimension on database
    EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH = 48
    # inquire size of MAX_NAME_LENGTH dimension on database
    EX_INQ_DB_MAX_USED_NAME_LENGTH = 49
    # inquire client-specified max size of returned names
    EX_INQ_MAX_READ_NAME_LENGTH = 50
    # inquire size of floating-point values stored on database
    EX_INQ_DB_FLOAT_SIZE = 51
    EX_INQ_ASSEMBLY = 60
    EX_INQ_BLOB = 61
    EX_INQ_INVALID = -1


class ex_type(Enum):
    EX_CHAR = 2     # NC_CHAR (from netcdf.h)
    EX_INTEGER = 4  # NC_INT (from netcdf.h)
    EX_DOUBLE = 6   # NC_DOUBLE (from netcdf.h)
    EX_INVALID = -1


def ex_type_map(type):
    """
    Map the exodus inquiry flags to an enum value.
    """
    return ex_type[type].value


def ex_inquiry_map(inquiry):
    """
    Map the exodus inquiry flags to an enum value.
    """
    return ex_inquiry[inquiry].value


def ex_obj_to_inq(objType):
    """
    Return the ex_inquiry string corresponding to the specified objType.
    This can be passed to the ex_inquiry_map() function to get the number of
    objects of the specified objType
    """
    entity_dictionary = {
        'EX_ASSEMBLY': 'EX_INQ_ASSEMBLY',
        'EX_BLOB': 'EX_INQ_BLOB',
        'EX_EDGE_BLOCK': 'EX_INQ_EDGE_BLK',
        'EX_FACE_BLOCK': 'EX_INQ_FACE_BLK',
        'EX_ELEM_BLOCK': 'EX_INQ_ELEM_BLK',
        'EX_NODE_SET': 'EX_INQ_NODE_SETS',
        'EX_EDGE_SET': 'EX_INQ_EDGE_SETS',
        'EX_FACE_SET': 'EX_INQ_FACE_SETS',
        'EX_ELEM_SET': 'EX_INQ_ELEM_SETS',
        'EX_SIDE_SET': 'EX_INQ_SIDE_SETS',
        'EX_NODE_MAP': 'EX_INQ_NODES',
        'EX_EDGE_MAP': 'EX_INQ_EDGE',
        'EX_FACE_MAP': 'EX_INQ_FACE',
        'EX_ELEM_MAP': 'EX_INQ_ELEM',
    }

    return entity_dictionary.get(objType, -1)


def ex_entity_type_to_objType(entity_type):
    """
    """
    entity_dictionary = {
        get_entity_type('EX_ASSEMBLY'): "assembly",
        get_entity_type('EX_BLOB'): "blob",
        get_entity_type('EX_EDGE_BLOCK'): "edge block",
        get_entity_type('EX_FACE_BLOCK'): "face block",
        get_entity_type('EX_ELEM_BLOCK'): "element block",
        get_entity_type('EX_NODE_SET'): "node set",
        get_entity_type('EX_EDGE_SET'): "edge set",
        get_entity_type('EX_FACE_SET'): "face set",
        get_entity_type('EX_ELEM_SET'): "element set",
        get_entity_type('EX_SIDE_SET'): "side set",
        get_entity_type('EX_NODE_MAP'): "node map",
        get_entity_type('EX_EDGE_MAP'): "edge map",
        get_entity_type('EX_FACE_MAP'): "face map",
        get_entity_type('EX_ELEM_MAP'): "element map"
    }

    return entity_dictionary.get(entity_type, 'EX_INVALID')


class ex_entity_type(Enum):
    """
    The ex_entity_type enum from the exodusII.h include file

    Attributes
    ----------
    EX_NODAL
         nodal \"block\" for variables
    EX_NODE_BLOCK
         alias for EX_NODAL
    EX_NODE_SET
         node set property code
    EX_EDGE_BLOCK
         edge block property code
    EX_EDGE_SET
         edge set property code
    EX_FACE_BLOCK
         face block property code
    EX_FACE_SET
         face set property code
    EX_ELEM_BLOCK
         element block property code
    EX_ELEM_SET
         face set property code
    EX_SIDE_SET
         side set property code
    EX_ELEM_MAP
         element map property code
    EX_NODE_MAP
         node map property code
    EX_EDGE_MAP
         edge map property code
    EX_FACE_MAP
         face map property code
    EX_GLOBAL
         global \"block\" for variables
    EX_COORDINATE
         kluge so some internal wrapper functions work
    EX_ASSEMBLY
         assemblies (collections of other entities)
    EX_BLOB
         blob (arbitrary sized binary object)
    EX_INVALID
         invalid
    """
    EX_NODAL = 14
    EX_NODE_BLOCK = 14
    EX_NODE_SET = 2
    EX_EDGE_BLOCK = 6
    EX_EDGE_SET = 7
    EX_FACE_BLOCK = 8
    EX_FACE_SET = 9
    EX_ELEM_BLOCK = 1
    EX_ELEM_SET = 10
    EX_SIDE_SET = 3
    EX_ELEM_MAP = 4
    EX_NODE_MAP = 5
    EX_EDGE_MAP = 11
    EX_FACE_MAP = 12
    EX_GLOBAL = 13
    EX_COORDINATE = 15
    EX_ASSEMBLY = 16
    EX_BLOB = 17
    EX_INVALID = -1


def get_entity_type(varType):
    """
    Map the exodus ex_entity_type flags to an integer value.
    """
    try:
        return ex_entity_type[varType].value
    except KeyError:
        return varType if (isinstance(varType, int)) else varType.value


# init params struct
class ex_init_params(ctypes.Structure):
    """
    Parameters defining the model dimension, note that many are optional.

    Attributes
    ----------
    num_dim : int
        number of model dimensions
    num_nodes : int
        number of model nodes
    num_edge : int
        number of model edges
    num_edge_blk : int
        number of model edge blocks
    num_face : int
        number of model faces
    num_face_blk : int
        number of model face blocks
    num_elem : int
        number of model elements
    num_elem_blk : int
        number of model element blocks
    num_node_sets : int
        number of model node sets
    num_edge_sets : int
        number of model edge sets
    num_face_sets : int
        number of model face sets
    num_side_sets : int
        number of model side sets
    num_elem_sets : int
        number of model elem sets
    num_node_maps : int
        number of model node maps
    num_edge_maps : int
        number of model edge maps
    num_face_maps : int
        number of model face maps
    num_elem_maps : int
        number of model elem maps
    num_assembly : int
        number of assemblies
    num_blob : int
        number of blobs
    """
    _fields_ = [("title", ctypes.c_char * (MAX_LINE_LENGTH + 1)),
                ("num_dim", ctypes.c_longlong),
                ("num_nodes", ctypes.c_longlong),
                ("num_edge", ctypes.c_longlong),
                ("num_edge_blk", ctypes.c_longlong),
                ("num_face", ctypes.c_longlong),
                ("num_face_blk", ctypes.c_longlong),
                ("num_elem", ctypes.c_longlong),
                ("num_elem_blk", ctypes.c_longlong),
                ("num_node_sets", ctypes.c_longlong),
                ("num_edge_sets", ctypes.c_longlong),
                ("num_face_sets", ctypes.c_longlong),
                ("num_side_sets", ctypes.c_longlong),
                ("num_elem_sets", ctypes.c_longlong),
                ("num_node_maps", ctypes.c_longlong),
                ("num_edge_maps", ctypes.c_longlong),
                ("num_face_maps", ctypes.c_longlong),
                ("num_elem_maps", ctypes.c_longlong),
                ("num_assembly", ctypes.c_longlong),
                ("num_blob", ctypes.c_longlong)]


class assembly:
    def __init__(self, name, id, type):
        self.name = name
        self.id = id
        if isinstance(type, str):
            type = ex_entity_type[type].value
        if isinstance(type, ex_entity_type):
            type = type.value
        self.type = type
        self.entity_list = []

    def __repr__(self):
        return "assembly(name=%r, type=%r, id=%r, members=%r)" % (self.name, self.type, self.id, self.entity_list)


class ex_assembly(ctypes.Structure):
    """
    Structure defining the assembly parameters.

    Attributes
    ----------
    id : int64_t
    name : char *
    type : ex_entity_type
    entity_count : int
    entity_list : void_int *
    """
    _fields_ = [("id", ctypes.c_longlong),
                ("name", ctypes.c_char_p),
                ("type", ctypes.c_int),
                ("entity_count", ctypes.c_int),
                ("entity_list", ctypes.POINTER(ctypes.c_longlong))]


def setup_ex_assembly(assembly):
    ctype_assem = ex_assembly(id=assembly.id, name=assembly.name.encode(), type=assembly.type)
    ctype_assem.entity_count = len(assembly.entity_list)
    ctype_assem.entity_list = (ctypes.c_longlong * ctype_assem.entity_count)(*assembly.entity_list)
    return ctype_assem


class blob(object):
    def __init__(self, name, id, num_entry):
        self.name = name
        self.id = id
        self.num_entry = num_entry

    def __repr__(self):
        return "blob(name=%r, id=%r, num_entry=%r)" % (self.name, self.id, self.num_entry)


class ex_blob(ctypes.Structure):
    """
    Structure defining the blob parameters.

    Attributes
    ----------
    id : int64_t
    name : char *
    num_entry : int64_t
    """
    _fields_ = [("id", ctypes.c_longlong),
                ("name", ctypes.c_char_p),
                ("num_entry", ctypes.c_longlong)]


class attribute:
    def __init__(self, name, type, id):
        self.name = name
        if isinstance(type, str):
            type = ex_entity_type[type].value
        if isinstance(type, ex_entity_type):
            type = type.value
        self.entity_type = type
        self.entity_id = id
        self.values = []

    def __repr__(self):
        return "attribute(name=%r, entity_type=%r, entity_id=%r, values=%r)" % (self.name, self.entity_type, self.entity_id, self.values)


class ex_attribute(ctypes.Structure):
    """
    Used for accessing underlying exodus library...

    Attributes
    ----------
    entity_type : ex_entity_type
    entity_id : int64_t
    name : char[257]
    type : ex_type
    value_count : int
    values : void*
    """
    _fields_ = [("entity_type", ctypes.c_int),
                ("entity_id", ctypes.c_longlong),
                ("name", ctypes.c_char * 257),
                ("type", ctypes.c_int),
                ("value_count", ctypes.c_int),
                ("values", ctypes.c_void_p)]


class exodus:
    """
    The exodus model abstraction
    """

    #
    # construction of a new exodus object
    #
    # --------------------------------------------------------------------

    def __init__(self, file, mode=None, array_type='ctype', title=None,
                 numDims=None, numNodes=None, numElems=None, numBlocks=None,
                 numNodeSets=None, numSideSets=None, numAssembly=None,
                 numBlob=None, init_params=None, io_size=0):
        """
        Open exodus database for data insertion/extraction.

        Parameters
        ----------
        file_name : string
           name of exodus file to open
        mode : string
          'r' for read, 'a' for append, 'w' for write, 'w+' for write+clobber existing file
        title : string
           database title
        array_type : string
           'ctype' for c-type arrays, 'numpy' for numpy arrays
        num_dims : int
           number of model dimensions ('w'/'w+' mode only)
        num_nodes : int
           number of model nodes ('w'/'w+' mode only)
        num_elems : int
           number of model elements ('w'/'w+' mode only)
        num_blocks : int
           number of model element blocks ('w'/'w+' mode only)
        num_ns : int
           number of model node sets ('w'/'w+' mode only)
        num_ss : int
           number of model side sets ('w'/'w+' mode only)

        init_params : ex_init_params
           see :py:func:`exodus.ex_init_params` for more info.

        Returns
        -------
        exo : exodus object
            the open exodus database


        >>> ex_pars = ex_init_params(num_dim=numDims, num_nodes=numNodes,
        ...                          num_elem=numElems, num_elem_blk=numElemBlocks, num_assembly=numAssembly)
        >>> exo = exodus(file_name, mode=mode, title=title,
        ...             array_type=array_type, init_params=ex_pars)
        >>> with exodus(file_name, mode=mode, title=title,\
        ...             array_type=array_type, init_params=ex_pars) as exo:
        ...     pass
        """
        global SHOW_BANNER
        if SHOW_BANNER:
            print(EXODUS_PY_COPYRIGHT)
            SHOW_BANNER = False

        if mode is None:
            mode = 'r'

        if array_type == 'numpy':
            # Import numpy to convert from c-type arrays to numpy arrays
            # (Numpy is imported here rather than at the module level, so that
            # the import only occurs if the user specifies a numpy array type.
            # This way, platforms without numpy installed can still import the
            # exodus.py module and just use ctype arrays.)
            import numpy as np
            self.np = np
            self.use_numpy = True
            # Warnings module is needed to suppress the invalid warning when
            # converting from c-type arrays to numpy arrays
            # http://stackoverflow.com/questions/4964101/pep-3118-warning-when-using-ctypes-array-as-numpy-array
            import warnings
            self.warnings = warnings
        else:
            self.use_numpy = False

        self.EXODUS_LIB = EXODUS_LIB
        self.fileName = str(file)
        self.basename = basename(file)
        self.modeChar = mode
        self.fileId = None
        self.__open(io_size=io_size)
        EXODUS_LIB.ex_set_max_name_length(self.fileId, MAX_NAME_LENGTH)
        if mode.lower() in ['w', 'w+']:
            if init_params is not None:
                self.init_params = init_params
                if title is not None:
                    self.init_params.title = title
                self.put_info_ext(self.init_params)
            else:
                if numNodeSets is None:
                    numNodeSets = 0
                if numSideSets is None:
                    numSideSets = 0
                if numNodes is None:
                    numNodes = 0
                if numElems is None:
                    numElems = 0
                if numBlocks is None:
                    numBlocks = 0

                info = [title, numDims, numNodes, numElems, numBlocks,
                        numNodeSets, numSideSets]
                if None not in info:
                    self.__ex_put_info(info)

            self.numTimes = ctypes.c_int(0)
        else:
            self.__ex_get_info()
            self.numTimes = ctypes.c_int(
                self.__ex_inquire_int(
                    ex_inquiry_map('EX_INQ_TIME')))

        self.coordsX = None
        self.coordsY = None
        self.coordsZ = None
        self.times = None

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()
        if not traceback:
            return True

    def summarize(self):
        """
        Outputs a summary of the exodus file data. Output is similar to::

           Database: base_ioshell_copy.e
           Title:  This is the title

           Number of spatial dimensions = 3                                                 Number of global variables     = 10
           Number of node blocks        = 1         Number of nodes              = 1,331    Number of nodal variables      =  2
           Number of element blocks     = 1         Number of elements           = 1,000    Number of element variables    =  5
           Number of node sets          = 3         Length of node list          =   363    Number of nodeset variables    =  4
           Number of element side sets  = 3         Length of element sides      =   300    Number of sideset variables    =  3
           Number of assemblies         = 4                                                 Number of assembly variables   = 10
           Number of blobs              = 0                                                 Number of blob     variables   =  0
           Number of time steps         = 5
        """

        sidesets = self.get_ids('EX_SIDE_SET')
        total_sides = sum(self.num_faces_in_side_set(sideset) for sideset in sidesets)
        nodesets = self.get_ids('EX_NODE_SET')
        total_ns_nodes = sum(self.num_nodes_in_node_set(nodeset) for nodeset in nodesets)

        num_glo_vars = self.get_variable_number('EX_GLOBAL')
        num_nod_vars = self.get_variable_number('EX_NODAL')
        num_ele_vars = self.get_variable_number('EX_ELEM_BLOCK')
        num_ns_vars = self.get_variable_number('EX_NODE_SET')
        num_ss_vars = self.get_variable_number('EX_SIDE_SET')
        num_assem_vars = self.get_reduction_variable_number('EX_ASSEMBLY')
        num_blob_vars = self.get_reduction_variable_number('EX_BLOB')

        print("\n Database: {0}\n"
              " Title:\t{17}\n\n"
              " Number of spatial dimensions = {1:3d}\t"
              "                                   {2:11s}\t"
              " Number of global variables     = {11:6d}\n"
              " Number of node blocks        = {5:3d}\t"
              " Number of nodes              = {3:10n}\t"
              " Number of nodal variables      = {12:6d}\n"
              " Number of element blocks     = {6:3n}\t"
              " Number of elements           = {4:10n}\t"
              " Number of element variables    = {13:6d}\n"
              " Number of node sets          = {7:3n}\t"
              " Length of node list          = {9:10n}\t"
              " Number of nodeset variables    = {14:6d}\n"
              " Number of element side sets  = {8:3n}\t"
              " Length of element sides      = {10:10n}\t"
              " Number of sideset variables    = {15:6d}\n"
              " Number of assemblies         = {18:3n}\t"
              "                                   {2:11s}\t"
              " Number of assembly red vars    = {19:6d}\n"
              " Number of blobs              = {20:3n}\t"
              "                                   {2:11s}\t"
              " Number of blob red vars        = {21:6d}\n"
              " Number of time steps         = {16:3n}\n"
              .format(self.fileName,
                      self.num_dimensions(), "",
                      self.num_nodes(),
                      self.num_elems(),
                      1,
                      self.num_blks(),
                      self.num_node_sets(),
                      self.num_side_sets(),
                      total_ns_nodes, total_sides,
                      num_glo_vars, num_nod_vars, num_ele_vars,
                      num_ns_vars, num_ss_vars, self.num_times(), self.title(),
                      self.num_assembly(), num_assem_vars,
                      self.num_blob(), num_blob_vars))

    def put_info_ext(self, p):
        """
        put initialization information into exodus file

        >>> e.put_info_ext(info_struct)
        """
        if len(p.title) > MAX_LINE_LENGTH:
            print(f'WARNING: Exodus title \"{p.title}\" exceeds maximum line length ({MAX_LINE_LENGTH}). It will be truncated.')

            p.title = p.title[-1 * MAX_LINE_LENGTH:]

        self.Title = ctypes.create_string_buffer(p.title, MAX_LINE_LENGTH + 1)
        self.numDim = ctypes.c_longlong(p.num_dim)
        self.numNodes = ctypes.c_longlong(p.num_nodes)
        self.numElem = ctypes.c_longlong(p.num_elem)
        self.numElemBlk = ctypes.c_longlong(p.num_elem_blk)
        self.numNodeSets = ctypes.c_longlong(p.num_node_sets)
        self.numSideSets = ctypes.c_longlong(p.num_side_sets)
        self.numAssembly = ctypes.c_longlong(p.num_assembly)

        EXODUS_LIB.ex_put_init_ext(self.fileId, ctypes.byref(p))
        return True

    def copy(self, fileName, include_transient=False, mode='a'):
        """
        Copies exodus database to file_name and returns an opened copy as a
        new exodus object. This object will need to be closed when it is done
        being used.

        >>> exo_copy = exo.copy(file_name)
        >>> exo_copy.close()

        Parameters
        ----------
        file_name : str
            name of exodus file to open

        Returns
        -------
        exo_copy : exodus object opened in append mode by default
        """
        i64Status = EXODUS_LIB.ex_int64_status(self.fileId)
        fileId = EXODUS_LIB.ex_create_int(fileName.encode('ascii'), EX_NOCLOBBER | i64Status,
                                          ctypes.byref(self.comp_ws),
                                          ctypes.byref(self.io_ws),
                                          EX_API_VERSION_NODOT)

        self.copy_file(fileId, include_transient)
        EXODUS_LIB.ex_close(fileId)

        return exodus(fileName, mode)

    def copy_file(self, file_id, include_transient=False):
        """
        Copies exodus database to the database pointed to by `fileId`
        Returns the passed in `file_id`.

        >>> with exodus.exodus(file_name, mode='w') as exofile:
        >>>     with exo.copy_file(exofile.fileId, include_transient=True) as exo_copy:
        >>>         exo_copy.close()

        Parameters
        ----------
        fileId : str
            name of exodus file to open
        include_transient: bool
            should the transient data in the original file also be copied to the output file
            or just the mesh (non-transient) portion.

        Returns
        -------
        file_id: The file_id of the copied to file

        """
        EXODUS_LIB.ex_copy(self.fileId, file_id)
        if include_transient:
            EXODUS_LIB.ex_copy_transient(self.fileId, file_id)

        return file_id

    def title(self):
        """
        get the database title

        >>> title = exo.title()

        Returns
        -------
        title : string
        """
        return self.Title.value.decode('utf8')

    def version_num(self):
        """
        get exodus version number used to create the database

        >>> version = exo.version_num()

        Returns
        -------
        version : string
            representation of version number
        """
        return "%1.2f" % self.version.value

    def put_info(self, Title, numDim, numNodes, numElem, numElemBlk,
                 numNodeSets, numSideSets):
        """
        Initialize static metadata for the database.

        >>> status = exo.put_info(title, num_dims, num_nodes, num_elems,
        ...                      num_blocks, num_ns, num_ss)

        Parameters
        ----------
        title : string
            database title
        num_dims : int
            number of model dimensions
        num_nodes : int
            number of model nodes
        num_elems : int
            number of model elements
        num_blocks : int
            number of model element blocks
        num_ns : int
            number of model node sets
        num_ss : int
            number of model side sets

        Returns
        -------
        status : bool
            True = successful execution

        """
        self.__ex_put_info([Title, numDim, numNodes, numElem,
                            numElemBlk, numNodeSets, numSideSets])
        return True

    def inquire(self, inquiry):
        """
        Inquire about various properties of the database

        Returns
        -------
        inq_res : int
        """
        return self.__ex_inquire_int(ex_inquiry_map(inquiry))

    def num_qa_records(self):
        """
        get the number of qa records

        >>> num_qa_recs = exo.num_qa_records()

        Returns
        -------
        num_qa_recs : int
        """
        return int(self.__ex_inquire_int(ex_inquiry_map('EX_INQ_QA')))

    def get_qa_records(self):
        """
        get a list of QA records where each QA record is a length-4 tuple of strings:
          1. the software name that accessed/modified the database
          2. the software descriptor, e.g. version
          3. additional software data
          4. time stamp

        >>> qa_recs = exo.get_qa_records()


        Returns
        -------
        qa_recs : list<tuple[4]<string>>
        """
        return self.__ex_get_qa()

    def put_qa_records(self, records):
        """
        store a list of QA records where each QA record is a length-4 tuple of strings:
          1. the software name that accessed/modified the database
          2. the software descriptor, e.g. version
          3. additional software data
          4. time stamp

        >>> status = exo.put_qa_records()

        Parameters
        ----------
        qa_recs : list<tuple[4]<string>>

        Returns
        ------
        status : bool
            True = successful execution
        """
        for rec in records:
            assert len(rec) == 4
            for recEntry in rec:
                assert len(str(recEntry).encode('ascii')) < MAX_STR_LENGTH
        return self.__ex_put_qa(records)

    def num_info_records(self):
        """
        get the number of info records

        >>> num_info_recs = exo.num_info_records()

        Returns
        -------
        num_info_recs : int
        """
        return int(self.__ex_inquire_int(ex_inquiry_map('EX_INQ_INFO')))

    def get_info_records(self):
        """
        get a list info records where each entry in the list is one info
        record, e.g. a line of an input deck

        >>> info_recs = exo.get_info_records()

        Returns
        -------
        info_recs : list<string>

        """
        return self.__ex_get_info_recs()

    def put_info_records(self, info):
        """
        store a list of info records where each entry in the list is
        one line of info, e.g. a line of an input deck

        >>> status = exo.put_info_records(info)

        Parameters
        ----------
        info_recs : list<tuple[4]<string>>

        Returns
        -------
        status : bool
            True = successful execution
        """
        for rec in info:
            if len(str(rec).encode('ascii')) > MAX_LINE_LENGTH:
                print("WARNING: max line length reached for one or more info records;")
                print(
                    "         info stored to exodus file is incomplete for these records")
                break
        return self.__ex_put_info_recs(info)

    # --------------------------------------------------------------------

    def get_sierra_input(self, inpFileName=None):
        """
        parse sierra input deck from the info records if inp_file_name
        is passed the deck is written to this file; otherwise a list
        of input deck file lines is returned

        >>> inp = exo.get_sierra_input(inpFileName=inp_file_name)

        Parameters
        ----------
        inp_file_name : string, optional
           Name of text file where info records corresponding to the Sierra input deck will be written

        Returns
        -------
        inp : list<string>
           lines if inp_file_name not provided; otherwise, an empty list

        """
        info_recs = self.__ex_get_info_recs()
        sierra_inp = []
        begin = False
        for rec in info_recs:
            vals = rec.split()
            if not begin:  # have not reached Sierra block
                if len(vals) >= 2 and vals[0].lower() == "begin" and vals[1].lower() == "sierra":
                    begin = True
            if begin:  # inside Sierra block
                sierra_inp.append(rec)
                if len(rec) > MAX_LINE_LENGTH:
                    print(
                        "WARNING: max line length reached for one or more input lines;")
                    print("         input data might be incomplete for these lines")
                    break
                if len(vals) >= 2 and vals[0].lower() == "end" and vals[1].lower() == "sierra":
                    break  # end of Sierra block

        if inpFileName:
            with open(inpFileName.encode('ascii'), "w") as fd:
                for fileLine in sierra_inp:
                    fd.write(fileLine + "\n")
            return []

        return sierra_inp

    #
    # time steps
    #
    # --------------------------------------------------------------------

    def num_times(self):
        """
        get the number of time steps

        >>> num_times = exo.num_times()

        Returns
        -------
        num_times : int
        """
        return self.numTimes.value

    # --------------------------------------------------------------------

    def get_times(self):
        """
        get the time values

        >>> time_vals = exo.get_times()

        Returns
        -------
            if array_type == 'ctype' :
              <list<ctypes.c_double>>  time_vals

            if array_type == 'numpy' :
              <np_array<double>>  time_vals
        """
        if self.numTimes.value == 0:
            self.times = []
        else:
            self.__ex_get_all_times()
        if self.use_numpy:
            self.times = ctype_to_numpy(self, self.times)
        return self.times

    # --------------------------------------------------------------------

    def put_time(self, step, value):
        """
        store a new time

        >>> exo.put_time(time_step, time_val)

        Parameters
        ----------
        time_step : int
            time step index (1-based)
        time_val : double
            time value for this step

        Returns
        -------
        status : bool
            True = successful execution
        """
        self.__ex_put_time(step, value)
        self.numTimes = ctypes.c_int(self.__ex_inquire_int(ex_inquiry_map('EX_INQ_TIME')))
        return True

    #
    # coordinate system
    #
    # --------------------------------------------------------------------

    def num_dimensions(self):
        """
        get the number of model spatial dimensions

        >>> num_dims = exo.num_dimensions()

        Returns
        -------
        num_dims : int
        """
        return self.numDim.value

    # --------------------------------------------------------------------

    def get_coord_names(self):
        """
        get a list of length exo.num_dimensions() that has the name
        of each model coordinate direction, e.g. ['x', 'y', 'z']

        >>> coord_names = exo.get_coord_names()

        Returns
        -------
        coord_names : list<string>
        """
        return self.__ex_get_coord_names()

    # --------------------------------------------------------------------

    def put_coord_names(self, names):
        """
        store a list of length exo.num_dimensions() that has the name
        of each model coordinate direction, e.g. ['x', 'y', 'z']

        >>> exo.put_coord_names()

        Parameters
        ----------
        coord_names : list<string>
        """
        self.__ex_put_coord_names(names)

    #
    # nodes
    #
    # --------------------------------------------------------------------

    def num_nodes(self):
        """
        get the number of nodes in the model

        >>> num_nodes = exo.num_nodes()

        Returns
        -------
        num_nodes : int
        """
        return self.numNodes.value

    # --------------------------------------------------------------------

    def get_coords(self):
        """
        get model coordinates of all nodes; for each coordinate
        direction, a length exo.num_nodes() list is returned

        >>> x_coords, y_coords, z_coords = exo.get_coords()

        Returns
        -------

            if array_type == 'ctype':
              <list<ctypes.c_double>>  x_coords  global x-direction coordinates
              <list<ctypes.c_double>>  y_coords  global y-direction coordinates
              <list<ctypes.c_double>>  z_coords  global z-direction coordinates

            if array_type == 'numpy':
              <np_array<double>>  x_coords  global x-direction coordinates
              <np_array<double>>  y_coords  global y-direction coordinates
              <np_array<double>>  z_coords  global z-direction coordinates
        """
        self.__ex_get_coord()
        if self.use_numpy:
            self.coordsX = ctype_to_numpy(self, self.coordsX)
            self.coordsY = ctype_to_numpy(self, self.coordsY)
            self.coordsZ = ctype_to_numpy(self, self.coordsZ)
        return self.coordsX, self.coordsY, self.coordsZ

    # --------------------------------------------------------------------

    def get_coord(self, i):
        """
        get model coordinates of a single node

        >>> x_coord, y_coord, z_coord = exo.get_coord(node_index)

        Parameters
        ----------
        node_index : int
            the 1-based node index (indexing is from 1 to exo.num_nodes())

        Returns
        -------
        x_coord : double
            global x-direction coordinate
        y_coord : double
            global y-direction coordinate
        z_coord : double
            global z-direction coordinate

        Note
        -----
        >>> x_coords, y_coords, z_coords = exo.get_coords()
        >>> x_coord = x_coords[node_index-1]
        >>> y_coord = y_coords[node_index-1]
        >>> z_coord = z_coords[node_index-1]
            ... is equivalent to ...
        >>> x_coord, y_coord, z_coord = exo.get_coords(node_index)

        """
        listX, listY, listZ = self.__ex_get_partial_coord(i, 1)
        return listX[0], listY[0], listZ[0]

    # --------------------------------------------------------------------

    def put_coords(self, xCoords, yCoords, zCoords):
        """
        store model coordinates of all nodes; for each coordinate
        direction, a length exo.num_nodes() list is input

        >>> status = exo.put_coords(x_coords, y_coords, z_coords)

        Parameters
        ----------
        x_coord : list<double>
            global x-direction coordinates
        y_coord : list<double>
            global y-direction coordinates
        z_coord : list<double>
            global z-direction coordinates

        Returns
        -------
        status : bool
            True = successful execution
        """
        self.__ex_put_coord(xCoords, yCoords, zCoords)
        return True

    # --------------------------------------------------------------------

    def get_node_num_map(self):
        """
        **DEPRECATED** use: `get_node_id_map()`

        get mapping of exodus node index to user- or application-
        defined node id; node_id_map is ordered the same as the nodal
        coordinate arrays returned by exo.get_coords() -- this ordering
        follows the exodus node *INDEX* order, a 1-based system going
        from 1 to exo.num_nodes(); a user or application can optionally
        use a separate node *ID* numbering system, so the node_id_map
        points to the node *ID* for each node *INDEX*

        >>> node_id_map = exo.get_node_num_map()

        Returns
        -------
        node_id_map : list<ctypes.c_int>
        """
        return self.__ex_get_node_num_map()

    # --------------------------------------------------------------------

    def put_node_id_map(self, id_map):
        """
        store mapping of exodus node index to user- or application-
        defined node id; node_id_map is ordered the same as the nodal
        coordinate arrays returned by exo.get_coords() -- this ordering
        follows the exodus node *INDEX* order, a 1-based system going
        from 1 to exo.num_nodes(); a user or application can optionally
        use a separate node *ID* numbering system, so the node_id_map
        points to the node *ID* for each node *INDEX*

        >>> status = exo.put_node_id_map(node_id_map)

        Parameters
        ----------
        node_id_map : list<int>

        Returns
        -------
        status : bool
            True = successful execution
        """
        return self.__ex_put_id_map('EX_NODE_MAP', id_map)

    # --------------------------------------------------------------------

    def get_node_variable_names(self):
        """
        get the list of nodal variable names in the model

        >>> nvar_names = exo.get_node_variable_names()

        Returns
        -------
        nvar_names : list<string>
        """
        if self.__ex_get_variable_param('EX_NODAL').value == 0:
            return []
        return self.__ex_get_variable_names('EX_NODAL')

    # --------------------------------------------------------------------

    def get_node_variable_number(self):
        """
        get the number of nodal variables in the model

        >>> num_nvars = exo.get_node_variable_number()

        Returns
        -------
        num_nvars : int
        """
        return self.__ex_get_variable_param('EX_NODAL').value

    # --------------------------------------------------------------------

    def set_node_variable_number(self, number):
        """
        update the number of nodal variables in the model

        >>> status = exo.set_node_variable_number(num_nvars)

        Parameters
        ----------
        num_nvars : int

        Returns
        -------
        status : bool
            True = successful execution
        """
        self.__ex_put_variable_param('EX_NODAL', number)
        return True

    # --------------------------------------------------------------------

    def put_node_variable_name(self, name, index):
        """
        add the name and index of a new nodal variable to the model;
        nodal variable indexing goes from 1 to exo.get_node_variable_number()

        >>> status = exo.put_node_variable_name(nvar_name, nvar_index)

        Parameters
        ----------
        nvar_name : string
            name of new nodal variable
        nvar_index : int
            1-based index of new nodal variable

        Returns
        -------
        status : bool
            True = successful execution

        Example
        -------
        this method is often called within the following sequence:

        >>> num_nvars = exo.get_node_variable_number()
        >>> new_nvar_index = num_nvars + 1
        >>> num_nvars += 1
        >>> exo.set_node_variable_number(num_nvars)
        >>> exo.put_node_variable_name("new_nvar_name", new_nvar_index)
        """
        NDvarNames = self.get_variable_names('EX_NODAL')
        if name in NDvarNames:
            print(f'WARNING: node variable \"{name}\" already exists.')
        if index > len(NDvarNames):
            raise Exception("ERROR: variable index out of range.")
        self.__ex_put_variable_name('EX_NODAL', index, name)
        return True

    # --------------------------------------------------------------------

    def get_node_variable_values(self, name, step):
        """
        get list of nodal variable values for a nodal variable name
        and time step

        >>> nvar_vals = exo.get_node_variable_values(nvar_name, time_step)

        Parameters
        ----------
        name : string
            name of nodal variable
        step : int
            1-based index of time step

        Returns
        -------

            if array_type == 'ctype':
              <list<ctypes.c_double>>  nvar_vals

            if array_type == 'numpy':
              <np_array<double>>  nvar_vals
        """
        names = self.get_variable_names('EX_NODAL')
        var_id = names.index(name) + 1
        numVals = self.num_nodes()
        values = self.__ex_get_var(step, 'EX_NODAL', var_id, 0, numVals)
        if self.use_numpy:
            values = ctype_to_numpy(self, values)
        return values

    # --------------------------------------------------------------------

    def get_partial_node_variable_values(self, name, step, start_index, num_nodes):
        """
        get partial list of nodal variable values for a nodal variable name
        and time step.  Start at node `node_index` (1-based) and return `num_nodes`
        from that point.

        >>> nvar_vals = exo.get_partial_node_variable_values(nvar_name, time_step, 10, 100)

        Parameters
        ----------
        nvar_name : string
             name of nodal variable
        time_step  : int
             1-based index of time step
        start_index : int
             1-based index of node to start returning data
        num_nodes : int
             number of nodes to return data for.

        Returns
        -------

            if array_type == 'ctype':
              <list<ctypes.c_double>>  nvar_vals

            if array_type == 'numpy':
              <np_array<double>>  nvar_vals
        """
        names = self.get_variable_names('EX_NODAL')
        var_id = names.index(name) + 1
        values = self.__ex_get_partial_var(step, 'EX_NODAL', var_id, 0, start_index, num_nodes)
        if self.use_numpy:
            values = ctype_to_numpy(self, values)
        return values

    # --------------------------------------------------------------------

    def put_node_variable_values(self, name, step, values):
        """
        store a list of nodal variable values for a nodal variable
        name and time step

        >>> status = exo.put_node_variable_values(nvar_name, time_step, nvar_vals)

        Parameters
        ----------
        nvar_name : string
             name of nodal variable
        time_step : int
             1-based index of time step
        nvar_vals : list<double>

        Returns
        -------
        status : bool
            True = successful execution
        """
        names = self.get_variable_names('EX_NODAL')
        var_id = names.index(name) + 1
        numVals = self.num_nodes()
        self.__ex_put_var(step, 'EX_NODAL', var_id, 0, numVals, values)
        return True

    #
    # elements
    #
    # --------------------------------------------------------------------

    def num_elems(self):
        """
        get the number of elements in the model

        >>> num_elems = exo.num_elems()

        Returns
        -------
        num_elems : int
        """
        return self.numElem.value

    # --------------------------------------------------------------------

    def get_num_map(self, mapType, idx):
        """
        get user-defined map of exodus element/node/edge/face index to user- or
        application- defined element/node/edge/face values. Map values are arbitrary integers

        >>> elem_num_map = exo.get_num_map('EX_ELEM_MAP', 1)

        Parameters
        ----------
        mapType   : ex_entity_type
             type of map being queried ('EX_ELEM_MAP', 'EX_NODE_MAP', 'EX_FACE_MAP', 'EX_EDGE_MAP')
        idx       : int
             which map to return (1-based).  Use `inquire(mapType)` to get number of maps stored on database.
        Returns
        -------

            if array_type == 'ctype':
              <list<int>>  num_map

            if array_type == 'numpy':
              <np_array<int>>  num_map


        >>> em_cnt = exo.inquire('EX_INQ_ELEM_MAP')
        >>> em     = exo.get_names('EX_ELEM_MAP')
        >>> map    = exo.get_num_map('EX_ELEM_MAP', 2)

        """
        return self.__ex_get_num_map(mapType, idx)

    # --------------------------------------------------------------------

    def put_num_map(self, mapType, idx, num_map):
        """
        put user-defined map of exodus element/node/edge/face index to user- or
        application- defined element/node/edge/face values. Map values are arbitrary integers


        Parameters
        ----------
            mapType   : ex_entity_type
                        type of map being written ('EX_ELEM_MAP', 'EX_NODE_MAP', 'EX_FACE_MAP', 'EX_EDGE_MAP')
            idx       : int
                        which map to write (1-based).  Use `put_map_param(node_map_cnt, elem_map_cnt)` prior to this
                        function to define number of maps on the database.
            elem_id_map : list<int>


        >>> exo.put_map_param(nm_cnt, em_cnt)
        >>> nm[0] = "My_Node_Map"
        >>> exo.put_names('EX_NODE_MAP', nm);
        >>> exo.put_num_map('EX_NODE_MAP', 1, scale_map)

        """
        return self.__ex_put_num_map(mapType, idx, num_map)

    def put_map_param(self, node_map_cnt, elem_map_cnt):
        """
        Define number of node and element maps that will be written to the database

        Parameters
        ----------
        node_map_cnt  : int
                        number of node maps
        elem_map_cnt  : int
                        number of element maps

        """
        return self.__ex_put_map_param(node_map_cnt, elem_map_cnt)

    # --------------------------------------------------------------------

    def get_id_map(self, mapType):
        """
        get mapping of exodus element/node/edge/face index to user- or
        application- defined element/node/edge/face id; id_map is ordered by the
        *INDEX* ordering, a 1-based system going from 1 to
        `exo.num_elem`, `exo.num_node`, used by exodus for storage and input/output
        of array data stored on the elements/nodes/edges/faces; a user or application
        can optionally use a separate *ID* numbering system,
        so the id_map points to the element/node/edge/face *ID* for each
        *INDEX*

        >>> elem_id_map = exo.get_id_map('EX_ELEM_MAP')

        Parameters
        ----------
        mapType   : ex_entity_type
                   type of map being queried ('EX_ELEM_MAP', 'EX_NODE_MAP', 'EX_FACE_MAP', 'EX_EDGE_MAP')

        Returns
        -------

            if array_type == 'ctype':
              <list<int>>  id_map

            if array_type == 'numpy':
              <np_array<int>>  id_map

        """
        return self.__ex_get_id_map(mapType)

    # --------------------------------------------------------------------

    def get_elem_id_map(self):
        """
        get mapping of exodus element index to user- or application-
        defined element id; elem_id_map is ordered by the element
        *INDEX* ordering, a 1-based system going from 1 to
        exo.num_elems(), used by exodus for storage and input/output
        of array data stored on the elements; a user or application
        can optionally use a separate element *ID* numbering system,
        so the elem_id_map points to the element *ID* for each
        element *INDEX*

        >>> elem_id_map = exo.get_elem_id_map()

        Returns
        -------

            if array_type == 'ctype':
              <list<int>>  elem_id_map

            if array_type == 'numpy':
              <np_array<int>>  elem_id_map
        """
        return self.__ex_get_id_map('EX_ELEM_MAP')

    # --------------------------------------------------------------------

    def put_elem_id_map(self, id_map):
        """
        store mapping of exodus element index to user- or application-
        defined element id; elem_id_map is ordered by the element
        *INDEX* ordering, a 1-based system going from 1 to
        exo.num_elems(), used by exodus for storage and input/output
        of array data stored on the elements; a user or application
        can optionally use a separate element *ID* numbering system,
        so the elem_id_map points to the element *ID* for each
        element *INDEX*

        >>> status = exo.put_elem_id_map(elem_id_map)

        Parameters
        ----------
            elem_id_map : list<int>

        Returns
        -------
        status : bool
            True = successful execution
        """
        return self.__ex_put_id_map('EX_ELEM_MAP', id_map)

    # --------------------------------------------------------------------

    def put_id_map(self, map_type, id_map):
        """
        store mapping of exodus entity index to user- or application-
        defined entity id; id_map is ordered by the entity
        *INDEX* ordering, a 1-based system going from 1 to
        exo.num_XXX(), used by exodus for storage and input/output
        of array data stored on the entities; a user or application
        can optionally use a separate entity *ID* numbering system,
        so the elem_id_map points to the entity *ID* for each
        entity *INDEX*

        >>> status = exo.put_id_map('EX_ELEM_MAP', elem_id_map)

        Parameters
        ----------
            map_type : ex_entity_type
                        type of map being queried ('EX_ELEM_MAP', 'EX_NODE_MAP', 'EX_FACE_MAP', 'EX_EDGE_MAP')
            elem_id_map : list<int>

        Returns
        -------
        status : bool
            True = successful execution
        """
        return self.__ex_put_id_map(map_type, id_map)

    # --------------------------------------------------------------------

    def get_block_id_map(self, obj_type, entity_id):
        """
        Gets the map of elements found in the given entity_id of an object
        of obj_type.

        *INDEX* ordering, a 1-based system going from 1 to
        number of elements in the elem_block, used by exodus for
        storage and input/output of array data stored on the elements;
        a user or application can optionally use a separate element *ID* numbering system,
        so the elem_id_map points to the element *ID* for each
        element *INDEX*

        >>> elem_block_id_map = exo.get_block_id_map("EX_ELEM_BLOCK", 100)

        Returns
        -------

            if array_type == 'ctype':
              <list<int>>  elem_id_map

            if array_type == 'numpy':
              <np_array<int>>  elem_id_map
        """
        return self.__ex_get_block_id_map(obj_type, entity_id)

    # --------------------------------------------------------------------

    def get_elem_num_map(self):
        """
        **DEPRECATED** use: `get_elem_id_map()`

        get mapping of exodus element index to user- or application-
        defined element id; elem_id_map is ordered by the element
        *INDEX* ordering, a 1-based system going from 1 to
        exo.num_elems(), used by exodus for storage and input/output
        of array data stored on the elements; a user or application
        can optionally use a separate element *ID* numbering system,
        so the elem_id_map points to the element *ID* for each
        element *INDEX*

        >>> elem_id_map = exo.get_elem_num_map()

        Returns
        -------
        elem_id_map : list<ctypes.c_int>
        """
        return self.__ex_get_elem_num_map()

    # --------------------------------------------------------------------

    def get_elem_order_map(self):
        """
        get mapping of exodus element index to application-defined
        optimal ordering; elem_order_map is ordered by the element
        index ordering used by exodus for storage and input/output
        of array data stored on the elements; a user or application
        can optionally use a separate element ordering, e.g. for
        optimal solver performance, so the elem_order_map points to
        the index used by the application for each exodus element
        index

        >>> elem_order_map = exo.get_elem_order_map()

        Returns
        -------

            if array_type == 'ctype':
              <list<int>>  elem_order_map

            if array_type == 'numpy':
              <np_array<int>>  elem_order_map
        """

        elemOrderMap = self.__ex_get_elem_order_map()
        if self.use_numpy:
            elemOrderMap = ctype_to_numpy(self, elemOrderMap)
        return elemOrderMap

    # Generic (objType) get/put/query...
    # --------------------------------------------------------------------

    # --------------------------------------------------------------------

    def get_node_id_map(self):
        """
        get mapping of exodus node index to user- or application-
        defined node id; node_id_map is ordered the same as the nodal
        coordinate arrays returned by exo.get_coords() -- this ordering
        follows the exodus node *INDEX* order, a 1-based system going
        from 1 to exo.num_nodes(); a user or application can optionally
        use a separate node *ID* numbering system, so the node_id_map
        points to the node *ID* for each node *INDEX*

        >>> node_id_map = exo.get_node_id_map()

        Returns
        -------

            if array_type == 'ctype':
              <list<int>>  node_id_map

            if array_type == 'numpy':
              <np_array<int>>  node_id_map
        """
        return self.__ex_get_id_map('EX_NODE_MAP')

    # --------------------------------------------------------------------

    def get_name(self, object_type, object_id):
        """
        get the name of the specified entity_type and entity

        >>> elem_blk_name = exo.get_name('EX_ELEM_BLOCK', elem_blk_id)

        Parameters
        ----------
        object_type : int
            block/set type
        object_id : ex_entity_type
            block/set *ID* (not *INDEX*)

        Returns
        -------
        name : string
        """
        return self.__ex_get_name(object_type, object_id)

    # --------------------------------------------------------------------

    def put_name(self, object_type, object_id, name):
        """
        put the name of the specified entity_type and entity

        >>> exo.put_name('EX_ELEM_BLOCK', elem_blk_id, block_name)

        Parameters
        ----------
        object_type : int
            block/set type
        object_id : ex_entity_id
            block/set *ID* (not *INDEX*)
        name : string
            block/set name

        Returns
        -------
        elem_blk_name : string
        """
        self.__ex_put_name(object_type, object_id, name)

    # --------------------------------------------------------------------

    def get_ids(self, objType):
        """
        get mapping of exodus block/set index to user- or application-
        defined block/set id; ids is ordered
        by the *INDEX* ordering, a 1-based system going from
        1 to number_set_or_block, used by exodus for storage
        and input/output of array data stored on the blocks/sets; a
        user or application can optionally use a separate block/set
        *ID* numbering system, so the ids array points to the
        block/set *ID* for each set *INDEX*

        >>> node_set_ids = exo.get_ids('EX_NODE_SET')

        Returns
        -------

            if array_type == 'ctype':
              <list<int>>  ids

            if array_type == 'numpy':
              <np_array<int>>  ids
        """
        ids = self.__ex_get_ids(objType)
        if self.use_numpy:
            ids = self.np.array(ids)
        return ids

    # --------------------------------------------------------------------

    def get_names(self, object_type):
        """
        get a list of all block/set names ordered by block/set *INDEX*;
        (see :py:func:`exodus.get_ids` for explanation of the
        difference between *ID* and *INDEX*)

        >>> blk_names = exo.get_names('EX_ELEM_BLOCK')

        Parameters
        ----------
        object_type : int
            block/set type

        Returns
        -------
        names : list<string>
        """
        return self.__ex_get_names(object_type)

    # --------------------------------------------------------------------

    def put_names(self, object_type, names):
        """
        store a list of all block/set names of the specified
        `object_type` ordered by *INDEX*;
        (see :py:func:`exodus.get_ids` for explanation of the
        difference between *ID* and *INDEX*)

        >>> exo.put_names('EX_ELEM_BLOCK', elem_blk_names)

        Parameters
        ----------
        object_type : int
        names : list<string>
        """

        self.__ex_put_names(object_type, names)

    # --------------------------------------------------------------------

    def get_reduction_variable_values(self, objType, id, step):
        """
        get list of reduction variable values for a specified entity type and
        id, and time step

        >>> evar_vals = exo.get_reduction_variable_values('EX_ELEM_BLOCK', elem_blk_id, time_step)

        Parameters
        ----------
        objType   : ex_entity_type
            type of object being queried
        id        : int
            entity *ID* (not *INDEX*)
        time_step : int
            1-based index of time step

        Returns
        -------

            if array_type == 'ctype':
              <list<ctypes.c_double>>  evar_vals

            if array_type == 'numpy':
              <np_array<double>>  evar_vals

        """
        numVals = self.get_reduction_variable_number(objType)
        values = self.__ex_get_reduction_vars(step, objType, id, numVals)
        if self.use_numpy:
            values = ctype_to_numpy(self, values)
        return values

    # --------------------------------------------------------------------

    def put_reduction_variable_values(self, objType, id, step, values):
        """
        store a list of 'objType' variable values for a specified entity,
        and time step

        >>> status = exo.put_redcution_variable_values('EX_ELEM_BLOCK', elem_blk_id,
        ...             time_step, evar_vals)

        Parameters
        ----------
        objType : ex_entity_type
            type of object begin queried
        id : ex_entity_id
            element block *ID* (not *INDEX*)
        step : int
            1-based index of time step
        values : list<double>

        Returns
        -------
        status : bool
            True = successful execution
        """
        numVals = self.get_reduction_variable_number(objType)
        self.__ex_put_reduction_vars(step, objType, id, numVals, values)
        return True

    # --------------------------------------------------------------------
    def get_variable_truth_table(self, objType, entId=None):
        """
        gets a truth table indicating which variables are defined for
        specified entity type; if entId is not passed, then a concatenated
        truth table for all entities is returned with variable index
        cycling faster than entity index

        >>> ssvar_truth_tab = exo.get_variable_truth_table('EX_SIDE_SET', sideSetID=side_set_id)

        Parameters
        ----------
        objType : ex_entity_type
            type of object begin queried
        entid : ex_entity_id, optional
            entity *ID* (not *INDEX*)

        Returns
        -------
        truth_tab : list<bool>
            True for variable defined in an entity, False otherwise
        """
        if entId is None:
            return self.__ex_get_truth_table(objType)
        else:
            return self.__ex_get_object_truth_vector(objType, entId)

    # --------------------------------------------------------------------

    def set_variable_truth_table(self, objType, table):
        """
        stores a truth table indicating which variables are defined for
        all sets/blocks of the specified `objType` and all variables; variable index cycles
        faster than entity index

        >>> status = exo.set_variable_truth_table('EX_NODE_SET', nsvar_truth_tab)

        Parameters
        ----------
        objType : ex_entity_type
            type of object begin queried
        table : list<bool>
            True for variable defined in a node set, False otherwise

        Returns
        -------
        status : bool
            True = successful execution
        """
        return self.__ex_put_truth_table(objType, table)

    # --------------------------------------------------------------------

    def get_variable_names(self, objType):
        """
        get the list of variable names in the model for the specified object type.

        >>> nar_names = exo.get_variable_names('EX_NODAL')

        Returns
        -------
        nvar_names : list<string>
        """
        if self.__ex_get_variable_param(objType).value == 0:
            return []
        return self.__ex_get_variable_names(objType)

    # --------------------------------------------------------------------

    def get_reduction_variable_names(self, objType):
        """
        get the list of reduction variable names in the model for the specified object type.

        >>> nar_names = exo.get_reduction_variable_names('EX_ASSEMBL"Y')

        Returns
        -------
        nvar_names : list<string>
        """
        if self.__ex_get_reduction_variable_param(objType).value == 0:
            return []
        return self.__ex_get_reduction_variable_names(objType)

    # --------------------------------------------------------------------

    def get_reduction_variable_name(self, objType, varId):
        """
        get a single reduction variable name in the model for the specified object type and index.

        >>> nar_name = exo.get_reduction_variable_name('EX_ASSEMBL"Y', 100)

        Returns
        -------
        nvar_name : string
        """
        if self.__ex_get_reduction_variable_param(objType).value == 0:
            return ""
        return self.__ex_get_reduction_variable_name(objType, varId)

    # --------------------------------------------------------------------

    def get_variable_number(self, objType):
        """
        get the number of variables of the specified type in the model

        >>> num_nvars = exo.get_variable_number('EX_NODAL')

        Returns
        -------
        num_nvars : int
        """
        return self.__ex_get_variable_param(objType).value

    # --------------------------------------------------------------------

    def get_reduction_variable_number(self, objType):
        """
        get the number of reduction variables of the specified type in the model

        >>> num_nvars = exo.get_reduction_variable_number('EX_ASSEMBLY')

        Returns
        -------
        num_nvars : int
        """
        return self.__ex_get_reduction_variable_param(objType).value

    # --------------------------------------------------------------------

    def set_variable_number(self, objType, number):
        """
        update the number of variables in the model

        >>> status = exo.set_variable_number('EX_NODAL', num_nvars)

        Parameters
        ----------
        objType : ex_entity_type
            type of object begin queried
        number : int

        Returns
        -------
        status : bool
            True = successful execution
        """
        self.__ex_put_variable_param(objType, number)
        return True

    # --------------------------------------------------------------------

    def set_reduction_variable_number(self, objType, number):
        """
        update the number of reduction variables in the model

        >>> status = exo.set_reduction_variable_number('EX_ASSEMBLY', num_nvars)

        Parameters
        ----------
        objType : ex_entity_type
            type of object begin queried
        number : int

        Returns
        -------
        status : bool
            True = successful execution
        """
        self.__ex_put_reduction_variable_param(objType, number)
        return True

    # --------------------------------------------------------------------

    def put_variable_name(self, objType, name, index):
        """
        add the name and index of a new variable to the model;
        variable indexing goes from 1 to exo.get_variable_number()

        >>> status = exo.put_variable_name('EX_NODAL', nvar_name, nvar_index)

        Parameters
        ----------
        objType : ex_entity_type
            type of object begin queried
        var_name : string
            name of new variable
        nvar_index : int
            1-based index of new nodal variable

        Returns
        -------
        status : bool
            True = successful execution

        Example
        -------
        this method is often called within the following sequence:

        >>> num_nvars = exo.get_variable_number('EX_NODAL')
        >>> new_nvar_index = num_nvars + 1
        >>> num_nvars += 1
        >>> exo.set_variable_number('EX_NODAL', num_nvars)
        >>> exo.put_variable_name('EX_NODAL', "new_nvar_name", new_nvar_index)
        """
        varNames = self.get_variable_names(objType)
        if name in varNames:
            print(f'WARNING: variable \"{name}\" already exists.')
        if index > len(varNames):
            raise Exception("ERROR: variable index out of range.")
        self.__ex_put_variable_name(objType, index, name)
        return True

    # --------------------------------------------------------------------

    def put_reduction_variable_name(self, objType, name, index):
        """
        add the name and index of a new reduction variable to the model;
        variable indexing goes from 1 to `get_reduction_variable_number()`

        >>> status = exo.put_reduction_variable_name('EX_ASSEMBLY', assemvar_name, assemvar_index)

        Parameters
        ----------
        objType : ex_entity_type
            type of object begin queried
        var_name : string
            name of new variable
        index : int
            1-based index of new variable

        Returns
        -------
        status : bool
            True = successful execution

        Example
        -------
        this method is often called within the following sequence:

        >>> num_assem_vars = exo.get_reduction_variable_number('EX_ASSEMBLY')
        >>> new_assem_var_index = num_assem_vars + 1
        >>> num_assem_vars += 1
        >>> exo.set_reduction_variable_number('EX_ASSEMBLY', num_assem_vars)
        >>> exo.put_reduction_variable_name('EX_ASSEMBLY', "new_assem_var_name", new_assem_var_index)
        """
        varNames = self.get_reduction_variable_names(objType)
        if name in varNames:
            print(f'WARNING: variable \"{name}\" already exists.')
        if index > len(varNames):
            raise Exception("ERROR: variable index out of range.")
        self.__ex_put_reduction_variable_name(objType, index, name)
        return True

    def get_entity_count(self, objType, entityId):
        """
        get number of nodes/elements in the specified set/block. Typically used internally
        by other user-callable functions.

        Parameters
        ----------
        objType   : ex_entity_type
            type of object being queried
        entityid : ex_entity_id
            id of the entity (block, set) *ID* (not *INDEX*)
        """

        numVals = 0
        if objType == "EX_NODAL":
            numVals = self.num_nodes()
        elif objType in ['EX_ELEM_BLOCK', 'EX_FACE_BLOCK', 'EX_EDGE_BLOCK']:
            (_elemType, numVals, _nodesPerElem, _numAttr) = self.__ex_get_block(objType, entityId)
        elif objType in ['EX_NODE_SET', 'EX_EDGE_SET', 'EX_FACE_SET', 'EX_SIDE_SET']:
            (numVals, _numDistFactInSet) = self.__ex_get_set_param(objType, entityId)

        return numVals

    def get_variable_values_time(self, objType, entityId, var_name, start_step, end_step):
        """
        get list of `objType` variable values for a specified object id
        block, variable name, and range of time steps

        >>> evar_vals = exo.get_variable_values_time('EX_ELEM_BLOCK', entity_id,
        ...                                            var_name, start_step, end_step)

        Parameters
        ----------
        objType   : ex_entity_type
            type of object being queried
        entityid : ex_entity_id
            id of the entity (block, set) *ID* (not *INDEX*)
        var_name : string
            name of variable
        start_step : int
            1-based index of time step
        end_step : int
            1-based index of time step

        Returns
        -------

            if array_type == 'ctype':
              <list<ctypes.c_double>>  evar_vals

            if array_type == 'numpy':
              <np_array<double>>  evar_vals
        """
        names = self.get_variable_names(objType)
        var_id = names.index(var_name) + 1
        values = self.__ex_get_var_time(objType, var_id, entityId, start_step, end_step)
        if self.use_numpy:
            values = ctype_to_numpy(self, values)
        return values

    def get_variable_values(self, objType, entityId, name, step):
        """
        get list of `objType` variable values for a specified object id
        block, variable name, and time step

        >>> evar_vals = exo.get_variable_values('EX_ELEM_BLOCK', elem_blk_id,
        ...                                            evar_name, time_step)

        Parameters
        ----------
        objType   : ex_entity_type
            type of object being queried
        entityid : ex_entity_id
            id of the entity (block, set) *ID* (not *INDEX*)
        name : string
            name of variable
        time_step : int
            1-based index of time step

        Returns
        -------

            if array_type == 'ctype':
              <list<ctypes.c_double>>  evar_vals

            if array_type == 'numpy':
              <np_array<double>>  evar_vals
        """
        names = self.get_variable_names(objType)
        var_id = names.index(name) + 1
        numVals = self.get_entity_count(objType, entityId)

        values = self.__ex_get_var(step, objType, var_id, entityId, numVals)
        if self.use_numpy:
            values = ctype_to_numpy(self, values)
        return values

    def put_variable_values(self, objType, entityId, name, step, values):
        """
        store a list of element variable values for a specified element
        block, element variable name, and time step

        >>> status = exo.put_variable_values('EX_ELEM_BLOCK', elem_blk_id,
        ...             evar_name, time_step, evar_vals)

        Parameters
        ----------
        entityid : ex_entity_id
           entity *ID* (not *INDEX*)
        name : string
           name of variable
        time_step : int
           1-based index of time step
        values : list<double>
           the variable values to be output

        Returns
        -------
        status : bool
            True = successful execution
        """
        names = self.get_variable_names(objType)
        var_id = names.index(name) + 1
        numVals = self.get_entity_count(objType, entityId)

        self.__ex_put_var(step, objType, var_id, entityId, numVals, values)
        return True

    def get_attribute_count(self, objType, objId):
        """
        IS THIS NEEDED, PYTHONIC WAY MAY BE TO JUST GET THEM...

        get the number of attributes on the specified entity

        >>> num_attribute = exo.get_attribute_count('EX_ASSEMBLY', 100)

        Parameters
        ----------
        objType   : ex_entity_type
            type of object being queried
        objId        : int
            entity *ID* (not *INDEX*)

        Returns
        -------
        num_attribute : int
        """
        return self.__ex_get_attribute_count(objType, objId)

    def get_attributes(self, objType, objId):
        """
        >>> attributes = exo.get_attributes('EX_ASSEMBLY', 100)

        Parameters
        ----------
        objType   : ex_entity_type
            type of object being queried
        objId        : int
            entity *ID* (not *INDEX*)

        Returns
        -------
        attributes : ex_attribute list
        """

        return self.__ex_get_attributes(objType, objId)

    def put_attribute(self, attribute):
        """
        >>> attribute = exodus.attribute('Scale', 'EX_ASSEMBLY', 100)
        >>> attribute.values = [1.1, 1.0, 1.2]
        >>> attributes = exo.put_attribute(attribute)

        Returns
        -------
        attributes : ex_attribute list
        """

        return self.__ex_put_attribute(attribute)

    def num_assembly(self):
        """
        get the number of assemblies in the model

        >>> num_assembly = exo.num_assembly()

        Returns
        -------
        num_assembly : int
        """
        return self.inquire('EX_INQ_ASSEMBLY')

    def get_assembly(self, object_id):
        """
        reads the assembly parameters and assembly data for one assembly
        """
        assem = ex_assembly(id=object_id)
        self.__ex_get_assembly(assem)
        assmbly = assembly(assem.name.decode('utf8'), assem.id, assem.type)
        for j in range(assem.entity_count):
            assmbly.entity_list.append(assem.entity_list[j])
        return assmbly

    def get_assemblies(self, object_ids):
        """
        reads the assembly parameters and assembly data for all assemblies
        with ids in object_ids
        """
        assemblies = [ex_assembly(id=object_id) for object_id in object_ids]
        assems = (ex_assembly * len(assemblies))(*assemblies)
        self.__ex_get_assemblies(assems)
        assembs = [assembly(assem.name.decode('utf8'), assem.id, assem.type) for assem in
                   assems]
        for i, a in enumerate(assems):
            for j in range(a.entity_count):
                assembs[i].entity_list.append(a.entity_list[j])
        return assembs

    def put_assembly(self, assembly):
        """
        writes the assembly parameters and assembly data for one assembly
        """
        self.__ex_put_assembly(assembly)

    def put_assemblies(self, assemblies):
        """
        writes the assembly parameters and assembly data for multiple assemblies
        """
        self.__ex_put_assemblies(assemblies)

    def num_blob(self):
        """
        get the number of blobs in the model

        >>> num_assembly = exo.num_blob()

        Returns
        -------
        num_blob : int
        """
        return self.numBlob.value

    def get_blob(self, object_id):
        """
        reads the blob parameters and blob data for one blob
        """
        blob = ex_blob(id=object_id)
        self.__ex_get_blob(blob)
        return blob

    def num_blks(self):
        """
        get the number of element blocks in the model

        >>> num_elem_blks = exo.num_blks()

        Returns
        -------
        num_elem_blks : int
        """
        return self.numElemBlk.value

    def get_elem_blk_ids(self):
        """
        get mapping of exodus element block index to user- or
        application-defined element block id; elem_blk_ids is ordered
        by the element block *INDEX* ordering, a 1-based system going
        from 1 to exo.num_blks(), used by exodus for storage
        and input/output of array data stored on the element blocks; a
        user or application can optionally use a separate element block
        *ID* numbering system, so the elem_blk_ids array points to the
        element block *ID* for each element block *INDEX*

        >>> elem_blk_ids = exo.get_elem_blk_ids()

        Returns
        -------

            if array_type == 'ctype':
              <list<int>>  elem_blk_ids

            if array_type == 'numpy':
              <np_array<int>>  elem_blk_ids
        """
        return self.get_ids('EX_ELEM_BLOCK')

    def get_elem_blk_name(self, object_id):
        """
        get the element block name

        >>> elem_blk_name = exo.get_elem_blk_name(elem_blk_id)

        Parameters
        ----------
        object_id : ex_entity_id
            element block *ID* (not *INDEX*)

        Returns
        -------
        elem_blk_name : string
        """
        return self.__ex_get_name('EX_ELEM_BLOCK', object_id)

    def put_elem_blk_name(self, object_id, name):
        """
        store the element block name

        >>> exo.put_elem_blk_name(elem_blk_id, elem_blk_name)

        Parameters
        ----------
        elem_blk_id : ex_entity_id
            element block *ID* (not *INDEX*)
        elem_blk_name : string
        """
        self.__ex_put_name('EX_ELEM_BLOCK', object_id, name)

    def get_elem_blk_names(self):
        """
        get a list of all element block names ordered by block *INDEX*;
        (see :py:func:`exodus.get_ids` for explanation of the
        difference between block *ID* and block *INDEX*)

        >>> elem_blk_names = exo.get_elem_blk_names()

        Returns
        -------
        elem_blk_names : list<string>
        """
        return self.__ex_get_names('EX_ELEM_BLOCK')

    def put_elem_blk_names(self, names):
        """
        store a list of all element block names ordered by block *INDEX*;
        (see :py:func:`exodus.get_ids` for explanation of the
        difference between block *ID* and block *INDEX*)

        >>> exo.put_elem_blk_names(elem_blk_names)

        Parameters
        ----------
        elem_blk_names : list<string>
        """
        self.__ex_put_names('EX_ELEM_BLOCK', names)

    def elem_blk_info(self, object_id):
        """
        get the element block info

        >>> elem_type, num_blk_elems, num_elem_nodes, num_elem_attrs
        ...       = exo.elem_blk_info(elem_blk_id)

        Parameters
        ----------
        elem_blk_id : ex_entity_id
            element block *ID* (not *INDEX*)

        Returns
        -------
        elem_type : string
            element type, e.g. 'HEX8'
        num_blk_elems : int
            number of elements in the block
        num_elem_nodes : int
            number of nodes per element
        num_elem_attrs : int
            number of attributes per element
        """
        (elemType, numElem, nodesPerElem, numAttr) = self.__ex_get_block('EX_ELEM_BLOCK', object_id)
        return elemType, numElem, nodesPerElem, numAttr

    def put_elem_blk_info(self, elem_blk_id, elem_type, num_blk_elems,
                          num_elem_nodes, num_elem_attrs):
        """
        store the element block *ID* and element block info

        >>> exo.put_elem_blk_info(elem_blk_id, elem_type, num_blk_elems,
        ...                      num_elem_nodes, num_elem_attrs)

        Parameters
        ----------
        elem_blk_id : ex_entity_id
            element block *ID* (not *INDEX*)
        elem_type : string
            element type (all caps), e.g. 'HEX8'
        num_blk_elems : int
            number of elements in the block
        num_elem_nodes : int
            number of nodes per element
        num_elem_attrs : int
            number of attributes per element
        """
        self.__ex_put_block('EX_ELEM_BLOCK', elem_blk_id, elem_type, num_blk_elems,
                            num_elem_nodes, num_elem_attrs)

    def put_concat_elem_blk(self, elem_blk_ids, elem_type, num_blk_elems,
                            num_elem_nodes, num_elem_attrs, defineMaps):
        """
        same as exo.put_elem_blk_info() but for all blocks at once

        >>> status = exo.put_concat_elem_blk(elem_blk_ids, elem_types,
        ...                                 num_blk_elems, num_elem_nodes, num_elem_attrs)

        Parameters
        ----------
        elem_blk_ids : list<int>
              element block *ID* (not *INDEX*) for each block
        elem_types   : list<string>
              element type for each block
        num_blk_elems : list<int>
              number of elements for each block
        num_elem_nodes : list<int>
              number of nodes per element for each block
        num_elem_attrs : list<int>
              number of attributes per element for each block

        Returns
        -------
        status : bool
            True = successful execution
        """
        self.__ex_put_concat_elem_blk(
            elem_blk_ids,
            elem_type,
            num_blk_elems,
            num_elem_nodes,
            num_elem_attrs,
            defineMaps)
        return True

    def get_elem_connectivity(self, object_id):
        """
        get the nodal connectivity, number of elements, and
        number of nodes per element for a single block

        >>> elem_conn, num_blk_elems, num_elem_nodes
        ...        = exo.get_elem_connectivity(elem_blk_id)

        Parameters
        ----------
        elem_blk_id : ex_entity_id
            element block *ID* (not *INDEX*)

        Returns
        -------
        elem_conn : <list<int>>  (if array_type == 'ctype')
        elem_conn : <np_array<int>>  (if array_type == 'numpy')
            ordered list of node *INDICES* that define the connectivity of each element
            in the block; the list cycles through all nodes of the first element, then
            all nodes of the second element, etc. (see :py:func:`exodus.get_id_map` for explanation
            of node *INDEX* versus node *ID*)
        num_blk_elems : int
            number of elements in the block
        num_elem_nodes : int
            number of nodes per element
        """
        (elem_block_connectivity, num_elem_this_blk,
         num_nodes_per_elem) = self.__ex_get_elem_conn(object_id)
        if self.use_numpy:
            elem_block_connectivity = ctype_to_numpy(
                self, elem_block_connectivity)
        return elem_block_connectivity, num_elem_this_blk, num_nodes_per_elem

    def put_elem_connectivity(self, object_id, connectivity):
        """
        store the nodal connectivity, number of elements, and
        number of nodes per element for a single block

        >>> exo.put_elem_connectivity(elem_blk_id, elem_conn)

        Parameters
        ----------
        elem_blk_id : ex_entity_id
            element block *ID* (not *INDEX*)
        connectivity : list<int>
            ordered list of node *INDICES* that
            define the connectivity of each
            element in the block; the list cycles
            through all nodes of the first element,
            then all nodes of the second element,
            etc.
            (see :py:func:`exodus.get_id_map` for explanation
            of node *INDEX* versus node *ID*)
        """
        _d1, numBlkElems, numNodesPerElem, _d2 = self.elem_blk_info(object_id)
        assert len(connectivity) == (numBlkElems * numNodesPerElem)
        self.__ex_put_elem_conn(object_id, connectivity)

    def get_elem_attr(self, elem_blk_id):
        """
        get all attributes for each element in a block

        >>> elem_attrs = exo.get_elem_attr(elem_blk_id)

        Parameters
        ----------
        elem_blk_id : ex_entity_id
            element block *ID* (not *INDEX*)

        Returns
        -------
        if array_type == 'ctype' : <list<double>> elem_attrs
        if array_type == 'numpy' : <np_array<double>> elem_attrs
            list of attribute values for all
            elements in the block; the list cycles
            through all attributes of the first
            element, then all attributes of the
            second element, etc. Attributes are
            ordered by the ordering of the names
            returned by :py:func:`exodus.get_attribute_names`
        """
        elem_attrs = self.__ex_get_elem_attr(elem_blk_id)
        if self.use_numpy:
            elem_attrs = ctype_to_numpy(self, elem_attrs)
        return elem_attrs

    def get_elem_attr_values(self, elem_blk_id, elem_attr_name):
        """
        get an attribute for each element in a block

        >>> elem_attrs = exo.get_elem_attr(elem_blk_id)

        Parameters
        ----------
        elem_blk_id : ex_entity_id
             element block *ID* (not *INDEX*)
        elem_attr_name : string
             element attribute name

        Returns
        -------
            if array_type == 'ctype': <list<double>>  values
            if array_type == 'numpy': <np_array<double>>  values
                array of values for the requested
                attribute.  Array has dimensions of
                1 x num_elem, where num_elem is the
                number of elements on the element block.
        """
        # Determine index of requested attribute in attribute list
        elem_attr_names = self.get_attribute_names('EX_ELEM_BLOCK', elem_blk_id)
        a_ndx = elem_attr_names.index(elem_attr_name)

        values = self.__ex_get_one_attr('EX_ELEM_BLOCK', elem_blk_id, a_ndx + 1)
        if self.use_numpy:
            values = ctype_to_numpy(self, values)
        return values

    def get_attr_values(self, objType, elem_blk_id, elem_attr_name):
        """
        get an attribute for each element in a block

        >>> elem_attrs = exo.get_elem_attr(elem_blk_id)

        Parameters
        ----------
        objType   : ex_entity_type
            type of object being queried
        elem_blk_id : ex_entity_id
             element block *ID* (not *INDEX*)
        elem_attr_name : string
             element attribute name

        Returns
        -------
            if array_type == 'ctype': <list<double>>  values
            if array_type == 'numpy': <np_array<double>>  values
                array of values for the requested
                attribute.  Array has dimensions of
                1 x num_elem, where num_elem is the
                number of elements on the element block.
        """
        # Determine index of requested attribute in attribute list
        elem_attr_names = self.get_attribute_names('EX_ELEM_BLOCK', elem_blk_id)
        a_ndx = elem_attr_names.index(elem_attr_name)

        values = self.__ex_get_one_attr('EX_ELEM_BLOCK', elem_blk_id, a_ndx + 1)
        if self.use_numpy:
            values = ctype_to_numpy(self, values)
        return values

    def put_elem_attr(self, elem_blk_id, elem_attrs):
        """
        store all attributes for each element in a block

        >>> exo.put_elem_attr(elem_blk_id, elem_attrs)

        Parameters
        ----------
        elem_blk_id : ex_entity_id
            element block *ID* (not *INDEX*)
        elem_attrs  : list<double>
            list of all attribute values for all
            elements in the block; the list
            cycles through all attributes of
            the first element, then all attributes
            of the second element, etc. Attributes
            are ordered by the ordering of the
            names returned by
            exo.get_attribute_names()
        """
        self.__ex_put_elem_attr(elem_blk_id, elem_attrs)

    def put_elem_attr_values(self, elem_blk_id, elem_attr_name, values):
        """
        store an attribute for each element in a block

        >>> exo.put_elem_attr_values(elem_blk_id, elem_attr_name, values)

        Parameters
        ----------
        elem_blk_id : ex_entity_id
            element block *ID* (not *INDEX*)
        elem_attr_name : string
            element attribute name
        values : list<double>
            list of values for a single attribute
            on a element block.  List dimensions
            should be 1 x N_elem, where N_elem is
            the number of elements on the element block.
        """
        # Determine index of requested attribute in attribute list
        elem_attr_names = self.get_attribute_names('EX_ELEM_BLOCK', elem_blk_id)
        a_ndx = elem_attr_names.index(elem_attr_name)
        self.__ex_put_one_attr('EX_ELEM_BLOCK', elem_blk_id, a_ndx + 1, values)

    def elem_type(self, object_id):
        """
        get the element type, e.g. "HEX8", for an element block

        >>> elem_type = exo.elem_type(elem_blk_id)

        Parameters
        ----------
        object_id : ex_entity_id
            element block *ID* (not *INDEX*)

        Returns
        -------
        elem_type : string
        """
        (elemType, _numElem, _nodesPerElem, _numAttr) = self.__ex_get_block('EX_ELEM_BLOCK', object_id)
        return elemType

    def num_attr(self, object_id):
        """
        get the number of attributes per element for an element block

        >>> num_elem_attrs = exo.num_attr(elem_blk_id)

        Parameters
        ----------
        object_id : ex_entity_id
            element block *ID* (not *INDEX*)

        Returns
        -------
        num_elem_attrs : int
        """
        (_elemType, _numElem, _nodesPerElem, numAttr) = self.__ex_get_block('EX_ELEM_BLOCK', object_id)
        return numAttr

    def num_elems_in_blk(self, object_id):
        """
        get the number of elements in an element block

        >>> num_blk_elems = exo.num_elems_in_blk(elem_blk_id)

        Parameters
        ----------
        object_id : ex_entity_id
            element block *ID* (not *INDEX*)

        Returns
        -------
        num_blk_elems : int
        """
        vals = self.get_entity_count('EX_ELEM_BLOCK', object_id)
        return vals

    def num_nodes_per_elem(self, object_id):
        """
        get the number of nodes per element for an element block

        >>> num_elem_nodes = exo.num_nodes_per_elem(elem_blk_id)

        Parameters
        ----------
        object_id : ex_entity_id
            element block *ID* (not *INDEX*)

        Returns
        -------
        num_elem_nodes : int
        """
        (_elemType, _numElem, nodesPerElem, _numAttr) = self.__ex_get_block('EX_ELEM_BLOCK', object_id)
        return nodesPerElem

    def get_element_variable_truth_table(self, entId=None):
        """
        See :py:func:`exodus.get_variable_truth_table`
        """
        return self.get_variable_truth_table('EX_ELEM_BLOCK', entId)

    def set_element_variable_truth_table(self, table):
        """
        See :py:func:`exodus.set_variable_truth_table`
        """
        return self.set_variable_truth_table('EX_ELEM_BLOCK', table)

    def get_element_variable_values(self, blockId, name, step):
        """
        get list of element variable values for a specified element
        block, element variable name, and time step

        >>> evar_vals = exo.get_element_variable_values(elem_blk_id,
        ...                                            evar_name, time_step)

        Parameters
        ----------
        elem_blk_id : ex_entity_id
            element block *ID* (not *INDEX*)
        evar_name : string
            name of element variable
        time_step : int
            1-based index of time step

        Returns
        -------

            if array_type == 'ctype':
              <list<ctypes.c_double>>  evar_vals

            if array_type == 'numpy':
              <np_array<double>>  evar_vals
        """
        return self.get_variable_values('EX_ELEM_BLOCK', blockId, name, step)

    def get_partial_element_variable_values(self, blockId, name, step, start_index, num_elements):
        """
        get list of element variable values for a specified element
        block, element variable name, and time step

        >>> evar_vals = exo.get_element_variable_values(elem_blk_id,
        ...                                            evar_name, time_step)

        Parameters
        ----------
        blockid : ex_entity_id
            element block *ID* (not *INDEX*)
        name : string
            name of element variable
        step : int
            1-based index of time step
        start_index: int
            1-based index of element in block to start returning data
        num_elements: int
            number of elements to return data for.

        Returns
        -------

            if array_type == 'ctype':
              <list<ctypes.c_double>>  evar_vals

            if array_type == 'numpy':
              <np_array<double>>  evar_vals
        """
        names = self.get_variable_names('EX_ELEM_BLOCK')
        var_id = names.index(name) + 1
        values = self.__ex_get_partial_var(step, 'EX_ELEM_BLOCK', var_id, blockId, start_index, num_elements)
        if self.use_numpy:
            values = ctype_to_numpy(self, values)
        return values

    # --------------------------------------------------------------------

    def put_element_variable_values(self, blockId, name, step, values):
        """
        store a list of element variable values for a specified element
        block, element variable name, and time step

        >>> status = exo.put_element_variable_values(elem_blk_id,
        ...             evar_name, time_step, evar_vals)

        Parameters
        ----------
        blockid : ex_entity_id
            element block *ID* (not *INDEX*)
        name : string
            name of element variable
        step : int
            1-based index of time step
        values : list<double>

        Returns
        -------
        status : bool
            True = successful execution
        """
        self.put_variable_values('EX_ELEM_BLOCK', blockId, name, step, values)
        return True

    # --------------------------------------------------------------------

    def get_element_variable_number(self):
        """
        get the number of element variables in the model

        >>> num_evars = exo.get_element_variable_number()

        Returns
        -------
        num_evars : int
        """
        return self.__ex_get_variable_param('EX_ELEM_BLOCK').value

    # --------------------------------------------------------------------

    def set_element_variable_number(self, number):
        """
        update the number of element variables in the model

        >>> status = exo.set_element_variable_number(num_evars)

        Parameters
        ----------
        number : int

        Returns
        -------
        status : bool
            True = successful execution
        """
        self.__ex_put_variable_param('EX_ELEM_BLOCK', number)
        return True

    # --------------------------------------------------------------------

    def get_element_variable_names(self):
        """
        get the list of element variable names in the model

        >>> evar_names = exo.get_element_variable_names()

        Returns
        -------
        evar_names : list<string>
        """
        if self.__ex_get_variable_param('EX_ELEM_BLOCK').value == 0:
            return []
        return self.__ex_get_variable_names('EX_ELEM_BLOCK')

    # --------------------------------------------------------------------

    def put_element_variable_name(self, name, index):
        """
        add the name and index of a new element variable to the model;
        element variable indexing goes from 1 to
        exo.get_element_variable_number()

        >>> status = exo.put_element_variable_name(evar_name, evar_index)

        Parameters
        ----------
            name : string
               name of new element variable
            index : int
               1-based index of new element variable

        Returns
        -------
        status : bool
            True = successful execution

        Example
        -------
        this method is often called within the following sequence:

        >>> num_evars = exo.get_element_variable_number()
        >>> new_evar_index = num_evars + 1
        >>> num_evars += 1
        >>> exo.set_element_variable_number(num_evars)
        >>> exo.put_element_variable_name("new_evar", new_evar_index)
        """
        EBvarNames = self.get_variable_names('EX_ELEM_BLOCK')
        if name in EBvarNames:
            print(f'WARNING: element variable \"{name}\" already exists.')
        if index > len(EBvarNames):
            print(("index", index, "len", len(EBvarNames)))
            raise Exception("ERROR: variable index out of range.")
        self.__ex_put_variable_name('EX_ELEM_BLOCK', index, name)
        return True

    # --------------------------------------------------------------------

    def get_attribute_names(self, objType, blkId):
        """
        get the list of attribute names for the specified block/set

        >>> attr_names = exo.get_attribute_names('EX_ELEM_BLOCK', elem_blk_id)

        Parameters
        ----------
        objType : ex_entity_type
            entity type
        blkid : ex_entity_id
            block/set *ID* (not *INDEX*)

        Returns
        -------
        attr_names : list<string>
        """
        names = self.__ex_get_attr_names(objType, blkId)
        return list(names)

    # --------------------------------------------------------------------

    def get_element_attribute_names(self, blkId):
        """
        get the list of element attribute names for a block

        >>> attr_names = exo.get_element_attribute_names(elem_blk_id)

        Parameters
        ----------
        blkid : ex_entity_id
            block/set *ID* (not *INDEX*)

        Returns
        -------
        attr_names : list<string>
        """
        names = self.__ex_get_attr_names('EX_ELEM_BLOCK', blkId)
        return list(names)

    # --------------------------------------------------------------------

    def put_attribute_names(self, objType, blkId, names):
        """
        store the list of element attribute names for a block

        >>> status = exo.put_attribute_names('EX_ELEM_BLOCK', elem_blk_id, attr_names)

        Parameters
        ----------
        objType:
            entity type
        blkid : ex_entity_id
            block/set  *ID* (not *INDEX*)
        attr_names : list<string>

        Returns
        -------
        status : bool
            True = successful execution
        """
        return self.__ex_put_attr_names(objType, blkId, names)

    # --------------------------------------------------------------------

    def put_element_attribute_names(self, blkId, names):
        """
        store the list of element attribute names for a block

        >>> status = exo.put_attribute_names(elem_blk_id, attr_names)

        Parameters
        ----------
        blkid : ex_entity_id
            block/set *ID* (not *INDEX*)
        names : list<string>

        Returns
        -------
        status : bool
            True = successful execution
        """
        return self.__ex_put_attr_names('EX_ELEM_BLOCK', blkId, names)

    # --------------------------------------------------------------------

    def get_element_property_names(self):
        """
        get the list of element property names for all element blocks
        in the model

        >>> eprop_names = exo.get_element_property_names()

        Returns
        -------
        eprop_names : list<string>
        """
        names = self.__ex_get_prop_names('EX_ELEM_BLOCK', 'EX_INQ_EB_PROP')
        return list(names)

    # --------------------------------------------------------------------

    def get_element_property_value(self, object_id, name):
        """
        get element property value (an integer) for a specified element
        block and element property name

        >>> eprop_val = exo.get_element_property_value(elem_blk_id, eprop_name)

        Parameters
        ----------
        elem_blk_id : ex_entity_id
            element block *ID* (not *INDEX*)
        name : string

        Returns
        -------
        eprop_val : int
        """
        propVal = self.__ex_get_prop('EX_ELEM_BLOCK', object_id, name)
        return int(propVal)

    # --------------------------------------------------------------------

    def put_element_property_value(self, object_id, name, value):
        """
        store an element property name and its integer value for an
        element block

        >>> status = exo.put_element_property_value(elem_blk_id,
        ...                                         eprop_name, eprop_val)


        Parameters
        ----------
        elem_blk_id : ex_entity_id
            element block *ID* (not *INDEX*)
        eprop_name : string

        eprop_val : int

        Returns
        -------
        status : bool
            True = successful execution
        """
        return self.__ex_put_prop('EX_ELEM_BLOCK', object_id, name, value)

    # --------------------------------------------------------------------

    #
    # nodesets
    #
    # --------------------------------------------------------------------

    def num_node_sets(self):
        """
        get the number of node sets in the model

        >>> num_node_sets = exo.num_node_sets()

        Returns
        -------
        num_node_sets : int
        """
        return self.numNodeSets.value

    # --------------------------------------------------------------------

    def get_node_set_ids(self):
        """
        get mapping of exodus node set index to user- or application-
        defined node set id; node_set_ids is ordered
        by the *INDEX* ordering, a 1-based system going from
        1 to exo.num_node_sets(), used by exodus for storage
        and input/output of array data stored on the node sets; a
        user or application can optionally use a separate node set
        *ID* numbering system, so the node_set_ids array points to the
        node set *ID* for each node set *INDEX*

        >>> node_set_ids = exo.get_ids('EX_NODE_SET')

        Returns
        -------

            if array_type == 'ctype':
              <list<int>>  node_set_ids

            if array_type == 'numpy':
              <np_array<int>>  node_set_ids
        """
        return self.get_ids('EX_NODE_SET')

    # --------------------------------------------------------------------

    def get_node_set_name(self, object_id):
        """
        get the name of a node set

        >>> node_set_name = exo.get_node_set_name(node_set_id)

        Parameters
        ----------
        node_set_id : ex_entity_id
          node set *ID* (not *INDEX*)

        Returns
        -------
        node_set_name : string
        """
        return self.__ex_get_name('EX_NODE_SET', object_id)

    # --------------------------------------------------------------------

    def put_node_set_name(self, object_id, name):
        """
        store the name of a node set

        >>> exo.put_node_set_name(node_set_id, node_set_name)

        Parameters
        ----------
            node_set_id : ex_entity_id
               node set *ID* (not *INDEX*)
            node_set_name : string
        """
        self.__ex_put_name('EX_NODE_SET', object_id, name)

    # --------------------------------------------------------------------

    def get_node_set_names(self):
        """
        get a list of all node set names ordered by node set *INDEX*;
        (see :py:func:`exodus.get_ids` for explanation of the
        difference between node set *ID* and node set *INDEX*)

        >>> node_set_names = exo.get_node_set_names()

        Returns
        -------
        node_set_names : list<string>
        """
        return self.__ex_get_names('EX_NODE_SET')

    # --------------------------------------------------------------------

    def put_node_set_names(self, names):
        """
        store a list of all node set names ordered by node set *INDEX*;
        (see :py:func:`exodus.get_ids` for explanation of the
        difference between node set *ID* and node set *INDEX*)

        >>> exo.put_node_set_names(node_set_names)

        Parameters
        ----------
            names : list<string>
        """
        self.__ex_put_names('EX_NODE_SET', names)

    # --------------------------------------------------------------------

    def num_nodes_in_node_set(self, object_id):
        """
        get the number of nodes in a node set

        >>> num_ns_nodes = exo.num_nodes_in_node_set(node_set_id)

        Parameters
        ----------
        node_set_id  : ex_entity_id
           node set *ID* (not *INDEX*)

        Returns
        -------
        num_ns_nodes : int
        """
        node_set_nodes = self.get_node_set_nodes(object_id)
        return len(node_set_nodes)

    # --------------------------------------------------------------------

    def get_node_set_nodes(self, object_id):
        """
        get the list of node *INDICES* in a node set
        (see :py:func:`exodus.get_id_map` for explanation of node *INDEX*
        versus node *ID*)

        >>> ns_nodes = exo.get_node_set_nodes(node_set_id)

        Parameters
        ----------
        node_set_id : ex_entity_id
           node set *ID* (not *INDEX*)

        Returns
        -------

            if array_type == 'ctype':
              <list<int>>  ns_nodes

            if array_type == 'numpy':
              <np_array<int>>  ns_nodes
        """
        node_set_ids = self.get_ids('EX_NODE_SET')
        assert object_id in node_set_ids
        node_set_nodes = self.__ex_get_node_set(object_id)
        node_set_nodes = list(node_set_nodes)
        if self.use_numpy:
            node_set_nodes = self.np.array(node_set_nodes)
        return node_set_nodes

    # --------------------------------------------------------------------

    def put_node_set(self, object_id, nodeSetNodes):
        """
        store a node set by its id and the list of node *INDICES* in
        the node set (see :py:func:`exodus.get_id_map` for explanation of node
        *INDEX* versus node *ID*)

        >>> exo.put_node_set(node_set_id, ns_nodes)

        Parameters
        ----------
        node_set_id : ex_entity_id
           node set *ID* (not *INDEX*)
        ns_nodes : list<int>
        """
        self.__ex_put_node_set(object_id, nodeSetNodes)

    # --------------------------------------------------------------------

    def get_node_set_dist_facts(self, object_id):
        """
        get the list of distribution factors for nodes in a node set

        >>> ns_dist_facts = exo.get_node_set_dist_facts(node_set_id)

        Parameters
        ----------
        node_set_id : ex_entity_id
           node set *ID* (not *INDEX*)

        Returns
        -------

            if array_type == 'ctype':
              <list<double>>  ns_dist_facts  a list of distribution factors,
                e.g. nodal 'weights'

            if array_type == 'numpy':
              <np_array<double>>  ns_dist_facts  a list of distribution
                factors, e.g. nodal
                'weights'
        """
        node_set_dfs = self.__ex_get_node_set_dist_fact(object_id)
        node_set_dfs = list(node_set_dfs)
        if self.use_numpy:
            node_set_dfs = self.np.array(node_set_dfs)
        return node_set_dfs

    # --------------------------------------------------------------------

    def put_node_set_dist_fact(self, object_id, nodeSetDistFact):
        """
        store the list of distribution factors for nodes in a node set

        >>> exo.put_node_set_dist_fact(node_set_id, ns_dist_facts)

        Parameters
        ----------
        object_id : ex_entity_id
            node set *ID* (not *INDEX*)
        nodeSetDistFact : list<double>
            a list of distribution factors, e.g. nodal 'weights'
        """
        self.__ex_put_node_set_dist_fact(object_id, nodeSetDistFact)

    # --------------------------------------------------------------------

    def get_node_set_variable_number(self):
        """
        get the number of node set variables in the model

        >>> num_nsvars = exo.get_node_set_variable_number()

        Returns
        -------
        num_nsvars : int
        """
        return self.__ex_get_variable_param('EX_NODE_SET').value

    # --------------------------------------------------------------------

    def set_node_set_variable_number(self, number):
        """
        update the number of node set variables in the model

        >>> status = exo.set_node_set_variable_number(num_nsvars)

        Parameters
        ----------
        number : int

        Returns
        -------
        status : bool
            True = successful execution
        """
        self.__ex_put_variable_param('EX_NODE_SET', number)
        return True

    # --------------------------------------------------------------------

    def get_node_set_variable_truth_table(self, entId=None):
        """
        See :py:func:`exodus.get_variable_truth_table`
        """
        return self.get_variable_truth_table('EX_NODE_SET', entId)

    # --------------------------------------------------------------------

    def set_node_set_variable_truth_table(self, table):
        """
        See :py:func:`exodus.set_variable_truth_table`
        """
        return self.set_variable_truth_table('EX_NODE_SET', table)

    # --------------------------------------------------------------------

    def get_node_set_variable_names(self):
        """
        get the list of node set variable names in the model

        >>> nsvar_names = exo.get_node_set_variable_names()

        Returns
        -------
        nsvar_names : list<string>
        """
        if self.__ex_get_variable_param('EX_NODE_SET').value == 0:
            return []
        return self.__ex_get_variable_names('EX_NODE_SET')

    # --------------------------------------------------------------------

    def put_node_set_variable_name(self, name, index):
        """
        add the name and index of a new node set variable to the model;
        node set variable indexing goes from 1 to
        exo.get_node_set_variable_number()

        >>> status = exo.put_node_set_variable_name(nsvar_name, nsvar_index)

        Parameters
        ----------
        name   : string
            name of new node set variable
        index  : int
            1-based index of new node set variable

        Returns
        -------
        status : bool
            True = successful execution

        Example
        -------
        this method is often called within the following sequence:

        >>> num_nsvars = exo.get_node_set_variable_number()
        >>> new_nsvar_index = num_nsvars + 1
        >>> num_nsvars += 1
        >>> exo.set_node_set_variable_number(num_nsvars)
        >>> exo.put_node_set_variable_name("new_nsvar", new_nsvar_index)
        """
        NSvarNames = self.get_variable_names('EX_NODE_SET')
        if name in NSvarNames:
            print(f'WARNING: Node set variable \"{name}\" already exists.')
        if index > len(NSvarNames):
            raise Exception("ERROR: variable index out of range.")
        self.__ex_put_variable_name('EX_NODE_SET', index, name)
        return True

    # --------------------------------------------------------------------

    def get_node_set_variable_values(self, object_id, name, step):
        """
        get list of node set variable values for a specified node
        set, node set variable name, and time step; the list has
        one variable value per node in the set

        >>> nsvar_vals =
        ...   exo.get_node_set_variable_values(node_set_id,
        ...    nsvar_name, time_step)

        Parameters
        ----------
        node_set_id : ex_entity_id
             node set *ID* (not *INDEX*)
        nsvar_name  : string
             name of node set variable
        time_step   : int
             1-based index of time step

        Returns
        -------

            if array_type == 'ctype':
              <list<ctypes.c_double>>  nsvar_vals

            if array_type == 'numpy':
              <np_array<double>>  nsvar_vals
        """
        return self.get_variable_values('EX_NODE_SET', object_id, name, step)

    # --------------------------------------------------------------------

    def get_partial_node_set_variable_values(self, object_id, name, step, start_index, num_nodes):
        """
        get list of node set variable values for a specified node
        set, node set variable name, and time step; the list has
        one variable value per node in the set

        >>> nsvar_vals =
        ...   exo.get_node_set_variable_values(node_set_id,
        ...    nsvar_name, time_step)

        Parameters
        ----------
        node_set_id : ex_entity_id
           node set *ID* (not *INDEX*)
        nsvar_name  : string
           name of node set variable
        time_step   : int
           1-based index of time step
        start_index : int
           1-based index of node to start returning data
        num_nodes   : int
           number of nodes to return data for.

        Returns
        -------

            if array_type == 'ctype':
              <list<ctypes.c_double>>  nsvar_vals

            if array_type == 'numpy':
              <np_array<double>>  nsvar_vals
        """
        names = self.get_variable_names('EX_NODE_SET')
        var_id = names.index(name) + 1
        values = self.__ex_get_partial_var(step, 'EX_NODE_SET', var_id, object_id, start_index, num_nodes)
        if self.use_numpy:
            values = ctype_to_numpy(self, values)
        return values

    # --------------------------------------------------------------------

    def put_node_set_variable_values(self, object_id, name, step, values):
        """
        store a list of node set variable values for a specified node
        set, node set variable name, and time step; the list has one
        variable value per node in the set

        >>> status =
        ... exo.put_node_set_variable_values(node_set_id,
        ...     nsvar_name, time_step, nsvar_vals)

        Parameters
        ----------
        node_set_id : ex_entity_id
           node set *ID* (not *INDEX*)
        nsvar_name  : string
           name of node set variable
        time_step   : int
           1-based index of time step
        nsvar_vals  : list<double>

        Returns
        -------
        status : bool
            True = successful execution
        """
        self.put_variable_values('EX_NODE_SET', object_id, name, step, values)
        return True

    # --------------------------------------------------------------------

    def get_all_node_set_params(self):
        """
        get total number of nodes and distribution factors (e.g. nodal
        'weights') combined among all node sets

        >>> tot_num_ns_nodes,
        ... tot_num_ns_dist_facts = exo.get_all_node_set_params()

        Returns
        -------
        tot_num_ns_nodes : int
        tot_num_ns_dist_facts : int
        """
        nodeSetIds = self.__ex_get_ids('EX_NODE_SET')
        totNumSetNodes, totNumSetDistFacts = 0, 0
        for nodeSetId in nodeSetIds:
            (numSetNodes, numSetDistFacts) = self.__ex_get_set_param('EX_NODE_SET', nodeSetId)
            totNumSetNodes += numSetNodes
            totNumSetDistFacts += numSetDistFacts
        return totNumSetNodes, totNumSetDistFacts

    # --------------------------------------------------------------------

    def get_set_params(self, object_type, object_id):
        """
        get number of entities and distribution factors (e.g. nodal
        'weights') in the specified set

        >>> num_ns_nodes, num_ns_dist_facts =
        ...     exo.get_set_params('EX_NODE_SET', node_set_id)

        Parameters
        ----------
        set_id : ex_entity_id
            set *ID* (not *INDEX*)

        Returns
        -------
        num_set_entities : int
        num_set_dist_facts : int
        """
        (numSetEntities, numSetDistFacts) = self.__ex_get_set_param(object_type, object_id)
        return numSetEntities, numSetDistFacts

    # --------------------------------------------------------------------

    def put_set_params(self, object_type, object_id, numSetEntity, numSetDistFacts=None):
        """
        initialize a new set of the specified type

        >>> exo.put_set_params('EX_NODE_SET', node_set_id,
        ...                 num_ns_nodes, num_ns_dist_facts)

        Parameters
        ----------
        set_id : ex_entity_id
            set *ID* (not *INDEX*)
        num_set_entity : int
            number of nodes/edges/faces/elements to be added to set
        num_dist_facts : int, optional
            number of distribution factors (e.g. nodal 'weights') --
            must be equal to zero or num_set_entity
        """
        if numSetDistFacts is None:
            numSetDistFacts = numSetEntity
        assert numSetDistFacts in (0, numSetEntity)
        self.__ex_put_set_param(object_type, object_id, numSetEntity, numSetDistFacts)

    # --------------------------------------------------------------------

    def get_node_set_params(self, object_id):
        """ See :py:func:`exodus.put_set_params` """

        (numSetNodes, numSetDistFacts) = self.__ex_get_set_param('EX_NODE_SET', object_id)
        return numSetNodes, numSetDistFacts

    # --------------------------------------------------------------------

    def put_node_set_params(self, object_id, numSetNodes, numSetDistFacts=None):
        """ See :py:func:`exodus.put_set_params` """
        if numSetDistFacts is None:
            numSetDistFacts = numSetNodes
        assert numSetDistFacts in (0, numSetNodes)
        self.__ex_put_set_param('EX_NODE_SET', object_id, numSetNodes, numSetDistFacts)

    # --------------------------------------------------------------------

    def get_node_set_property_names(self):
        """
        get the list of node set property names for all node sets in
        the model

        >>> nsprop_names = exo.get_node_set_property_names()

        Returns
        -------
        nsprop_names : list<string>
        """
        names = self.__ex_get_prop_names('EX_NODE_SET', 'EX_INQ_NS_PROP')
        return list(names)

    # --------------------------------------------------------------------

    def get_node_set_property_value(self, object_id, name):
        """
        get node set property value (an integer) for a specified node
        set and node set property name

        >>> nsprop_val = exo.get_node_set_property_value(node_set_id, nsprop_name)

        Parameters
        ----------
        node_set_id : ex_entity_id
          node set *ID* (not *INDEX*)
        name : string

        Returns
        -------
        nsprop_val : int
        """
        propVal = self.__ex_get_prop('EX_NODE_SET', object_id, name)
        return int(propVal)

    # --------------------------------------------------------------------

    def put_node_set_property_value(self, object_id, name, value):
        """
        store a node set property name and its integer value for a
        node set

        >>> status = exo.put_node_set_property_value(node_set_id,
        ...                   nsprop_name, nsprop_val)

        Parameters
        ----------
        node_set_id : ex_entity_id
            node set *ID* (not *INDEX*)
        name : string
        value : int

        Returns
        -------
        status : bool
            True = successful execution
        """
        return self.__ex_put_prop('EX_NODE_SET', object_id, name, value)

    #
    # sidesets
    #
    # --------------------------------------------------------------------

    def num_side_sets(self):
        """
        get the number of side sets in the model

        >>> num_side_sets = exo.num_side_sets()

        Returns
        -------
        num_side_sets : int
        """
        return self.numSideSets.value

    # --------------------------------------------------------------------

    def get_side_set_ids(self):
        """
        get mapping of exodus side set index to user- or application-
        defined side set id; side_set_ids is ordered
        by the *INDEX* ordering, a 1-based system going from
        1 to exo.num_side_sets(), used by exodus for storage
        and input/output of array data stored on the side sets; a
        user or application can optionally use a separate side set
        *ID* numbering system, so the side_set_ids array points to the
        side set *ID* for each side set *INDEX*

        >>> side_set_ids = exo.get_ids('EX_SIDE_SET')

        Returns
        -------

            if array_type == 'ctype':
              <list<int>>  side_set_ids

            if array_type == 'numpy':
              <np_array<int>>  side_set_ids
        """
        return self.get_ids('EX_SIDE_SET')

    # --------------------------------------------------------------------

    def get_side_set_name(self, object_id):
        """
        get the name of a side set

        >>> side_set_name = exo.get_side_set_name(side_set_id)

        Parameters
        ----------
        object_id : ex_entity_id
           side set *ID* (not *INDEX*)

        Returns
        -------
        side_set_name : string
        """
        return self.__ex_get_name('EX_SIDE_SET', object_id)

    # --------------------------------------------------------------------

    def put_side_set_name(self, object_id, name):
        """
        store the name of a side set

        >>> exo.put_side_set_name(side_set_id, side_set_name)

        Parameters
        ----------
        side_set_id : ex_entity_id
           side set *ID* (not *INDEX*)
        side_set_name : string
        """
        self.__ex_put_name('EX_SIDE_SET', object_id, name)

    # --------------------------------------------------------------------

    def get_side_set_names(self):
        """
        get a list of all side set names ordered by side set *INDEX*;
        (see :py:func:`exodus.get_ids` for explanation of the
        difference between side set *ID* and side set *INDEX*)

        >>> side_set_names = exo.get_side_set_names()

        Returns
        -------
        side_set_names : list<string>
        """
        return self.__ex_get_names('EX_SIDE_SET')

    # --------------------------------------------------------------------

    def put_side_set_names(self, names):
        """
        store a list of all side set names ordered by side set *INDEX*;
        (see :py:func:`exodus.get_ids` for explanation of the
        difference between side set *ID* and side set *INDEX*)

        >>> exo.put_side_set_names(side_set_names)

        Parameters
        ----------
        side_set_names : list<string>
        """
        self.__ex_put_names('EX_SIDE_SET', names)

    # --------------------------------------------------------------------

    def num_faces_in_side_set(self, object_id):
        """
        get the number of faces in a side set

        >>> num_ss_faces = exo.num_faces_in_side_set(side_set_id)

        Parameters
        ----------
        side_set_id : ex_entity_id
            side set *ID* (not *INDEX*)

        Returns
        -------
        num_ss_faces : int
        """
        ssids = self.get_ids('EX_SIDE_SET')
        if object_id not in ssids:
            print("WARNING: queried side set ID does not exist in database")
            return 0
        (num_side_in_set, _num_dist_fact_in_set) = self.__ex_get_set_param('EX_SIDE_SET', object_id)
        return num_side_in_set

    # --------------------------------------------------------------------

    def get_all_side_set_params(self):
        """
        get total number of sides, nodes, and distribution factors
        (e.g. nodal 'weights') combined among all side sets

        >>> tot_num_ss_sides, tot_num_ss_nodes, tot_num_ss_dist_facts =
        ...          exo.get_all_side_set_params()

        Returns
        -------
        tot_num_ss_sides : int
        tot_num_ss_nodes : int
        tot_num_ss_dist_facts : int

        Note
        -----
        The number of nodes (and distribution factors) in a side set is
        the sum of all face nodes.  A single node can be counted more
        than once, i.e. once for each face it belongs to in the side set.
        """
        ids = self.__ex_get_ids('EX_SIDE_SET')
        totNumSetSides, totNumSetDistFacts = 0, 0  # totNumSetDistFacts = totNumSetNodes
        for sideSetId in ids:
            (numSetSides, numSetDistFacts) = self.__ex_get_set_param('EX_SIDE_SET', sideSetId)
            totNumSetSides += numSetSides
            totNumSetDistFacts += numSetDistFacts
        totNumSetNodes = totNumSetDistFacts
        return totNumSetSides, totNumSetNodes, totNumSetDistFacts

    # --------------------------------------------------------------------

    def get_side_set_params(self, object_id):
        """
        get number of sides and nodal distribution factors (e.g. nodal
        'weights') in a side set

        >>> num_ss_sides, num_ss_dist_facts = exo.get_side_set_params(side_set_id)

        Parameters
        ----------
        side_set_id : ex_entity_id
             side set *ID* (not *INDEX*)

        Returns
        -------
        num_ss_sides : int
        num_ss_dist_facts : int

        Note
        -----
        The number of nodes (and distribution factors) in a side set is
        the sum of all face nodes.  A single node can be counted more
        than once, i.e. once for each face it belongs to in the side set.
        """
        (numSetSides, numSetDistFacts) = self.__ex_get_set_param('EX_SIDE_SET', object_id)
        return numSetSides, numSetDistFacts

    # --------------------------------------------------------------------

    def put_side_set_params(self, object_id, numSetSides, numSetDistFacts):
        """
        initialize a new side set

        >>> exo.put_side_set_params(side_set_id, num_ss_sides, num_ss_dist_facts)

        Parameters
        ----------
        side_set_id : ex_entity_id
            side set *ID* (not *INDEX*)
        num_ss_sides : int
            number of sides to be added to set
        num_ss_dist_facts : int
            number of nodal distribution factors (e.g. nodal 'weights')

        Note
        -----
        The number of nodes (and distribution factors) in a side set is
        the sum of all face nodes.  A single node can be counted more
        than once, i.e. once for each face it belongs to in the side set.
        """
        self.__ex_put_set_param('EX_SIDE_SET', object_id, numSetSides, numSetDistFacts)

    # --------------------------------------------------------------------

    def get_side_set(self, object_id):
        """
        get the lists of element and side indices in a side set; the
        two lists correspond: together, ss_elems[i] and ss_sides[i]
        define the face of an element

        >>> ss_elems, ss_sides = exo.get_side_set(side_set_id)

        Parameters
        ----------
        side_set_id : ex_entity_id
            side set *ID* (not *INDEX*)

        Returns
        -------

            if array_type == 'ctype':
              <list<int>>  ss_elems
              <list<int>>  ss_sides

            if array_type == 'numpy':
              <np_array<int>>  ss_elems
              <np_array<int>>  ss_sides
        """
        (side_set_elem_list, side_set_side_list) = self.__ex_get_side_set(object_id)
        if self.use_numpy:
            side_set_elem_list = ctype_to_numpy(self, side_set_elem_list)
            side_set_side_list = ctype_to_numpy(self, side_set_side_list)
        return side_set_elem_list, side_set_side_list

    # --------------------------------------------------------------------

    def put_side_set(self, object_id, sideSetElements, sideSetSides):
        """
        store a side set by its id and the lists of element and side
        indices in the side set; the two lists correspond: together,
        ss_elems[i] and ss_sides[i] define the face of an element

        >>> exo.put_side_set(side_set_id, ss_elems, ss_sides)

        Parameters
        ----------
        side_set_id : ex_entity_id
            side set *ID* (not *INDEX*)
        ss_elems : list<int>
        ss_sides : list<int>
        """
        self.__ex_put_side_set(object_id, sideSetElements, sideSetSides)

    # --------------------------------------------------------------------

    def get_side_set_dist_fact(self, object_id):
        """
        get the list of distribution factors for nodes in a side set

        >>> ss_dist_facts = exo.get_side_set_dist_fact(side_set_id)

        Parameters
        ----------
        side_set_id : ex_entity_id
             side set *ID* (not *INDEX*)

        Returns
        -------

            if array_type == 'ctype':
              <list<double>>  ss_dist_facts  a list of distribution factors,
                e.g. nodal 'weights'

            if array_type == 'numpy':
              <np_array<double>>  ss_dist_facts  a list of distribution
                factors, e.g. nodal
                'weights'

        Note
        -----
        The number of nodes (and distribution factors) in a side set is
        the sum of all face nodes.  A single node can be counted more
        than once, i.e. once for each face it belongs to in the side set.
        """
        side_set_dfs = list(self.__ex_get_side_set_dist_fact(object_id))
        if self.use_numpy:
            side_set_dfs = self.np.array(side_set_dfs)
        return side_set_dfs

    # --------------------------------------------------------------------

    def put_side_set_dist_fact(self, object_id, sideSetDistFact):
        """
        store the list of distribution factors for nodes in a side set

        >>> exo.put_side_set_dist_fact(side_set_id, ss_dist_facts)

        Parameters
        ----------
        object_id : ex_entity_id
            node set *ID* (not *INDEX*)
        sideSetDistFact : list<double>
            a list of distribution factors, e.g. nodal 'weights'

        Note
        -----
        The number of nodes (and distribution factors) in a side set is
        the sum of all face nodes.  A single node can be counted more
        than once, i.e. once for each face it belongs to in the side set.
        """
        self.__ex_put_side_set_dist_fact(object_id, sideSetDistFact)

    # --------------------------------------------------------------------

    def get_side_set_node_list(self, object_id):
        """
        get two lists:
         1. number of nodes for each side in the set
         2. concatenation of the nodes for each side in the set

        >>> ss_num_nodes_per_side, ss_nodes = exo.get_side_set_node_list(side_set_id)

        Parameters
        ----------
        side_set_id : ex_entity_id
             side set *ID* (not *INDEX*)

        Returns
        -------

            if array_type == 'ctype':
              <list<int>>  ss_num_side_nodes
              <list<int>>  ss_nodes

            if array_type == 'numpy':
              <np_array<int>>  ss_num_side_nodes
              <np_array<int>>  ss_nodes

        Note
        -----
        The number of nodes (and distribution factors) in a side set is
        the sum of all face nodes.  A single node can be counted more
        than once, i.e. once for each face it belongs to in the side set.
        """
        (side_set_node_cnt_list,
         side_set_node_list) = self.__ex_get_side_set_node_list(object_id)
        if self.use_numpy:
            side_set_node_cnt_list = ctype_to_numpy(
                self, side_set_node_cnt_list)
            side_set_node_list = ctype_to_numpy(self, side_set_node_list)
        return side_set_node_cnt_list, side_set_node_list

    # --------------------------------------------------------------------

    def get_side_set_variable_truth_table(self, entId=None):
        """
        See :py:func:`exodus.get_variable_truth_table`
        """
        return self.get_variable_truth_table('EX_SIDE_SET', entId)

    # --------------------------------------------------------------------

    def set_side_set_variable_truth_table(self, table):
        """
        See :py:func:`exodus.set_variable_truth_table`
        """
        return self.set_variable_truth_table('EX_SIDE_SET', table)

    # --------------------------------------------------------------------

    def get_side_set_variable_number(self):
        """
        get the number of side set variables in the model

        >>> num_ssvars = exo.get_side_set_variable_number()

        Returns
        -------
        num_ssvars : int
        """
        return self.__ex_get_variable_param('EX_SIDE_SET').value

    # --------------------------------------------------------------------

    def set_side_set_variable_number(self, number):
        """
        update the number of side set variables in the model

        >>> status = exo.set_side_set_variable_number(num_ssvars)

        Parameters
        ----------
        number : int

        Returns
        -------
        status : bool
            True = successful execution
        """
        self.__ex_put_variable_param('EX_SIDE_SET', number)
        return True

    # --------------------------------------------------------------------

    def get_side_set_variable_names(self):
        """
        get the list of side set variable names in the model

        >>> ssvar_names = exo.get_side_set_variable_names()

        Returns
        -------
        ssvar_names : list<string>
        """
        if self.__ex_get_variable_param('EX_SIDE_SET').value == 0:
            return []
        return self.__ex_get_variable_names('EX_SIDE_SET')

    # --------------------------------------------------------------------

    def put_side_set_variable_name(self, name, index):
        """
        add the name and index of a new side set variable to the model;
        side set variable indexing goes from 1 to
        exo.get_side_set_variable_number()

        >>> status = exo.put_side_set_variable_name(ssvar_name, ssvar_index)

        Parameters
        ----------
        name : string
           name of new side set variable
        index : int
           1-based index of new side set variable

        Returns
        -------
        status : bool
            True = successful execution

        Example
        -------
        this method is often called within the following sequence:

        >>> num_ssvars = exo.get_side_set_variable_number()
        >>> new_ssvar_index = num_ssvars + 1
        >>> num_ssvars += 1
        >>> exo.set_side_set_variable_number(num_ssvars)
        >>> exo.put_side_set_variable_name("new_ssvar", new_ssvar_index)

        """
        SSvarNames = self.get_variable_names('EX_SIDE_SET')
        if name in SSvarNames:
            print(f'WARNING: Side set variable \"{name}\" already exists.')
        if index > len(SSvarNames):
            raise Exception("ERROR: variable index out of range.")
        self.__ex_put_variable_name('EX_SIDE_SET', index, name)
        return True

    # --------------------------------------------------------------------

    def get_side_set_variable_values(self, object_id, name, step):
        """
        get list of side set variable values for a specified side
        set, side set variable name, and time step; the list has
        one variable value per side in the set

        >>> ssvar_vals = exo.get_side_set_variable_values(side_set_id,
        ...    ssvar_name, time_step)

        Parameters
        ----------
        side_set_id : ex_entity_id
            side set *ID* (not *INDEX*)
        ssvar_name  : string
            name of side set variable
        time_step   : int
            1-based index of time step

        Returns
        -------

            if array_type == 'ctype':
              <list<ctypes.c_double>>  ssvar_vals

            if array_type == 'numpy':
              <np_array<double>>  ssvar_vals
        """
        return self.get_variable_values('EX_SIDE_SET', object_id, name, step)

    # --------------------------------------------------------------------

    def get_partial_side_set_variable_values(self, object_id, name, step, start_index, num_sides):
        """
        get list of side set variable values for a specified side
        set, side set variable name, and time step; the list has
        one variable value per side in the set

        >>> ssvar_vals = exo.get_side_set_variable_values(side_set_id,
        ...    ssvar_name, time_step)

        Parameters
        ----------
        object_id : ex_entity_id
            side set *ID* (not *INDEX*)
        name  : string
            name of side set variable
        step   : int
            1-based index of time step
        start_index : int
            1-based index of side to start returning data
        num_sides : int
            number of sides to return data for.

        Returns
        -------

            if array_type == 'ctype':
              <list<ctypes.c_double>>  ssvar_vals

            if array_type == 'numpy':
              <np_array<double>>  ssvar_vals
        """
        names = self.get_variable_names('EX_SIDE_SET')
        var_id = names.index(name) + 1
        values = self.__ex_get_partial_var(step, 'EX_SIDE_SET', var_id, object_id, start_index, num_sides)
        if self.use_numpy:
            values = ctype_to_numpy(self, values)
        return values

    # --------------------------------------------------------------------

    def put_side_set_variable_values(self, object_id, name, step, values):
        """
        store a list of side set variable values for a specified side
        set, side set variable name, and time step; the list has one
        variable value per side in the set

        >>> status = exo.put_side_set_variable_values(side_set_id,
        ...              ssvar_name, time_step, ssvar_vals)

        Parameters
        ----------
        object_id  : ex_entity_id
            side set *ID* (not *INDEX*)
        name   : string
            name of side set variable
        step    : int
            1-based index of time step
        values   : list<double>

        Returns
        -------
        status : bool
            True = successful execution
        """
        self.put_variable_values('EX_SIDE_SET', object_id, name, step, values)
        return True

    # --------------------------------------------------------------------

    def get_side_set_property_names(self):
        """
        get the list of side set property names for all side sets in
        the model

        >>> ssprop_names = exo.get_side_set_property_names()

        Returns
        -------
        ssprop_names : list<string>
        """
        names = self.__ex_get_prop_names('EX_SIDE_SET', 'EX_INQ_SS_PROP')
        return list(names)

    # --------------------------------------------------------------------

    def get_side_set_property_value(self, object_id, name):
        """
        get side set property value (an integer) for a specified side
        set and side set property name

        >>> ssprop_val = exo.get_side_set_property_value(side_set_id, ssprop_name)

        Parameters
        ----------
        object_id  : ex_entity_id
            side set *ID* (not *INDEX*)
        name   : string
            name of side set property

        Returns
        -------
        ssprop_val : int
        """
        propVal = self.__ex_get_prop('EX_SIDE_SET', object_id, name)
        return int(propVal)

    # --------------------------------------------------------------------

    def put_side_set_property_value(self, object_id, name, value):
        """
        store a side set property name and its integer value for a
        side set

        >>> status = exo.put_side_set_property_value(side_set_id,
        ...               ssprop_name, ssprop_val)

        Parameters
        ----------
        object_id  : ex_entity_id
            side set *ID* (not *INDEX*)
        name   : string
            name of side set property
        value : int

        Returns
        -------
        status : bool
            True = successful execution
        """
        return self.__ex_put_prop('EX_SIDE_SET', object_id, name, value)

    #
    # global variables
    #
    # --------------------------------------------------------------------

    def get_global_variable_number(self):
        """
        get the number of global variables in the model

        >>> num_gvars = exo.get_global_variable_number()

        Returns
        -------
        num_gvars : int
        """
        return self.__ex_get_variable_param('EX_GLOBAL').value

    # --------------------------------------------------------------------

    def set_global_variable_number(self, number):
        """
        update the number of global variables in the model

        >>> status = exo.set_global_variable_number(num_gvars)

        Parameters
        ----------
        number : int

        Returns
        -------
        status : bool
            True = successful execution
        """
        self.__ex_put_variable_param('EX_GLOBAL', number)
        return True

    # --------------------------------------------------------------------

    def get_global_variable_names(self):
        """
        get the list of global variable names in the model

        >>> gvar_names = exo.get_global_variable_names()

        Returns
        -------
        gvar_names : list<string>
        """
        if self.get_variable_number('EX_GLOBAL') == 0:
            return []
        return self.__ex_get_variable_names('EX_GLOBAL')

    # --------------------------------------------------------------------

    def put_global_variable_name(self, name, index):
        """
        add the name and index of a new global variable to the model;
        global variable indexing goes from 1 to
        exo.get_global_variable_number()

        >>> status = exo.put_global_variable_name(gvar_name, gvar_index)

        Parameters
        ----------
        name : string
           name of new global variable
        index : int
           1-based index of new global variable

        Returns
        -------
        status : bool
            True = successful execution

        Example
        -------
        this method is often called within the following sequence:

        >>> num_gvars = exo.get_global_variable_number()
        >>> new_gvar_index = num_gvars + 1
        >>> num_gvars += 1
        >>> exo.set_global_variable_number(num_gvars)
        >>> exo.put_global_variable_name("new_gvar", new_gvar_index)
        """
        GlobVarNames = self.get_variable_names('EX_GLOBAL')
        if name in GlobVarNames:
            print(f'WARNING: Global variable \"{name}\" already exists.')
        if index > len(GlobVarNames):
            print(("index", index, "len", len(GlobVarNames)))
            raise Exception("ERROR: variable index out of range.")
        self.__ex_put_variable_name('EX_GLOBAL', index, name)
        return True

    # --------------------------------------------------------------------

    def get_global_variable_value(self, name, step):
        """
        get a global variable value for a specified global variable
        name and time step

        >>> gvar_val = exo.get_global_variable_value(gvar_name, time_step)

        Parameters
        ----------
        name : string
           name of global variable
        step : int
           1-based index of time step

        Returns
        -------
        gvar_val : double
        """
        names = self.get_variable_names('EX_GLOBAL')
        var_id = names.index(name)
        num = self.__ex_get_variable_param('EX_GLOBAL')
        gvalues = self.__ex_get_var(step, 'EX_GLOBAL', 0, 1, num.value)
        return gvalues[var_id]

    # --------------------------------------------------------------------

    def get_all_global_variable_values(self, step):
        """
        get all global variable values (one for each global variable
        name, and in the order given by exo.get_global_variable_names())
        at a specified time step

        >>> gvar_vals = exo.get_all_global_variable_values(time_step)

        Parameters
        ----------
        step : int
           1-based index of time step

        Returns
        -------

            if array_type == 'ctype':
              <list<double>>  gvar_vals

            if array_type == 'numpy':
              <np_array<double>>  gvar_vals
        """
        num = self.__ex_get_variable_param('EX_GLOBAL')
        gvalues = self.__ex_get_var(step, 'EX_GLOBAL', 0, 1, num.value)
        values = [gvalues[i] for i in range(num.value)]
        if self.use_numpy:
            values = self.np.array(values)
        return values

    # --------------------------------------------------------------------

    def put_global_variable_value(self, name, step, value):
        """
        store a global variable value for a specified global variable
        name and time step

        >>> status = exo.put_global_variable_value(gvar_name, time_step, gvar_val)

        Parameters
        ----------
        name : string
            name of global variable
        step : int
            1-based index of time step
        value  : double

        Returns
        -------
        status : bool
            True = successful execution
        """
        # we must write all values at once, not individually
        names = self.get_variable_names('EX_GLOBAL')
        # get all values
        numVals = self.get_variable_number('EX_GLOBAL')
        values = (ctypes.c_double * numVals)()
        for i in range(numVals):
            values[i] = ctypes.c_double(
                self.get_global_variable_value(
                    names[i], step))
        # adjust one of them
        values[names.index(name)] = ctypes.c_double(value)
        # write them all
        EXODUS_LIB.ex_put_var(self.fileId,
                              ctypes.c_int(step),
                              ctypes.c_int(get_entity_type('EX_GLOBAL')),
                              ctypes.c_int(1),
                              ctypes.c_int(0),
                              ctypes.c_int(numVals),
                              values)
        return True

    # --------------------------------------------------------------------

    def put_all_global_variable_values(self, step, values):
        """
        store all global variable values (one for each global variable
        name, and in the order given by exo.get_global_variable_names())
        at a specified time step

        >>> status = exo.put_all_global_variable_values(time_step, gvar_vals)

        Parameters
        ----------
        step : int
           1-based index of time step
        values : list<double>

        Returns
        -------
        status : bool
            True = successful execution
        """
        numVals = self.get_variable_number('EX_GLOBAL')
        gvalues = (ctypes.c_double * numVals)()
        for i in range(numVals):
            gvalues[i] = ctypes.c_double(values[i])
        EXODUS_LIB.ex_put_var(self.fileId,
                              ctypes.c_int(step),
                              ctypes.c_int(get_entity_type('EX_GLOBAL')),
                              ctypes.c_int(1),
                              ctypes.c_int(0),
                              ctypes.c_int(numVals),
                              gvalues)
        return True

    # --------------------------------------------------------------------

    def get_global_variable_values(self, name):
        """
        get global variable values over all time steps for one global
        variable name

        >>> gvar_vals = exo.get_global_variable_values(gvar_name)

        Parameters
        ----------
        name : string
           name of global variable

        Returns
        -------

            if array_type == 'ctype':
              <list<double>>  gvar_vals

            if array_type == 'numpy':
              <np_array<double>>  gvar_vals
        """
        names = self.get_variable_names('EX_GLOBAL')
        var_id = names.index(name)
        num = self.__ex_get_variable_param('EX_GLOBAL')
        values = []
        for i in range(self.numTimes.value):
            gvalues = self.__ex_get_var(i + 1, 'EX_GLOBAL', 0, 1, num.value)
            values.append(gvalues[var_id])
        if self.use_numpy:
            values = self.np.array(values)
        return values

    # --------------------------------------------------------------------

    def put_polyhedra_elem_blk(self, blkID,
                               num_elems_this_blk,
                               num_faces,
                               num_attr_per_elem):
        """
        put in an element block with polyhedral elements

        >>> status = exo.put_polyhedra_elem_blk(blkID, num_elems_this_blk,
        ...                                     num_faces, num_attr_per_elem)

        Parameters
        ----------
        blkID : ex_entity_id
            id of the block to be added
        num_elems_this_blk : int
        num_faces  : int
            total number of faces in this block
        num_attr_per_elem : int

        Returns
        -------
        status : bool
            True = successful execution
        """

        ebType = ctypes.c_int(get_entity_type('EX_ELEM_BLOCK'))
        EXODUS_LIB.ex_put_block(self.fileId, ebType, ctypes.c_longlong(blkID),
                                ctypes.create_string_buffer(b"NFACED"),
                                ctypes.c_longlong(num_elems_this_blk),
                                ctypes.c_longlong(0),
                                ctypes.c_longlong(0),
                                ctypes.c_longlong(num_faces),
                                ctypes.c_longlong(num_attr_per_elem))
        return True

    # --------------------------------------------------------------------

    def put_polyhedra_face_blk(self, blkID,
                               num_faces_this_blk,
                               num_nodes,
                               num_attr_per_face):
        """
        put in a block of faces

        >>> status = exo.put_polyhedra_face_blk(blkID, num_faces_this_blk,
        ...                                     num_nodes, num_attr_per_face)

        Parameters
        ----------
        blkID : ex_entity_id
            id of the block to be added
        num_faces_this_blk : int
        num_nodes  : int
            total number of nodes in this block
        num_attr_per_face : int

        Returns
        -------
        status : bool
            True = successful execution
        """
        fbType = ctypes.c_int(get_entity_type('EX_FACE_BLOCK'))
        EXODUS_LIB.ex_put_block(self.fileId, fbType, ctypes.c_longlong(blkID),
                                ctypes.create_string_buffer(b"NSIDED"),
                                ctypes.c_longlong(num_faces_this_blk),
                                ctypes.c_longlong(num_nodes),
                                ctypes.c_longlong(0),
                                ctypes.c_longlong(0),
                                ctypes.c_longlong(num_attr_per_face))
        return True

    # --------------------------------------------------------------------

    def put_face_count_per_polyhedra(self, blkID, entityCounts):
        """
        put in a count of faces in for each polyhedra in an elem block

        >>> status = exo.put_face_count_per_polyhedra(blkID, entityCounts)

        Parameters
        ----------
        blkID : ex_entity_id
            id of the block to be added

        if array_type == 'ctype':
            <list<int>>  entityCounts

        if array_type == 'numpy':
              <np_array<int>>  entityCounts

        Returns
        -------
        status : bool
            True = successful execution
        """
        ebType = ctypes.c_int(get_entity_type('EX_ELEM_BLOCK'))
        entity_counts = (ctypes.c_int * len(entityCounts))()
        entity_counts[:] = entityCounts
        EXODUS_LIB.ex_put_entity_count_per_polyhedra(
            self.fileId, ebType, ctypes.c_longlong(blkID), entity_counts)
        return True

    # --------------------------------------------------------------------

    def put_node_count_per_face(self, blkID, entityCounts):
        """
        put in a count of nodes in for each face in a polygonal face block

        >>> status = exo.put_node_count_per_face(blkID, entityCounts)

        Parameters
        ----------
        blkID : ex_entity_id
            id of the block to be added

        if array_type == 'ctype':
            <list<int>>  entityCounts

        if array_type == 'numpy':
            <np_array<int>>  entityCounts

        Returns
        -------
        status : bool
            True = successful execution
        """
        ebType = ctypes.c_int(get_entity_type('EX_FACE_BLOCK'))
        entity_counts = (ctypes.c_int * len(entityCounts))()
        entity_counts[:] = entityCounts
        EXODUS_LIB.ex_put_entity_count_per_polyhedra(
            self.fileId, ebType, ctypes.c_longlong(blkID), entity_counts)
        return True

    # --------------------------------------------------------------------

    def put_elem_face_conn(self, blkId, elemFaceConn):
        """
        put in connectivity information from elems to faces

        >>> status = exo.put_elem_face_conn(blkID, elemFaceConn)

        Parameters
        ----------
        blkID : ex_entity_id
            id of the block to be added

        if array_type == 'ctype':
            <list<int>>  elemFaceConn  (raveled/flat list)

        if array_type == 'numpy':
            <np_array<int>>  elemFaceConn  (raveled/flat array)

        Returns
        -------
        status : bool
            True = successful execution
        """
        ebType = ctypes.c_int(get_entity_type('EX_ELEM_BLOCK'))
        elem_face_conn = (ctypes.c_int * len(elemFaceConn))()
        elem_face_conn[:] = elemFaceConn
        EXODUS_LIB.ex_put_conn(self.fileId, ebType, ctypes.c_longlong(blkId),
                               None, None, elem_face_conn)
        return True

    # --------------------------------------------------------------------

    def put_face_node_conn(self, blkId, faceNodeConn):
        """
        put in connectivity information from faces to nodes

        >>> status = exo.put_face_node_conn(blkID, faceNodeConn)

        Parameters
        ----------
        blkID : ex_entity_id
            id of the block to be added

        if array_type == 'ctype':
            <list<int>>  faceNodeConn  (raveled/flat list)

        if array_type == 'numpy':
            <np_array<int>>  faceNodeConn  (raveled/flat array)

        Returns
        -------
        status : bool
            True = successful execution
        """
        ebType = ctypes.c_int(get_entity_type('EX_FACE_BLOCK'))
        node_conn = (ctypes.c_int * len(faceNodeConn))()
        node_conn[:] = faceNodeConn
        EXODUS_LIB.ex_put_conn(self.fileId, ebType, ctypes.c_longlong(blkId),
                               node_conn, None, None)
        return True

    # --------------------------------------------------------------------

    def close(self):
        """
        close the exodus file

        >>> exo.close()

        Note
        -----
        Can only be called once for an exodus object, and once called
        all methods for that object become inoperable
        """
        print(f"Closing exodus file: {self.fileName}")
        errorInt = EXODUS_LIB.ex_close(self.fileId)
        if errorInt != 0:
            raise Exception(
                "ERROR: Closing file " + self.fileName + " had problems.")

    # --------------------------------------------------------------------
    #
    # Private Exodus API calls
    #
    # --------------------------------------------------------------------

    def __open(self, io_size=0):
        print(f"Opening exodus file: {self.fileName}")
        self.mode = EX_READ
        if self.modeChar.lower() == "a":
            self.mode = EX_WRITE
        if self.modeChar.lower() == "w+":
            self.mode = EX_CLOBBER

        if self.modeChar.lower() in [
                "a", "r"] and not os.path.isfile(self.fileName):
            raise Exception(
                "ERROR: Cannot open " + self.fileName + " for read. Does not exist.")
        elif self.modeChar.lower() == "w" and os.path.isfile(self.fileName):
            raise Exception(f"ERROR: Cowardly not opening {self.fileName} for write. File already exists.")

        elif self.modeChar.lower() not in ["a", "r", "w", "w+"]:
            raise Exception(
                "ERROR: File open mode " + self.modeChar + " unrecognized.")

        self.comp_ws = ctypes.c_int(8)
        self.io_ws = ctypes.c_int(io_size)
        self.version = ctypes.c_float(0.0)
        if self.modeChar.lower() in ["a", "r"]:  # open existing file
            self.fileId = EXODUS_LIB.ex_open_int(self.fileName.encode('ascii'), self.mode,
                                                 ctypes.byref(self.comp_ws),
                                                 ctypes.byref(self.io_ws),
                                                 ctypes.byref(self.version),
                                                 EX_API_VERSION_NODOT)
        else:  # create file
            if io_size == 0:
                io_size = 8
                self.io_ws = ctypes.c_int(io_size)
            self.__create()

    # --------------------------------------------------------------------

    def __create(self):
        self.fileId = EXODUS_LIB.ex_create_int(self.fileName.encode('ascii'), self.mode,
                                               ctypes.byref(self.comp_ws),
                                               ctypes.byref(self.io_ws),
                                               EX_API_VERSION_NODOT)

    # --------------------------------------------------------------------

    def __ex_get_info(self):
        self.Title = ctypes.create_string_buffer(MAX_LINE_LENGTH + 1)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API:
            self.numDim = ctypes.c_longlong(0)
            self.numNodes = ctypes.c_longlong(0)
            self.numElem = ctypes.c_longlong(0)
            self.numElemBlk = ctypes.c_longlong(0)
            self.numNodeSets = ctypes.c_longlong(0)
            self.numSideSets = ctypes.c_longlong(0)
            self.numAssembly = ctypes.c_longlong(0)
            self.numBlob = ctypes.c_longlong(0)
        else:
            self.numDim = ctypes.c_int(0)
            self.numNodes = ctypes.c_int(0)
            self.numElem = ctypes.c_int(0)
            self.numElemBlk = ctypes.c_int(0)
            self.numNodeSets = ctypes.c_int(0)
            self.numSideSets = ctypes.c_int(0)
            self.numAssembly = ctypes.c_int(0)
            self.numBlob = ctypes.c_int(0)
        EXODUS_LIB.ex_get_init(
            self.fileId, self.Title,
            ctypes.byref(self.numDim),
            ctypes.byref(self.numNodes),
            ctypes.byref(self.numElem),
            ctypes.byref(self.numElemBlk),
            ctypes.byref(self.numNodeSets),
            ctypes.byref(self.numSideSets))

    # --------------------------------------------------------------------

    def __ex_put_info(self, info):
        self.Title = ctypes.create_string_buffer(info[0].encode('ascii'), MAX_LINE_LENGTH + 1)
        self.numDim = ctypes.c_longlong(info[1])
        self.numNodes = ctypes.c_longlong(info[2])
        self.numElem = ctypes.c_longlong(info[3])
        self.numElemBlk = ctypes.c_longlong(info[4])
        self.numNodeSets = ctypes.c_longlong(info[5])
        self.numSideSets = ctypes.c_longlong(info[6])
        EXODUS_LIB.ex_put_init(
            self.fileId,
            self.Title,
            self.numDim,
            self.numNodes,
            self.numElem,
            self.numElemBlk,
            self.numNodeSets,
            self.numSideSets)
        self.version = self.__ex_inquire_float(ex_inquiry_map('EX_INQ_DB_VERS'))

    # --------------------------------------------------------------------

    def __ex_put_concat_elem_blk(self, elemBlkIDs, elemType, numElemThisBlk,
                                 numNodesPerElem, numAttr, defineMaps):
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API:
            elem_blk_ids = (ctypes.c_longlong * len(elemBlkIDs))()
            elem_blk_ids[:] = elemBlkIDs
            num_elem_this_blk = (ctypes.c_longlong * len(elemBlkIDs))()
            num_elem_this_blk[:] = numElemThisBlk
            num_nodes_per_elem = (ctypes.c_longlong * len(elemBlkIDs))()
            num_nodes_per_elem[:] = numNodesPerElem
            num_attr = (ctypes.c_longlong * len(elemBlkIDs))()
        else:
            elem_blk_ids = (ctypes.c_int * len(elemBlkIDs))()
            elem_blk_ids[:] = elemBlkIDs
            num_elem_this_blk = (ctypes.c_int * len(elemBlkIDs))()
            num_elem_this_blk[:] = numElemThisBlk
            num_nodes_per_elem = (ctypes.c_int * len(elemBlkIDs))()
            num_nodes_per_elem[:] = numNodesPerElem
            num_attr = (ctypes.c_int * len(elemBlkIDs))()
        num_attr[:] = numAttr
        elem_type = (ctypes.c_char_p * len(elemBlkIDs))()
        elem_type[:] = elemType
        define_maps = ctypes.c_int(defineMaps)
        EXODUS_LIB.ex_put_concat_elem_block(
            self.fileId,
            elem_blk_ids,
            elem_type,
            num_elem_this_blk,
            num_nodes_per_elem,
            num_attr,
            define_maps)

    # --------------------------------------------------------------------

    def __ex_get_qa(self):
        num_qa_recs = ctypes.c_int(self.__ex_inquire_int(ex_inquiry_map('EX_INQ_QA')))
        qa_rec_ptrs = ((ctypes.POINTER(ctypes.c_char * (MAX_STR_LENGTH + 1)) * 4) * num_qa_recs.value)()
        for i in range(num_qa_recs.value):
            for j in range(4):
                qa_rec_ptrs[i][j] = ctypes.pointer(
                    ctypes.create_string_buffer(MAX_STR_LENGTH + 1))
        if num_qa_recs.value:
            EXODUS_LIB.ex_get_qa(self.fileId, ctypes.byref(qa_rec_ptrs))
        qa_recs = []
        for qara in qa_rec_ptrs:
            qa_rec_list = [ptr.contents.value.decode("utf8") for ptr in qara]
            qa_rec_tuple = tuple(qa_rec_list)
            assert len(qa_rec_tuple) == 4
            qa_recs.append(qa_rec_tuple)
        return qa_recs

    # --------------------------------------------------------------------

    def __ex_put_qa(self, qaRecs):
        num_qa_recs = ctypes.c_int(len(qaRecs))
        qa_rec_ptrs = ((ctypes.POINTER(ctypes.c_char * (MAX_STR_LENGTH + 1)) * 4) * num_qa_recs.value)()
        for i in range(num_qa_recs.value):
            for j in range(4):
                qa_rec_ptrs[i][j] = ctypes.pointer(ctypes.create_string_buffer(
                    str(qaRecs[i][j]).encode('ascii'), MAX_STR_LENGTH + 1))
        EXODUS_LIB.ex_put_qa(self.fileId, num_qa_recs, ctypes.byref(qa_rec_ptrs))
        return True

    # --------------------------------------------------------------------

    def _ex_get_info_recs_quietly(self):
        num_infos = ctypes.c_int(self.__ex_inquire_int(ex_inquiry_map('EX_INQ_INFO')))
        info_ptrs = (ctypes.POINTER(ctypes.c_char * (MAX_LINE_LENGTH + 1)) * num_infos.value)()
        for i in range(num_infos.value):
            info_ptrs[i] = ctypes.pointer(ctypes.create_string_buffer(MAX_LINE_LENGTH + 1))
        if num_infos.value:
            EXODUS_LIB.ex_get_info(self.fileId, ctypes.byref(info_ptrs))
        return [irp.contents.value.decode("utf8") for irp in info_ptrs]

    # --------------------------------------------------------------------

    def __ex_get_info_recs(self):
        num_infos = ctypes.c_int(self.__ex_inquire_int(ex_inquiry_map('EX_INQ_INFO')))
        info_ptrs = (ctypes.POINTER(ctypes.c_char * (MAX_LINE_LENGTH + 1)) * num_infos.value)()
        for i in range(num_infos.value):
            info_ptrs[i] = ctypes.pointer(ctypes.create_string_buffer(MAX_LINE_LENGTH + 1))
        EXODUS_LIB.ex_get_info(self.fileId, ctypes.byref(info_ptrs))
        info_recs = [irp.contents.value.decode("utf8") for irp in info_ptrs]
        for rec in info_recs:
            if len(rec) > MAX_LINE_LENGTH:
                print("WARNING: max line length reached for one or more info records;")
                print("         info might be incomplete for these records")
                break
        return info_recs

    # --------------------------------------------------------------------

    def __ex_put_info_recs(self, infoRecs):
        num_infos = ctypes.c_int(len(infoRecs))
        info_ptrs = (ctypes.POINTER(ctypes.c_char * (MAX_LINE_LENGTH + 1)) * num_infos.value)()
        for i in range(num_infos.value):
            info_ptrs[i] = ctypes.pointer(ctypes.create_string_buffer(
                str(infoRecs[i]).encode('ascii'), MAX_LINE_LENGTH + 1))
        EXODUS_LIB.ex_put_info(self.fileId, num_infos, ctypes.byref(info_ptrs))
        return True

    # --------------------------------------------------------------------

    def __ex_inquire_float(self, inq_id):
        dummy_char = ctypes.create_string_buffer(MAX_LINE_LENGTH + 1)
        ret_float = ctypes.c_float(0.0)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_INQ_INT64_API:
            dummy_int = ctypes.c_longlong(0)
        else:
            dummy_int = ctypes.c_int(0)
        val = EXODUS_LIB.ex_inquire(
            self.fileId,
            inq_id,
            ctypes.byref(dummy_int),
            ctypes.byref(ret_float),
            dummy_char)
        if val < 0:
            raise Exception(
                "ERROR: ex_inquire(" + str(inq_id) + ") failed on " + self.fileName)
        return ret_float

    # --------------------------------------------------------------------

    def __ex_inquire_int(self, inq_id):
        val = EXODUS_LIB.ex_inquire_int(self.fileId, inq_id)
        if val < 0:
            raise Exception(
                "ERROR: ex_inquire_int(" + str(inq_id) + ") failed on " + self.fileName)
        return val

    # --------------------------------------------------------------------

    def __ex_get_coord_names(self):
        coord_name_ptrs = (
            ctypes.POINTER(ctypes.c_char * (MAX_NAME_LENGTH + 1)) * self.numDim.value)()
        for i in range(self.numDim.value):
            coord_name_ptrs[i] = ctypes.pointer(
                ctypes.create_string_buffer(
                    MAX_NAME_LENGTH + 1))
        EXODUS_LIB.ex_get_coord_names(self.fileId, ctypes.byref(coord_name_ptrs))
        return [cnp.contents.value.decode('utf8') for cnp in coord_name_ptrs]

    # --------------------------------------------------------------------

    def __ex_put_coord_names(self, names):
        coord_name_ptrs = (
            ctypes.POINTER(ctypes.c_char * (MAX_NAME_LENGTH + 1)) * self.numDim.value)()
        assert len(names) == self.numDim.value
        for i in range(self.numDim.value):
            coord_name_ptrs[i] = ctypes.pointer(
                ctypes.create_string_buffer(
                    names[i].encode('ascii'), MAX_NAME_LENGTH + 1))
        EXODUS_LIB.ex_put_coord_names(self.fileId, ctypes.byref(coord_name_ptrs))

    # --------------------------------------------------------------------

    def __ex_get_all_times(self):
        self.times = (ctypes.c_double * self.numTimes.value)()
        EXODUS_LIB.ex_get_all_times(self.fileId, ctypes.byref(self.times))

    # --------------------------------------------------------------------

    def __ex_get_time(self, timeStep):
        time_step = ctypes.c_int(timeStep)
        time_val = ctypes.c_double(0.0)
        EXODUS_LIB.ex_get_time(self.fileId, time_step, ctypes.byref(time_val))
        return time_val.value()

    # --------------------------------------------------------------------

    def __ex_put_time(self, timeStep, timeVal):
        time_step = ctypes.c_int(timeStep)
        time_val = ctypes.c_double(timeVal)
        EXODUS_LIB.ex_put_time(self.fileId, time_step, ctypes.byref(time_val))
        return True

    # --------------------------------------------------------------------

    def __ex_get_name(self, objType, objId):
        obj_type = ctypes.c_int(get_entity_type(objType))
        obj_id = ctypes.c_longlong(objId)
        obj_name = ctypes.create_string_buffer(MAX_NAME_LENGTH + 1)
        EXODUS_LIB.ex_get_name(self.fileId, obj_type, obj_id, ctypes.byref(obj_name))
        return obj_name.value.decode('utf8')

    # --------------------------------------------------------------------

    def __ex_put_name(self, objType, objId, objName):
        obj_type = ctypes.c_int(get_entity_type(objType))
        obj_id = ctypes.c_longlong(objId)
        obj_name = ctypes.create_string_buffer(objName.encode('ascii'), MAX_NAME_LENGTH + 1)
        EXODUS_LIB.ex_put_name(self.fileId, obj_type, obj_id, obj_name)

    # --------------------------------------------------------------------

    def __ex_get_names(self, objType):
        inqType = ex_inquiry_map(ex_obj_to_inq(objType))
        num_objs = ctypes.c_int(self.__ex_inquire_int(inqType)).value
        obj_name_ptrs = (ctypes.POINTER(ctypes.c_char * (MAX_NAME_LENGTH + 1)) * num_objs)()
        for i in range(num_objs):
            obj_type = ctypes.c_int(get_entity_type(objType))
            obj_name_ptrs[i] = ctypes.pointer(
                ctypes.create_string_buffer(
                    MAX_NAME_LENGTH + 1))

        EXODUS_LIB.ex_get_names(self.fileId, obj_type, ctypes.byref(obj_name_ptrs))
        return [onp.contents.value.decode('utf8') for onp in obj_name_ptrs]

    def __ex_put_names(self, objType, objNames):
        inqType = ex_inquiry_map(ex_obj_to_inq(objType))
        numObjs = ctypes.c_int(self.__ex_inquire_int(inqType)).value
        assert numObjs == len(objNames)
        obj_name_ptrs = (ctypes.POINTER(ctypes.c_char * (MAX_NAME_LENGTH + 1)) * numObjs)()
        obj_type = ctypes.c_int(get_entity_type(objType))
        for i in range(numObjs):
            obj_name_ptrs[i] = ctypes.pointer(
                ctypes.create_string_buffer(
                    objNames[i].encode('ascii'), MAX_NAME_LENGTH + 1))
        EXODUS_LIB.ex_put_names(self.fileId, obj_type, ctypes.byref(obj_name_ptrs))

    def __ex_get_ids(self, objType):
        inqType = ex_inquiry_map(ex_obj_to_inq(objType))
        numObjs = ctypes.c_int(self.__ex_inquire_int(inqType)).value
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API:
            ids = (ctypes.c_longlong * numObjs)()
        else:
            ids = (ctypes.c_int * numObjs)()
        if numObjs > 0:
            obj_type = ctypes.c_int(get_entity_type(objType))
            EXODUS_LIB.ex_get_ids(self.fileId, obj_type, ctypes.byref(ids))
        return ids

    def __ex_get_assembly(self, assem_struct):
        EXODUS_LIB.ex_get_assembly(self.fileId, ctypes.byref(assem_struct))
        ptr = ctypes.create_string_buffer(MAX_NAME_LENGTH + 1)
        assem_struct.name = ctypes.cast(ptr, ctypes.c_char_p)
        eptr = (ctypes.c_longlong * assem_struct.entity_count)()
        assem_struct.entity_list = eptr
        EXODUS_LIB.ex_get_assembly(self.fileId, ctypes.byref(assem_struct))

    def __ex_get_assemblies(self, assem_list):
        EXODUS_LIB.ex_get_assemblies(self.fileId, assem_list)
        for assem_struct in assem_list:
            ptr = ctypes.create_string_buffer(MAX_NAME_LENGTH + 1)
            assem_struct.name = ctypes.cast(ptr, ctypes.c_char_p)
            eptr = (ctypes.c_longlong * assem_struct.entity_count)()
            assem_struct.entity_list = eptr
        EXODUS_LIB.ex_get_assemblies(self.fileId, assem_list)

    def __ex_get_blob(self, blob_struct):
        EXODUS_LIB.ex_get_blob(self.fileId, ctypes.byref(blob_struct))
        ptr = ctypes.create_string_buffer(MAX_NAME_LENGTH + 1)
        blob_struct.name = ctypes.cast(ptr, ctypes.c_char_p)
        EXODUS_LIB.ex_get_blob(self.fileId, ctypes.byref(blob_struct))

    def __ex_put_assembly(self, assembly):
        assem = setup_ex_assembly(assembly)
        EXODUS_LIB.ex_put_assembly(self.fileId, assem)

    def __ex_put_assemblies(self, assemblies):
        assembly_list = []
        for assembly in assemblies:
            assem = setup_ex_assembly(assembly)
            assembly_list.append(assem)
        assems = (ex_assembly * len(assemblies))(*assembly_list)

        EXODUS_LIB.ex_put_assemblies(self.fileId, len(assembly_list), assems)

    def __ex_get_attribute_count(self, objType, objId):
        # Get attribute count...
        obj_type = ctypes.c_int(get_entity_type(objType))
        obj_id = ctypes.c_longlong(objId)
        return EXODUS_LIB.ex_get_attribute_count(self.fileId, obj_type, obj_id)

    def __ex_get_attributes(self, objType, objId):
        # Get attribute count...
        obj_type = ctypes.c_int(get_entity_type(objType))
        obj_id = ctypes.c_longlong(objId)
        att_count = EXODUS_LIB.ex_get_attribute_count(self.fileId, obj_type, obj_id)

        attributes = {}
        if att_count > 0:
            att = (ex_attribute * att_count)()
            EXODUS_LIB.ex_get_attribute_param(self.fileId, obj_type, obj_id, ctypes.byref(att))
            for i in range(att_count):
                EXODUS_LIB.ex_get_attribute(self.fileId, ctypes.byref(att[i]))
                tmp_att = attribute(att[i].name.decode('utf8'), att[i].entity_type, att[i].entity_id)

                if (att[i].type == 2):
                    vals = ctypes.cast(att[i].values, ctypes.POINTER(ctypes.c_char))
                    tmp = [vals[j] for j in range(att[i].value_count - 1)]
                    tmp_att.values = b''.join(tmp).decode('utf8')

                if (att[i].type == 4):
                    vals = ctypes.cast(att[i].values, ctypes.POINTER(ctypes.c_int))
                    for j in range(att[i].value_count):
                        tmp_att.values.append(vals[j])

                if (att[i].type == 6):
                    vals = ctypes.cast(att[i].values, ctypes.POINTER(ctypes.c_double))
                    for j in range(att[i].value_count):
                        tmp_att.values.append(vals[j])

                attributes[att[i].name.decode('utf8')] = tmp_att

        return attributes

    def __ex_put_attribute(self, attribute):
        att_id = ctypes.c_longlong(attribute.entity_id)
        att = ex_attribute(entity_id=att_id)
        att.name = attribute.name.encode('ascii')
        att.entity_type = ctypes.c_int(get_entity_type(attribute.entity_type))
        att.value_count = len(attribute.values)

        if (isinstance(attribute.values[0], int)):
            eptr = (ctypes.c_int * len(attribute.values))()
            for i in range(len(attribute.values)):
                eptr[i] = ctypes.c_int(attribute.values[i])
            att.values = ctypes.cast(eptr, ctypes.c_void_p)
            att.type = 4

        elif (isinstance(attribute.values[0], float)):
            eptr = (ctypes.c_double * len(attribute.values))()
            for i in range(len(attribute.values)):
                eptr[i] = ctypes.c_double(attribute.values[i])
            att.values = ctypes.cast(eptr, ctypes.c_void_p)
            att.type = 6

        elif (isinstance(attribute.values[0], str)):
            eptr = (ctypes.c_char * (len(attribute.values) + 1))()
            eptr = attribute.values[0].encode('ascii')
            att.values = ctypes.cast(eptr, ctypes.c_void_p)
            att.type = 2

        EXODUS_LIB.ex_put_attribute(self.fileId, att)

    def __ex_get_node_set(self, nodeSetId):
        node_set_id = ctypes.c_longlong(nodeSetId)
        num_node_set_nodes = self.__ex_get_set_param('EX_NODE_SET', nodeSetId)[0]
        if num_node_set_nodes == 0:
            return []
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API:
            set_nodes = (ctypes.c_longlong * num_node_set_nodes)()
        else:
            set_nodes = (ctypes.c_int * num_node_set_nodes)()
        EXODUS_LIB.ex_get_set(self.fileId, ctypes.c_int(get_entity_type('EX_NODE_SET')), node_set_id, ctypes.byref(set_nodes), None)
        return set_nodes

    def __ex_put_node_set(self, nodeSetId, nodeSetNodes):
        node_set_id = ctypes.c_longlong(nodeSetId)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API:
            node_set_nodes = (ctypes.c_longlong * len(nodeSetNodes))()
            for i, node_set_node in enumerate(nodeSetNodes):
                node_set_nodes[i] = ctypes.c_longlong(node_set_node)
        else:
            node_set_nodes = (ctypes.c_int * len(nodeSetNodes))()
            for i, node_set_node in enumerate(nodeSetNodes):
                node_set_nodes[i] = ctypes.c_int(node_set_node)
        EXODUS_LIB.ex_put_set(self.fileId, ctypes.c_int(get_entity_type('EX_NODE_SET')), node_set_id, node_set_nodes, None)

    def __ex_get_node_set_dist_fact(self, nodeSetId):
        node_set_id = ctypes.c_longlong(nodeSetId)
        num_node_set_nodes = self.__ex_get_set_param('EX_NODE_SET', nodeSetId)[0]
        set_dfs = (ctypes.c_double * num_node_set_nodes)()
        EXODUS_LIB.ex_get_node_set_dist_fact(
            self.fileId, node_set_id, ctypes.byref(set_dfs))
        return set_dfs

    def __ex_put_node_set_dist_fact(self, nodeSetId, nodeSetDistFact):
        node_set_id = ctypes.c_longlong(nodeSetId)
        node_set_dist_fact = (ctypes.c_double * len(nodeSetDistFact))()
        for i, dist_fact in enumerate(nodeSetDistFact):
            node_set_dist_fact[i] = ctypes.c_double(dist_fact)
        EXODUS_LIB.ex_put_node_set_dist_fact(
            self.fileId, node_set_id, node_set_dist_fact)

    def __ex_get_object_truth_vector(self, objType, entId):
        obj_type = ctypes.c_int(get_entity_type(objType))
        entity_id = ctypes.c_longlong(entId)
        variable_count = self.__ex_get_variable_param(objType)
        truth_table = (ctypes.c_int * (variable_count.value))()

        EXODUS_LIB.ex_get_object_truth_vector(self.fileId, obj_type,
                                              entity_id, variable_count,
                                              ctypes.byref(truth_table))
        truthTab = []
        for val in truth_table:
            if val:
                truthTab.append(True)
            else:
                truthTab.append(False)
        return truthTab

    def __ex_get_truth_table(self, objType):
        inqType = ex_inquiry_map(ex_obj_to_inq(objType))
        num_objs = ctypes.c_int(self.__ex_inquire_int(inqType)).value

        obj_type = ctypes.c_int(get_entity_type(objType))
        variable_count = self.__ex_get_variable_param(objType)

        truth_table = (ctypes.c_int * (num_objs * variable_count.value))()
        EXODUS_LIB.ex_get_truth_table(self.fileId, obj_type,
                                      num_objs, variable_count,
                                      ctypes.byref(truth_table))
        truthTab = []
        for val in truth_table:
            if val:
                truthTab.append(True)
            else:
                truthTab.append(False)
        return truthTab

    def __ex_put_truth_table(self, objType, truthTab):
        inqType = ex_inquiry_map(ex_obj_to_inq(objType))
        num_objs = ctypes.c_int(self.__ex_inquire_int(inqType)).value

        obj_type = ctypes.c_int(get_entity_type(objType))
        num_vars = self.__ex_get_variable_param(objType).value

        assert len(truthTab) == (num_objs * num_vars)

        truth_tab = (ctypes.c_int * (num_objs * num_vars))()
        for i, boolVal in enumerate(truthTab):
            truth_tab[i] = ctypes.c_int(1) if boolVal else ctypes.c_int(0)
        EXODUS_LIB.ex_put_truth_table(
            self.fileId, obj_type, num_objs, num_vars, truth_tab)
        return True

    def __ex_get_coord(self):
        self.coordsX = (ctypes.c_double * self.numNodes.value)()
        self.coordsY = (ctypes.c_double * self.numNodes.value)()
        self.coordsZ = (ctypes.c_double * self.numNodes.value)()
        EXODUS_LIB.ex_get_coord(
            self.fileId,
            ctypes.byref(self.coordsX),
            ctypes.byref(self.coordsY),
            ctypes.byref(self.coordsZ))

    def __ex_put_coord(self, xCoords, yCoords, zCoords):
        self.coordsX = (ctypes.c_double * self.numNodes.value)()
        self.coordsY = (ctypes.c_double * self.numNodes.value)()
        self.coordsZ = (ctypes.c_double * self.numNodes.value)()
        for i in range(self.numNodes.value):
            self.coordsX[i] = xCoords[i]
            self.coordsY[i] = yCoords[i]
            self.coordsZ[i] = zCoords[i]
        EXODUS_LIB.ex_put_coord(
            self.fileId,
            ctypes.byref(self.coordsX),
            ctypes.byref(self.coordsY),
            ctypes.byref(self.coordsZ))

    def __ex_get_partial_coord(self, startNodeId, numNodes):
        start_node_num = ctypes.c_longlong(startNodeId)
        num_nodes = ctypes.c_longlong(numNodes)
        coordsX = (ctypes.c_double * numNodes)()
        coordsY = (ctypes.c_double * numNodes)()
        coordsZ = (ctypes.c_double * numNodes)()
        EXODUS_LIB.ex_get_partial_coord(
            self.fileId,
            start_node_num,
            num_nodes,
            ctypes.byref(coordsX),
            ctypes.byref(coordsY),
            ctypes.byref(coordsZ))
        return list(coordsX), list(coordsY), list(coordsZ)

    def __ex_get_id_map(self, objType):
        inqType = ex_obj_to_inq(objType)
        obj_type = ctypes.c_int(get_entity_type(objType))
        inq_type = ctypes.c_int(ex_inquiry_map(inqType))
        num_objs = ctypes.c_int(self.__ex_inquire_int(inq_type))
        numObjs = num_objs.value
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API:
            id_map = (ctypes.c_longlong * numObjs)()
        else:
            id_map = (ctypes.c_int * numObjs)()
        EXODUS_LIB.ex_get_id_map(self.fileId, obj_type, ctypes.byref(id_map))
        idMap = [id_map[i] for i in range(numObjs)]
        if self.use_numpy:
            idMap = self.np.array(idMap)
        return idMap

    def __ex_get_num_map(self, objType, idx):
        inqType = ex_obj_to_inq(objType)
        map_id = ctypes.c_longlong(idx)
        obj_type = ctypes.c_int(get_entity_type(objType))
        inq_type = ctypes.c_int(ex_inquiry_map(inqType))
        num_objs = ctypes.c_int(self.__ex_inquire_int(inq_type))
        numObjs = num_objs.value
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API:
            id_map = (ctypes.c_longlong * numObjs)()
        else:
            id_map = (ctypes.c_int * numObjs)()
        EXODUS_LIB.ex_get_num_map(self.fileId, obj_type, map_id, ctypes.byref(id_map))
        idMap = [id_map[i] for i in range(numObjs)]
        if self.use_numpy:
            idMap = self.np.array(idMap)
        return idMap

    def __ex_put_map_param(self, nodeMapCnt, elemMapCnt):
        node_map_cnt = ctypes.c_int(nodeMapCnt)
        elem_map_cnt = ctypes.c_int(elemMapCnt)
        errorInt = EXODUS_LIB.ex_put_map_param(
            self.fileId, node_map_cnt, elem_map_cnt)
        if errorInt != 0:
            print(("ERROR code =", errorInt))
            raise Exception(
                "ERROR: ex_put_map_param had problems.")
        return True

    def __ex_put_num_map(self, objType, idx, numMap):
        inqType = ex_obj_to_inq(objType)
        map_id = ctypes.c_longlong(idx)
        obj_type = ctypes.c_int(get_entity_type(objType))
        inq_type = ctypes.c_int(ex_inquiry_map(inqType))
        num_objs = ctypes.c_int(self.__ex_inquire_int(inq_type))
        numObjs = num_objs.value
        assert numObjs == len(numMap)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API:
            num_map = (ctypes.c_longlong * numObjs)()
            for i in range(numObjs):
                num_map[i] = ctypes.c_longlong(numMap[i])
        else:
            num_map = (ctypes.c_int * numObjs)()
            for i in range(numObjs):
                num_map[i] = ctypes.c_int(numMap[i])
        EXODUS_LIB.ex_put_num_map(self.fileId, obj_type, map_id, ctypes.byref(num_map))
        return True

    def __ex_get_block_id_map(self, obj_type, id):
        obj_type = ctypes.c_int(get_entity_type(obj_type))
        entity_id = ctypes.c_longlong(id)
        _, numObjs, _, _ = self.__ex_get_block('EX_ELEM_BLOCK', id)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API:
            id_map = (ctypes.c_longlong * numObjs)()
        else:
            id_map = (ctypes.c_int * numObjs)()
        EXODUS_LIB.ex_get_block_id_map(self.fileId, obj_type, entity_id, id_map)
        if self.use_numpy:
            id_map = ctype_to_numpy(self, id_map)
            return id_map
        else:
            idMap = [id_map[i] for i in range(numObjs)]
            return idMap

    def __ex_put_id_map(self, objType, idMap):
        inqType = ex_obj_to_inq(objType)
        obj_type = ctypes.c_int(get_entity_type(objType))
        inq_type = ctypes.c_int(ex_inquiry_map(inqType))
        num_objs = ctypes.c_int(self.__ex_inquire_int(inq_type))
        numObjs = num_objs.value
        assert numObjs == len(idMap)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API:
            id_map = (ctypes.c_longlong * numObjs)()
            for i in range(numObjs):
                id_map[i] = ctypes.c_longlong(idMap[i])
        else:
            id_map = (ctypes.c_int * numObjs)()
            for i in range(numObjs):
                id_map[i] = ctypes.c_int(idMap[i])
        EXODUS_LIB.ex_put_id_map(self.fileId, obj_type, ctypes.byref(id_map))
        return True

    def __ex_get_elem_num_map(self):
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_MAPS_INT64_API:
            elemNumMap = (ctypes.c_longlong * self.numElem.value)()
        else:
            elemNumMap = (ctypes.c_int * self.numElem.value)()
        EXODUS_LIB.ex_get_elem_num_map(self.fileId, ctypes.byref(elemNumMap))
        return elemNumMap

    def __ex_get_node_num_map(self):
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_MAPS_INT64_API:
            nodeNumMap = (ctypes.c_longlong * self.numNodes.value)()
        else:
            nodeNumMap = (ctypes.c_int * self.numNodes.value)()
        EXODUS_LIB.ex_get_node_num_map(self.fileId, ctypes.byref(nodeNumMap))
        return nodeNumMap

    def __ex_get_elem_order_map(self):
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_MAPS_INT64_API:
            elemOrderMap = (ctypes.c_longlong * self.numElem.value)()
        else:
            elemOrderMap = (ctypes.c_int * self.numElem.value)()
        EXODUS_LIB.ex_get_map(self.fileId, ctypes.byref(elemOrderMap))
        return elemOrderMap

    def __ex_get_block(self, object_type, object_id):
        obj_type = ctypes.c_int(get_entity_type(object_type))
        block_id = ctypes.c_longlong(object_id)
        blk_type = ctypes.create_string_buffer(MAX_STR_LENGTH + 1)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API:
            num_elem_this_blk = ctypes.c_longlong(0)
            num_nodes_per_elem = ctypes.c_longlong(0)
            num_edges_per_elem = ctypes.c_longlong(0)
            num_faces_per_elem = ctypes.c_longlong(0)
            num_attr = ctypes.c_longlong(0)
        else:
            num_elem_this_blk = ctypes.c_int(0)
            num_nodes_per_elem = ctypes.c_int(0)
            num_edges_per_elem = ctypes.c_int(0)
            num_faces_per_elem = ctypes.c_int(0)
            num_attr = ctypes.c_int(0)
        EXODUS_LIB.ex_get_block(
            self.fileId,
            obj_type,
            block_id,
            blk_type,
            ctypes.byref(num_elem_this_blk),
            ctypes.byref(num_nodes_per_elem),
            ctypes.byref(num_edges_per_elem),
            ctypes.byref(num_faces_per_elem),
            ctypes.byref(num_attr))
        return blk_type.value, int(num_elem_this_blk.value), int(num_nodes_per_elem.value), int(num_attr.value)

    def __ex_put_block(
            self,
            object_type,
            object_id,
            eType,
            numElems,
            numNodesPerElem,
            numAttrsPerElem):
        obj_type = ctypes.c_int(get_entity_type(object_type))
        block_id = ctypes.c_longlong(object_id)
        if isinstance(eType, str):
            eType = eType.encode('ascii')
        elem_type = ctypes.create_string_buffer(eType.upper(), MAX_NAME_LENGTH + 1)
        num_elem_this_blk = ctypes.c_longlong(numElems)
        num_nodes_per_elem = ctypes.c_longlong(numNodesPerElem)
        num_edges_per_elem = ctypes.c_longlong(0)
        num_faces_per_elem = ctypes.c_longlong(0)
        num_attr = ctypes.c_longlong(numAttrsPerElem)
        EXODUS_LIB.ex_put_block(self.fileId, obj_type, block_id, elem_type,
                                num_elem_this_blk, num_nodes_per_elem,
                                num_edges_per_elem, num_faces_per_elem, num_attr)

    def __ex_get_elem_conn(self, object_id):
        (_elem_type, num_elem_this_blk, num_nodes_per_elem,
         _num_attr) = self.__ex_get_block('EX_ELEM_BLOCK', object_id)
        elem_block_id = ctypes.c_longlong(object_id)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API:
            elem_block_connectivity = (
                ctypes.c_longlong * (num_elem_this_blk * num_nodes_per_elem))()
        else:
            elem_block_connectivity = (
                ctypes.c_int * (num_elem_this_blk * num_nodes_per_elem))()
        EXODUS_LIB.ex_get_conn(
            self.fileId,
            ctypes.c_int(get_entity_type('EX_ELEM_BLOCK')),
            elem_block_id,
            ctypes.byref(elem_block_connectivity), None, None)
        return elem_block_connectivity, num_elem_this_blk, num_nodes_per_elem

    def __ex_put_elem_conn(self, object_id, connectivity):
        (_elem_type, num_elem_this_blk, num_nodes_per_elem,
         _num_attr) = self.__ex_get_block('EX_ELEM_BLOCK', object_id)
        elem_block_id = ctypes.c_longlong(object_id)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API:
            elem_block_connectivity = (
                ctypes.c_longlong * (num_elem_this_blk * num_nodes_per_elem))()
            for i in range(num_elem_this_blk * num_nodes_per_elem):
                elem_block_connectivity[i] = ctypes.c_longlong(connectivity[i])
        else:
            elem_block_connectivity = (
                ctypes.c_int * (num_elem_this_blk * num_nodes_per_elem))()
            for i in range(num_elem_this_blk * num_nodes_per_elem):
                elem_block_connectivity[i] = ctypes.c_int(connectivity[i])
        EXODUS_LIB.ex_put_conn(
            self.fileId,
            ctypes.c_int(get_entity_type('EX_ELEM_BLOCK')),
            elem_block_id,
            elem_block_connectivity, None, None)

    def __ex_put_one_attr(self, objType, elemBlkID, attrIndx, Attr):
        elem_blk_id = ctypes.c_longlong(elemBlkID)
        obj_type = ctypes.c_int(get_entity_type(objType))
        attr_index = ctypes.c_longlong(attrIndx)
        attrib = (ctypes.c_double * len(Attr))()
        for i, attr in enumerate(Attr):
            attrib[i] = attr
        EXODUS_LIB.ex_put_one_attr(
            self.fileId,
            obj_type,
            elem_blk_id,
            attr_index,
            attrib)

    def __ex_get_one_attr(self, objType, entityId, attrIndx):
        entity_id = ctypes.c_longlong(entityId)
        obj_type = ctypes.c_int(get_entity_type(objType))
        attr_index = ctypes.c_longlong(attrIndx)

        numVals = self.get_entity_count(objType, entityId)
        attrib = (ctypes.c_double * numVals)()

        EXODUS_LIB.ex_get_one_attr(
            self.fileId,
            obj_type,
            entity_id,
            attr_index,
            ctypes.byref(attrib))
        return attrib

    def __ex_put_elem_attr(self, elemBlkID, Attr):
        elem_blk_id = ctypes.c_longlong(elemBlkID)
        attrib = (ctypes.c_double * len(Attr))()
        for i, attr in enumerate(Attr):
            attrib[i] = ctypes.c_double(attr)
        EXODUS_LIB.ex_put_attr(
            self.fileId,
            ctypes.c_int(get_entity_type('EX_ELEM_BLOCK')),
            elem_blk_id,
            attrib)

    def __ex_get_elem_attr(self, elemBlkID):
        elem_blk_id = ctypes.c_longlong(elemBlkID)
        numAttrThisBlk = self.num_attr(elemBlkID)
        numElemsThisBlk = self.get_entity_count('EX_ELEM_BLOCK', elemBlkID)
        totalAttr = numAttrThisBlk * numElemsThisBlk
        attrib = (ctypes.c_double * totalAttr)()
        EXODUS_LIB.ex_get_attr(
            self.fileId,
            ctypes.c_int(get_entity_type('EX_ELEM_BLOCK')),
            elem_blk_id,
            ctypes.byref(attrib))
        return attrib

    def __ex_get_variable_param(self, varType):
        var_type = ctypes.c_int(get_entity_type(varType))
        num_vars = ctypes.c_int()
        EXODUS_LIB.ex_get_variable_param(
            self.fileId, var_type, ctypes.byref(num_vars))
        return num_vars

    def __ex_get_variable_names(self, varType):
        num_vars = self.__ex_get_variable_param(varType)
        var_name_ptrs = (
            ctypes.POINTER(ctypes.c_char * (MAX_NAME_LENGTH + 1)) * num_vars.value)()

        for i in range(num_vars.value):
            var_name_ptrs[i] = ctypes.pointer(
                ctypes.create_string_buffer(
                    MAX_NAME_LENGTH + 1))

        var_type = ctypes.c_int(get_entity_type(varType))
        EXODUS_LIB.ex_get_variable_names(
            self.fileId,
            var_type,
            num_vars,
            ctypes.byref(var_name_ptrs))
        return [vnp.contents.value.decode('utf8') for vnp in var_name_ptrs]

    def __ex_get_var(self, timeStep, varType, varId, blkId, numValues):
        step = ctypes.c_int(timeStep)
        var_type = ctypes.c_int(get_entity_type(varType))
        var_id = ctypes.c_int(varId)
        block_id = ctypes.c_longlong(blkId)
        num_values = ctypes.c_longlong(numValues)
        var_vals = (ctypes.c_double * num_values.value)()
        EXODUS_LIB.ex_get_var(
            self.fileId,
            step,
            var_type,
            var_id,
            block_id,
            num_values,
            var_vals)
        return var_vals

    def __ex_get_var_time(self, varType, varId, entityID, start_step, end_step):
        s_step = ctypes.c_int(start_step)
        e_step = ctypes.c_int(end_step)
        var_type = ctypes.c_int(get_entity_type(varType))
        var_id = ctypes.c_int(varId)
        entity_id = ctypes.c_longlong(entityID)
        num_steps = end_step - start_step + 1
        var_vals = (ctypes.c_double * num_steps)()
        EXODUS_LIB.ex_get_var_time(self.fileId, var_type, var_id, entity_id, s_step, e_step, var_vals)
        return var_vals

    def __ex_get_partial_var(self, timeStep, varType, varId, blkId, startIndex, numValues):
        step = ctypes.c_int(timeStep)
        var_type = ctypes.c_int(get_entity_type(varType))
        var_id = ctypes.c_int(varId)
        block_id = ctypes.c_longlong(blkId)
        start_index = ctypes.c_longlong(startIndex)
        num_values = ctypes.c_longlong(numValues)
        var_vals = (ctypes.c_double * num_values.value)()
        EXODUS_LIB.ex_get_partial_var(
            self.fileId,
            step,
            var_type,
            var_id,
            block_id,
            start_index,
            num_values,
            var_vals)
        return var_vals

    def __ex_put_var(self, timeStep, varType, varId, blkId, numValues, values):
        step = ctypes.c_int(timeStep)
        var_type = ctypes.c_int(get_entity_type(varType))
        var_id = ctypes.c_int(varId)
        block_id = ctypes.c_longlong(blkId)
        num_values = ctypes.c_longlong(numValues)
        var_vals = (ctypes.c_double * num_values.value)()
        for i in range(num_values.value):
            var_vals[i] = values[i]
        EXODUS_LIB.ex_put_var(
            self.fileId,
            step,
            var_type,
            var_id,
            block_id,
            num_values,
            var_vals)
        return True

    def __ex_put_reduction_variable_param(self, varType, numVars):
        num_vars = ctypes.c_int(numVars)
        current_num = self.__ex_get_reduction_variable_param(varType)
        if current_num.value == num_vars.value:
            # print "value already set"
            return True

        var_type = ctypes.c_int(get_entity_type(varType))
        errorInt = EXODUS_LIB.ex_put_reduction_variable_param(
            self.fileId, var_type, num_vars)
        if errorInt != 0:
            print(("ERROR code =", errorInt))
            raise Exception(
                "ERROR: ex_put_reduction_variable_param had problems."
                " This can only be called once per varType.")
        return True

    def __ex_get_reduction_variable_param(self, varType):
        var_type = ctypes.c_int(get_entity_type(varType))
        num_vars = ctypes.c_int()
        EXODUS_LIB.ex_get_reduction_variable_param(
            self.fileId, var_type, ctypes.byref(num_vars))
        return num_vars

    def __ex_get_reduction_variable_name(self, varType, varId):
        var_type = ctypes.c_int(get_entity_type(varType))
        var_id = ctypes.c_int(varId)
        name = ctypes.create_string_buffer(MAX_NAME_LENGTH + 1)
        EXODUS_LIB.ex_get_reduction_variable_name(self.fileId, var_type, var_id, name)
        return name.value.decode("utf8")

    def __ex_put_reduction_variable_name(self, varType, varId, varName):
        var_type = ctypes.c_int(get_entity_type(varType))
        var_id = ctypes.c_int(varId)
        name = ctypes.create_string_buffer(varName.encode('ascii'), MAX_NAME_LENGTH + 1)
        EXODUS_LIB.ex_put_reduction_variable_name(self.fileId, var_type, var_id, name)
        return True

    def __ex_get_reduction_variable_names(self, varType):
        num_vars = self.__ex_get_reduction_variable_param(varType)
        var_name_ptrs = (
            ctypes.POINTER(ctypes.c_char * (MAX_NAME_LENGTH + 1)) * num_vars.value)()

        for i in range(num_vars.value):
            var_name_ptrs[i] = ctypes.pointer(
                ctypes.create_string_buffer(
                    MAX_NAME_LENGTH + 1))

        var_type = ctypes.c_int(get_entity_type(varType))
        EXODUS_LIB.ex_get_reduction_variable_names(
            self.fileId,
            var_type,
            num_vars,
            ctypes.byref(var_name_ptrs))
        return [vnp.contents.value.decode('utf8') for vnp in var_name_ptrs]

    def __ex_get_reduction_vars(self, timeStep, varType, blkId, numValues):
        step = ctypes.c_int(timeStep)
        var_type = ctypes.c_int(get_entity_type(varType))
        block_id = ctypes.c_longlong(blkId)
        num_values = ctypes.c_longlong(numValues)
        var_vals = (ctypes.c_double * num_values.value)()
        EXODUS_LIB.ex_get_reduction_vars(
            self.fileId,
            step,
            var_type,
            block_id,
            num_values,
            var_vals)
        return var_vals

    def __ex_put_reduction_vars(self, timeStep, varType, blkId, numValues, values):
        step = ctypes.c_int(timeStep)
        var_type = ctypes.c_int(get_entity_type(varType))
        block_id = ctypes.c_longlong(blkId)
        num_values = ctypes.c_longlong(numValues)
        var_vals = (ctypes.c_double * num_values.value)()
        for i in range(num_values.value):
            var_vals[i] = values[i]
        EXODUS_LIB.ex_put_reduction_vars(
            self.fileId,
            step,
            var_type,
            block_id,
            num_values,
            var_vals)
        return True

    def __ex_get_side_set_node_list_len(self, object_id):
        side_set_id = ctypes.c_longlong(object_id)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API:
            side_set_node_list_len = ctypes.c_longlong(0)
        else:
            side_set_node_list_len = ctypes.c_int(0)
        EXODUS_LIB.ex_get_side_set_node_list_len(
            self.fileId, side_set_id, ctypes.byref(side_set_node_list_len))
        return side_set_node_list_len

    def __ex_get_set_param(self, objType, object_id):
        object_type = ctypes.c_int(get_entity_type(objType))
        side_set_id = ctypes.c_longlong(object_id)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API:
            num_side_in_set = ctypes.c_longlong(0)
            num_dist_fact_in_set = ctypes.c_longlong(0)
        else:
            num_side_in_set = ctypes.c_int(0)
            num_dist_fact_in_set = ctypes.c_int(0)
        EXODUS_LIB.ex_get_set_param(
            self.fileId,
            object_type,
            side_set_id,
            ctypes.byref(num_side_in_set),
            ctypes.byref(num_dist_fact_in_set))
        return int(num_side_in_set.value), int(num_dist_fact_in_set.value)

    def __ex_put_set_param(self, objType, object_id, numSides, numDistFacts):
        object_type = ctypes.c_int(get_entity_type(objType))
        side_set_id = ctypes.c_longlong(object_id)
        num_side_in_set = ctypes.c_longlong(numSides)
        num_dist_fact_in_set = ctypes.c_longlong(numDistFacts)
        EXODUS_LIB.ex_put_set_param(
            self.fileId,
            object_type,
            side_set_id,
            num_side_in_set,
            num_dist_fact_in_set)
        return True

    def __ex_get_side_set(self, sideSetId):
        side_set_id = ctypes.c_longlong(sideSetId)
        (num_side_in_set, _num_dist_fact_in_set) = self.__ex_get_set_param('EX_SIDE_SET', sideSetId)
        if num_side_in_set == 0:
            return [], []
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API:
            side_set_elem_list = (ctypes.c_longlong * num_side_in_set)()
            side_set_side_list = (ctypes.c_longlong * num_side_in_set)()
        else:
            side_set_elem_list = (ctypes.c_int * num_side_in_set)()
            side_set_side_list = (ctypes.c_int * num_side_in_set)()
        EXODUS_LIB.ex_get_set(
            self.fileId,
            ctypes.c_int(get_entity_type('EX_SIDE_SET')),
            side_set_id,
            ctypes.byref(side_set_elem_list),
            ctypes.byref(side_set_side_list))
        return side_set_elem_list, side_set_side_list

    def __ex_put_side_set(self, object_id, sideSetElements, sideSetSides):
        side_set_id = ctypes.c_longlong(object_id)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API:
            side_set_elem_list = (ctypes.c_longlong * len(sideSetElements))()
            side_set_side_list = (ctypes.c_longlong * len(sideSetSides))()
            for i, sse in enumerate(sideSetElements):
                side_set_elem_list[i] = ctypes.c_longlong(sse)
                side_set_side_list[i] = ctypes.c_longlong(sideSetSides[i])
        else:
            side_set_elem_list = (ctypes.c_int * len(sideSetElements))()
            side_set_side_list = (ctypes.c_int * len(sideSetSides))()
            for i, sse in enumerate(sideSetElements):
                side_set_elem_list[i] = ctypes.c_int(sse)
                side_set_side_list[i] = ctypes.c_int(sideSetSides[i])
        EXODUS_LIB.ex_put_set(
            self.fileId,
            ctypes.c_int(get_entity_type('EX_SIDE_SET')),
            side_set_id,
            side_set_elem_list,
            side_set_side_list)
        return True

    def __ex_get_side_set_dist_fact(self, sideSetId):
        side_set_id = ctypes.c_longlong(sideSetId)
        side_set_node_list_len = self.__ex_get_side_set_node_list_len(
            sideSetId)
        set_dfs = (ctypes.c_double * side_set_node_list_len.value)()
        EXODUS_LIB.ex_get_side_set_dist_fact(
            self.fileId, side_set_id, ctypes.byref(set_dfs))
        return set_dfs

    def __ex_put_side_set_dist_fact(self, sideSetId, sideSetDistFact):
        side_set_id = ctypes.c_longlong(sideSetId)
        side_set_dist_fact = (ctypes.c_double * len(sideSetDistFact))()
        for i, df in enumerate(sideSetDistFact):
            side_set_dist_fact[i] = ctypes.c_double(df)
        EXODUS_LIB.ex_put_side_set_dist_fact(
            self.fileId, side_set_id, side_set_dist_fact)

    def __ex_get_side_set_node_list(self, object_id):
        side_set_id = ctypes.c_longlong(object_id)
        side_set_node_list_len = self.__ex_get_side_set_node_list_len(object_id)
        (num_side_in_set, _num_dist_fact_in_set) = self.__ex_get_set_param('EX_SIDE_SET', object_id)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API:
            side_set_node_cnt_list = (ctypes.c_longlong * num_side_in_set)()
            side_set_node_list = (ctypes.c_longlong * side_set_node_list_len.value)()
        else:
            side_set_node_cnt_list = (ctypes.c_int * num_side_in_set)()
            side_set_node_list = (ctypes.c_int * side_set_node_list_len.value)()
        EXODUS_LIB.ex_get_side_set_node_list(self.fileId, side_set_id,
                                             ctypes.byref(side_set_node_cnt_list),
                                             ctypes.byref(side_set_node_list))
        return side_set_node_cnt_list, side_set_node_list

    def __ex_put_variable_param(self, varType, numVars):
        num_vars = ctypes.c_int(numVars)
        current_num = self.__ex_get_variable_param(varType)
        if current_num.value == num_vars.value:
            # print "value already set"
            return True

        var_type = ctypes.c_int(get_entity_type(varType))
        errorInt = EXODUS_LIB.ex_put_variable_param(
            self.fileId, var_type, num_vars)
        if errorInt != 0:
            print(("ERROR code =", errorInt))
            raise Exception(
                "ERROR: ex_put_variable_param had problems."
                " This can only be called once per varType.")
        return True

    def __ex_get_variable_name(self, varType, varId):
        var_type = ctypes.c_int(varType)
        var_id = ctypes.c_int(varId)
        name = ctypes.create_string_buffer(MAX_NAME_LENGTH + 1)
        EXODUS_LIB.ex_get_variable_name(self.fileId, var_type, var_id, name)
        return name.decode('utf8')

    def __ex_put_variable_name(self, varType, varId, varName):
        var_type = ctypes.c_int(get_entity_type(varType))
        var_id = ctypes.c_int(varId)
        name = ctypes.create_string_buffer(varName.encode('ascii'), MAX_NAME_LENGTH + 1)
        EXODUS_LIB.ex_put_variable_name(self.fileId, var_type, var_id, name)
        return True

    def __ex_get_attr_names(self, objType, blkId):
        obj_type = ctypes.c_int(get_entity_type(objType))
        object_id = ctypes.c_longlong(blkId)
        num_attr = ctypes.c_int(self.num_attr(blkId))
        len_name = self.__ex_inquire_int(ex_inquiry_map('EX_INQ_MAX_READ_NAME_LENGTH'))
        attr_name_ptrs = (ctypes.POINTER(ctypes.c_char * (len_name + 1)) * num_attr.value)()
        for i in range(num_attr.value):
            attr_name_ptrs[i] = ctypes.pointer(ctypes.create_string_buffer(len_name + 1))
        EXODUS_LIB.ex_get_attr_names(
            self.fileId, obj_type, object_id, ctypes.byref(attr_name_ptrs))
        return [cnp.contents.value.decode('utf8') for cnp in attr_name_ptrs]

    def __ex_put_attr_names(self, objType, blkId, varNames):
        obj_type = ctypes.c_int(get_entity_type(objType))
        object_id = ctypes.c_int(blkId)
        num_attr = ctypes.c_int(self.num_attr(blkId))
        len_name = self.__ex_inquire_int(ex_inquiry_map('EX_INQ_MAX_READ_NAME_LENGTH'))
        attr_name_ptrs = (ctypes.POINTER(ctypes.c_char * (len_name + 1)) * num_attr.value)()
        assert len(varNames) == num_attr.value
        for i in range(num_attr.value):
            attr_name_ptrs[i] = ctypes.pointer(
                ctypes.create_string_buffer(
                    varNames[i].encode('ascii'), len_name + 1))
        EXODUS_LIB.ex_put_attr_names(
            self.fileId, obj_type, object_id, ctypes.byref(attr_name_ptrs))
        return True

    def __ex_get_prop_names(self, varType, inqType):
        var_type = ctypes.c_int(get_entity_type(varType))
        num_props = ctypes.c_int(self.__ex_inquire_int(ex_inquiry_map(inqType)))
        prop_name_ptrs = (
            ctypes.POINTER(ctypes.c_char * (MAX_STR_LENGTH + 1)) * num_props.value)()
        for i in range(num_props.value):
            prop_name_ptrs[i] = ctypes.pointer(
                ctypes.create_string_buffer(
                    MAX_STR_LENGTH + 1))
        EXODUS_LIB.ex_get_prop_names(
            self.fileId, var_type, ctypes.byref(prop_name_ptrs))
        return [cnp.contents.value.decode('utf8') for cnp in prop_name_ptrs]

    def __ex_get_prop(self, objType, objId, propName):
        obj_type = ctypes.c_int(get_entity_type(objType))
        obj_id = ctypes.c_longlong(objId)
        prop_name = ctypes.create_string_buffer(propName.encode('ascii'), MAX_STR_LENGTH + 1)
        if EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API:
            prop_val = ctypes.c_longlong(0)
        else:
            prop_val = ctypes.c_int(0)
        EXODUS_LIB.ex_get_prop(
            self.fileId,
            obj_type,
            obj_id,
            ctypes.byref(prop_name),
            ctypes.byref(prop_val))
        return prop_val.value

    def __ex_put_prop(self, objType, objId, propName, propVal):
        obj_type = ctypes.c_int(get_entity_type(objType))
        obj_id = ctypes.c_longlong(objId)
        prop_name = ctypes.create_string_buffer(propName.encode('ascii'), MAX_STR_LENGTH + 1)
        prop_val = ctypes.c_longlong(propVal)
        EXODUS_LIB.ex_put_prop(
            self.fileId,
            obj_type,
            obj_id,
            ctypes.byref(prop_name),
            prop_val)
        return True

    def __ex_update(self):
        EXODUS_LIB.ex_update(self.fileId)
        return True

# --------------------------------------------------------------------
# Utility Functions
# --------------------------------------------------------------------


def collectElemConnectivity(exodusHandle, connectivity):
    """
      This function generates a list of lists that represent the element connectivity.

    >>> with exodus("file.g", "r") as exodusHandle:
    >>>     connectivity = []
    >>>     collectElemConnectivity(exodusHandle, connectivity)
    """

    if not isinstance(connectivity, list):
        raise Exception(
            "ERROR: connectivity is not a list in call to collectElemConnectivity().")
    if connectivity:
        raise Exception(
            "ERROR: connectivity is not empty in call to collectElemConnectivity().")

    blockIds = exodusHandle.get_ids('EX_ELEM_BLOCK')
    for blId in blockIds:
        (elem_block_conn, num_elem, num_nodes) = exodusHandle.get_elem_connectivity(blId)
        for k in range(num_elem):
            i = k * num_nodes
            j = i + num_nodes
            local_elem_conn = elem_block_conn[i:j]
            connectivity.append(local_elem_conn)


def collectLocalNodeToLocalElems(
        exodusHandle,
        connectivity,
        localNodeToLocalElems):
    """
      This function generates a list of lists to go from local node id
      to local elem id.

    >>> connectivity = [] ## If this is not empty it will assume it is already filled.
    >>> localNodeToLocalElems = []
    >>> with exodus("file.g", "r") as exodusHandle:
    >>>     collectLocalNodeToLocalElems(exodusHandle, connectivity, localNodeToLocalElems)
    """

    if not isinstance(connectivity, list):
        raise Exception(
            "ERROR: connectivity is not a list in call to collectLocalNodeToLocalElems().")
    if not isinstance(localNodeToLocalElems, list):
        raise Exception(
            "ERROR: localNodeToLocalElems is not a list in call to collectLocalNodeToLocalElems().")
    if localNodeToLocalElems:
        raise Exception(
            "ERROR: localNodeToLocalElems is not empty in call to collectLocalNodeToLocalElems().")

    if not connectivity:
        collectElemConnectivity(exodusHandle, connectivity)

    numNodes = exodusHandle.num_nodes()
    for _i in range(numNodes + 1):
        localNodeToLocalElems.append([])

    localElemId = 0
    for local_elem_conn in connectivity:
        for n in local_elem_conn:
            localNodeToLocalElems[n].append(localElemId)
        localElemId = localElemId + 1


def collectLocalElemToLocalElems(
        exodusHandle,
        connectivity,
        localNodeToLocalElems,
        localElemToLocalElems):
    """
      This function generates a list of lists to go from local elem id
      to connected local elem ids.

    >>> connectivity = [] ## If this is not empty it will assume it is already filled.
    >>> localNodeToLocalElems = [] ## If this is not empty it will assume it is already filled.
    >>> localElemToLocalElems = []
    >>> with exodus("file.g", "r") as exodusHandle:
    >>>     collectLocalElemToLocalElems(exodusHandle, connectivity, localNodeToLocalElems,
    ...                              localElemToLocalElems)
    """

    if not isinstance(connectivity, list):
        raise Exception(
            "ERROR: connectivity is not a list in call to collectLocalElemToLocalElems().")
    if not isinstance(localNodeToLocalElems, list):
        raise Exception(
            "ERROR: localNodeToLocalElems is not a list in call to collectLocalElemToLocalElems().")
    if not isinstance(localElemToLocalElems, list):
        raise Exception(
            "ERROR: localElemToLocalElems is not a list in call to collectLocalElemToLocalElems().")
    if localElemToLocalElems:
        raise Exception(
            "ERROR: localElemToLocalElems is not empty in call to collectLocalElemToLocalElems().")

    if not connectivity:
        collectElemConnectivity(exodusHandle, connectivity)
    if not localNodeToLocalElems:
        collectLocalNodeToLocalElems(
            exodusHandle, connectivity, localNodeToLocalElems)

    numElems = exodusHandle.num_elems()
    for _i in range(numElems):
        localElemToLocalElems.append([])
    for localElemId in range(numElems):
        nodeList = list(connectivity[localElemId])
        newConnectedElems = []
        for n in nodeList:
            newConnectedElems.extend(iter(localNodeToLocalElems[n]))
        localElemToLocalElems[localElemId] = list(set(newConnectedElems))


def copy_mesh(fromFileName, toFileName, exoFromObj=None, additionalElementAttributes=None, array_type='ctype'):
    """
    Copies the mesh data from an existing exodus database to a new exodus
    database.

    Parameters
    ----------
    fromFileName : string
        File name of the exodus mesh to be copied
    toFileName : string
        File name of the new exodus mesh
    exoFromObj : exodus object, optional
        Exodus object to be copied from.  If an exodus object is supplied, the
        fromFileName string will be ignored.
    additionalElementAttributes : list
        list of element attribute names to add to all blocks or tuples
        ( name, blkIds ) where name is the element attribute to add and blkIds is
        a list of blkIds to add it to.
    array_type : 'ctype' | 'numpy'
        Specifies whether arrays will be imported and copied as ctype or numpy
        arrays.  (This option should make no difference to the user, but it can
        be used by developers to test whether the commands within this function
        handle both array types correctly.)

    Returns
    -------
    exo_to : exodus object
        New exodus mesh

    Note
    -----
    This function also allows one to add new element attributes during the copy
    process.  The number of element attributes is permanently set when the
    block is created, meaning new element attributes can only be added to an
    existing mesh by copying it to a new mesh.  The values of the element
    attributes are set to their defaults so that the user can populate them
    later.
    """
    if additionalElementAttributes is None:
        additionalElementAttributes = []
    # If the user did not supply a exodus object to copy from, attempt to read an
    # exodus database with the name "fromFileName"
    if exoFromObj is None:
        exoFrom = exodus(fromFileName, "r", array_type=array_type)
    else:
        exoFrom = exoFromObj

    if os.path.isfile(toFileName):
        raise Exception(
            "ERROR: ",
            toFileName,
            " file already exists cowardly exiting instead of overwriting in call to copy_mesh().")

    title = exoFrom.title().encode('ascii')
    ex_pars = ex_init_params(num_dim=exoFrom.num_dimensions(),
                             num_nodes=exoFrom.num_nodes(),
                             num_elem=exoFrom.num_elems(),
                             num_elem_blk=exoFrom.num_blks(),
                             num_node_sets=exoFrom.num_node_sets(),
                             num_side_sets=exoFrom.num_side_sets(),
                             num_assembly=exoFrom.num_assembly(),
                             num_blob=exoFrom.num_blob())

    exo_to = exodus(toFileName, mode="w", array_type=array_type,
                    title=title, init_params=ex_pars)

    debugPrint = False
    if debugPrint:
        print("Transfer QA records")
    qaRecords = exoFrom.get_qa_records()
    exo_to.put_qa_records(qaRecords)

    if debugPrint:
        print("Transfer Nodal Coordinates and Names")
    exo_to.put_coord_names(exoFrom.get_coord_names())
    (xCoords, yCoords, zCoords) = exoFrom.get_coords()
    exo_to.put_coords(xCoords, yCoords, zCoords)

    if debugPrint:
        print("Transfer Node Id Map")
    nodeIdMap = exoFrom.get_node_id_map()
    exo_to.put_node_id_map(nodeIdMap)

    if debugPrint:
        print("Construct mapping from block ID to element attribute data")
    # The exodus library does not provide a way to add only new element
    # attributes, so we must collect both the new and the old element
    # attributes
    e_attr_names = {}
    e_attr_vals = {}
    # Collect the old element attribute names and the number of elements in each
    # block
    blk_ids = exoFrom.get_ids('EX_ELEM_BLOCK')
    blk_num_elem = {}
    for blk_id in blk_ids:
        (elemType, numElem, nodesPerElem, numAttr) = exoFrom.elem_blk_info(blk_id)
        e_attr_names[blk_id] = []
        e_attr_vals[blk_id] = []
        if numAttr > 0:
            e_attr_names[blk_id].extend(
                exoFrom.get_attribute_names('EX_ELEM_BLOCK', blk_id))
            e_attr_vals[blk_id].extend(exoFrom.get_elem_attr(blk_id))
        blk_num_elem[blk_id] = numElem
    # Collect the new element attribute names
    # (The new names are mapped from "attribute name" to "list of block IDs that
    # contain that attribute".  We need to have them be mapped as "block ID" to
    # "list of attribute names contained in that block".)
    for item in additionalElementAttributes:
        if isinstance(item, tuple):
            e_attr_name = item[0]
            e_attr_blk_ids = item[1]
        elif isinstance(item, str):
            e_attr_name = item
            e_attr_blk_ids = blk_ids
        else:
            print(f"Warning additional element attribute item {item} is not right type to add.")

            print("should be a string or tuple, skipping")
        for blk_id in e_attr_blk_ids:
            if blk_id in blk_ids:
                e_attr_names[blk_id].append(e_attr_name)
                # Concatenate all element attribute values into a single big list,
                # because that is format required by exo.put_elem_attr().
                e_attr_vals[blk_id].extend([0.0] * blk_num_elem[blk_id])

    if debugPrint:
        print("Transfer Element Data")
    blkIds = exoFrom.get_ids('EX_ELEM_BLOCK')
    for blkId in blkIds:
        (elemType, numElem, nodesPerElem, _oldnumAttr) = exoFrom.elem_blk_info(blkId)
        numAttr = len(e_attr_names[blkId])
        exo_to.put_elem_blk_info(blkId, elemType, numElem, nodesPerElem, numAttr)
        (connectivity, numElem, nodesPerElem) = exoFrom.get_elem_connectivity(blkId)
        exo_to.put_elem_connectivity(blkId, connectivity)
        if numAttr > 0:
            exo_to.put_attribute_names('EX_ELEM_BLOCK', blkId, e_attr_names[blkId])
            exo_to.put_elem_attr(blkId, e_attr_vals[blkId])
        elemProps = exoFrom.get_element_property_names()
        for elemProp in elemProps:
            propVal = exoFrom.get_element_property_value(blkId, elemProp)
            if elemProp == "ID" and propVal == blkId:
                continue
            else:
                exo_to.put_element_property_value(blkId, elemProp, propVal)
        blockName = exoFrom.get_name('EX_ELEM_BLOCK', blkId)
        exo_to.put_name('EX_ELEM_BLOCK', blkId, blockName)

    if debugPrint:
        print("Transfer Element Id Map")
    elemIdMap = exoFrom.get_elem_id_map()
    exo_to.put_elem_id_map(elemIdMap)

    if debugPrint:
        print("Transfer Node Sets")
    if exoFrom.num_node_sets() > 0:
        nodeSetProps = exoFrom.get_node_set_property_names()
        nodeSetIds = exoFrom.get_ids('EX_NODE_SET')
        for nsId in nodeSetIds:
            (numSetNodes, numSetDistFacts) = exoFrom.get_set_params('EX_NODE_SET', nsId)
            exo_to.put_node_set_params(nsId, numSetNodes, numSetDistFacts)
            nsNodes = exoFrom.get_node_set_nodes(nsId)
            exo_to.put_node_set(nsId, nsNodes)
            if numSetDistFacts > 0:
                nsDF = exoFrom.get_node_set_dist_facts(nsId)
                exo_to.put_node_set_dist_fact(nsId, nsDF)
            nodeSetName = exoFrom.get_name('EX_NODE_SET', nsId)
            exo_to.put_name('EX_NODE_SET', nsId, nodeSetName)
            for nodeSetProp in nodeSetProps:
                propVal = exoFrom.get_node_set_property_value(
                    nsId, nodeSetProp)
                if nodeSetProp == "ID" and propVal == nsId:
                    continue
                else:
                    exo_to.put_node_set_property_value(
                        nsId, nodeSetProp, propVal)

    if debugPrint:
        print("Transfer Side Sets")
    if exoFrom.num_side_sets() > 0:
        sideSetProps = exoFrom.get_side_set_property_names()
        sideSetIds = exoFrom.get_ids('EX_SIDE_SET')
        for ssId in sideSetIds:
            (numSetSides, numSetDistFacts) = exoFrom.get_set_params('EX_SIDE_SET', ssId)
            exo_to.put_side_set_params(ssId, numSetSides, numSetDistFacts)
            (elemList, sideList) = exoFrom.get_side_set(ssId)
            exo_to.put_side_set(ssId, elemList, sideList)
            if numSetDistFacts > 0:
                ssDF = exoFrom.get_side_set_dist_fact(ssId)
                exo_to.put_side_set_dist_fact(ssId, ssDF)
            sideSetName = exoFrom.get_name('EX_SIDE_SET', ssId)
            exo_to.put_name('EX_SIDE_SET', ssId, sideSetName)
            for sideSetProp in sideSetProps:
                propVal = exoFrom.get_side_set_property_value(
                    ssId, sideSetProp)
                if sideSetProp == "ID" and propVal == ssId:
                    continue
                else:
                    exo_to.put_side_set_property_value(
                        ssId, sideSetProp, propVal)

    # If the user did not supply an exodus object to copy from, then close the
    # database.
    if exoFromObj is None:
        exoFrom.close()

    return exo_to


def transfer_variables(exoFrom, exo_to, array_type='ctype', additionalGlobalVariables=None, additionalNodalVariables=None, additionalElementVariables=None, additionalNodeSetVariables=None, additionalSideSetVariables=None):
    """
    This function transfers variables from `exoFrom` to `exo_to` and allows
    additional variables to be added with `additionalGlobalVariables`,
    `additionalNodalVariables`, and `additionalElementVariables`.  Additional
    variables values are set to their defaults so that the user can populate
    them later.

    Parameters
    ----------
    exoFrom : exodus database object
        exodus object to transfer from

    exo_to : exodus database object
        exodus object to transfer to

    additionalGlobalVariables : list
        list of global variable names to add.

    additionalNodalVariables : list
        list of nodal variable names to add.

    additionalElementVariables : list
        should be a list of element variable names to add to all blocks or
        tuples `(name, blkIds)` where `name` is the element variable to add
        and `blkIds` is a list of block ids to add it to.
    """
    if additionalGlobalVariables is None:
        additionalGlobalVariables = []
    if additionalNodalVariables is None:
        additionalNodalVariables = []
    if additionalElementVariables is None:
        additionalElementVariables = []
    if additionalNodeSetVariables is None:
        additionalNodeSetVariables = []
    if additionalSideSetVariables is None:
        additionalSideSetVariables = []
    # IDEA: It may make sense to make transfer_variables() strictly transfer
    # variables, and use add_variables() to add new variables.  Alternatively,
    # add_variables() could be called within transfer_variables() to add
    # new variables to the exo_to database.  Either way, we should minimize
    # duplicate code if possible.

    debugPrint = False

    if not isinstance(additionalGlobalVariables, list):
        raise Exception("ERROR: additionalGlobalVariables is not a list.")
    if not isinstance(additionalNodalVariables, list):
        raise Exception("ERROR: additionalNodalVariables is not a list.")
    if not isinstance(additionalElementVariables, list):
        raise Exception("ERROR: additionalElementVariables is not a list.")
    if not isinstance(additionalNodeSetVariables, list):
        raise Exception("ERROR: additionalNodeSetVariables is not a list.")
    if not isinstance(additionalSideSetVariables, list):
        raise Exception("ERROR: additionalSideSetVariables is not a list.")

    if debugPrint:
        print("Transfer Info records")
    numInfoRecs = exoFrom.num_info_records()
    if numInfoRecs > 0:
        infoRecs = exoFrom.get_info_records()
        exo_to.put_info_records(infoRecs)
    if debugPrint:
        print("Transfer time values")

    nSteps = exoFrom.num_times()
    if nSteps == 0:
        return exo_to

    timeVals = exoFrom.get_times()
    for step in range(nSteps):
        exo_to.put_time(step + 1, timeVals[step])

    if debugPrint:
        print("Add Global Variables")
    nNewGlobalVars = len(additionalGlobalVariables)
    nGlobalVars = exoFrom.get_variable_number('EX_GLOBAL') + nNewGlobalVars
    defaultNewVarVals = [0.0 for _ in range(nNewGlobalVars)]
    if nGlobalVars > 0:
        exo_to.set_variable_number('EX_GLOBAL', nGlobalVars)
        gVarNames = exoFrom.get_variable_names('EX_GLOBAL')
        gVarNames.extend(additionalGlobalVariables)
        for nameIndex in range(nGlobalVars):
            globalVarName = gVarNames[nameIndex]
            exo_to.put_variable_name('EX_GLOBAL', globalVarName, nameIndex + 1)
        for step in range(nSteps):
            gValues = exoFrom.get_all_global_variable_values(step + 1)
            if array_type == 'numpy':
                gValues = exo_to.np.append(gValues, defaultNewVarVals)
            else:
                gValues.extend(defaultNewVarVals)
            exo_to.put_all_global_variable_values(step + 1, gValues)

    if debugPrint:
        print("Add Nodal Variables")
    nNewNodalVars = len(additionalNodalVariables)
    nOrigNodalVars = exoFrom.get_variable_number('EX_NODAL')
    nNodalVars = nOrigNodalVars + nNewNodalVars
    if nNodalVars > 0:
        exo_to.set_variable_number('EX_NODAL', nNodalVars)
        nVarNames = exoFrom.get_variable_names('EX_NODAL')
        nVarNames.extend(additionalNodalVariables)
        for nameIndex in range(nNodalVars):
            nodalVarName = nVarNames[nameIndex]
            exo_to.put_variable_name('EX_NODAL', nodalVarName, nameIndex + 1)
            if nameIndex < nOrigNodalVars:
                for step in range(nSteps):
                    nValues = exoFrom.get_node_variable_values(
                        nodalVarName, step + 1)
                    exo_to.put_node_variable_values(
                        nodalVarName, step + 1, nValues)

    internal_transfer_variables(exoFrom, exo_to, 'EX_ELEM_BLOCK', additionalElementVariables, debugPrint)
    internal_transfer_variables(exoFrom, exo_to, 'EX_NODE_SET', additionalNodeSetVariables, debugPrint)
    internal_transfer_variables(exoFrom, exo_to, 'EX_SIDE_SET', additionalSideSetVariables, debugPrint)
    return exo_to


def internal_transfer_variables(exoFrom, exo_to, obj_type, additionalVariables, debugPrint):
    """ Internal support function for :py:func:`exodus.transfer_variables` """
    if debugPrint:
        print("Construct Truth Table for additionalVariables")
    blkIds = exoFrom.get_ids(obj_type)
    numBlks = len(blkIds)
    newVariableNames = []
    newVariableBlocks = []
    for item in additionalVariables:
        if isinstance(item, tuple):
            newVariableNames.append(item[0])
            inBlks = [blkId for blkId in item[1] if blkId in blkIds]
            newVariableBlocks.append(inBlks)
        elif isinstance(item, str):
            newVariableNames.append(item)
            newVariableBlocks.append(blkIds)
        else:
            print(("Warning additionalVariable item ",
                   item, " is not right type to add."))
            print("should be a string or tuple, skipping")

    nSteps = exoFrom.num_times()
    if debugPrint:
        print("Add Variables")
    nNewVars = len(newVariableNames)
    nOrigVars = exoFrom.get_variable_number(obj_type)
    nVars = nOrigVars + nNewVars
    if nVars > 0:
        exo_to.set_variable_number(obj_type, nVars)
        origVarNames = exoFrom.get_variable_names(obj_type)
        varNames = []
        varNames.extend(origVarNames)
        varNames.extend(newVariableNames)

        truthTable = []
        if nOrigVars > 0:
            truthTable = exoFrom.get_variable_truth_table(obj_type)
        if nNewVars > 0:
            newTruth = []
            for j in range(numBlks):

                for k in range(nOrigVars):
                    index = j * nOrigVars + k
                    newTruth.append(truthTable[index])
                for m in range(nNewVars):
                    if blkIds[j] in newVariableBlocks[m]:
                        newTruth.append(True)
                    else:
                        newTruth.append(False)
            truthTable = newTruth
        exo_to.set_variable_truth_table(obj_type, truthTable)
        for nameIndex in range(nVars):
            varName = varNames[nameIndex]
            exo_to.put_variable_name(obj_type, varName, nameIndex + 1)
        truthIndex = 0
        for blkId in blkIds:
            for origVarName in origVarNames:
                if truthTable[truthIndex]:
                    for step in range(nSteps):
                        eValues = exoFrom.get_variable_values(obj_type, blkId,
                                                              origVarName, step + 1)
                        exo_to.put_variable_values(obj_type, blkId,
                                                   origVarName, step + 1, eValues)
                truthIndex = truthIndex + 1
            truthIndex = truthIndex + nNewVars


def add_variables(exo, global_vars=None, nodal_vars=None, element_vars=None, node_set_vars=None, side_set_vars=None):
    """
    This function adds variables to the exodus object.  The values of the
    variables are set to their defaults so that the user can populate them later.

    Parameters
    ----------
    exo : exodus database object
        Exodus database that variables will be added to.
    global_vars : list
        global variable names to add.
    nodal_vars : list
        nodal variable names to add.
    element_vars : list
        list of element variable names to add to all blocks or tuples
        ( name, blkIds ) where name is the element variable to add and blkIds is
        a list of blkIds to add it to.
    node_set_vars : list
        list of node set variable names to add to all sets or tuples
        ( name, setIds ) where name is the node set variable to add and `setIds` is
        a list of `setIds` to add it to.
    side_set_vars : list
        list of side set variable names to add to all sets or tuples
        ( name, setIds ) where name is the side set variable to add and `setIds` is
        a list of `setIds` to add it to.

    Returns
    -------
    exo : exodus database object
        Exodus database with variables added to it.  (The values of the variables
        are set to their defaults so that the user can populate them later.)

    Note
    -----
    This function does not allow one to add element attributes to an exodus
    database.  See :py:func:`exodus.copy_mesh` function for that capability.
    """
    if global_vars is None:
        global_vars = []
    if nodal_vars is None:
        nodal_vars = []
    if element_vars is None:
        element_vars = []
    if node_set_vars is None:
        node_set_vars = []
    if side_set_vars is None:
        side_set_vars = []
    debugPrint = False

    if not isinstance(global_vars, list):
        raise Exception("ERROR: global_vars is not a list.")
    if not isinstance(nodal_vars, list):
        raise Exception("ERROR: nodal_vars is not a list.")
    if not isinstance(element_vars, list):
        raise Exception("ERROR: element_vars is not a list.")
    if not isinstance(node_set_vars, list):
        raise Exception("ERROR: node_set_vars is not a list.")
    if not isinstance(side_set_vars, list):
        raise Exception("ERROR: side_set_vars is not a list.")

    if exo.modeChar == 'r':
        raise Exception(
            "ERROR: variables cannot be added to an exodus object in read only mode")

    if debugPrint:
        print("Add Global Variables")
    n_new_vars = len(global_vars)
    n_old_vars = exo.get_variable_number('EX_GLOBAL')
    n_vars = n_old_vars + n_new_vars
    default_vals = [0.0] * n_new_vars
    if n_new_vars > 0:
        exo.set_variable_number('EX_GLOBAL', n_vars)
        for i, var_name in enumerate(global_vars):
            exo.put_variable_name('EX_GLOBAL', var_name, n_old_vars + i + 1)
        # One might wish to put all the values for a given global variable in the
        # database at once, but exo.put_global_variable_value() ends up loading
        # all the global variables for a given step and then putting them all back
        # in, so we might as well just use
        # exo.put_all_global_variable_values().
        nSteps = exo.num_times()
        for step in range(nSteps):
            gValues = exo.get_all_global_variable_values(step + 1)
            if exo.use_numpy:
                gValues = exo.np.append(gValues, default_vals)
            else:
                gValues.extend(default_vals)
            exo.put_all_global_variable_values(step + 1, gValues)

    if debugPrint:
        print("Add Nodal Variables")
    n_new_vars = len(nodal_vars)
    n_old_vars = exo.get_variable_number('EX_NODAL')
    n_vars = n_old_vars + n_new_vars
    if n_new_vars > 0:
        exo.set_variable_number('EX_NODAL', n_vars)
        for i, var_name in enumerate(nodal_vars):
            exo.put_variable_name('EX_NODAL', var_name, i + n_old_vars + 1)

    internal_add_variables(exo, 'EX_ELEM_BLOCK', element_vars, debugPrint)
    internal_add_variables(exo, 'EX_NODE_SET', node_set_vars, debugPrint)
    internal_add_variables(exo, 'EX_SIDE_SET', side_set_vars, debugPrint)

    return exo


def internal_add_variables(exo, obj_type, entvars, debugPrint):
    """ Internal support function for :py:func:`exodus.add_variables` """

    if len(entvars) == 0:
        return

    if debugPrint:
        print("Construct Truth Table for additional variables")

    new_var_blks = []
    blk_ids = exo.get_ids(obj_type)
    new_var_names = []
    for item in entvars:
        if isinstance(item, tuple):
            new_var_names.append(item[0])
            in_blks = [blk_id for blk_id in item[1] if blk_id in blk_ids]
            new_var_blks.append(in_blks)
        elif isinstance(item, str):
            new_var_names.append(item)
            new_var_blks.append(blk_ids)
        else:
            print(f"Warning additional variable item {item} is not right type to add.")
            print("should be a string or tuple, skipping")

    if debugPrint:
        print("Add Variables")
    n_new_vars = len(new_var_names)
    if n_new_vars > 0:
        n_old_vars = exo.get_variable_number(obj_type)
        n_vars = n_old_vars + n_new_vars
        exo.set_variable_number(obj_type, n_vars)
        old_truth_table = []
        if n_old_vars > 0:
            old_truth_table = exo.get_variable_truth_table(obj_type)
        truth_table = []
        n_blks = len(blk_ids)
        for j in range(n_blks):
            for k in range(n_old_vars):
                ndx = j * n_old_vars + k
                truth_table.append(old_truth_table[ndx])
            for m in range(n_new_vars):
                if blk_ids[j] in new_var_blks[m]:
                    truth_table.append(True)
                else:
                    truth_table.append(False)
        exo.set_variable_truth_table(obj_type, truth_table)
        for i, var_name in enumerate(new_var_names):
            exo.put_variable_name(obj_type, var_name, n_old_vars + i + 1)


def copyTransfer(fromFileName, toFileName, array_type='ctype', additionalGlobalVariables=None, additionalNodalVariables=None, additionalElementVariables=None, additionalNodeSetVariables=None, additionalSideSetVariables=None, additionalElementAttributes=None):
    """
    This function creates an exodus file `toFileName` and copies
    everything from exodus file `fromFileName` returning a file handle
    to `toFileName`.

    Additional space is allocated for `additionalGlobalVariables`,
    `additionalNodalVariables` and `additionalElementVariables` if
    specified.

    Parameters
    ----------
    fromFileName : string
        File name of the exodus mesh to be copied
    toFileName : string
        File name of the new exodus mesh
    array_type : 'ctype' | 'numpy'
        Specifies whether arrays will be imported and copied as ctype or numpy
        arrays.  (This option should make no difference to the user, but it can
        be used by developers to test whether the commands within this function
        handle both array types correctly.)
    additionalGlobalVariables : list
        list of global variable names to add.
    additionalNodalVariables: list
        list of nodal variable names to add.
    additionalElementVariables: list
        should be a list of element variable
        names to add to all blocks or tuples (name, blkIds) where name is
        the element variable to add and blkIds is a list of blkIds to add
        it to.
    additionalElementAttributes: list
        list of element attribute names to
        add to all blocks or tuples ( name, blkIds ) where name is the
        element attribute to add and blkIds is a list of blkIds to add it
        to.


    >>> fromFileName = "input.e"
    >>> toFileName = "output.e"
    >>> addGlobalVariables = [] ## Do not add any new global variables
    >>> ## Add node_dummy1 and node_dummy2 as new node variables
    >>> addNodeVariables = ["node_dummy1", "node_dummy2"]
    >>> ## Add elem_dummy1 on blkIds 1, 2, 3 and elem_dummy2 on all blocks
    >>> addElementVariables = [ ("elem_dummy1", [1, 2, 3]), "elem_dummy2" ]
    >>> ## Add elem_attr_dummy1 on blkIds 1,2,3 and elem_attr_dummy2 on all blocks
    >>> addElementAttributes = [ ("elem_attr_dummy1",[1,2,3]), "elem_attr_dummy2" ]
    >>>
    >>> toFileHandle = copyTransfer(fromFileName,toFileName,addGlobalVariables,addNodeVariables,
    ...                             addElementVariables,addElementAttributes)
    ...
    >>> ## Fill in new variables
    ...
    >>> toFileHandle.close()
    """

    if additionalGlobalVariables is None:
        additionalGlobalVariables = []
    if additionalNodalVariables is None:
        additionalNodalVariables = []
    if additionalElementVariables is None:
        additionalElementVariables = []
    if additionalNodeSetVariables is None:
        additionalNodeSetVariables = []
    if additionalSideSetVariables is None:
        additionalSideSetVariables = []
    if additionalElementAttributes is None:
        additionalElementAttributes = []
    with exodus(fromFileName, "r", array_type=array_type) as exoFrom:
        exo_to = copy_mesh(fromFileName, toFileName, exoFromObj=exoFrom,
                           additionalElementAttributes=additionalElementAttributes,
                           array_type=array_type)

        exo_to = transfer_variables(
            exoFrom,
            exo_to,
            additionalGlobalVariables=additionalGlobalVariables,
            additionalNodalVariables=additionalNodalVariables,
            additionalElementVariables=additionalElementVariables,
            additionalNodeSetVariables=additionalNodeSetVariables,
            additionalSideSetVariables=additionalSideSetVariables,
            array_type=array_type)

    return exo_to


def ctype_to_numpy(exo, c_array):
    """
    Converts a c-type array into a numpy array

    Parameters
    ----------
    exo : exodus object
        exodus database object initialized with the option `array_type = 'numpy'`
    c_array : c-type array
        c-type array to be converted into a numpy array

    Returns
    -------
    np_array : numpy array
        Numpy array converted from c-type array
    """
    # ctypes currently produce invalid PEP 3118 type codes, which causes numpy
    # to issue a warning.  This is a bug and can be ignored.
    # http://stackoverflow.com/questions/4964101/pep-3118-warning-when-using-ctypes-array-as-numpy-array
    if not c_array:
        return exo.np.array([])

    with exo.warnings.catch_warnings():
        exo.warnings.simplefilter('ignore')
        np_array = exo.np.ctypeslib.as_array(c_array)
    return np_array
