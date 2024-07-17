"""
Exomerge is a lightweight Python interface for manipulating ExodusII files.

Copyright(C) 1999-2020, 2023, 2024 National Technology & Engineering Solutions
of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
NTESS, the U.S. Government retains certain rights in this software.

See packages/seacas/LICENSE for details


Author: Tim Kostka (tdkostk@sandia.gov)
Created: May 3, 2012

Simple example:
>>> import exomerge
>>> model = exomerge.import_model('results.e')
>>> model.delete_element_block(1)
>>> model.export_model('most_results.e')

Documentation can be accessed online through Python's pydoc utility:
$ pydoc exomerge
$ pydoc exomerge.ExodusModel.delete_element_block

Documentation and sample scripts can be found in the following report.
"Exomerge User's Manual: A lightweight Python interface for manipulating
Exodus files" (SAND2013-0725)

For additional information or to provide feedback or feature requests, please
contact the author.

This is the python-2 version
"""

# import standard modules
import string
import textwrap
import sys
import traceback
import os
import inspect
import re
import datetime
import itertools
import math
import struct
import bisect
import colorsys
import difflib
import operator

if sys.version_info[0] >= 3:
    raise Exception("Python-2 version. If using python-3, try `import exomerge3 as exomerge`")

# import exodus module
# (exodus.py should be in the same directory as this file)
import exodus2 as exodus

# informal version number of this module
__version__ = 8.6
VERSION = __version__

# contact person for issues
CONTACT = 'Tim Kostka <tdkostk@sandia.gov>'

# show the banner on first use
SHOW_BANNER = True

# if true, will crash if warnings are generated
EXIT_ON_WARNING = False

# if true, will suppress output from exodus.py
# this will not suppress the banner output
SUPPRESS_EXODUS_OUTPUT = True

# a list of deprecated or renamed functions
# When the user calls one of these, it will issue a warning message
# saying it is deprecated.  If it has simply been renamed, we will call that
# function.  If it has been deleted, we will error out.
# This is a dict where the key is the deprecated function name and the value
# is the renamed function, or None.
# e.g. DEPRECATED_FUNCTIONS['build_hex8_cube'] = 'create_hex8_cube'
DEPRECATED_FUNCTIONS = dict()


class DummyFile(object):
    """Dummy class used to suppress output to stdout."""

    def write(self, x):
        """Ignore the write command."""
        pass


def import_model(filename, *args, **kwargs):
    """
    Load information from an ExodusII file.

    This function is a wrapper around 'ExodusModel.import_model(...)' and is
    provided for convenience.  Internally, this is equivalent to executing the
    following two statements.

    >>> model = ExodusModel()
    >>> model.import_model(...)

    See 'ExodusModel.import_model' for additional information.

    """
    model = ExodusModel()
    model.import_model(filename, *args, **kwargs)
    return model


class ExodusModel(object):
    """This class holds all information stored within an ExodusII file."""

    # define translation from two faces to a 3D cohesive element
    COHESIVE_FORMULA = dict()
    COHESIVE_FORMULA['quad4'] = ['hex8', (0, 1, 2, 3, 4, 5, 6, 7)]
    COHESIVE_FORMULA['tri3'] = ['wedge6', (0, 1, 2, 3, 4, 5)]
    COHESIVE_FORMULA['tri6'] = ['wedge12', (0, 1, 2, 6, 7, 8,
                                            3, 4, 5, 9, 10, 11)]

    # define calculation for element length/area/volume
    VOLUME_FORMULA = dict()
    VOLUME_FORMULA['line2'] = [1.0, (0, 1)]
    VOLUME_FORMULA['line3'] = VOLUME_FORMULA['line2']
    VOLUME_FORMULA['tri3'] = [0.5, (0, 1), (0, 2)]
    VOLUME_FORMULA['tri6'] = VOLUME_FORMULA['tri3']
    VOLUME_FORMULA['quad4'] = [0.5, (0, 2), (1, 3)]
    VOLUME_FORMULA['quad6'] = VOLUME_FORMULA['quad4']
    VOLUME_FORMULA['quad8'] = VOLUME_FORMULA['quad4']
    VOLUME_FORMULA['tet4'] = [1.0 / 6.0, (0, 1), (0, 2), (0, 3)]
    VOLUME_FORMULA['tet10'] = VOLUME_FORMULA['tet4']
    VOLUME_FORMULA['wedge6'] = [0.5,
                                ((0, 3), (1, 4)),
                                ((0, 3), (2, 5)),
                                ((0, 1, 2), (3, 4, 5))]
    VOLUME_FORMULA['wedge12'] = VOLUME_FORMULA['wedge6']
    VOLUME_FORMULA['wedge15'] = VOLUME_FORMULA['wedge6']
    VOLUME_FORMULA['wedge16'] = VOLUME_FORMULA['wedge6']
    VOLUME_FORMULA['hex8'] = [1.0,
                              ((0, 3, 4, 7), (1, 2, 5, 6)),
                              ((0, 1, 4, 5), (2, 3, 6, 7)),
                              ((0, 1, 2, 3), (4, 5, 6, 7))]
    VOLUME_FORMULA['hex20'] = VOLUME_FORMULA['hex8']

    # define mapping of the external faces of an element
    FACE_MAPPING = dict()
    FACE_MAPPING['hex8'] = [('quad4', (0, 1, 5, 4)),
                            ('quad4', (1, 2, 6, 5)),
                            ('quad4', (2, 3, 7, 6)),
                            ('quad4', (3, 0, 4, 7)),
                            ('quad4', (0, 3, 2, 1)),
                            ('quad4', (4, 5, 6, 7))]
    FACE_MAPPING['hex20'] = [('quad8', (0, 1, 5, 4, 8, 13, 16, 12)),
                             ('quad8', (1, 2, 6, 5, 9, 14, 17, 13)),
                             ('quad8', (2, 3, 7, 6, 10, 15, 18, 14)),
                             ('quad8', (3, 0, 4, 7, 11, 12, 19, 15)),
                             ('quad8', (0, 3, 2, 1, 11, 10, 9, 8)),
                             ('quad8', (4, 5, 6, 7, 16, 17, 18, 19))]
    FACE_MAPPING['tet4'] = [('tri3', (0, 1, 3)),
                            ('tri3', (1, 2, 3)),
                            ('tri3', (2, 0, 3)),
                            ('tri3', (0, 2, 1))]
    FACE_MAPPING['tet10'] = [('tri6', (0, 1, 3, 4, 8, 7)),
                             ('tri6', (1, 2, 3, 5, 9, 8)),
                             ('tri6', (2, 0, 3, 6, 7, 9)),
                             ('tri6', (0, 2, 1, 6, 5, 4))]
    FACE_MAPPING['wedge6'] = [('quad4', (0, 1, 4, 3)),
                              ('quad4', (1, 2, 5, 4)),
                              ('quad4', (2, 0, 3, 5)),
                              ('tri3', (0, 2, 1)),
                              ('tri3', (3, 4, 5))]
    FACE_MAPPING['wedge12'] = [('quad6', (0, 1, 4, 3, 6, 9)),
                               ('quad6', (1, 2, 5, 4, 7, 10)),
                               ('quad6', (2, 0, 3, 5, 8, 11)),
                               ('tri6', (0, 2, 1, 8, 7, 6)),
                               ('tri6', (3, 4, 5, 11, 10, 9))]
    FACE_MAPPING['wedge15'] = [('quad8', (0, 1, 4, 3, 6, 10, 12, 9)),
                               ('quad8', (1, 2, 5, 4, 7, 11, 13, 10)),
                               ('quad8', (2, 0, 3, 5, 8, 9, 14, 11)),
                               ('tri6', (0, 2, 1, 8, 7, 6)),
                               ('tri6', (3, 4, 5, 12, 13, 14))]
    FACE_MAPPING['wedge16'] = FACE_MAPPING['wedge15']
    FACE_MAPPING['tri3'] = [('line2', (0, 1)),
                            ('line2', (1, 2)),
                            ('line2', (2, 0))]
    FACE_MAPPING['tri6'] = [('line3', (0, 1, 3)),
                            ('line3', (1, 2, 4)),
                            ('line3', (2, 0, 5))]
    FACE_MAPPING['quad4'] = [('line2', (0, 1)),
                             ('line2', (1, 2)),
                             ('line2', (2, 3)),
                             ('line2', (3, 0))]
    FACE_MAPPING['quad6'] = [('line3', (0, 1, 4)),
                             ('line2', (1, 2)),
                             ('line3', (2, 5, 3)),
                             ('line2', (3, 0))]
    FACE_MAPPING['quad8'] = [('line3', (0, 1, 4)),
                             ('line3', (1, 2, 5)),
                             ('line3', (2, 3, 6)),
                             ('line3', (3, 0, 7))]
    FACE_MAPPING['line2'] = [('point', tuple([0])),
                             ('point', tuple([1]))]
    FACE_MAPPING['line3'] = FACE_MAPPING['line2']
    FACE_MAPPING['point'] = []

    # define all standard element types
    STANDARD_ELEMENT_TYPES = set(FACE_MAPPING.keys())

    # define a connectivity permutation which inverts the element
    # for 2d elements, this flips the element normal
    # for 1d elements, this changes the direction
    INVERTED_CONNECTIVITY = dict()
    INVERTED_CONNECTIVITY['hex8'] = (4, 5, 6, 7, 0, 1, 2, 3)
    INVERTED_CONNECTIVITY['hex20'] = (4, 5, 6, 7, 0, 1, 2, 3, 16, 17, 18, 19,
                                      8, 9, 10, 11, 12, 13, 14, 15)
    INVERTED_CONNECTIVITY['tet4'] = (0, 2, 1, 3)
    INVERTED_CONNECTIVITY['tet10'] = (0, 2, 1, 3, 6, 5, 4, 7, 9, 8)
    INVERTED_CONNECTIVITY['wedge6'] = (3, 4, 5, 0, 1, 2)
    INVERTED_CONNECTIVITY['wedge12'] = (3, 4, 5, 0, 1, 2, 9, 10, 11, 6, 7, 8)
    INVERTED_CONNECTIVITY['wedge15'] = (3, 4, 5, 0, 1, 2, 12, 13, 14, 6, 7, 8,
                                        9, 10, 11)
    INVERTED_CONNECTIVITY['wedge16'] = (3, 4, 5, 0, 1, 2, 12, 13, 14, 6, 7, 8,
                                        9, 10, 11, 15)
    INVERTED_CONNECTIVITY['tri3'] = (0, 2, 1)
    INVERTED_CONNECTIVITY['tri6'] = (0, 2, 1, 5, 4, 3)
    INVERTED_CONNECTIVITY['quad4'] = (0, 3, 2, 1)
    INVERTED_CONNECTIVITY['quad6'] = (3, 2, 1, 0, 5, 4)
    INVERTED_CONNECTIVITY['quad8'] = (0, 3, 2, 1, 7, 6, 5, 4)
    INVERTED_CONNECTIVITY['line2'] = (1, 0)
    INVERTED_CONNECTIVITY['line3'] = (1, 0, 2)
    INVERTED_CONNECTIVITY['point'] = tuple([0])

    # define a connectivity permutation which rotates the given 2d face
    ROTATED_CONNECTIVITY = dict()
    ROTATED_CONNECTIVITY['quad4'] = (3, 0, 1, 2)
    ROTATED_CONNECTIVITY['quad8'] = (3, 0, 1, 2, 7, 4, 5, 6)
    ROTATED_CONNECTIVITY['tri3'] = (2, 0, 1)
    ROTATED_CONNECTIVITY['tri6'] = (2, 0, 1, 5, 3, 4)

    # define the topological dimension of each element
    DIMENSION = dict()
    DIMENSION['point'] = 0
    DIMENSION['line2'] = 1
    DIMENSION['line3'] = 1
    DIMENSION['tri3'] = 2
    DIMENSION['tri6'] = 2
    DIMENSION['quad4'] = 2
    DIMENSION['quad6'] = 2
    DIMENSION['quad8'] = 2
    DIMENSION['hex8'] = 3
    DIMENSION['hex20'] = 3
    DIMENSION['tet4'] = 3
    DIMENSION['tet10'] = 3
    DIMENSION['wedge6'] = 3
    DIMENSION['wedge12'] = 3
    DIMENSION['wedge15'] = 3
    DIMENSION['wedge16'] = 3

    # define the number of nodes per element
    NODES_PER_ELEMENT = dict((key, len(value))
                             for key, value in INVERTED_CONNECTIVITY.items())

    # define how to triangulate faces of 3D elements
    TRIANGULATED_FACES = dict()
    TRIANGULATED_FACES['tri3'] = [(0, 1, 2)]
    TRIANGULATED_FACES['tri6'] = [(0, 3, 5),
                                  (1, 4, 3),
                                  (2, 5, 4),
                                  (3, 4, 5)]
    TRIANGULATED_FACES['quad4'] = [(0, 1, (0, 1, 2, 3)),
                                   (1, 2, (0, 1, 2, 3)),
                                   (2, 3, (0, 1, 2, 3)),
                                   (3, 0, (0, 1, 2, 3))]
    TRIANGULATED_FACES['quad6'] = [(0, 4, (4, 5)),
                                   (4, 1, (4, 5)),
                                   (1, 2, (4, 5)),
                                   (2, 5, (4, 5)),
                                   (5, 3, (4, 5)),
                                   (3, 0, (4, 5))]
    TRIANGULATED_FACES['quad8'] = [(4, 7, 0),
                                   (5, 4, 1),
                                   (6, 5, 2),
                                   (7, 6, 3),
                                   (7, 4, (4, 5, 6, 7)),
                                   (4, 5, (4, 5, 6, 7)),
                                   (5, 6, (4, 5, 6, 7)),
                                   (6, 7, (4, 5, 6, 7))]

    # define formulas for converting between element types
    ELEMENT_CONVERSIONS = dict()
    ELEMENT_CONVERSIONS['hex8'] = dict()
    ELEMENT_CONVERSIONS['hex8']['hex20'] = [[0, 1, 2, 3, 4, 5, 6, 7,
                                             (0, 1), (1, 2), (2, 3), (0, 3),
                                             (0, 4), (1, 5), (2, 6), (3, 7),
                                             (4, 5), (5, 6), (6, 7), (4, 7)]]
    ELEMENT_CONVERSIONS['hex20'] = dict()
    ELEMENT_CONVERSIONS['hex20']['hex8'] = [range(8)]
    ELEMENT_CONVERSIONS['tet4'] = dict()
    ELEMENT_CONVERSIONS['tet4']['tet10'] = [[0, 1, 2, 3,
                                             (0, 1), (1, 2), (0, 2),
                                             (0, 3), (1, 3), (2, 3)]]
    ELEMENT_CONVERSIONS['tet10'] = dict()
    ELEMENT_CONVERSIONS['tet10']['tet4'] = [range(4)]
    ELEMENT_CONVERSIONS['wedge6'] = dict()
    ELEMENT_CONVERSIONS['wedge6']['wedge15'] = [[0, 1, 2, 3, 4, 5,
                                                 (0, 1), (1, 2), (0, 2),
                                                 (0, 3), (1, 4), (2, 5),
                                                 (3, 4), (4, 5), (3, 5)]]
    ELEMENT_CONVERSIONS['wedge6']['wedge16'] = [[0, 1, 2, 3, 4, 5,
                                                 (0, 1), (1, 2), (0, 2),
                                                 (0, 3), (1, 4), (2, 5),
                                                 (3, 4), (4, 5), (3, 5),
                                                 (0, 1, 2, 3, 4, 5)]]
    ELEMENT_CONVERSIONS['wedge12'] = dict()
    ELEMENT_CONVERSIONS['wedge12']['wedge6'] = [range(6)]
    ELEMENT_CONVERSIONS['wedge12']['wedge15'] = [[0, 1, 2, 3, 4, 5, 6, 7, 8,
                                                  (0, 3), (1, 4), (2, 5),
                                                  9, 10, 11]]
    ELEMENT_CONVERSIONS['wedge15'] = dict()
    ELEMENT_CONVERSIONS['wedge15']['wedge6'] = [range(6)]
    ELEMENT_CONVERSIONS['wedge15']['wedge16'] = [[0, 1, 2, 3, 4, 5,
                                                  6, 7, 8, 9, 10, 11,
                                                  12, 13, 14,
                                                  (0, 1, 2, 3, 4, 5)]]
    ELEMENT_CONVERSIONS['wedge16'] = dict()
    ELEMENT_CONVERSIONS['wedge16']['wedge6'] = [range(6)]
    ELEMENT_CONVERSIONS['wedge16']['wedge15'] = [range(15)]

    # define the order of each element
    ELEMENT_ORDER = dict()
    ELEMENT_ORDER['hex8'] = 1
    ELEMENT_ORDER['hex20'] = 2
    ELEMENT_ORDER['tet4'] = 1
    ELEMENT_ORDER['tet10'] = 2
    ELEMENT_ORDER['wedge6'] = 1
    ELEMENT_ORDER['wedge12'] = 1
    ELEMENT_ORDER['wedge15'] = 2
    ELEMENT_ORDER['wedge16'] = 2
    ELEMENT_ORDER['tri3'] = 1
    ELEMENT_ORDER['tri6'] = 2
    ELEMENT_ORDER['quad4'] = 1
    ELEMENT_ORDER['quad6'] = 1
    ELEMENT_ORDER['quad8'] = 2
    ELEMENT_ORDER['line2'] = 1
    ELEMENT_ORDER['line3'] = 2
    ELEMENT_ORDER['point'] = 1

    # define components of multi-component fields
    MULTI_COMPONENT_FIELD_SUBSCRIPTS = dict()
    MULTI_COMPONENT_FIELD_SUBSCRIPTS['vector'] = ('x', 'y', 'z')
    MULTI_COMPONENT_FIELD_SUBSCRIPTS['symmetric_3x3_tensor'] = (
        'xx', 'yy', 'zz', 'xy', 'yz', 'zx')
    MULTI_COMPONENT_FIELD_SUBSCRIPTS['full_3x3_tensor'] = (
        'xx', 'yy', 'zz', 'xy', 'yz', 'zx', 'yx', 'zy', 'xz')
    ALL_MULTI_COMPONENT_FIELD_SUBSCRIPTS = set(
        itertools.chain(*MULTI_COMPONENT_FIELD_SUBSCRIPTS.values()))

    def __init__(self):
        """Initialize the model."""
        # (only) the first time this module is used, show an info banner
        global SHOW_BANNER
        if SHOW_BANNER:
            print('\n\nYou are using Exomerge v%s -- A lightweight Python '
                  'interface for manipulating\nExodusII files. (Python-2 version)'
                  % (VERSION))
            # print out the exodus banner after formatting it to fit within
            # 79 characters
            exodus_text = exodus.EXODUS_PY_COPYRIGHT
            exodus_text = exodus_text.strip().replace('\n', ' ')
            exodus_text = exodus_text.replace('. ', '.  ')
            exodus_text = exodus_text.replace('.   ', '.  ')
            exodus_text = textwrap.fill(exodus_text, width=79)
            print('\n%s\n' % exodus_text)
            SHOW_BANNER = False
        # list of [x, y, z] coordinates for each node
        self.nodes = []
        # self.node_fields[name] = values
        # with values[timestep_index][node_index]
        self.node_fields = {}
        # self.global_variables[name] = values
        # with values[timestep_index]
        self.global_variables = {}
        # self.element_blocks[id] = [name, info, connectivity, fields]
        # with name = string name of element block or '' if unnamed
        # with info = ['hex8', elements, 8, 0]
        # with connectivity a shallow list of connectivity
        # with fields['eqps'] = values
        # with values[timestep][local_element_index]
        self.element_blocks = {}
        # self.side_sets[side_set_id] = [name, members, fields]
        # with name = string name of side set or '' if unnamed
        # with members = list of (element_block_id, element_index, side_index)
        # with fields[name] = values
        # with values[timestep_index][local_member_index]
        self.side_sets = {}
        # self.node_sets[node_set_id] = [name, members, fields]
        # with name = string name of node set or '' if unnamed
        # with members = list of node indices
        # with fields[name] = values
        # with values[timestep_index][local_member_index]
        self.node_sets = {}
        # list of timesteps
        self.timesteps = []
        # list of info records
        self.info_records = []
        # list of qa records
        self.qa_records = []
        # title of the database
        self.title = ''

    def __getattr__(self, name):
        """
        Try to find the given attribute.

        This is a special Python method which gets called if the attribute
        cannot be found.  We use to to suggest similar names in the case the
        user made a typo.

        We also use this to alias the plural forms of the delete functions.

        """
        # for special methods, use default behavior
        if name[:2] == '__':
            raise AttributeError
        # get non-special function names
        names = [x
                 for x in dir(self.__class__)
                 if x[:2] != '__']
        # if the name appears to be singular, search for the plural version
        if not name.endswith('s'):
            trial = name + 's'
            if trial in names:
                return getattr(self, trial)
        # if the name appears to be plural, search for the singular version
        if name.endswith('s'):
            trial = name[:-1]
            if not trial.endswith('s') and trial in names:
                return getattr(self, trial)
        # we can't find the name, so find names close to it to offer as
        # suggestions
        # filter names by closeness
        names = dict((x, difflib.SequenceMatcher(None, x, name).ratio())
                     for x in names)
        sorted_names = sorted(names.iteritems(),
                              key=operator.itemgetter(1),
                              reverse=True)
        # get the closest 5 matching function names
        closest_names = [x[0] for x in sorted_names[0:5]]
        # in an interactive shell, just issue a warning, else
        # issue an error an exit
        if self._is_interactive():
            call = self._warning
        else:
            call = self._error
        call('Function not found.',
             'The function "%s" does not exist.  Perhaps you meant to '
             'call one of the following functions:\n\n%s'
             % (name, '\n'.join(closest_names)))

    @staticmethod
    def _is_interactive():
        """Return True if we're in an interactive shell."""
        import __main__ as main
        return not hasattr(main, '__file__')

    @staticmethod
    def _assert(condition):
        """
        Raise a 'ValueError' if the given condition is not true.

        Unlike the assert builtin function, this checks the condition
        regardless of the state of '__debug__'.

        """
        if not condition:
            raise ValueError

    @staticmethod
    def _cubic_interpolation(x, x0, x1, x2, x3):
        """
        Return proportions using the cubic interpolation formula.

        Find the proportions of 'y0', 'y1', 'y2', 'y3' to take to find 'y(x)'
        for 'x1 <= x <= x2'.  This requires 'x0 < x1 < x2 < x3'.

        Example:
        >>> model._cubic_interpolation(1.71, 0.0, 1.0, 2.0, 3.0)
        [-0.0298555, 0.2766165, 0.8263335, -0.0730945]

        """
        assert x0 < x1 and x1 < x2 and x2 < x3
        assert x1 <= x and x <= x2
        values = []
        # proportion of y0
        values.append(((x - x1) * (x - x2) ** 2) /
                      ((x0 - x2) * (x1 - x2) ** 2))
        # proportion of y1
        values.append(-(((x - x2) * (-(x * x1 * (x1 + 3 * x2)) - x1 *
                                     (x1 ** 2 - 4 * x1 * x2 + x2 ** 2) +
                                     x ** 2 * (x1 + x2 - 2 * x3) + x2 *
                                     (-3 * x1 + x2) * x3 + x * (3 * x1 + x2) *
                                     x3)) / ((x1 - x2) ** 3 * (x1 - x3))))
        # proportion of y2
        values.append(((x - x1) * (x0 * x1 * (x1 - 3 * x2) + x ** 2 *
                                   (-2 * x0 + x1 + x2) - x2 *
                                   (x1 ** 2 - 4 * x1 * x2 + x2 ** 2) + x *
                                   (-(x2 * (3 * x1 + x2)) + x0 *
                                    (x1 + 3 * x2)))) /
                      ((x0 - x2) * (-x1 + x2) ** 3))
        # proportion of y3
        values.append(((x - x1) ** 2 * (x - x2)) /
                      ((x1 - x2) ** 2 * (x3 - x1)))
        return values

    def _new_element_block_id(self):
        """Return an element block id which is not used in the model."""
        id_ = 1
        while self.element_block_exists(id_):
            id_ += 1
        return id_

    def _new_side_set_id(self):
        """Return a side set id which is not used in the model."""
        id_ = 1
        while self.side_set_exists(id_):
            id_ += 1
        return id_

    def _new_node_set_id(self):
        """Return a node set id which is not used in the model."""
        id_ = 1
        while self.node_set_exists(id_):
            id_ += 1
        return id_

    def _new_element_field_name(self, quantity=1):
        """
        Return an element field name which is not used in the model.

        If 'quantity=1', this will return a single string.  Else, it will
        return a list of strings.
        """
        id_ = 1
        names = []
        for _ in range(quantity):
            name = 'temp%d' % id_
            while name in self.get_element_field_names():
                id_ += 1
                name = 'temp%d' % id_
            names.append(name)
            id_ += 1
        if quantity == 1:
            return names[0]
        else:
            return names

    def _new_node_field_name(self):
        """Return a node field name which is not used in the model."""
        id_ = 1
        name = 'temp%d' % id_
        while name in self.get_node_field_names():
            id_ += 1
            name = 'temp%d' % id_
        return name

    def _new_side_set_field_name(self):
        """Return a side set field name which is not used in the model."""
        id_ = 1
        name = 'temp%d' % id_
        while name in self.get_side_set_field_names():
            id_ += 1
            name = 'temp%d' % id_
        return name

    def _new_node_set_field_name(self):
        """Return a node set field name which is not used in the model."""
        id_ = 1
        name = 'temp%d' % id_
        while name in self.get_node_set_field_names():
            id_ += 1
            name = 'temp%d' % id_
        return name

    def _delete_elements(self, element_block_id, element_indices):
        """
        Delete the specified elements from the given element block.

        This will also delete all references to those elements in element
        fields, side sets, and side set fields.  This will also delete nodes
        which were used by those elements which are no longer used.

        Element indices are local to that element block.

        """
        # validate input
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        self._input_check(element_indices, [list, 0, int])
        # get element block information
        element_count = self.get_element_count(element_block_id)
        # if not deleting any indices, we're done
        if not element_indices:
            return
        # error if values are outside valid range
        max_element_index = max(element_indices)
        if max_element_index >= element_count:
            self._error('Invalid element indices',
                        'The element index list contains invalid entries.  '
                        'There are only %d elements in this element block and '
                        'the list contains a reference to index %d.'
                        % (max_element_index))
        # select valid elements
        indices_to_delete = set()
        invalid_indices = []
        duplicate_indices = []
        for x in element_indices:
            if x >= 0 and x < element_count:
                if x not in indices_to_delete:
                    indices_to_delete.add(x)
                else:
                    duplicate_indices.append(x)
            else:
                invalid_indices.append(x)
        # warn about duplicate indices
        if duplicate_indices:
            self._warning('Duplicate element indices',
                          'The element index list contains duplicate '
                          'indices.  Elements can only be deleted once so '
                          'duplicates will be ignored.  There were %d '
                          'duplicates.'
                          % (len(duplicate_indices)))
        del duplicate_indices
        # error if invalid indices encountered
        if invalid_indices:
            self._error('Invalid element indices',
                        'The element index list contains invalid entries.  '
                        'There were a total of %d invalid indices.'
                        % len(invalid_indices))
        del invalid_indices
        # create list of nodes used by elements to delete
        nodes_per_element = self.get_nodes_per_element(element_block_id)
        connectivity = self.get_connectivity(element_block_id)
        used_nodes = []
        for element_index in indices_to_delete:
            used_nodes.extend(
                connectivity[element_index * nodes_per_element:
                             (element_index + 1) * nodes_per_element])
        used_nodes = set(used_nodes)
        # remove references to these elements within each side set
        for id_ in self.get_side_set_ids():
            # look for members to delete
            members = self.get_side_set_members(id_)
            members_to_delete = [x
                                 for x in members
                                 if (x[0] == element_block_id and
                                     x[1] in indices_to_delete)]
            # delete those members (if any exist)
            if members_to_delete:
                self._delete_side_set_members(id_, members_to_delete)
        # create translation list
        element_count = self.get_element_count(element_block_id)
        remaining_indices = sorted(set(range(element_count)) -
                                   set(indices_to_delete))
        # create element map
        # old element i is new element new_index[i]
        new_index = [False] * element_count
        for i, x in enumerate(remaining_indices):
            new_index[x] = i
        # delete elements from element fields
        del indices_to_delete
        fields = self._get_element_block_fields(element_block_id)
        for this_field in fields.values():
            this_field[:] = [[values[x] for x in remaining_indices]
                             for values in this_field]
        # delete elements from the block
        new_connectivity = []
        for element_index in remaining_indices:
            new_connectivity.extend(
                connectivity[element_index * nodes_per_element:
                             (element_index + 1) * nodes_per_element])
        # set the new connectivity
        self.get_connectivity(element_block_id)[:] = new_connectivity
        # set the number of elements
        self.element_blocks[element_block_id][1][1] = len(remaining_indices)
        # change side set element numbering
        for side_set_id in self.get_side_set_ids():
            members = self.get_side_set_members(side_set_id)
            new_members = []
            for member in members:
                if member[0] == element_block_id:
                    new_members.append((member[0],
                                        new_index[member[1]],
                                        member[2]))
                else:
                    new_members.append(member)
            members[:] = new_members
        # make sure we didn't mess anything up
        # find nodes which are not used by any element
        abandoned_nodes = set(self._get_unreferenced_nodes())
        # find nodes to delete
        nodes_to_delete = sorted(abandoned_nodes & used_nodes)
        # delete those nodes
        self._delete_nodes(nodes_to_delete)

    def threshold_element_blocks(self,
                                 expression,
                                 element_block_ids='all',
                                 timestep='last',
                                 new_element_block_id=None):
        """
        Delete elements which do not satisfy a condition.

        The expression can contain element field values.

        If 'new_element_block_id' is given, all elements which don't pass
        the threshold will be moved into that element block.  Else, they are
        deleted.

        Example:
        >>> model.threshold_element_block('all', 'eqps >= 0.01')

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        [timestep] = self._format_id_list(
            [timestep],
            self.get_timesteps(),
            'timestep',
            single=True)
        if new_element_block_id is not None:
            if self.element_block_exists(new_element_block_id):
                self._exists_error(new_element_block_id, 'element block')
        # create a temporary field name on these elements
        name = self._new_element_field_name()
        # evaluate the expression
        self.calculate_element_field('%s = %s' % (name, expression),
                                     element_block_ids=element_block_ids,
                                     timesteps=[timestep])
        # go through each element block and delete elements
        for id_ in element_block_ids:
            values = self.get_element_field_values(name,
                                                   element_block_id=id_,
                                                   timestep=timestep)
            # if we're saving elements, duplicate the block, and delete the
            # ones we're keeping
            if new_element_block_id is not None:
                indices_to_keep = [i
                                   for i, x in enumerate(values)
                                   if x]
                temporary_element_block_id = self._new_element_block_id()
                self.duplicate_element_block(id_,
                                             temporary_element_block_id,
                                             duplicate_nodes=False)
                self._delete_elements(temporary_element_block_id,
                                      indices_to_keep)
                if temporary_element_block_id != new_element_block_id:
                    if self.element_block_exists(new_element_block_id):
                        self.combine_element_blocks(
                            [temporary_element_block_id,
                             new_element_block_id],
                            new_element_block_id)
                    else:
                        self.rename_element_block(temporary_element_block_id,
                                                  new_element_block_id)
            # delete elements which don't meet the threshold
            indices_to_delete = [i
                                 for i, x in enumerate(values)
                                 if not x]
            self._delete_elements(id_, indices_to_delete)
        # delete the temporary element field
        self.delete_element_field(name, element_block_ids=element_block_ids)
        if new_element_block_id is not None:
            self.delete_element_field(name,
                                      element_block_ids=new_element_block_id)

    def _get_list_ids(self, thing):
        """
        Return a list of list ids found within a variable.

        This handles lists, dictionaries, and nested combinations of these.
        Other objects are ignored.

        This method is used to check aliasing of lists.

        """
        ids = []
        if isinstance(thing, list):
            ids.append(id(thing))
            if isinstance(thing, list) or isinstance(thing[0], dict):
                for other_thing in thing:
                    ids.extend(self._get_list_ids(other_thing))
            return ids
        elif isinstance(thing, dict):
            for other_thing in thing.values():
                ids.extend(self._get_list_ids(other_thing))
        return ids

    def _verify(self, allow_aliased_lists=True):
        """
        Verify model information is valid and arrays are appropriately sized.

        If there is a problem somewhere, this method will throw an error
        message and exit.

        """
        timestep_count = len(self.timesteps)
        try:
            # verify self.node_fields
            node_count = len(self.nodes)
            for all_values in self.node_fields.values():
                self._assert(len(all_values) == timestep_count)
                for values in all_values:
                    self._assert(len(values) == node_count)
            # verify self.element_blocks
            for _, info, connectivity, fields in self.element_blocks.values():
                self._assert(len(info) == 4)
                self._assert(len(connectivity) == info[1] * info[2])
                if connectivity:
                    self._assert(min(connectivity) >= 0)
                    self._assert(max(connectivity) < node_count)
                for _, all_values in fields.items():
                    self._assert(len(all_values) == timestep_count)
                    for values in all_values:
                        self._assert(len(values) == info[1])
            # verify self.node_sets
            for _, members, fields in self.node_sets.values():
                member_count = len(members)
                if members:
                    self._assert(min(members) >= 0)
                    self._assert(max(members) < node_count)
                for _, all_values in fields.items():
                    self._assert(len(all_values) == timestep_count)
                    for values in all_values:
                        self._assert(len(values) == member_count)
            # verify self.side_sets
            element_count = dict(
                (id_, info[1])
                for id_, (_, info, _, _) in self.element_blocks.items())
            for _, members, fields in self.side_sets.values():
                member_count = len(members)
                if members:
                    self._assert(min(x[1] for x in members) >= 0)
                for id_, element_index, _ in members:
                    self._assert(element_index < element_count[id_])
                for _, all_values in fields.items():
                    self._assert(len(all_values) == timestep_count)
                    for values in all_values:
                        self._assert(len(values) == member_count)
            # verify self.global_variables
            for _, values in self.global_variables.items():
                self._assert(len(values) == timestep_count)
            # optionally verify that no lists are aliases of one another
            if not allow_aliased_lists:
                # hold list ids for later use
                list_ids = []
                for name, thing in inspect.getmembers(self):
                    if not name[0].isalpha() or not name[0].islower():
                        continue
                    if isinstance(thing, dict) or isinstance(thing, list):
                        list_ids.extend(self._get_list_ids(thing))
                # ensure none are copies
                unique_list_ids = list(set(list_ids))
                self._assert(len(unique_list_ids) == len(list_ids))
        except ValueError:
            # get the code line the assertion failed on
            _, _, trace = sys.exc_info()
            error_line = traceback.extract_tb(trace)[0][3]
            if error_line.startswith('self._assert('):
                error_line = error_line[13:-1]
            error_line_number = traceback.extract_tb(trace)[0][1]
            # write out the error/warning message
            self._error(
                'Corrupted database information.',
                'The finite element model information failed a validity '
                'check.  The following assertion failed:\n'
                '\n'
                '[%d] %s\n'
                '\n'
                'This resulted from one of the two scenarios:\n'
                '* You changed information in the database '
                'inconsistently.\n'
                '* There is a bug within the code.\n'
                '\n'
                'If you believe it is the latter, please contact support '
                'with as much information as necessary.\n'
                '* Support contact: %s'
                % (error_line_number, error_line, CONTACT))
            exit(1)

    @staticmethod
    def _print_message(messages, width=79):
        """
        Print out a nicely formatted message which fits in a specified width.

        Example:
        >>> exomerge.ExodusModel._print_message(
        ...     [("ERROR: ", "I'm sorry, Dave"),
        ...      ("Message: ", "I'm afraid I can't let you do that.")],
        ...     width=30)

        ***********************************
        *                                 *
        *  ERROR: I'm sorry, Dave.        *
        *                                 *
        *  Message: I'm afraid I can't    *
        *           let you do that.      *
        *                                 *
        ***********************************

        """
        # format messages
        if not isinstance(messages, list):
            messages = [messages]
        for i in xrange(len(messages)):
            if not isinstance(messages[i], tuple):
                messages[i] = ('', messages[i])
        # store the message width
        max_header_width = max(len(x[0]) for x in messages)
        width = max(width, max_header_width + 20 + 6)
        # print the top frame after a blank line
        print('')
        print('*' * width)
        print('*' + (' ' * (width - 2)) + '*')
        # process each message
        # *  Header: Text             *
        for message in messages:
            header = message[0]
            # format text into broken lines
            text = ''
            for line in message[1].split('\n'):
                text += textwrap.fill(
                    line,
                    width=width - 6 - len(header))
                text += '\n'
            # remove trailing newlines
            while text and text[-1] == '\n':
                text = text[:-1]
            # break into a list
            text = text.split('\n')
            # process each line
            for i in xrange(len(text)):
                if i == 0:
                    line = '%s%s' % (header, text[i])
                else:
                    line = '%s%s' % (' ' * len(header), text[i])
                if len(line) < width - 6:
                    line += ' ' * (width - 6 - len(line))
                print('*  %s  *' % (line))
            print('*' + (' ' * (width - 2)) + '*')
        # print the bottom frame
        print('*' * width)

    @staticmethod
    def _get_trace(omitted_frames=1):
        """Return a compressed stack trace."""
        last_file = ''
        message = ''
        for frame in inspect.stack()[omitted_frames:]:
            this_file = frame[1].split('/')[-1]
            if this_file != last_file:
                message += this_file + ':'
                last_file = this_file
            else:
                message += ' ' * len(last_file) + ' '
            message += frame[3] + ':' + str(frame[2])
            message += '\n'
        return message

    def _bug(self, short='Possible bug', detailed='', exit_code=1):
        """
        Print out an error message possibly caused by a bug and exit.

        As opposed to the _error function, this also outputs support contact
        information.

        """
        messages = [('  ERROR: ', short)]
        if detailed:
            messages.append(('Message: ', detailed))
        messages.append(('Support: ',
                         'This error may have been caused by a bug in the '
                         'code.  If you believe this is the case, please '
                         'contact support and provide enough information '
                         'to duplicate the issue.\n\nSupport contact: %s'
                         % CONTACT))
        messages.append(('  Trace: ', self._get_trace()))
        self._print_message(messages)
        exit(exit_code)

    def _error(self, short='Unspecified error', detailed='', exit_code=1):
        """Print out an error message and exit."""
        messages = [('  ERROR: ', short)]
        if detailed:
            messages.append(('Message: ', detailed))
        messages.append(('  Trace: ', self._get_trace()))
        self._print_message(messages)
        exit(exit_code)

    def _warning(self, short='Unspecified warning', detailed=''):
        """Print out a warning message, but don't exit."""
        messages = [('WARNING: ', short)]
        if detailed:
            messages.append(('Message: ', detailed))
        if EXIT_ON_WARNING:
            messages.append(('  Trace: ', self._get_trace()))
        self._print_message(messages)
        if EXIT_ON_WARNING:
            exit(1)

    @staticmethod
    def _get_qa_record():
        """
        Return a new qa record for this program.

        QA records are stored within the ExodusII file.

        """
        now = datetime.datetime.today()
        return ('exomerge',
                VERSION,
                now.strftime('%Y/%m/%d'),
                now.strftime('%H:%M:%S'))

    def _ensure_acceptable_id(self, this_id, acceptable_ids, entity=''):
        """Ensure ID is acceptable, or error out."""
        if this_id not in acceptable_ids:
            if entity == "":
                self._error("Reference to undefined entity.")
            else:
                id_list = string.join(
                    [str(x) + ',' for x in acceptable_ids])[:-1]
                self._error(
                    "Reference to undefined %s." % (entity),
                    "A reference to %s \"%s\" was encountered but is "
                    "not defined.  There are %d defined %ss: %s" %
                    (entity, str(this_id), len(acceptable_ids),
                     entity, id_list))

    @staticmethod
    def _get_matching_strings(string_pattern, list_of_strings):
        """
        Return list of string which match the pattern.

        The pattern can contain a single wildcard character (*) which
        represents any number of characters.  No other special regex
        characters are allowed.

        Example:
        >>> model._get_matching_strings('f*x',
        ...                             ['fox', 'faux', 'far', 'fx', 'foxes'])
        ['fox', 'faux', 'fx']

        """
        # convert pattern to a regex expression
        regex_pattern = '^' + string_pattern + '$'
        regex_pattern = regex_pattern.replace('.', r'\.').replace('*', '.*')
        matches = [x
                   for x in list_of_strings
                   if re.match(regex_pattern, x)]
        return matches

    @staticmethod
    def _plural(name):
        """
        Return the likely plural form of the given word.

        Example:
        >>> model._plural('rose')
        roses
        >>> model._plural('bunny')
        bunnies

        """
        if name[-1] == 'y':
            return name[:-1] + 'ies'
        else:
            return name + 's'

    def _format_element_block_id_list(self, id_list, *args, **kwargs):
        """Return a validated list of element block IDs."""
        object = self.element_blocks
        entity = 'element block'
        # get element block translations
        tuples = [(value[0], key)
                  for key, value in object.items()
                  if value[0]]
        # ensure element block names are unique
        if len(set(x[0] for x in tuples)) != len(tuples):
            self._warning('Duplicate %s names' % entity,
                          'At least two %s have the same name.  '
                          'This may cause problems.' % self._plural(entity))
        return self._format_id_list(id_list,
                                    sorted(object.keys()),
                                    entity=entity,
                                    *args,
                                    translator=dict(tuples),
                                    **kwargs)

    def _format_side_set_id_list(self, id_list, *args, **kwargs):
        """Return a validated list of side set IDs."""
        object = self.side_sets
        entity = 'side set'
        # get element block translations
        tuples = [(value[0], key)
                  for key, value in object.items()
                  if value[0]]
        # ensure element block names are unique
        if len(set(x[0] for x in tuples)) != len(tuples):
            self._warning('Duplicate %s names' % entity,
                          'At least two %s have the same name.  '
                          'This may cause problems.' % self._plural(entity))
        return self._format_id_list(id_list,
                                    sorted(object.keys()),
                                    entity=entity,
                                    *args,
                                    translator=dict(tuples),
                                    **kwargs)

    def _format_node_set_id_list(self, id_list, *args, **kwargs):
        """Return a validated list of node set IDs."""
        object = self.node_sets
        entity = 'node set'
        # get element block translations
        tuples = [(value[0], key)
                  for key, value in object.items()
                  if value[0]]
        # ensure element block names are unique
        if len(set(x[0] for x in tuples)) != len(tuples):
            self._warning('Duplicate %s names' % entity,
                          'At least two %s have the same name.  '
                          'This may cause problems.' % self._plural(entity))
        return self._format_id_list(id_list,
                                    sorted(object.keys()),
                                    entity=entity,
                                    *args,
                                    translator=dict(tuples),
                                    **kwargs)

    def _format_id_list(self,
                        id_list,
                        acceptable_ids,
                        entity='',
                        single=False,
                        empty_list_okay=True,
                        translator={}):
        """
        Return a validated list of IDs.

        If 'single' is 'True', this will error out if zero or more than one id
        is returned.

        If 'empty_list_okay=False', this will error out if the list is
        empty.

        If 'translator' is given, it must be a dictionary which converts
        strings to ids.

        """
        new_id_list = []
        if not isinstance(id_list, list):
            id_list = [id_list]
        for id_ in id_list:
            if id_ == 'none':
                pass
            elif id_ == 'all':
                new_id_list.extend(acceptable_ids)
            elif id_ == 'first':
                if not acceptable_ids:
                    self._error(
                        'Empty entity list.',
                        'The "first" identifier was used to refer to a %s '
                        'but no %ss exist.' % (entity, entity))
                new_id_list.append(acceptable_ids[0])
            elif id_ == 'last':
                if not acceptable_ids:
                    self._error(
                        'Empty entity list.',
                        'The "last" identifier was used to refer to a %s '
                        'but no %ss exist.' % (entity, entity))
                new_id_list.append(acceptable_ids[-1])
            elif id_ == 'auto':
                if len(acceptable_ids) != 1:
                    self._error(
                        'Ambiguous reference.',
                        'The "auto" identifier was used to refer to a %s '
                        'but exactly one %s does not exist.  There are %d '
                        'defined %ss.' % (entity, entity, len(acceptable_ids),
                                          entity))
                new_id_list.append(acceptable_ids[0])
            elif id_ == 'last_if_any':
                if acceptable_ids:
                    new_id_list.append(acceptable_ids[-1])
            elif id_ == 'first_if_any':
                if acceptable_ids:
                    new_id_list.append(acceptable_ids[-1])
            else:
                new_id_list.append(id_)
        id_list = new_id_list
        # now apply the translation list
        if translator:
            id_list = [translator[x] if x in translator else x
                       for x in id_list]
        # for ids with a wildcard, find all that match
        if id_list and isinstance(id_list[0], basestring):
            for i in reversed(xrange(len(id_list))):
                this_id = id_list[i]
                if '*' in this_id:
                    del id_list[i]
                    for new_id in reversed(self._get_matching_strings(
                            this_id, acceptable_ids)):
                        id_list.insert(i, new_id)
        # do not allow duplicates
        # and ensure all given ids are valid
        found = set()
        unique_ids = list()
        for id_ in id_list:
            if id_ not in found:
                found.add(id_)
                if id_ in acceptable_ids:
                    unique_ids.append(id_)
                else:
                    id_list = ', '.join(
                        [str(x) for x in acceptable_ids])
                    self._warning(
                        'Reference to undefined %s.' % (entity),
                        'A reference to %s \"%s\" was encountered but is '
                        'not defined.  This %s will not be included in the '
                        'operation.\n\nThere are %d defined %s: %s' %
                        (entity, id_,
                         entity,
                         len(acceptable_ids),
                         self._plural(entity),
                         id_list))
        # return the unique list
        id_list = unique_ids
        # if only a single item was requested, ensure it exists
        if single and len(id_list) != 1:
            if len(id_list) == 0:
                self._error('No %s specified.' % entity,
                            'A single %s is required but none '
                            'were specified.' % entity)
            self._error(
                'Multiple %ss specified.' % entity,
                'A single %s is required but %d %ss were specified: %s'
                % (entity, len(id_list), entity,
                   ', '.join([str(x) for x in id_list])))
        # if list is empty and that's a problem, error out
        if not empty_list_okay and len(id_list) == 0:
            self._error('No %ss specified.' % entity,
                        'We require at least one %s to be specified for this '
                        'opration but none were specified.\n\nThere are %d '
                        'defined %s: %s' %
                        (entity,
                         len(acceptable_ids),
                         self._plural(entity),
                         ', '.join([str(x) for x in acceptable_ids])))
        return id_list

    def _get_standard_element_type(self, element_type, warning=True):
        """
        Return a standardized string for the given element type.

        Within an ExodusII file, the element type is stored as a string.
        Because of this, there is no standardized way to store this field
        and different codes can and do use different strings for the same
        element type.  This translation function interprets this field and
        returns a standardized name if possible.

        If an element type is not recognized and 'warning=True', a warning
        message will be issued.

        Example:
        >>> model._get_standard_element_type('HEX')
        'hex8'

        """
        element_type = element_type.lower()
        translation_list = {'hex': 'hex8',
                            'tet': 'tet4',
                            'tetra': 'tet4',
                            'tetra10': 'tet10',
                            'tri': 'tri3',
                            'triangle': 'tri3',
                            'quad': 'quad4',
                            'shell4': 'quad4',
                            'wedge': 'wedge6'}
        if element_type in translation_list:
            element_type = translation_list[element_type]
        if warning and element_type not in self.STANDARD_ELEMENT_TYPES:
            self._warning('Nonstandard element type',
                          'The element type "%s" is not a standard '
                          'element type.  This may cause issues with '
                          'handling element blocks with elements of '
                          'this type.' %
                          (element_type))
        return element_type

    def _get_default_field_value(self, name):
        """Return the default value for a field with the given name."""
        if self._is_displacement_field_prefix(name[:-2]):
            return 0.0
        return float('nan')

    @staticmethod
    def _is_displacement_field_prefix(name):
        """Return 'True' if name is a displacement field prefix."""
        default_name = 'displacement'
        if len(name) >= 3 and name.lower() == default_name[:len(name)]:
            return True
        return False

    def _get_displacement_field_prefix(self):
        """
        Return the prefix of the displacement field.

        If no displacement field exists, this will return the default prefix.

        """
        prefix = 'disp'
        for node_field_name in self.get_node_field_names():
            # if it doesn't end in "_x", it's not a prefix
            if len(node_field_name) < 3 or node_field_name[-2] != '_':
                continue
            # check against acceptable names
            this_prefix = node_field_name[:-2]
            if self._is_displacement_field_prefix(this_prefix):
                prefix = this_prefix
        return prefix

    def _get_element_faces(self, element_block_ids):
        """
        Return a list of all element faces.

        Note that internal element faces will be duplicated in the list.

        A list is returned with members of the form:
        * '(element_block_id, element_index, face_index)'

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids,
            empty_list_okay=False)
        members = []
        for id_ in element_block_ids:
            mapping = self._get_face_mapping_from_id(id_)
            element_count = self.get_element_count(id_)
            for i in xrange(element_count):
                members.extend([(id_, i, f)
                                for f in range(len(mapping))])
        return members

    def _get_block_info(self, element_block_id):
        """Return the info of an element block."""
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        return self.element_blocks[element_block_id][1]

    def _get_external_element_faces(self, element_block_ids='all'):
        """
        Return a list of external element faces.

        External faces are element faces which are shared by no other elements.

        A list is returned with members of the form:
        * '(element_block_id, element_index, face_index)'

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        external_faces = dict()
        for id_ in element_block_ids:
            info = self._get_block_info(id_)
            connectivity = self.get_connectivity(id_)
            face_mapping = self._get_face_mapping_from_id(id_)
            for element_index in range(info[1]):
                for face_index, (_, face) in enumerate(face_mapping):
                    sorted_nodes = tuple(sorted([
                        connectivity[element_index * info[2] + x]
                        for x in face]))
                    if sorted_nodes in external_faces:
                        external_faces[sorted_nodes] = None
                    else:
                        this_face = (id_, element_index, face_index)
                        external_faces[sorted_nodes] = this_face
        return [value
                for value in external_faces.values()
                if value is not None]

    @staticmethod
    def _mean(values):
        """Return the mean of a list of values."""
        return sum(values) / float(len(values))

    @staticmethod
    def _order_element_faces_by_block(members):
        """
        Sort element faces by element block id.

        This takes in a list of members of the form:
        * '(element_block_id, element_index, face_index)'

        This outputs a dictionary with the form
        * 'output[element_block_id] = list of (element_index, face_index)'

        """
        members_by_block = dict()
        for element_block_id, element_index, face_index in members:
            if element_block_id not in members_by_block:
                members_by_block[element_block_id] = []
            members_by_block[element_block_id].append(
                (element_index, face_index))
        return members_by_block

    def _get_dimension(self, element_type):
        """Return the dimension of the given element type."""
        element_type = self._get_standard_element_type(element_type)
        return self.DIMENSION[element_type]

    def _get_triangulated_faces(self, element_type):
        """Return the rule for how to triangulate the given 2D face."""
        element_type = self._get_standard_element_type(element_type)
        if element_type not in self.TRIANGULATED_FACES:
            self._bug('Invalid element type.',
                      'We do not know how to subdivide an element of type '
                      '"%s" into triangles.' % element_type)
        return self.TRIANGULATED_FACES[element_type]

    def _convert_side_set_to_triangle_block(self,
                                            side_set_ids,
                                            new_element_block_id='auto'):
        """
        Create a new 'tri3' element block composed of the given side set.

        The side set must contain faces of three dimensional elements.

        """
        side_set_ids = self._format_side_set_id_list(side_set_ids)
        if new_element_block_id == 'auto':
            new_element_block_id = self._new_element_block_id()
        new_nodes = []
        # find members of all the given side sets by block
        members_by_block = self._order_element_faces_by_block(
            itertools.chain(*[self.get_side_set_members(x)
                              for x in side_set_ids]))
        # ensure all blocks are 3D
        for id_ in members_by_block.keys():
            element_type = self._get_element_type(id_)
            dimension = self._get_dimension(element_type)
            if dimension != 3:
                self._error('Incompatible element dimension.',
                            'We expect the side set to contain faces of 3D '
                            'elements.  However, the side set contains at '
                            'least one face from element block %d which '
                            'contains %dD elements of type "%s".' %
                            (id_, dimension, element_type))
        # Now we need to create the triangles from each member of the side set.
        faces = dict()
        for element_block_id, members in members_by_block.items():
            connectivity = self.get_connectivity(element_block_id)
            nodes_per_element = self.get_nodes_per_element(element_block_id)
            face_mapping = self._get_face_mapping_from_id(element_block_id)
            for element_index, face_index in members:
                face_element_map = face_mapping[face_index]
                local_node = connectivity[
                    element_index * nodes_per_element:
                    (element_index + 1) * nodes_per_element]
                local_node = tuple(local_node[x] for x in face_element_map[1])
                if not face_element_map[0] in faces:
                    faces[face_element_map[0]] = [local_node]
                else:
                    faces[face_element_map[0]].append(local_node)
        # now triangulate each face
        connectivity = []
        next_node_index = len(self.nodes)
        new_nodes = []
        for element_type, connectivity_list in faces.items():
            triangulated_faces = self._get_triangulated_faces(element_type)
            for local_face_indices in connectivity_list:
                for face_indices in triangulated_faces:
                    for index in face_indices:
                        if isinstance(index, int):
                            connectivity.append(local_face_indices[index])
                            continue
                        averaged_nodes = tuple(local_face_indices[x]
                                               for x in index)
                        if not new_nodes or new_nodes[-1] != averaged_nodes:
                            new_nodes.append(averaged_nodes)
                            next_node_index += 1
                        connectivity.append(next_node_index - 1)
        # create the new nodes
        self._create_averaged_nodes(new_nodes, [])
        # create the new element block
        self.create_element_block(new_element_block_id,
                                  ['tri3', len(connectivity) / 3, 3, 0],
                                  connectivity)

    @staticmethod
    def _sign(x):
        """
        Return the sign of a value.

        This returns 1 if the value positive, -1 if negative, and 0 otherwise.

        """
        if x > 0:
            return 1
        elif x < 0:
            return -1
        return 0

    def _partition_triangle_block_from_node_field(self,
                                                  element_block_id,
                                                  new_element_block_id,
                                                  node_field_name,
                                                  timestep,
                                                  interval_list):
        """
        Subdivide a 'tri3' element block based on the given field intervals.

        The target element block will have a element field named 'interval'
        corresponding to the given interval the element falls into.

        This is a helper function for generating nice looking WRL files.

        """
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        [node_field_name] = self._format_id_list(
            [node_field_name],
            self.get_node_field_names(),
            'node field',
            single=True)
        [timestep] = self._format_id_list(
            [timestep],
            self.get_timesteps(),
            'timestep',
            single=True)
        # verify block is actually a tri3 block
        if self._get_element_type(element_block_id) != 'tri3':
            self._error('Invalid element type',
                        'This function can only be used to divide element '
                        'blocks composed of "tri3" elements.')
        timestep_index = self.timesteps.index(timestep)
        connectivity = self.get_connectivity(element_block_id)
        element_count = self.get_element_count(element_block_id)
        new_connectivity = []
        element_interval_values = []
        triangles = [tuple(connectivity[x * 3: (x + 1) * 3])
                     for x in xrange(element_count)]
        for index, upper_bound in enumerate(interval_list):
            # hold new vertices we have to create
            # edge = (x, y)
            # new_vertex[edge] = mid_edge_vertex_index
            # we do this to keep the mesh watertight
            new_vertex = dict()
            new_nodes = []
            next_node_index = len(self.nodes)
            values = self.node_fields[node_field_name][timestep_index]
            for tri in triangles:
                for d in xrange(3):
                    d2 = (d + 1) % 3
                    edge = (tri[d], tri[d2])
                    if (self._sign(values[tri[d]] - upper_bound) *
                            self._sign(values[tri[d2]] - upper_bound) == -1):
                        if edge not in new_vertex:
                            opposite_edge = tuple(reversed(edge))
                            new_vertex[edge] = next_node_index
                            new_vertex[opposite_edge] = next_node_index
                            next_node_index += 1
                            phi = (upper_bound - values[tri[d]]) / (
                                values[tri[d2]] - values[tri[d]])
                            new_node = ((tri[d], 1.0 - phi), (tri[d2], phi))
                            new_nodes.append(new_node)
            self._create_averaged_nodes(new_nodes, [])
            # now form new triangles
            new_triangles = []
            this_partition = []
            values = self.node_fields[node_field_name][timestep_index]
            for tri in triangles:
                if (values[tri[0]] <= upper_bound and
                        values[tri[1]] <= upper_bound and
                        values[tri[2]] <= upper_bound):
                    this_partition.append(tri)
                elif (values[tri[0]] >= upper_bound and
                      values[tri[1]] >= upper_bound and
                      values[tri[2]] >= upper_bound):
                    new_triangles.append(tri)
                else:
                    # find vertex below the bound
                    d = 0
                    while (values[tri[d]] >= upper_bound or
                           values[tri[(d + 1) % 3]] < upper_bound):
                        d += 1
                    tri = tuple(tri[(d + x) % 3]
                                for x in xrange(3))
                    case = tuple(self._sign(values[tri[x]] - upper_bound)
                                 for x in xrange(3))
                    if case == (-1, 1, -1):
                        m1 = new_vertex[(tri[0], tri[1])]
                        m2 = new_vertex[(tri[1], tri[2])]
                        # below triangles
                        this_partition.append((tri[0], m1, tri[2]))
                        this_partition.append((m1, m2, tri[2]))
                        # above triangles
                        new_triangles.append((m1, tri[1], m2))
                    elif case == (-1, 1, 1):
                        m1 = new_vertex[(tri[0], tri[1])]
                        m2 = new_vertex[(tri[2], tri[0])]
                        # below triangles
                        this_partition.append((tri[0], m1, m2))
                        # above triangles
                        new_triangles.append((m1, tri[1], m2))
                        new_triangles.append((tri[1], tri[2], m2))
                    elif case == (-1, 0, 1):
                        m1 = new_vertex[(tri[2], tri[0])]
                        # below triangles
                        this_partition.append((tri[0], tri[1], m1))
                        # above triangles
                        new_triangles.append((tri[1], tri[2], m1))
                    elif case == (-1, 1, 0):
                        m1 = new_vertex[(tri[0], tri[1])]
                        # below triangles
                        this_partition.append((tri[0], m1, tri[2]))
                        # above triangles
                        new_triangles.append((m1, tri[1], tri[2]))
                    else:
                        self._bug('Unknown case')
            triangles = new_triangles
            new_connectivity.extend(itertools.chain(*this_partition))
            element_interval_values.extend(
                [float(index)] * len(this_partition))
        # add rest of triangle to last partition
        new_connectivity.extend(itertools.chain(*triangles))
        element_interval_values.extend(
            [float(len(interval_list))] * len(triangles))
        self.create_element_block(new_element_block_id,
                                  ['tri3', len(new_connectivity) / 3, 3, 0],
                                  new_connectivity)
        self.create_element_field('interval', new_element_block_id, 0.0)
        fields = self._get_element_block_fields(new_element_block_id)
        field = fields['interval']
        for index in xrange(len(self.timesteps)):
            field[index] = list(element_interval_values)

    def _get_coordinates_at_time(self, timestep):
        """
        Return the node coordinates list at the given timestep.

        This includes the displacement if it exists.

        """
        timestep = self._format_id_list(
            timestep,
            self.get_timesteps(),
            'timestep')
        if len(timestep) > 1:
            self._error('Ambiguous timestep.',
                        'More than one timestep was specified.  We expected '
                        'one or zero timesteps.')
        if len(timestep) == 0:
            return [tuple(x) for x in self.nodes]
        timestep_index = self.timesteps.index(timestep[0])
        displacement_values = [x[timestep_index]
                               for x in self._get_displacement_field_values()]
        return [(x + dx, y + dy, z + dz)
                for (x, y, z), dx, dy, dz
                in itertools.izip(self.nodes,
                                  *displacement_values)]

    def export_wrl_model(self,
                         filename,
                         node_field_name,
                         element_block_ids='all',
                         timestep='last',
                         field_range='auto',
                         intervals=9,
                         colorspace='rgb',
                         displacement_timestep='auto',
                         export_exodus_copy=True):
        """
        Export the exterior of the model to a colored WRL file.

        The WRL file format is used by 3D printing software.

        Example:
        >>> model.export_wrl_model('colored_eqps_model.wrl', 'eqps')

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        [node_field_name] = self._format_id_list(
            [node_field_name],
            self.get_node_field_names(),
            'node field',
            single=True)
        [timestep] = self._format_id_list(
            [timestep],
            self.get_timesteps(),
            'timestep',
            single=True)
        timestep_index = self.timesteps.index(timestep)
        if displacement_timestep == 'auto':
            if self.displacement_field_exists():
                displacement_timestep = timestep
            else:
                displacement_timestep = 'none'
        else:
            [displacement_timestep] = self._format_id_list(
                [displacement_timestep],
                self.get_timesteps(),
                'timestep',
                single=True)
        # all element blocks should either be 3D elements, or 'tri3' elements
        element_types = set(self._get_element_type(x)
                            for x in element_block_ids)
        element_dimensions = set(self._get_dimension(x)
                                 for x in element_types)
        # if we have all 3D elements, create a skin
        if element_dimensions == set([3]):
            triangle_element_block_id = self._new_element_block_id()
            external_side_set_id = self._new_side_set_id()
            self.create_side_set(
                external_side_set_id,
                self._get_external_element_faces(element_block_ids))
            self._convert_side_set_to_triangle_block(external_side_set_id,
                                                     triangle_element_block_id)
            self.delete_side_set(external_side_set_id)
        elif element_dimensions == set([2]) and element_types == set(['tri3']):
            # merge all blocks into another block
            ids_to_merge = []
            for id_ in element_block_ids:
                ids_to_merge.append(self._new_element_block_id())
                self.duplicate_element_block(id_, ids_to_merge[-1])
            triangle_element_block_id = self._new_element_block_id()
            self.combine_element_blocks(ids_to_merge,
                                        triangle_element_block_id)
            # these blocks better form a watertight enclosure
            external_sides = self._get_external_element_faces(
                triangle_element_block_id)
            if external_sides:
                self._warning(
                    'Enclosure is not watertight.',
                    'The "tri3" element blocks passed to this function do not '
                    'form a watertight enclosure.  An output file will still '
                    'be created, but it is unlikely to be useful as a WRL '
                    'model.  We found %d unmatched sides.' %
                    len(external_sides))
        else:
            self._error('Invalid input.',
                        'We expect a list of 3D element blocks or a list of '
                        '2D element blocks composed of triangles.')
        # now that we have a triangle block, subdivide it
        subdivided_element_block_id = self._new_element_block_id()
        # set up intervals
        if field_range == 'auto':
            indices = self.get_nodes_in_element_block(
                triangle_element_block_id)
            values = self.node_fields[node_field_name][timestep_index]
            this_values = [values[x] for x in indices]
            field_range = [min(this_values), max(this_values)]
        if not field_range[0] < field_range[1]:
            self._error('Invalid field range.',
                        'The given field range [%g - %g] was not valid.' %
                        (field_range[0], field_range[1]))
        if isinstance(intervals, int):
            intervals = [(1 - phi) * field_range[0] + phi * field_range[1]
                         for phi in [x / float(intervals)
                                     for x in xrange(1, intervals)]]
        # partition that block based on the node field
        self._partition_triangle_block_from_node_field(
            triangle_element_block_id,
            subdivided_element_block_id,
            node_field_name,
            timestep,
            intervals)
        # save model as an ExodusII file so we can view it
        if export_exodus_copy:
            temp_nodes = self.nodes
            self.nodes = self._get_coordinates_at_time(displacement_timestep)
            exodus_filename = filename.rsplit('.', 1)[0] + '_wrl.e'
            self.export_model(exodus_filename,
                              element_block_ids=subdivided_element_block_id,
                              element_field_names='interval',
                              node_field_names='none',
                              node_set_ids='none',
                              side_set_ids='none',
                              global_variable_names='none',
                              timesteps=[timestep])
            self.nodes = temp_nodes
        # save triangles to the WRL file
        with open(filename, 'w') as output:
            output.write('#VRML V2.0 utf8\n')
            output.write('Transform {\n')
            output.write('  translation 0.0 0.0 0.0\n')
            output.write('  children Shape {\n')
            output.write('    appearance Appearance {material Material {} }\n')
            output.write('    geometry IndexedFaceSet {\n')
            output.write('      coord Coordinate {\n')
            output.write('        point [\n')
            for position in self._get_coordinates_at_time(
                    displacement_timestep):
                output.write(','.join(str(x) for x in position))
                output.write(',')
            output.write('\n        ]\n')
            output.write('      }\n')
            output.write('      coordIndex [\n')
            connectivity = self.get_connectivity(subdivided_element_block_id)
            element_count = self.get_element_count(subdivided_element_block_id)
            for element_index in xrange(element_count):
                output.write('%d,%d,%d,-1,' %
                             (connectivity[element_index * 3 + 0],
                              connectivity[element_index * 3 + 1],
                              connectivity[element_index * 3 + 2]))
            output.write('\n      ]\n')
            output.write('      color Color {\n')
            output.write('        color [\n')
            if colorspace == 'rgb':
                colorspace = [
                    colorsys.hls_to_rgb(
                        2.0 / 3 - 2.0 / 3 * x / len(intervals),
                        0.5,
                        1)
                    for x in xrange(len(intervals) + 1)]
            if len(colorspace) != len(intervals) + 1:
                self._error('Unrecognized colorspace.',
                            'The given colorspace was not recognized.  We '
                            'expected a string such as "rgb" or a list of '
                            'length %d with RGB triplets.  Instead, we found '
                            '"%s".' % (len(intervals) + 1, str(colorspace)))
            for color in colorspace:
                output.write('          %s,\n' %
                             ' '.join(str(x) for x in color))
            output.write('        ]\n')
            output.write('      }\n')
            output.write('      colorIndex [\n')
            output.write('        ')
            values = self.get_element_field_values(
                'interval',
                subdivided_element_block_id,
                self.timesteps[timestep_index])
            output.write(','.join(str(int(x)) for x in values))
            output.write('\n')
            output.write('      ]\n')
            output.write('      colorPerVertex FALSE\n')
            output.write('    }\n')
            output.write('  }\n')
            output.write('}\n')
        # delete the temporary element blocks
        self.delete_element_block([triangle_element_block_id,
                                   subdivided_element_block_id])

    def export_stl_file(self,
                        filename,
                        element_block_ids='all',
                        displacement_timestep='auto'):
        """
        Export the exterior of the model to an STL file.

        By default, if timesteps exist and a displacement field exists, the
        displacements at the last timestep will be applied.

        Example:
        >>> model.export_stl_file('mesh_surface.stl')

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        if displacement_timestep == 'auto':
            if self.timesteps and self.displacement_field_exists():
                displacement_timestep = 'last'
            else:
                displacement_timestep = 'none'
        # create a new element block composed of triangles
        triangle_element_block_id = self._new_element_block_id()
        external_side_set_id = self._new_side_set_id()
        self.create_side_set(
            external_side_set_id,
            self._get_external_element_faces(element_block_ids))
        self._convert_side_set_to_triangle_block(external_side_set_id,
                                                 triangle_element_block_id)
        # export that set
        displacement_timestep = self._format_id_list(
            displacement_timestep,
            self.get_timesteps(),
            'timestep')
        # get coordinates at the time we care about
        c = self._get_coordinates_at_time(displacement_timestep)
        # create the STL file
        connectivity = self.get_connectivity(triangle_element_block_id)
        with open(filename, 'wb') as output:
            output.write(' ' * 80)
            output.write(struct.pack('<l', len(connectivity) / 3))
            for element_index in xrange(len(connectivity) / 3):
                tri = connectivity[element_index * 3:(element_index + 1) * 3]
                normal = (
                    c[tri[1]][1] * c[tri[2]][2] - c[tri[2]][1] * c[tri[1]][2],
                    c[tri[1]][2] * c[tri[2]][0] - c[tri[2]][2] * c[tri[1]][0],
                    c[tri[1]][0] * c[tri[2]][1] - c[tri[2]][0] * c[tri[1]][1])
                scale = (normal[0] ** 2 + normal[1] ** 2 +
                         normal[2] ** 2) ** 0.5
                if scale > 0.0:
                    normal = tuple(x / scale for x in normal)
                output.write(struct.pack('<3f', *normal))
                for vertex in tri:
                    output.write(struct.pack('<3f', *c[vertex]))
                output.write(struct.pack('<h', 0))
        # delete the temporary block we created
        self.delete_element_block(triangle_element_block_id)

    def export(self, filename, *args, **kwargs):
        """
        Export a model based on the filename extension.

        The following extensions will call the appropriate functions:
        * 'WRL --> export_wrl_model'
        * 'STL --> export_stl_file'
        * 'E, G, EXO --> export_model'

        Arguments are passed through to the export function.

        Examples:
        >>> model.export('result.e')
        >>> model.export('result.stl')
        >>> model.export('result.wrl', 'eqps')

        """
        exporters = dict()
        exporters['wrl'] = self.export_wrl_model
        exporters['stl'] = self.export_stl_file
        exporters['e'] = self.export_model
        exporters['g'] = self.export_model
        exporters['exo'] = self.export_model
        extension = filename.rsplit('.', 1)[-1].lower()
        if extension not in exporters:
            self._error('Unrecognized file extension.',
                        'The filename extension "%s" was not recognized.  The '
                        'list of recognized extensions is : %s' %
                        (extension, ', '.join(exporters.keys())))
        exporters[extension](filename, *args, **kwargs)

    def _error_evaluating_expression(self, expression, var):
        """Throw an error saying we could not evaluate the given expression."""
        self._error('Invalid expression',
                    'An error occurred while trying to evaluate the given '
                    'expression.  It is likely that this expression is '
                    'ill-formed, or that the variables it uses do not '
                    'exist:\n\n%s\n\nDefined variables: %s'
                    % (expression,
                       ', '.join(sorted(var.keys()))))

    @staticmethod
    def _transform_eval_expression(expression, variable_names):
        """Transform a string expression into one usable by eval."""
        eval_expression = ' %s ' % expression
        for name in variable_names:
            eval_expression = re.sub('([^a-zA-Z0-9_])%s([^a-zA-Z0-9_])' % name,
                                     r"\1var['%s']\2" % name,
                                     eval_expression)
        # add "math." to supported mathematical functions
        # Note: abs is a built-in function, it is not math.abs
        for name in ['atan', 'sinh', 'cosh', 'tanh', 'exp',
                     'sqrt', 'sin', 'cos', 'tan']:
            eval_expression = re.sub('([^a-zA-Z0-9_])%s([^a-zA-Z0-9_])' % name,
                                     r'\1math.%s\2' % name,
                                     eval_expression)
        # convert powers to from "^" to "**" format
        eval_expression = re.sub(r'\^', '**', eval_expression)
        return eval_expression.strip()

    def get_element_block_dimension(self, element_block_id):
        """
        Return the dimension of elements within this element block.

        If this cannot be determined, return -1.
        """
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        element_type = self._get_element_type(element_block_id)
        element_type = self._get_standard_element_type(element_type,
                                                       warning=False)
        if element_type in self.DIMENSION:
            return self.DIMENSION[element_type]
        else:
            return -1

    def get_nodes_per_element(self, element_block_id):
        """Return the nodes per element in the given element block."""
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        return self.element_blocks[element_block_id][1][2]

    def _get_element_type(self, element_block_id):
        """
        Return the element type for the given element block.

        This is a string value stored within the ExodusII file.  There is no
        standard naming convention and different programs can call the same
        element type by different names.

        """
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        return self.element_blocks[element_block_id][1][0]

    def _calculate_element_field_extreme(self,
                                         element_field_names,
                                         element_block_ids,
                                         function,
                                         function_name,
                                         calculate_block_id=False,
                                         calculate_location=False):
        """Store an element field extreme as a global variable."""
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        element_field_names = self._format_id_list(
            element_field_names,
            self.get_element_field_names(element_block_ids),
            'element field')
        for element_field_name in element_field_names:
            for element_block_id in element_block_ids:
                if not self.element_field_exists(element_field_name,
                                                 element_block_ids):
                    self._missing_on_entity_error(
                        element_field_name,
                        'element field',
                        element_block_id,
                        'element block')
        for element_field_name in element_field_names:
            name = element_field_name + '_' + function_name
            self.create_global_variable(name)
            if calculate_location:
                self.create_global_variable(name + '_x')
                self.create_global_variable(name + '_y')
                self.create_global_variable(name + '_z')
            if calculate_block_id:
                self.create_global_variable(name + '_block_id')
            for index in xrange(len(self.timesteps)):
                extreme = None
                extreme_block_id = None
                extreme_node_indices = None
                for element_block_id in element_block_ids:
                    fields = self._get_element_block_fields(element_block_id)
                    values = fields[element_field_name][index]
                    if not values:
                        continue
                    this_extreme = function(values)
                    if (extreme is None or
                            function([extreme, this_extreme]) == this_extreme):
                        extreme = this_extreme
                        extreme_block_id = element_block_id
                        connectivity = self.get_connectivity(element_block_id)
                        nodes_per_element = self.get_nodes_per_element(
                            element_block_id)
                        element_index = values.index(extreme)
                        extreme_node_indices = connectivity[
                            element_index * nodes_per_element:
                            (element_index + 1) * nodes_per_element]
                if extreme is None:
                    self._warning('No values encountered',
                                  'The element field extreme cannot be '
                                  'calculated because no values are present.')
                    continue
                self.global_variables[name][index] = extreme
                if calculate_block_id:
                    self.global_variables[
                        name + '_block_id'][index] = float(extreme_block_id)
                if calculate_location:
                    coords = self._get_coordinates_at_time(
                        self.timesteps[index])
                    locations = [coords[x] for x in extreme_node_indices]
                    centroid = [self._mean([x[i] for x in locations])
                                for i in xrange(3)]
                    self.global_variables[name + '_x'][index] = centroid[0]
                    self.global_variables[name + '_y'][index] = centroid[1]
                    self.global_variables[name + '_z'][index] = centroid[2]

    def calculate_element_field_maximum(self,
                                        element_field_names,
                                        element_block_ids='all',
                                        calculate_location=False,
                                        calculate_block_id=False):
        """
        Store an element field maximum as a global variable.

        If 'calculate_block_id=True', also store the element block id
        containing the element with this value.

        If 'calculate_location=True', also store the location of the centroid
        of the element with this value.

        Examples:
        >>> model.calculate_element_field_maximum('eqps')
        >>> model.global_variables['eqps_max']

        """
        self._calculate_element_field_extreme(
            element_field_names,
            element_block_ids,
            max,
            'max',
            calculate_location=calculate_location,
            calculate_block_id=calculate_block_id)

    def calculate_element_field_minimum(self,
                                        element_field_names,
                                        element_block_ids='all',
                                        calculate_location=False,
                                        calculate_block_id=False):
        """
        Store an element field minimum as a global variable.

        If 'calculate_block_id=True', also store the element block id
        containing the element with this value.

        If 'calculate_location=True', also store the location of the centroid
        of the element with this value.

        Examples:
        >>> model.calculate_element_field_minimum('eqps')
        >>> model.global_variables['eqps_min']

        """
        self._calculate_element_field_extreme(
            element_field_names,
            element_block_ids,
            min,
            'min',
            calculate_location=calculate_location,
            calculate_block_id=calculate_block_id)

    def _calculate_node_field_extreme(self,
                                      node_field_names,
                                      function,
                                      function_name,
                                      element_block_ids='auto',
                                      calculate_location=False):
        """Store a node field extreme as a global variable."""
        node_field_names = self._format_id_list(
            node_field_names,
            self.get_node_field_names(),
            'node field')
        if element_block_ids == 'auto':
            node_indices = range(len(self.nodes))
        else:
            node_indices = self.get_nodes_in_element_block(element_block_ids)
        if not node_indices:
            self._error('No node values',
                        'No nodes were specified, so we cannot calculate the '
                        'extreme value.')
        for node_field_name in node_field_names:
            name = node_field_name + '_' + function_name
            self.create_global_variable(name)
            if calculate_location:
                self.create_global_variable(name + '_x')
                self.create_global_variable(name + '_y')
                self.create_global_variable(name + '_z')
            for index in xrange(len(self.timesteps)):
                values = self.node_fields[node_field_name][index]
                values = [values[x] for x in node_indices]
                extreme = function(values)
                self.global_variables[name][index] = extreme
                if calculate_location:
                    coords = self._get_coordinates_at_time(
                        self.timesteps[index])
                    node_index = node_indices[values.index(extreme)]
                    location = coords[node_index]
                    self.global_variables[name + '_x'][index] = location[0]
                    self.global_variables[name + '_y'][index] = location[1]
                    self.global_variables[name + '_z'][index] = location[2]

    def calculate_node_field_maximum(self,
                                     node_field_names,
                                     element_block_ids='auto',
                                     calculate_location=False):
        """
        Store a node field maximum as a global variable.

        If 'calculate_location=True', also store the location of the node with
        this value.

        Examples:
        >>> model.calculate_node_field_maximum('temp')
        >>> model.global_variables['temp_max']

        """
        self._calculate_node_field_extreme(
            node_field_names,
            max,
            'max',
            element_block_ids=element_block_ids,
            calculate_location=calculate_location)

    def calculate_node_field_minimum(self,
                                     node_field_names,
                                     element_block_ids='auto',
                                     calculate_location=False):
        """
        Store a node field minimum as a global variable.

        If 'calculate_location=True', also store the location of the node with
        this value.

        Examples:
        >>> model.calculate_node_field_maximum('temp')
        >>> model.global_variables['temp_min']

        """
        self._calculate_node_field_extreme(
            node_field_names,
            min,
            'min',
            element_block_ids=element_block_ids,
            calculate_location=calculate_location)

    def output_global_variables(self,
                                filename=None,
                                global_variable_names='all',
                                timesteps='all'):
        """
        Output global variables in CSV format.

        If 'filename=None', output will be sent to stdout.  Else, it will be
        output to the given filename.

        By default, information is output for all timesteps and all global
        variables.  This may be changed by specifying which timesteps and
        variables are output.

        Examples:
        >>> model.output_global_variables('variables.csv')
        >>> model.output_global_variables('variables.csv', timesteps='last')

        """
        global_variable_names = self._format_id_list(
            global_variable_names,
            self.get_global_variable_names(),
            'global variable')
        timesteps = self._format_id_list(
            timesteps,
            self.get_timesteps(),
            'timestep')
        # output header
        names = ['time']
        names.extend(global_variable_names)
        output = ', '.join(names) + '\n'
        # output values
        for timestep in timesteps:
            index = self.timesteps.index(timestep)
            values = [timestep]
            for name in global_variable_names:
                values.append(self.global_variables[name][index])
            output += ', '.join(str(x) for x in values) + '\n'
        # output to file or screen
        if filename:
            open(filename, 'w').write(output)
        else:
            if output:
                output = output[:-1]
            print output

    def calculate_global_variable(self, expression):
        """
        Store a global variable calculated from the given expression.

        The expression may include the following variables:
        * 'time' to refer to the current time
        * global variables (by name)

        Examples:
        >>> model.calculate_global_variable('time_squared = time ^ 2')
        >>> model.calculate_global_variable('total = potential + kinetic')

        """
        if '=' not in expression:
            self._error('Invalid expression',
                        'A "=" sign must be present in the expression but '
                        'was not found.\n\nExpression: %s' % expression)
        (name, expression) = expression.split('=', 1)
        new_name = name.strip()
        self.create_global_variable(new_name)
        # create list of variable names and modify them in the expression
        variable_names = set(['time'])
        variable_names.update(self.get_global_variable_names())
        eval_expression = self._transform_eval_expression(expression,
                                                          variable_names)
        function = eval('lambda var: ' + eval_expression)
        var = dict()
        try:
            for index, time in enumerate(self.timesteps):
                var['time'] = time
                for name, values in self.global_variables.items():
                    var[name] = values[index]
                value = float(function(var))
                self.global_variables[new_name][index] = value
        except (SyntaxError, NameError):
            self._error_evaluating_expression(
                "%s = %s" % (new_name, eval_expression),
                var)

    def calculate_node_field(self, expression):
        """
        Store a node field calculated from the given expression.

        The expression may include the following variables:
        * 'time' to refer to the current time
        * global variables (by name)
        * node fields (by name)
        * model coordinates ('X', 'Y', and 'Z')

        Example:
        >>> model.calculate_node_field('temp_C = temp_K - 273.15')

        """
        if '=' not in expression:
            self._error('Invalid expression',
                        'A "=" sign must be present in the expression but '
                        'was not found.\n\nExpression: %s' % expression)
        (name, expression) = expression.split('=', 1)
        new_name = name.strip()
        self.create_node_field(new_name)
        new_values = self.node_fields[new_name]
        # create list of variable names and modify them in the expression
        variable_names = set(['time'])
        variable_names.update(['X', 'Y', 'Z'])
        variable_names.update(self.get_global_variable_names())
        variable_names.update(self.get_node_field_names())
        eval_expression = self._transform_eval_expression(expression,
                                                          variable_names)
        var = dict()
        function = eval('lambda var: ' + eval_expression)
        try:
            for time_index, time in enumerate(self.timesteps):
                # set time
                var['time'] = time
                # set global variables
                for name, values in self.global_variables.items():
                    var[name] = values[time_index]
                # go through each node
                for node_index in xrange(len(self.nodes)):
                    # set coordinate values
                    var['X'] = self.nodes[node_index][0]
                    var['Y'] = self.nodes[node_index][1]
                    var['Z'] = self.nodes[node_index][2]
                    # set node field values
                    for name, values in self.node_fields.items():
                        var[name] = values[time_index][node_index]
                    value = float(function(var))
                    new_values[time_index][node_index] = value
        except (SyntaxError, NameError):
            self._error_evaluating_expression(
                "%s = %s" % (new_name, eval_expression),
                var)

    def calculate_node_set_field(self, expression, node_set_ids='all'):
        """
        Store a node set field calculated from the given expression.

        The expression may include the following variables:
        * 'time' to refer to the current time
        * global variables (by name)
        * node fields (by name)
        * node set fields (by name)
        * model coordinates ('X', 'Y', and 'Z')

        Example:
        >>> model.calculate_node_set_field('temp_C = temp_K - 273.15')

        """
        node_set_ids = self._format_node_set_id_list(node_set_ids)
        if '=' not in expression:
            self._error('Invalid expression',
                        'A "=" sign must be present in the expression but '
                        'was not found.\n\nExpression: %s' % expression)
        (name, expression) = expression.split('=', 1)
        new_name = name.strip()
        # for each node set
        for node_set_id in node_set_ids:
            self.create_node_set_field(new_name, node_set_id)
            members = self.get_node_set_members(node_set_id)
            fields = self._get_node_set_fields(node_set_id)
            new_values = fields[new_name]
            # create list of variable names and modify them in the expression
            variable_names = set(['time'])
            variable_names.update(['X', 'Y', 'Z'])
            variable_names.update(self.get_global_variable_names())
            variable_names.update(self.get_node_field_names())
            variable_names.update(fields.keys())
            eval_expression = self._transform_eval_expression(expression,
                                                              variable_names)
            function = eval('lambda var: ' + eval_expression)
            var = dict()
            try:
                for time_index, time in enumerate(self.timesteps):
                    # set time
                    var['time'] = time
                    # set global variables
                    for name, values in self.global_variables.items():
                        var[name] = values[time_index]
                    # go through each node
                    for member_index, node_index in enumerate(members):
                        # set coordinate values
                        var['X'] = self.nodes[node_index][0]
                        var['Y'] = self.nodes[node_index][1]
                        var['Z'] = self.nodes[node_index][2]
                        # set node fields
                        for name, values in self.node_fields.items():
                            var[name] = values[time_index][node_index]
                        # set node set fields
                        for name, values in fields.items():
                            var[name] = values[time_index][member_index]
                        value = float(function(var))
                        new_values[time_index][member_index] = value
            except (SyntaxError, NameError):
                self._error_evaluating_expression(
                    "%s = %s" % (new_name, eval_expression),
                    var)

    def calculate_side_set_field(self, expression, side_set_ids='all'):
        """
        Store a side set field calculated from the given expression.

        The expression may include the following variables:
        * 'time' to refer to the current time
        * global variables (by name)

        Example:
        >>> model.calculate_side_set_field('temp_C = temp_K - 273.15')

        """
        side_set_ids = self._format_side_set_id_list(side_set_ids)
        if '=' not in expression:
            self._error('Invalid expression',
                        'A "=" sign must be present in the expression but '
                        'was not found.\n\nExpression: %s' % expression)
        (name, expression) = expression.split('=', 1)
        new_name = name.strip()
        # for each side set
        for side_set_id in side_set_ids:
            self.create_side_set_field(new_name, side_set_id)
            members = self.get_side_set_members(side_set_id)
            fields = self._get_side_set_fields(side_set_id)
            new_values = fields[new_name]
            # create list of variable names and modify them in the expression
            variable_names = set(['time'])
            variable_names.update(self.get_global_variable_names())
            variable_names.update(fields.keys())
            eval_expression = self._transform_eval_expression(expression,
                                                              variable_names)
            function = eval('lambda var: ' + eval_expression)
            var = dict()
            try:
                for time_index, time in enumerate(self.timesteps):
                    # set time
                    var['time'] = time
                    # set global variables
                    for name, values in self.global_variables.items():
                        var[name] = values[time_index]
                    # go through each face
                    for member_index in xrange(len(members)):
                        # set side set fields
                        for name, values in fields.items():
                            var[name] = values[time_index][member_index]
                        value = float(function(var))
                        new_values[time_index][member_index] = value
            except (SyntaxError, NameError):
                self._error_evaluating_expression(
                    "%s = %s" % (new_name, eval_expression),
                    var)

    def calculate_element_field(self,
                                expression,
                                element_block_ids='all',
                                timesteps='all'):
        """
        Store an element field calculated from the given expression.

        The expression may include the following variables:
        * 'time' to refer to the current time
        * global variables (by name)
        * element fields (by name)

        Example:
        >>> model.calculate_element_field('pressure = (stress_xx + '
        ...                               'stress_yy + stress_zz) / -3')

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        if '=' not in expression:
            self._error('Invalid expression',
                        'A "=" sign must be present in the expression but '
                        'was not found.\n\nExpression: %s' % expression)
        (name, expression) = expression.split('=', 1)
        new_name = name.strip()
        # for each element block
        for element_block_id in element_block_ids:
            self.create_element_field(new_name, element_block_id)
            fields = self._get_element_block_fields(element_block_id)
            element_count = self.get_element_count(element_block_id)
            new_values = fields[new_name]
            # create list of variable names and modify them in the expression
            variable_names = set(['time'])
            variable_names.update(self.get_global_variable_names())
            variable_names.update(fields.keys())
            eval_expression = self._transform_eval_expression(expression,
                                                              variable_names)
            function = eval('lambda var: ' + eval_expression)
            var = dict()
            try:
                for time_index, time in enumerate(self.timesteps):
                    # set time
                    var['time'] = time
                    # set global variables
                    for name, values in self.global_variables.items():
                        var[name] = values[time_index]
                    # go through element
                    for index in xrange(element_count):
                        # set element fields
                        for name, values in fields.items():
                            var[name] = values[time_index][index]
                        value = float(function(var))
                        new_values[time_index][index] = value
            except (SyntaxError, NameError):
                self._error_evaluating_expression(
                    "%s = %s" % (new_name, eval_expression),
                    var)

    def to_lowercase(self):
        """
        Convert the names of all entities within the model to lowercase.

        This affects global variables, node fields, element fields, side set
        fields, and node set fields.

        Example:
        >>> model.to_lowercase()

        """
        for name in self.get_global_variable_names():
            if name != name.lower():
                self.rename_global_variable(name, name.lower())
        for name in self.get_node_field_names():
            if name != name.lower():
                self.rename_node_field(name, name.lower())
        for id_ in self.get_element_block_ids():
            for name in self.get_element_field_names(id_):
                if name != name.lower():
                    self.rename_element_field(name, name.lower(), id_)
        for id_ in self.get_side_set_ids():
            for name in self.get_side_set_field_names(id_):
                if name != name.lower():
                    self.rename_side_set_field(name, name.lower(), id_)
        for id_ in self.get_node_set_ids():
            for name in self.get_node_set_field_names(id_):
                if name != name.lower():
                    self.rename_node_set_field(name, name.lower(), id_)

    def _translate_element_blocks(self, element_block_translation_list):
        """
        Convert element blocks to a new element type.

        >>> model._translate_element_blocks([(1, 'hex20'), (2, 'tet10')])

        """
        # validity check the inputs
        valid_element_block_translation_list = []
        for (id_, new_element_type) in element_block_translation_list:
            # check if element block exists
            if not self.element_block_exists(id_):
                self._missing_warning(id_, 'element block')
                continue
            old_element_type = self._get_standard_element_type(
                self._get_element_type(id_))
            # check if any schemes exist for this element block
            if old_element_type not in self.ELEMENT_CONVERSIONS:
                self._error('No element conversion schemes',
                            'There are no element conversion schemes to '
                            'convert elements of type "%s" into any other '
                            'element type.'
                            % (old_element_type))
                continue
            # check if the given scheme exists for this element block
            conversions = self.ELEMENT_CONVERSIONS[old_element_type]
            if new_element_type not in conversions:
                self._error('Invalid target element type',
                            'There are no element conversion schemes to '
                            'convert elements of type "%s" into elements of '
                            'type "%s".  Conversions are available to the '
                            'following element types: "%s".'
                            % (old_element_type,
                               new_element_type,
                               '", "'.join(conversions.keys())))
                continue
            # all is good, keep this one on the list
            valid_element_block_translation_list.append([id_,
                                                         new_element_type])
        # if we have nothing to do, just exit
        if not valid_element_block_translation_list:
            return
        # now we need to find the new nodes to create, if any
        averaged_nodes = set()
        for (id_, new_element_type) in valid_element_block_translation_list:
            old_element_type = self._get_standard_element_type(
                self._get_element_type(id_))
            scheme = self.ELEMENT_CONVERSIONS[old_element_type]
            scheme = scheme[new_element_type]
            # find new nodes which need created
            averaged_node_list = set()
            for formula in scheme:
                for new_node in formula:
                    # if it's just an integer, we don't need a new node created
                    if isinstance(new_node, int):
                        continue
                    # new node
                    averaged_node_list.add(new_node)
            # now loop through all elements and create a list of averaged
            # nodes to create
            connectivity = self.get_connectivity(id_)
            nodes_per_element = self.get_nodes_per_element(id_)
            # iterate over each element
            for local_node in [connectivity[x:x + nodes_per_element]
                               for x in range(0,
                                              len(connectivity),
                                              nodes_per_element)]:
                # iterate over each new averaged node formula
                for formula in averaged_node_list:
                    node_list = tuple(sorted(local_node[x] for x in formula))
                    averaged_nodes.add(node_list)
        # create the averaged node mapping
        averaged_nodes = sorted(averaged_nodes)
        averaged_node_map = dict((x, len(self.nodes) + index)
                                 for index, x in enumerate(averaged_nodes))
        # create new averaged nodes
        self._create_averaged_nodes(averaged_nodes, [])
        # now create the new element blocks
        for (element_block_id,
             new_element_type) in valid_element_block_translation_list:
            old_element_type = self._get_standard_element_type(
                self._get_element_type(element_block_id))
            scheme = self.ELEMENT_CONVERSIONS[old_element_type]
            scheme = scheme[new_element_type]
            # rename some things
            connectivity = self.get_connectivity(element_block_id)
            nodes_per_element = self.get_nodes_per_element(element_block_id)
            # create the connectivity for the new element block
            new_connectivity = []
            for local_node in [connectivity[x:x + nodes_per_element]
                               for x in range(0,
                                              len(connectivity),
                                              nodes_per_element)]:
                for new_element in scheme:
                    for new_node in new_element:
                        if isinstance(new_node, int):
                            new_connectivity.append(local_node[new_node])
                        else:
                            node_list = tuple(sorted(local_node[x]
                                                     for x in new_node))
                            new_connectivity.append(
                                averaged_node_map[node_list])
            # create a temporary element block
            temporary_element_block_id = self._new_element_block_id()
            new_nodes_per_element = self.NODES_PER_ELEMENT[new_element_type]
            self.create_element_block(
                temporary_element_block_id,
                [new_element_type,
                 len(new_connectivity) / new_nodes_per_element,
                 new_nodes_per_element,
                 0],
                new_connectivity)
            temporary_fields = self._get_element_block_fields(
                temporary_element_block_id)
            element_multiplier = len(scheme)
            # transfer element values
            fields = self._get_element_block_fields(element_block_id)
            for field_name in self.get_element_field_names(element_block_id):
                values = fields[field_name]
                new_values = [list(x) for x in values]
                new_values = [[x
                               for x in these_values
                               for _ in range(element_multiplier)]
                              for these_values in new_values]
                temporary_fields[field_name] = new_values
            # for each face in the old scheme, find all of its nodes
            old_face_mapping = self._get_face_mapping(old_element_type)
            old_face_nodes = [set(x) for _, x in old_face_mapping]
            new_face_mapping = self._get_face_mapping(new_element_type)
            face_translation = [[] for _ in range(len(old_face_nodes))]
            for new_element_index, new_element in enumerate(scheme):
                for new_face_index, (_, face_members) in enumerate(
                        new_face_mapping):
                    # find all nodes used by this new face
                    used_nodes = set()
                    for face_member in face_members:
                        if isinstance(new_element[face_member], int):
                            used_nodes.add(new_element[face_member])
                        else:
                            used_nodes.update(new_element[face_member])
                    # see if these are a subset of nodes in the old set
                    for (old_face_index,
                         old_members) in enumerate(old_face_nodes):
                        if used_nodes <= old_members:
                            face_translation[old_face_index].append(
                                (new_element_index, new_face_index))
            # update self.side_sets
            for id_ in self.get_side_set_ids():
                members = self.get_side_set_members(id_)
                fields = self._get_side_set_fields(id_)
                old_members = list(members)
                for member_index, (block_id,
                                   element_index,
                                   face_index) in enumerate(old_members):
                    if block_id == element_block_id:
                        # add some new faces
                        for new_faces in face_translation[face_index]:
                            members.append(
                                (temporary_element_block_id,
                                 element_index * len(scheme) + new_faces[0],
                                 new_faces[1]))
                        # add values for the new faces
                        for all_values in fields.values():
                            for values in all_values:
                                for _ in face_translation[face_index]:
                                    values.append(values[member_index])
            # delete the old block
            self.delete_element_block(element_block_id)
            # rename the temporary element block
            self.rename_element_block(temporary_element_block_id,
                                      element_block_id)

    def convert_element_blocks(self, element_block_ids, new_element_type):
        """Convert elements within a block to a new element type."""
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        # form the conversion list
        conversion_list = [(id_, new_element_type)
                           for id_ in element_block_ids]
        # convert each block
        self._translate_element_blocks(conversion_list)

    def _change_element_order(self, element_block_ids, target_order):
        """
        Convert elements to the given order.

        This will attempt to find the best element to convert to given the
        conversion schemes involved.  If more than one element option exists
        for a given element type, this will choose the option that produces
        the fewest nodes.

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        # find all element types we need to convert
        element_types = set()
        for element_block_id in element_block_ids:
            element_types.add(self._get_element_type(element_block_id))
        element_types = sorted(self._get_standard_element_type(x)
                               for x in element_types)
        # ignore non-standard element types
        element_types = [x
                         for x in element_types
                         if x in self.STANDARD_ELEMENT_TYPES]
        # find all elements of the target order
        target_element_types = [x
                                for x in self.STANDARD_ELEMENT_TYPES
                                if self.ELEMENT_ORDER[x] == target_order]
        # for each element type, find the optimal target type
        element_type_map = dict()
        for element_type in element_types:
            # if no conversion schemes, error out
            if element_type not in self.ELEMENT_CONVERSIONS:
                # if element type is already the desired order, no warning
                if self.ELEMENT_ORDER[element_type] == target_order:
                    continue
                self._error('Unable to convert',
                            'There is no valid scheme to convert elements '
                            'of type "%s" to an element type of order '
                            '%d.'
                            % (element_type, target_order))
            scheme = self.ELEMENT_CONVERSIONS[element_type]
            possible_types = set(scheme.keys()).intersection(
                target_element_types)
            if not possible_types:
                # if element type is already the desired order, no warning
                if self.ELEMENT_ORDER[element_type] == target_order:
                    continue
                self._error('Unable to convert',
                            'There is no valid scheme to convert elements '
                            'of type "%s" to an element type of order '
                            '%d.'
                            % (element_type, target_order))
            # out of the given options, choose the one that first creates the
            # least number of elements, and secondly the least number of nodes
            ranking = [(len(scheme[x]), len(scheme[x]) * len(scheme[x][0]), x)
                       for x in possible_types]
            ranking = sorted(ranking)
            element_type_map[element_type] = ranking[0][-1]
        # create the conversion list
        conversion_list = []
        for element_block_id in element_block_ids:
            element_type = self._get_element_type(element_block_id)
            element_type = self._get_standard_element_type(element_type)
            if element_type in element_type_map:
                conversion_list.append(
                    (element_block_id, element_type_map[element_type]))
        # now translate the elements
        self._translate_element_blocks(conversion_list)

    def make_elements_linear(self, element_block_ids='all'):
        """
        Convert elements in one or more element blocks to a linear type.

        This will attempt to find the best element to convert to given the
        conversion schemes involved.  If more than one element option exists
        for a given element type, this will choose the option that produces
        the fewest nodes.

        """
        self._change_element_order(element_block_ids, 1)

    def make_elements_quadratic(self, element_block_ids='all'):
        """
        Convert elements in one or more element blocks to a quadratic type.

        This will attempt to find the best element to convert to given the
        conversion schemes involved.  If more than one element option exists
        for a given element type, this will choose the option that produces
        the fewest nodes.

        """
        self._change_element_order(element_block_ids, 2)

    def _translate_element_type(self,
                                element_block_id,
                                new_element_type,
                                scheme):
        """Convert elements within a block to a new type."""
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        old_element_type = self._get_standard_element_type(
            self._get_element_type(element_block_id))
        new_element_type = self._get_standard_element_type(new_element_type)
        old_nodes_per_element = self.NODES_PER_ELEMENT[old_element_type]
        new_nodes_per_element = self.NODES_PER_ELEMENT[new_element_type]
        element_multiplier = len(scheme)
        # Make a list of nodes which need duplicated and nodes which need
        # averaged.  We will assign indices to these after we determine how
        # many there are of each.
        duplicate_nodes = set()
        averaged_nodes = set()
        connectivity = self.get_connectivity(element_block_id)
        element_count = self.get_element_count(element_block_id)
        for element_index in xrange(element_count):
            local_node = connectivity[
                element_index * old_nodes_per_element:
                (element_index + 1) * old_nodes_per_element]
            for new_element in scheme:
                for new_node in new_element:
                    if isinstance(new_node, int):
                        duplicate_nodes.add(local_node[new_node])
                    else:
                        averaged_nodes.add(
                            tuple(local_node[x] for x in new_node))
        # create new nodes
        next_node_index = len(self.nodes)
        duplicate_nodes = sorted(duplicate_nodes)
        averaged_nodes = sorted(averaged_nodes)
        self._duplicate_nodes(duplicate_nodes, [])
        self._create_averaged_nodes(averaged_nodes, [])
        # assign node indices
        duplicate_nodes = dict((x, next_node_index + index)
                               for index, x in enumerate(duplicate_nodes))
        next_node_index += len(duplicate_nodes)
        averaged_nodes = dict((x, next_node_index + index)
                              for index, x in enumerate(averaged_nodes))
        # create the connectivity for the new element block
        new_connectivity = []
        for element_index in xrange(element_count):
            local_node = connectivity[
                element_index * old_nodes_per_element:
                (element_index + 1) * old_nodes_per_element]
            for new_element in scheme:
                for new_node in new_element:
                    if isinstance(new_node, int):
                        new_connectivity.append(
                            duplicate_nodes[local_node[new_node]])
                    else:
                        new_connectivity.append(
                            averaged_nodes[tuple(local_node[x]
                                                 for x in new_node)])
        # create the new block
        temporary_element_block_id = self._new_element_block_id()
        self.create_element_block(
            temporary_element_block_id,
            [new_element_type,
             len(new_connectivity) / new_nodes_per_element,
             new_nodes_per_element,
             0],
            new_connectivity)
        temporary_element_fields = self._get_element_block_fields(
            temporary_element_block_id)
        # transfer element values
        fields = self._get_element_block_fields(element_block_id)
        for field_name in self.get_element_field_names(element_block_id):
            values = fields[field_name]
            new_values = [list(x) for x in values]
            new_values = [[x
                           for x in these_values
                           for _ in range(element_multiplier)]
                          for these_values in new_values]
            temporary_element_fields[field_name] = new_values
        # for each face in the old scheme, find all of its nodes
        old_face_mapping = self._get_face_mapping(old_element_type)
        old_face_nodes = [set(x) for _, x in old_face_mapping]
        new_face_mapping = self._get_face_mapping(new_element_type)
        face_translation = [[] for _ in xrange(len(old_face_nodes))]
        for new_element_index, new_element in enumerate(scheme):
            for new_face_index, (_, face_members) in enumerate(
                    new_face_mapping):
                # find all nodes used by this new face
                used_nodes = set()
                for face_member in face_members:
                    if isinstance(new_element[face_member], int):
                        used_nodes.add(new_element[face_member])
                    else:
                        used_nodes.update(new_element[face_member])
                # see if these are a subset of nodes in the old set
                for old_face_index, old_members in enumerate(old_face_nodes):
                    if used_nodes <= old_members:
                        face_translation[old_face_index].append(
                            (new_element_index, new_face_index))
        # update self.side_sets
        for side_set_id in self.get_side_set_ids():
            members = self.get_side_set_members(side_set_id)
            fields = self._get_side_set_fields(side_set_id)
            old_members = list(members)
            for member_index, (block_id,
                               element_index,
                               face_index) in enumerate(old_members):
                if block_id == element_block_id:
                    # add some new faces
                    for new_faces in face_translation[face_index]:
                        members.append(
                            (temporary_element_block_id,
                             element_index * len(scheme) + new_faces[0],
                             new_faces[1]))
                    # add values for the new faces
                    for all_values in fields.values():
                        for values in all_values:
                            for _ in face_translation[face_index]:
                                values.append(values[member_index])
        # delete the old block and rename the new one
        self.delete_element_block(element_block_id)
        self.rename_element_block(temporary_element_block_id, element_block_id)

    def convert_hex8_block_to_tet4_block(self,
                                         element_block_id,
                                         scheme='hex24tet'):
        """
        Convert a block of 'hex8' elements to a block of 'tet4' elements.

        Side sets are updated accordingly.

        Node sets are not modified.

        Currently, only the 'hex24tet' scheme is implemented, which creates 24
        'tet4' element for eeach 'hex8' element in the original mesh.

        Example:
        >>> model.convert_hex8_block_to_tet4_block(1)

        """
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        # ensure source block is actually hex8
        source_element_type = self._get_standard_element_type(
            self._get_element_type(element_block_id))
        if source_element_type != 'hex8':
            self._error('Incompatible element type.',
                        'We were expecting an element block composed of '
                        '"hex8" element but instead encountered one with "%s" '
                        'elements.' % source_element_type)
        # get the chosen scheme
        if scheme != 'hex24tet':
            self._error('Unsupported scheme.',
                        'The scheme "%s" was not recognized.' % scheme)
        # create the scheme
        scheme = []
        for _, face in self._get_face_mapping('hex8'):
            for index in xrange(4):
                scheme.append(tuple([face[(index + 1) % 4],
                                     face[index],
                                     tuple(sorted(face)),
                                     tuple(range(8))]))
        # apply it
        self._translate_element_type(element_block_id, 'tet4', scheme)

    def displacement_field_exists(self):
        """
        Return 'True' if the displacement field exists.

        Example:
        >>> model.displacement_field_exists()

        """
        prefix = self._get_displacement_field_prefix()
        return (self.node_field_exists(prefix + '_x') and
                self.node_field_exists(prefix + '_y') and
                self.node_field_exists(prefix + '_z'))

    def create_displacement_field(self):
        """
        Create the displacement field if it doesn't already exist.

        Example:
        >>> model.create_displacement_field()

        """
        prefix = self._get_displacement_field_prefix()
        for component in ['x', 'y', 'z']:
            this_name = prefix + '_' + component
            if not self.node_field_exists(this_name):
                self.create_node_field(this_name)

    def _get_displacement_field_values(self):
        """
        Return a list of the displacement field names.

        This will create them if they do not already exist.

        """
        self.create_displacement_field()
        prefix = self._get_displacement_field_prefix()
        displacement_field_list = []
        for component in ['x', 'y', 'z']:
            this_name = prefix + '_' + component
            displacement_field_list.append(self.node_fields[this_name])
        return displacement_field_list

    def _apply_displacements(self, timestep='last', scale_factor=1.0):
        """
        Apply the displacements to each node.

        If 'dx, dy, dz' are the components of the displacement fields, this
        applies the following transformation.
        * 'X --> X + dx'
        * 'Y --> Y + dy'
        * 'Z --> Z + dz'

        The nodal coordinates are modified but the displacement field itself
        is unchanged.

        """
        timestep = self._format_id_list(
            timestep,
            self.get_timesteps(),
            'timestep')
        if len(timestep) > 1:
            self._error('Ambiguous timestep.',
                        'More than one timestep was specified.  We expected '
                        'one or zero timesteps.')
        if len(timestep) == 0:
            return [tuple(x) for x in self.nodes]
        timestep_index = self.timesteps.index(timestep[0])
        displacement_values = [x[timestep_index]
                               for x in self._get_displacement_field_values()]
        new_nodes = [(x + dx * scale_factor,
                      y + dy * scale_factor,
                      z + dz * scale_factor)
                     for (x, y, z), dx, dy, dz
                     in itertools.izip(self.nodes,
                                       *displacement_values)]
        self.nodes = new_nodes

    def _get_local_index(self, this_id, id_list, entity='entity'):
        """
        Return the local index corresponding to the given id.

        If an entity is not present, throw an error.

        Example:
        >>> model._get_local_index(10, [10, 20, 30])
        0
        >>> model._get_local_index('first', [10, 20, 30])
        10

        """
        if this_id == 'first':
            if not id_list:
                self._error('Undefined %s reference.' % entity,
                            'A reference to the first %s was encountered but '
                            'no %ss are defined.' % (entity, entity))
            return 0
        if this_id == 'last':
            if not id_list:
                self._error('Undefined %s reference.' % entity,
                            'A reference to the last %s was encountered but '
                            'no %ss are defined.' % (entity, entity))
            return len(id_list) - 1
        if this_id not in id_list:
            entity_list = ', '.join([str(x) for x in id_list])
            self._error('Reference to undefined %s.' % entity,
                        'A reference to %s "%s" was encountered but is not '
                        'defined in the model.  There are %d defined %ss: %s'
                        % (entity, str(this_id), len(id_list), entity,
                           entity_list))
        return id_list.index(this_id)

    def get_element_field_values(self,
                                 element_field_name,
                                 element_block_id='auto',
                                 timestep='last'):
        """
        Return the list of element field values.

        The actual list of values is returned, so any modifications to it will
        be stored in the model.

        Examples:
        >>> model.get_element_field_values('strain', element_block_id=1)
        >>> model.get_element_field_values('strain', timestep=2.0)
        >>> model.get_element_field_values('strain',
        ...                                element_block_id=5,
        ...                                timestep='last')

        """
        [element_field_name] = self._format_id_list(
            [element_field_name],
            self.get_element_field_names(),
            'element field',
            single=True)
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        [timestep] = self._format_id_list(
            [timestep],
            self.get_timesteps(),
            'timestep',
            single=True)
        timestep_index = self._get_internal_timestep_index(timestep)
        if not self.element_field_exists(element_field_name,
                                         element_block_id):
            self._missing_on_entity_error(
                element_field_name,
                'element field',
                element_block_id,
                'element block')
        fields = self._get_element_block_fields(element_block_id)
        return fields[element_field_name][timestep_index]

    def get_side_set_field_values(self,
                                  side_set_field_name,
                                  side_set_id='auto',
                                  timestep='last'):
        """
        Return the list of side set field values.

        The actual list of values is returned, so any modifications to it will
        be stored in the model.

        Examples:
        >>> model.get_side_set_field_values('contact_pressure', side_set_id=1)
        >>> model.get_side_set_field_values('contact_pressure', timestep=2.0)
        >>> model.get_side_set_field_values('contact_pressure',
        ...                                 side_set_id=5,
        ...                                 timestep='last')

        """
        [side_set_id] = self._format_side_set_id_list(
            [side_set_id],
            single=True)
        [side_set_field_name] = self._format_id_list(
            [side_set_field_name],
            self.get_side_set_field_names(),
            'side set field',
            single=True)
        [timestep] = self._format_id_list(
            [timestep],
            self.get_timesteps(),
            'timestep',
            single=True)
        timestep_index = self._get_internal_timestep_index(timestep)
        if not self.side_set_field_exists(side_set_field_name,
                                          side_set_id):
            self._missing_on_entity_error(
                side_set_field_name,
                'side set field',
                side_set_id,
                'side set')
        fields = self._get_side_set_fields(side_set_id)
        return fields[side_set_field_name][timestep_index]

    def get_node_set_field_values(self,
                                  node_set_field_name,
                                  node_set_id='auto',
                                  timestep='last'):
        """
        Return the list of node set field values.

        The actual list of values is returned, so any modifications to it will
        be stored in the model.

        Examples:
        >>> model.get_node_set_field_values('contact_pressure', node_set_id=1)
        >>> model.get_node_set_field_values('contact_pressure', timestep=2.0)
        >>> model.get_node_set_field_values('contact_pressure',
        ...                                 node_set_id=5,
        ...                                 timestep='last')

        """
        [node_set_id] = self._format_node_set_id_list(
            [node_set_id],
            single=True)
        [node_set_field_name] = self._format_id_list(
            [node_set_field_name],
            self.get_node_set_field_names(),
            'node set field',
            single=True)
        [timestep] = self._format_id_list(
            [timestep],
            self.get_timesteps(),
            'timestep',
            single=True)
        timestep_index = self._get_internal_timestep_index(timestep)
        if not self.node_set_field_exists(node_set_field_name,
                                          node_set_id):
            self._missing_on_entity_error(
                node_set_field_name,
                'node set field',
                node_set_id,
                'node set')
        fields = self._get_node_set_fields(node_set_id)
        return fields[node_set_field_name][timestep_index]

    def _is_standard_element_type(self, element_type):
        """Return 'True' if 'element_type' is a standard element type."""
        name = self._get_standard_element_type(element_type,
                                               warning=False)
        return name in self.STANDARD_ELEMENT_TYPES

    def get_element_block_ids(self):
        """
        Return a list of all element block ids.

        Example:
        >>> model.get_element_block_ids()

        """
        return sorted(self.element_blocks.keys())

    def get_element_block_name(self, element_block_id):
        """
        Return the given element block name.

        Example:
        >>> model.get_element_block_names()

        """
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        return self.element_blocks[element_block_id][0]

    def get_all_element_block_names(self):
        """
        Return a list of all element block names.

        Example:
        >>> model.get_all_element_block_names()

        """
        return sorted(x[0]
                      for x in self.element_blocks.values()
                      if x[0])

    def _get_standard_element_block_ids(self):
        """
        Return a list of element block ids which use standard element types.

        Example:
        >>> model._get_standard_element_block_ids()

        """
        return [x
                for x in self.get_element_block_ids()
                if self._is_standard_element_type(
                    self._get_element_type(x))]

    def get_node_set_ids(self):
        """
        Return a list of all node set ids.

        Example:
        >>> model.get_node_set_ids()

        """
        return sorted(self.node_sets.keys())

    def get_side_set_ids(self):
        """
        Return a list of all side set ids.

        Example:
        >>> model.get_side_set_ids()

        """
        return sorted(self.side_sets.keys())

    def get_node_field_names(self):
        """
        Return a list of all node field names.

        Example:
        >>> model.get_node_field_names()

        """
        return self._sort_field_names(list(self.node_fields.keys()))

    def get_node_set_field_names(self, node_set_ids='all'):
        """
        Return a list of all node set field names.

        By default, this returns node set fields which are defined on any
        node set.  To return fields which are defined on a particular
        node set, pass a node set id.

        Examples:
        >>> model.get_node_set_field_names()
        >>> model.get_node_set_field_names(1)

        """
        node_set_ids = self._format_node_set_id_list(node_set_ids)
        names = set()
        for id_ in node_set_ids:
            names.update(self.node_sets[id_][-1].keys())
        return self._sort_field_names(list(names))

    def get_side_set_field_names(self, side_set_ids='all'):
        """
        Return a list of all side set field names.

        By default, this returns side set fields which are defined on any
        side set.  To return fields which are defined on a particular
        side set, pass a side set id.

        Examples:
        >>> model.get_side_set_field_names()
        >>> model.get_side_set_field_names(1)

        """
        side_set_ids = self._format_side_set_id_list(side_set_ids)
        names = set()
        for id_ in side_set_ids:
            names.update(self._get_side_set_fields(id_).keys())
        return self._sort_field_names(list(names))

    def get_element_field_names(self, element_block_ids='all'):
        """
        Return a list of all element field names.

        By default, this returns element fields which are defined on any
        element block.  To return fields which are defined on a particular
        element block, pass an element block id.

        Examples:
        >>> model.get_element_field_names()
        >>> model.get_element_field_names(1)

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        names = set()
        for id_ in element_block_ids:
            names.update(self.element_blocks[id_][-1].keys())
        return self._sort_field_names(list(names))

    def get_global_variable_names(self):
        """
        Return a list of all global variable names.

        Example:
        >>> model.get_global_variable_names()

        """
        return self._sort_field_names(list(self.global_variables.keys()))

    def get_timesteps(self):
        """
        Return a list of all timesteps.

        Example:
        >>> model.get_timesteps()

        """
        return sorted(self.timesteps)

    def _get_internal_timestep_index(self, timestep):
        """Return the local timestep index."""
        [timestep] = self._format_id_list(
            [timestep],
            self.get_timesteps(),
            'timestep',
            single=True)
        return self.timesteps.index(timestep)

    def _create_element_field_truth_table(self,
                                          element_block_ids,
                                          field_names):
        """
        Return the element field truth table.

        The element field truth table defines which element fields are defined
        for which element blocks according to the given order of element block
        ids and field names passed in.  This is needed within the ExodusII
        file format.

        """
        # go through each case and set values to True if they exist
        truth_table = [
            field_name in self._get_element_block_fields(element_block_id)
            for element_block_id in element_block_ids
            for field_name in field_names]
        return truth_table

    def _create_side_set_field_truth_table(self, side_set_ids, field_names):
        """
        Return the side set field truth table.

        The side set field truth table defines which side set fields are
        defined for which side sets according to the given order of side set
        ids and field names passed in.  This is needed within the ExodusII
        file format.

        """
        # go through each case and set values to True if they exist
        truth_table = [
            self.side_set_field_exists(field_name, side_set_id)
            for side_set_id in side_set_ids
            for field_name in field_names]
        return truth_table

    def _create_node_set_field_truth_table(self, node_set_ids, field_names):
        """
        Return the node set field truth table.

        The node set field truth table defines which node set fields are
        defined for which node sets according to the given order of node set
        ids and field names passed in.  This is needed within the ExodusII
        file format.

        """
        # go through each case and set values to True if they exist
        truth_table = [
            self.node_set_field_exists(field_name, node_set_id)
            for node_set_id in node_set_ids
            for field_name in field_names]
        return truth_table

    def get_connectivity(self, element_block_id='auto'):
        """
        Return the connectivity list of an element block.

        Example:
        >>> model.get_connectivity(1)

        """
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        return self.element_blocks[element_block_id][2]

    def get_element_block_connectivity(self, element_block_id='auto'):
        """Alias for 'get_connectivity()'."""
        return self.get_connectivity(element_block_id)

    def get_nodes_in_element_block(self, element_block_ids):
        """
        Return a list of all node indices used in the given element blocks.

        Examples:
        >>> model.get_nodes_in_element_block(1)
        >>> model.get_nodes_in_element_block([1, 3])

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        node_list = set()
        for id_ in element_block_ids:
            connectivity = self.get_connectivity(id_)
            node_list.update(connectivity)
        return sorted(node_list)

    def _rotate_nodes(self,
                      axis,
                      angle_in_degrees,
                      node_indices='all',
                      adjust_displacement_field='auto'):
        """Rotate nodes about an axis by the given angle."""
        if adjust_displacement_field == 'auto':
            adjust_displacement_field = self.displacement_field_exists()
        # Create rotation matrix.
        # x --> R * x
        scale = math.sqrt(sum(x * x for x in axis))
        (ux, uy, uz) = [float(x) / scale for x in axis]
        theta = float(angle_in_degrees) * math.pi / 180
        cost = math.cos(theta)
        sint = math.sin(theta)
        # if angle_in_degrees is a multiple of 90 degrees, then make sin and
        # cos exactly 0, 1, or -1 to avoid roundoff errors.
        if angle_in_degrees % 90 == 0:
            sint = math.floor(sint + 0.5)
            cost = math.floor(cost + 0.5)
        rxx = cost + ux * ux * (1 - cost)
        rxy = ux * uy * (1 - cost) - uz * sint
        rxz = ux * uz * (1 - cost) + uy * sint
        ryx = uy * ux * (1 - cost) + uz * sint
        ryy = cost + uy * uy * (1 - cost)
        ryz = uy * uz * (1 - cost) - ux * sint
        rzx = uz * ux * (1 - cost) - uy * sint
        rzy = uz * uy * (1 - cost) + ux * sint
        rzz = cost + uz * uz * (1 - cost)
        # Rotate nodes.
        if node_indices == 'all':
            self.nodes = [[rxx * x + rxy * y + rxz * z,
                           ryx * x + ryy * y + ryz * z,
                           rzx * x + rzy * y + rzz * z]
                          for x, y, z in self.nodes]
        else:
            for index in node_indices:
                n = self.nodes[index]
                self.nodes[index] = [rxx * n[0] + rxy * n[1] + rxz * n[2],
                                     ryx * n[0] + ryy * n[1] + ryz * n[2],
                                     rzx * n[0] + rzy * n[1] + rzz * n[2]]
        # Rotate the displacement field.
        if adjust_displacement_field:
            (disp_x, disp_y, disp_z) = self._get_displacement_field_values()
            for timestep_index in xrange(len(self.timesteps)):
                if node_indices == 'all':
                    new_disp_x = [
                        rxx * x + rxy * y + rxz * z
                        for x, y, z in itertools.izip(disp_x[timestep_index],
                                                      disp_y[timestep_index],
                                                      disp_z[timestep_index])]
                    new_disp_y = [
                        ryx * x + ryy * y + ryz * z
                        for x, y, z in itertools.izip(disp_x[timestep_index],
                                                      disp_y[timestep_index],
                                                      disp_z[timestep_index])]
                    new_disp_z = [
                        rzx * x + rzy * y + rzz * z
                        for x, y, z in itertools.izip(disp_x[timestep_index],
                                                      disp_y[timestep_index],
                                                      disp_z[timestep_index])]
                    disp_x[timestep_index] = new_disp_x
                    disp_y[timestep_index] = new_disp_y
                    disp_z[timestep_index] = new_disp_z
                else:
                    for index in node_indices:
                        x = disp_x[timestep_index][index]
                        y = disp_y[timestep_index][index]
                        z = disp_z[timestep_index][index]
                        (x, y, z) = [rxx * x + rxy * y + rxz * z,
                                     ryx * x + ryy * y + ryz * z,
                                     rzx * x + rzy * y + rzz * z]
                        disp_x[timestep_index][index] = x
                        disp_y[timestep_index][index] = y
                        disp_z[timestep_index][index] = z

    def _translate_nodes(self, offset, node_indices='all'):
        """Translate nodes by the given offset."""
        (dx, dy, dz) = [float(x) for x in offset]
        if node_indices == 'all':
            self.nodes = [[x + dx, y + dy, z + dz]
                          for x, y, z in self.nodes]
        else:
            for index in node_indices:
                self.nodes[index][0] += dx
                self.nodes[index][1] += dy
                self.nodes[index][2] += dz

    def _scale_nodes(self,
                     scale_factor,
                     node_indices='all',
                     adjust_displacement_field='auto'):
        """Scale nodes in the list by the given scale factor."""
        scale_factor = float(scale_factor)
        if adjust_displacement_field == 'auto':
            adjust_displacement_field = self.displacement_field_exists()
        # Scale the nodal coordinates.
        if node_indices == 'all':
            self.nodes = [[x * scale_factor for x in n]
                          for n in self.nodes]
        else:
            for index in node_indices:
                self.nodes[index] = [x * scale_factor
                                     for x in self.nodes[index]]
        # Scale the displacement field.
        if adjust_displacement_field:
            for all_values in self._get_displacement_field_values():
                for values in all_values:
                    if node_indices == 'all':
                        values[:] = [x * scale_factor for x in values]
                    else:
                        for index in node_indices:
                            values[index] *= scale_factor

    def rotate_geometry(self,
                        axis,
                        angle_in_degrees,
                        adjust_displacement_field='auto'):
        """
        Rotate the model about an axis by the given angle.

        The rotation axis includes the origin and points in the direction of
        the 'axis' parameter.

        Example:
        >>> model.rotate_geometry([1, 0, 0], 90)

        """
        self._rotate_nodes(axis,
                           angle_in_degrees,
                           adjust_displacement_field=adjust_displacement_field)

    def translate_geometry(self, offset):
        """
        Translate the model by the given offset.

        Example:
        >>> model.translate_geometry([1, 2, 3])

        """
        self._translate_nodes(offset)

    def scale_geometry(self, scale_factor, adjust_displacement_field='auto'):
        """
        Scale the model by the given factor.

        By default, if it exists, the displacement field will also be scaled
        accordingly.

        Example:
        >>> model.scale_geometry(0.0254)

        """
        self._scale_nodes(scale_factor,
                          adjust_displacement_field=adjust_displacement_field)

    def _ensure_no_shared_nodes(self, element_block_ids):
        """
        Ensure no nodes are shared outside the given element blocks.

        If nodes are shared between the given element blocks and all
        other element block then an error message is output.

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        affected_nodes = self.get_nodes_in_element_block(element_block_ids)
        nodes_in_other_blocks = self.get_nodes_in_element_block(
            list(set(self.get_element_block_ids()) - set(element_block_ids)))
        shared_nodes = sorted(set.intersection(set(affected_nodes),
                                               set(nodes_in_other_blocks)))
        if shared_nodes:
            max_nodes_to_display = 20
            node_list = ', '.join(
                [str(x) for x in shared_nodes[:max_nodes_to_display]])
            if len(shared_nodes) > max_nodes_to_display:
                node_list += ', ...'
            self._error('Unable to operate on merged nodes.',
                        'You are attempting to operate on some element blocks '
                        'while keeping others unaffected.  Because some nodes '
                        'are shared between the two groups, this cannot be '
                        'done.  Use unmerge_element_blocks() to unmerge the '
                        'blocks if that is desired.\n'
                        '\n'
                        'There are %d shared nodes: %s' %
                        (len(shared_nodes), node_list))

    def translate_element_blocks(self,
                                 element_block_ids,
                                 offset,
                                 check_for_merged_nodes=True):
        """
        Translate the specified element blocks by the given offset.

        Example:
        >>> model.translate_element_blocks(1, [1.0, 2.0, 3.0])

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        if check_for_merged_nodes:
            self._ensure_no_shared_nodes(element_block_ids)
        affected_nodes = self.get_nodes_in_element_block(element_block_ids)
        self._translate_nodes(offset, affected_nodes)

    def reflect_element_blocks(self,
                               element_block_ids,
                               point,
                               normal,
                               check_for_merged_nodes=True,
                               adjust_displacement_field='auto'):
        """
        Reflect the specified element blocks about the given plane.

        Since an element becomes inverted when it is reflected across a plane,
        this function will also uninvert the elements.

        Example:
        >>> model.reflect_element_blocks(1, [0, 0, 0], [1, 0, 0])

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        if adjust_displacement_field == 'auto':
            adjust_displacement_field = self.displacement_field_exists()
        point = [float(x) for x in point]
        scale = math.sqrt(sum(x * x for x in normal))
        normal = [float(x) / scale for x in normal]
        if check_for_merged_nodes:
            self._ensure_no_shared_nodes(element_block_ids)
        affected_nodes = self.get_nodes_in_element_block(element_block_ids)
        for index in affected_nodes:
            distance = sum((a - p) * n
                           for a, p, n
                           in zip(self.nodes[index], point, normal))
            self.nodes[index] = [a - 2 * distance * n
                                 for a, n
                                 in zip(self.nodes[index], normal)]
        # adjust displacement fields
        if adjust_displacement_field:
            displacement_fields = self._get_displacement_field_values()
            for timestep_index in xrange(len(self.timesteps)):
                for index in affected_nodes:
                    dx = displacement_fields[0][timestep_index][index]
                    dy = displacement_fields[1][timestep_index][index]
                    dz = displacement_fields[2][timestep_index][index]
                    distance = dx * normal[0] + dy * normal[1] + dz * normal[2]
                    dx -= 2 * distance * normal[0]
                    dy -= 2 * distance * normal[1]
                    dz -= 2 * distance * normal[2]
                    displacement_fields[0][timestep_index][index] = dx
                    displacement_fields[1][timestep_index][index] = dy
                    displacement_fields[2][timestep_index][index] = dz
        # uninvert elements
        self._invert_element_blocks(element_block_ids)

    def scale_element_blocks(self,
                             element_block_ids,
                             scale_factor,
                             check_for_merged_nodes=True,
                             adjust_displacement_field='auto'):
        """
        Scale all nodes in the given element blocks by the given amount.

        By default, if a displacement field exists, this will also scale the
        displacement field.

        Example:
        >>> model.scale_element_blocks(1, 0.0254)

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        if adjust_displacement_field == 'auto':
            adjust_displacement_field = self.displacement_field_exists()
        if check_for_merged_nodes:
            self._ensure_no_shared_nodes(element_block_ids)
        affected_nodes = self.get_nodes_in_element_block(element_block_ids)
        self._scale_nodes(scale_factor,
                          affected_nodes,
                          adjust_displacement_field=adjust_displacement_field)

    def rotate_element_blocks(self,
                              element_block_ids,
                              axis,
                              angle_in_degrees,
                              check_for_merged_nodes=True,
                              adjust_displacement_field='auto'):
        """
        Rotate all nodes in the given element blocks by the given amount.

        By default, if a displacement field exists, this will also rotate the
        displacement field.

        The rotation axis includes the origin and points in the direction of
        the 'axis' parameter.

        Example:
        >>> model.rotate_element_blocks(1, [1, 0, 0], 90)

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        if adjust_displacement_field == 'auto':
            adjust_displacement_field = self.displacement_field_exists()
        if check_for_merged_nodes:
            self._ensure_no_shared_nodes(element_block_ids)
        affected_nodes = self.get_nodes_in_element_block(element_block_ids)
        self._rotate_nodes(axis,
                           angle_in_degrees,
                           affected_nodes,
                           adjust_displacement_field=adjust_displacement_field)

    def displace_element_blocks(self,
                                element_block_ids,
                                offset,
                                timesteps='all',
                                check_for_merged_nodes=True):
        """
        Displace all nodes in the given element blocks.

        This function operates on the displacement field rather than the model
        coordinate field.  To operate on the model coordinates, use
        'translate_element_blocks()'.

        Examples:
        >>> model.displace_element_blocks(1, [1.0, 2.0, 3.0])
        >>> model.displace_element_blocks('all', [1.0, 2.0, 3.0])

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        timesteps = self._format_id_list(
            timesteps,
            self.get_timesteps(),
            'timestep')
        if check_for_merged_nodes:
            self._ensure_no_shared_nodes(element_block_ids)
        timestep_indices = [self.timesteps.index(x) for x in timesteps]
        offset = [float(x) for x in offset]
        # if no timesteps exist, issue an error
        if not self.timesteps:
            self._error('No timesteps defined',
                        'A displacement field cannot exist if no timesteps '
                        'are defined.  To avoid this error, create a timestep '
                        'before calling this function.')
        # get displacement field indices
        displacement_fields = self._get_displacement_field_values()
        # get affected nodes
        node_list = self.get_nodes_in_element_block(element_block_ids)
        # translate nodes
        for timestep_index in timestep_indices:
            for dimension in xrange(3):
                values = displacement_fields[dimension][timestep_index]
                for index in node_list:
                    values[index] += offset[dimension]

    def _rename_field_on_entity(self,
                                field_type,
                                field_name,
                                new_field_name,
                                field_name_list_function,
                                entity_type,
                                entity_list,
                                entity_list_function,
                                entity_objects):
        """Rename a field defined on an entity."""
        entity_list = self._format_id_list(
            entity_list,
            entity_list_function(),
            entity_type)
        [field_name] = self._format_id_list(
            [field_name],
            field_name_list_function(entity_list),
            field_type,
            single=True)
        for id_ in entity_list:
            fields = entity_objects[id_][-1]
            if field_name not in fields:
                self._warning(field_type + ' not defined.',
                              'The given %s "%s" is not '
                              'defined on %s %s.  It cannot '
                              'be renamed.'
                              % (field_type,
                                 field_name,
                                 entity_type,
                                 str(id_)))
                continue
            if field_name == new_field_name:
                continue
            if new_field_name in fields:
                self._exists_on_entity_warning(
                    new_field_name,
                    field_type,
                    id_,
                    entity_type)
            fields[new_field_name] = fields[field_name]
            del fields[field_name]

    def rename_element_field(self,
                             element_field_name,
                             new_element_field_name,
                             element_block_ids='all'):
        """
        Rename an element field.

        Example:
        >>> model.rename_element_field('p', 'pressure')

        """
        self._rename_field_on_entity('element field',
                                     element_field_name,
                                     new_element_field_name,
                                     self.get_element_field_names,
                                     'element block',
                                     element_block_ids,
                                     self.get_element_block_ids,
                                     self.element_blocks)

    def rename_side_set_field(self,
                              side_set_field_name,
                              new_side_set_field_name,
                              side_set_ids='all'):
        """
        Rename a side set field.

        Example:
        >>> model.rename_side_set_field('cp', 'contact_pressure')

        """
        self._rename_field_on_entity('side set field',
                                     side_set_field_name,
                                     new_side_set_field_name,
                                     self.get_side_set_field_names,
                                     'side set',
                                     side_set_ids,
                                     self.get_side_set_ids,
                                     self.side_sets)

    def get_side_set_name(self, side_set_id):
        """Return the name of a side set."""
        [side_set_id] = self._format_side_set_id_list([side_set_id],
                                                      single=True)
        return self.side_sets[side_set_id][0]

    def get_all_side_set_names(self):
        """Return a list of all side set names."""
        return sorted(x[0]
                      for x in self.side_sets.values()
                      if x[0])

    def get_side_set_members(self, side_set_id):
        """Return the members of a side set."""
        [side_set_id] = self._format_side_set_id_list([side_set_id],
                                                      single=True)
        return self.side_sets[side_set_id][1]

    def _get_side_set_fields(self, side_set_id):
        """Return the dictionary of side set fields."""
        [side_set_id] = self._format_side_set_id_list([side_set_id],
                                                      single=True)
        return self.side_sets[side_set_id][-1]

    def rename_node_set_field(self,
                              node_set_field_name,
                              new_node_set_field_name,
                              node_set_ids='all'):
        """
        Rename a node set field.

        Example:
        >>> model.rename_node_set_field('cp', 'contact_pressure')

        """
        self._rename_field_on_entity('node set field',
                                     node_set_field_name,
                                     new_node_set_field_name,
                                     self.get_node_set_field_names,
                                     'node set',
                                     node_set_ids,
                                     self.get_node_set_ids,
                                     self.node_sets)

    def _rename_entity(self,
                       entity_type,
                       entity_name,
                       new_entity_name,
                       entity_name_list_function,
                       entity_objects):
        """Rename an entity."""
        [entity_name] = self._format_id_list(
            [entity_name],
            entity_name_list_function(),
            entity_type,
            single=True)
        if new_entity_name == entity_name:
            return
        if new_entity_name in entity_name_list_function():
            self._exists_warning(new_entity_name, entity_type)
        entity_objects[new_entity_name] = entity_objects[entity_name]
        del entity_objects[entity_name]

    def rename_node_field(self, node_field_name, new_node_field_name):
        """
        Rename a node field.

        Example:
        >>> model.rename_node_field('temp', 'temperature')

        """
        self._rename_entity('node field',
                            node_field_name,
                            new_node_field_name,
                            self.get_node_field_names,
                            self.node_fields)

    def rename_global_variable(self,
                               global_variable_name,
                               new_global_variable_name):
        """
        Rename a global variable.

        Example:
        >>> model.rename_global_variable('ke', 'kinetic_energy')

        """
        self._rename_entity('global variable',
                            global_variable_name,
                            new_global_variable_name,
                            self.get_global_variable_names,
                            self.global_variables)

    def rename_element_block(self, element_block_id, new_element_block_id):
        """
        Change an element block id or name.

        This function can be used to change either the element block id or
        name.  If 'new_element_block_id' is an integer, it will change the id.
        If it is a string, it will change the name.

        Example:
        >>> model.rename_element_block(1, 100)
        >>> model.rename_element_block(1, 'block_1')

        """
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        # if we're just changing the name
        if type(new_element_block_id) is str:
            # if the same name already, just exit
            if (self.element_blocks[element_block_id][0] ==
                    new_element_block_id):
                return
            # if the name already exists, issue a warning
            if self.element_block_exists(new_element_block_id):
                self._exists_warning('"' + new_element_block_id + '"',
                                     'element block')
            # rename it
            self.element_blocks[
                element_block_id][0] = new_element_block_id
            return
        assert type(new_element_block_id) is int
        # rename the block
        self._rename_entity('element block',
                            element_block_id,
                            new_element_block_id,
                            self.get_element_block_ids,
                            self.element_blocks)
        # adjust side sets
        for side_set_id in self.get_side_set_ids():
            members = self.get_side_set_members(side_set_id)
            new_members = []
            for member in members:
                if member[0] == element_block_id:
                    member = list(member)
                    member[0] = new_element_block_id
                    member = tuple(member)
                new_members.append(member)
            members[:] = new_members

    def rename_node_set(self, node_set_id, new_node_set_id):
        """
        Change a node set id or name.

        This function can be used to change either the node set id or name.
        If 'new_node_set_id' is an integer, it will change the id. If it is a
        string, it will change the name.

        Example:
        >>> model.rename_node_set(1, 100)
        >>> model.rename_node_set(1, 'node_group_1')

        """
        [node_set_id] = self._format_node_set_id_list(
            [node_set_id],
            single=True)
        # if we're just changing the name:
        if type(new_node_set_id) is str:
            # if the same name already, just exit
            if self.get_node_set_name(node_set_id) == new_node_set_id:
                return
            # if the name already exists, issue a warning
            if self.node_set_exists(new_node_set_id):
                self._exists_warning('"' + new_node_set_id + '"',
                                     'node set')
            # rename it
            self.node_sets[node_set_id][0] = new_node_set_id
            return
        # rename it
        self._rename_entity('node set',
                            node_set_id,
                            new_node_set_id,
                            self.get_node_set_ids,
                            self.node_sets)

    def rename_side_set(self, side_set_id, new_side_set_id):
        """
        Change a side set id or name.

        This function can be used to change either the side set id or name.
        If 'new_side_set_id' is an integer, it will change the id. If it is a
        string, it will change the name.

        Example:
        >>> model.rename_side_set(1, 100)
        >>> model.rename_side_set(1, 'surface_1')

        """
        [side_set_id] = self._format_side_set_id_list(
            [side_set_id],
            single=True)
        # if we're just changing the name:
        if type(new_side_set_id) is str:
            # if the same name already, just exit
            if self.get_side_set_name(side_set_id) == new_side_set_id:
                return
            # if the name already exists, issue a warning
            if self.side_set_exists(new_side_set_id):
                self._exists_warning('"' + new_side_set_id + '"',
                                     'side set')
            # rename it
            self.side_sets[side_set_id][0] = new_side_set_id
            return
        self._rename_entity('side set',
                            side_set_id,
                            new_side_set_id,
                            self.get_side_set_ids,
                            self.side_sets)

    def _empty_field_warning(self):
        """Issue a warning if no timesteps exist."""
        assert len(self.get_timesteps()) == 0
        self._warning('Creating an empty field',
                      'No timesteps are defined.  Because fields are only '
                      'defined at each timestep, creating a field when no '
                      'timesteps are defined will allow let the field to '
                      'have any values.  By default, empty fields are not '
                      'exported.\n'
                      '\n'
                      'To avoid this warning, create a timestep before '
                      'creating a field.')

    def create_element_field(self,
                             element_field_name,
                             element_block_ids='all',
                             value='auto'):
        """
        Create an element field on one or more element blocks.

        A default value can be passed.  If no default value is given, 0.0 will
        be used if the element field appears to be a displacement field.
        Otherwise, NaN will be used.

        Example:
        >>> model.create_element_field('eqps', value=0.0)

        """
        # warn if no timesteps are present
        if len(self.get_timesteps()) == 0:
            self._empty_field_warning()
        # get list of block ids
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        # get the value to assign
        if value == 'auto':
            value = self._get_default_field_value(element_field_name)
        # for each block
        for element_block_id in element_block_ids:
            # if it exists, no need to do anything
            if self.element_field_exists(element_field_name,
                                         element_block_id):
                self._exists_on_entity_warning(
                    element_field_name,
                    'element field',
                    element_block_id,
                    'element block')
            # get the element block info
            element_count = self.get_element_count(element_block_id)
            # create an empty field
            values = []
            for _ in xrange(len(self.timesteps)):
                values.append([value] * element_count)
            # add the field
            fields = self._get_element_block_fields(element_block_id)
            fields[element_field_name] = values

    def create_node_field(self, node_field_name, value='auto'):
        """
        Create a node field and assign it a default value.

        A default value can be passed.  If no default value is given, 0.0 will
        be used if the element field appears to be a displacement field.
        Otherwise, NaN will be used.

        Example:
        >>> model.create_node_field('temperature', 298.15)

        """
        # issue warning if no timesteps exist
        if len(self.get_timesteps()) == 0:
            self._empty_field_warning()
        # if it exists, no need to do anything
        if self.node_field_exists(node_field_name):
            self._exists_warning(node_field_name, 'node field')
            return
        # get the value
        if value == 'auto':
            value = self._get_default_field_value(node_field_name)
        # create the new field
        new_field_values = []
        for _ in xrange(len(self.timesteps)):
            new_field_values.append([value] * len(self.nodes))
        self.node_fields[node_field_name] = new_field_values

    def create_node_set_field(self,
                              node_set_field_name,
                              node_set_ids='all',
                              value='auto'):
        """
        Create a node set field on the given node sets.

        A default value can be passed.  If no default value is given, 0.0 will
        be used if the element field appears to be a displacement field.
        Otherwise, NaN will be used.

        Example:
        >>> model.create_node_set_field('temperature', 13, 298.15)

        """
        node_set_ids = self._format_node_set_id_list(node_set_ids)
        if value == 'auto':
            value = self._get_default_field_value(node_set_field_name)
        for node_set_id in node_set_ids:
            members = self.get_node_set_members(node_set_id)
            fields = self._get_node_set_fields(node_set_id)
            member_count = len(members)
            # if it exists, warn about overwriting
            if self.node_set_field_exists(node_set_field_name, [node_set_id]):
                self._exists_on_entity_warning(
                    node_set_field_name,
                    'node set field',
                    node_set_id,
                    'node set')
            # create the new field
            new_field_values = []
            for _ in xrange(len(self.timesteps)):
                new_field_values.append([value] * member_count)
            fields[node_set_field_name] = new_field_values

    def create_side_set_field(self,
                              side_set_field_name,
                              side_set_ids='all',
                              value='auto'):
        """
        Create a side set field.

        By default the field will be created on all side sets.  To create a
        side set field on a particular field, pass in 'side_set_ids'.

        To set the value of the field, pass in 'value'.  By default this is
        0 for displacement fields and NaN for all other fields.

        Example:
        >>> model.create_side_set_field('temperature', 13, 298.15)

        """
        side_set_ids = self._format_side_set_id_list(side_set_ids)
        if value == 'auto':
            value = self._get_default_field_value(side_set_field_name)
        for side_set_id in side_set_ids:
            members = self.get_side_set_members(side_set_id)
            fields = self._get_side_set_fields(side_set_id)
            member_count = len(members)
            # if it exists, warn about overwriting
            if self.side_set_field_exists(side_set_field_name, side_set_id):
                self._exists_on_entity_warning(
                    side_set_field_name,
                    'side set field',
                    side_set_id,
                    'side set')
            # create the new field
            new_field_values = []
            for _ in xrange(len(self.timesteps)):
                new_field_values.append([value] * member_count)
            fields[side_set_field_name] = new_field_values

    def create_node_set_from_side_set(self, node_set_id, side_set_id):
        """
        Create a node set from a side set.

        The created node set will contain all nodes used by faces in the given
        side set.

        Example:
        >>> model.create_node_set_from_side_set(1, 1)

        """
        node_set_members = self.get_nodes_in_side_set(side_set_id)
        self.create_node_set(node_set_id, node_set_members)

    @staticmethod
    def _remove_duplicates(item_list, other_list=None, preserve_order=False):
        """
        Return a list with duplicates removed.

        If 'other_list' is used, this will also remove items which appear in
        that list.

        To preserve the relative order of items, set 'preserve_order=True'.

        Example:
        >>> _remove_duplicates([1, 3, 2, 3])
        [1, 3, 2]
        >>> _remove_duplicates([1, 3, 2, 3], [3])
        [1, 2]

        """
        if not preserve_order:
            if not other_list:
                return list(set(item_list))
            return list(set(item_list) - set(other_list))
        if other_list:
            appeared = set(other_list)
        else:
            appeared = set()
        unique_list = []
        for item in item_list:
            if item not in appeared:
                unique_list.append(item)
                appeared.add(item)
        return unique_list

    def create_node_set(self, node_set_id, node_set_members=None):
        """
        Create a node side from the list of node indices.

        Node sets are unnamed when created.  To name them, use the
        'rename_node_set' function.

        Example:
        >>> model.create_node_set(1, [0, 1, 2, 3])

        """
        # if it exists, warn that it will be overwritten
        if self.node_set_exists(node_set_id):
            self._exists_warning(node_set_id, 'node set')
        # reformat members if necessary
        if not node_set_members:
            node_set_members = []
        # do not allow duplicates
        new_members = self._remove_duplicates(
            node_set_members,
            preserve_order=True)
        if len(new_members) != len(node_set_members):
            self._warning('Duplicate nodes in set.',
                          'The node set member list contains multiple copies '
                          'of some nodes.  These nodes will only be added '
                          'once.')
        # add the new set
        self.node_sets[node_set_id] = ['', new_members, {}]

    def add_nodes_to_node_set(self, node_set_id, new_node_set_members):
        """
        Add nodes to an existing node set.

        Example:
        >>> model.add_nodes_to_node_set(1, [4, 5, 6, 7])

        """
        if not self.node_set_exists(node_set_id):
            self._warning('Node set does not exist.'
                          'The specified node set "%s% does not exist.  It '
                          'will be created.' % node_set_id)
            self.create_node_set(node_set_id, new_node_set_members)
            return
        members = self.get_node_set_members(node_set_id)
        fields = self._get_node_set_fields(node_set_id)
        # Remove duplicate nodes.
        new_nodes = self._remove_duplicates(
            new_node_set_members,
            other_list=members,
            preserve_order=True)
        if len(new_nodes) != len(new_node_set_members):
            self._warning('Duplicates nodes in set',
                          'The node set already contains some of the nodes '
                          'in the list.  These members will not be '
                          'duplicated.')
        # Add the members.
        members.extend(new_nodes)
        # Add new values to each field.
        for name, all_values in fields.items():
            value = self._get_default_field_value(name)
            for values in all_values:
                values.extend([value] * len(new_nodes))

    def _transform_expression(self, expression):
        """
        Transform the expression into one that evaluates to 0 when true.

        For example, 'x=y' becomes '(x)-(y)'.

        """
        # create transformation dictionary by precedence order
        transforms = []
        transforms.append(('||', 'min(abs(L), abs(R))'))
        transforms.append(('&&', 'max(abs(L), abs(R))'))
        transforms.append(('>=', '((L) - (R)) - abs((L) - (R))'))
        transforms.append(('>', '((L) - (R)) - abs((L) - (R))'))
        transforms.append(('<=', '((R) - (L)) - abs((R) - (L))'))
        transforms.append(('<', '((R) - (L)) - abs((R) - (L))'))
        transforms.append(('==', 'abs((L) - (R))'))
        transforms.append(('=', 'abs((L) - (R))'))
        # replace occurrences of each transform
        for separator, transform in transforms:
            while separator in expression:
                # ensure parenthesis count is identical
                if expression.count('(') != expression.count(')'):
                    self._error('Unbalances parenthesis.',
                                'We cannot transform the given expression '
                                'becuase it contains unbalanced '
                                'parenthesis:\n%s' % expression)
                # get parenthesis depth on each character
                next_depth = 0
                depth = []
                for c in expression:
                    if c == '(':
                        depth.append(next_depth)
                        next_depth += 1
                    elif c == ')':
                        next_depth -= 1
                        depth.append(next_depth)
                    else:
                        depth.append(next_depth)
                # find extents of the inner expression
                separator_index = expression.index(separator)
                left_index = separator_index
                right_index = separator_index + len(separator)
                while left_index > 0 and (
                        depth[left_index - 1] > depth[separator_index] or
                        (depth[left_index - 1] == depth[separator_index] and
                         expression[left_index - 1] != ',')):
                    left_index -= 1
                while right_index < len(expression) and (
                        depth[right_index] > depth[separator_index] or
                        (depth[right_index] == depth[separator_index] and
                         expression[right_index] != ',')):
                    right_index += 1
                new_expression = expression[:left_index]
                left_expression = expression[
                    left_index:separator_index].strip()
                right_expression = expression[
                    separator_index + len(separator):right_index].strip()
                for c in transform:
                    if c == 'L':
                        new_expression += left_expression
                    elif c == 'R':
                        new_expression += right_expression
                    else:
                        new_expression += c
                new_expression += expression[right_index:]
                expression = new_expression
        # return the result
        return expression

    def create_side_set_from_expression(self,
                                        side_set_id,
                                        expression,
                                        element_block_ids='all',
                                        tolerance='auto',
                                        timesteps='last_if_any',
                                        zero_member_warning=True):
        """
        Create a side set from faces which satisfy an expression.

        Only external element faces of the given element blocks are checked.

        For example, if the model had a symmetry plane at 'Y = 0', the
        expression 'Y == 0' would select all element faces on this plane.

        Example:
        >>> model.create_side_set_from_expression(1, 'Y == 0')

        """
        timesteps = self._format_id_list(
            timesteps,
            self.get_timesteps(),
            'timestep')
        if len(timesteps) > 2:
            self._error('Too many timesteps specified.',
                        'We were expecting zero or one timesteps but instead '
                        'found %d specified.' % len(timesteps))
        if tolerance == 'auto':
            tolerance = self.get_length_scale() * 1e-6
        if self.side_set_exists(side_set_id):
            self._exists_warning(side_set_id, 'side set')
        external_faces_by_block = self._order_element_faces_by_block(
            self._get_external_element_faces(element_block_ids))
        members = []
        # create list of variable names and modify them in the expression
        variable_names = set(['X', 'Y', 'Z'])
        if timesteps:
            timestep_index = self._get_internal_timestep_index(timesteps[0])
            variable_names.update(['time'])
            variable_names.update(self.get_global_variable_names())
            variable_names.update(self.get_node_field_names())
        expression = self._transform_expression(expression)
        eval_expression = self._transform_eval_expression(expression,
                                                          variable_names)
        var = dict()
        function = eval('lambda var: ' + eval_expression)
        if timesteps:
            var['time'] = timesteps[0]
            for name, values in self.global_variables.items():
                var[name] = values[timestep_index]
        try:
            for id_, external_faces in external_faces_by_block.items():
                connectivity = self.get_connectivity(id_)
                face_mapping = self._get_face_mapping(
                    self._get_element_type(id_))
                nodes_per_element = self.get_nodes_per_element(id_)
                for element_index, face_index in external_faces:
                    node_indices = [
                        connectivity[element_index * nodes_per_element + x]
                        for x in face_mapping[face_index][1]]
                    all_passed = True
                    for node_index in node_indices:
                        # set coordinate values
                        var['X'] = self.nodes[node_index][0]
                        var['Y'] = self.nodes[node_index][1]
                        var['Z'] = self.nodes[node_index][2]
                        # set node field values
                        if timesteps:
                            for name, values in self.node_fields.items():
                                var[name] = values[timestep_index][node_index]
                        value = float(function(var))
                        if not abs(value) <= tolerance:
                            all_passed = False
                            break
                    if all_passed:
                        members.append(
                            (id_, element_index, face_index))
        except (SyntaxError, NameError):
            self._error_evaluating_expression(
                "%s" % (eval_expression),
                var)
        # create the side set
        self.create_side_set(side_set_id, members)
        if zero_member_warning and not members:
            self._warning('Empty side set.',
                          'No external element faces satisfied the given '
                          'expression.  An empty side set was created.')

    def create_side_set(self, side_set_id, side_set_members=None):
        """
        Create a side set from the given element faces.

        If the side set already exists, the faces will be added.

        side_set_members should be a list of tuples of the form:
        * '(element_block_id, local_element_index, element_side_index)'

        Side sets are unnamed when created.  To name them, use the
        'rename_side_set' function.

        Example:
        >>> model.create_side_set(1, [(1, 0, 1), (1, 0, 2)])

        """
        if not side_set_members:
            side_set_members = []
        # ensure set set id doesn't exist
        if self.side_set_exists(side_set_id):
            self._exists_warning(side_set_id, 'side set')
        unique_members = self._remove_duplicates(
            side_set_members,
            preserve_order=True)
        if len(unique_members) != len(side_set_members):
            self._warning('Duplicate faces in set.',
                          'The face set member list contains multiple copies '
                          'of some faces.  These faces will only be added '
                          'once.')
        self.side_sets[side_set_id] = ['', unique_members, {}]

    def add_faces_to_side_set(self, side_set_id, new_side_set_members):
        """
        Add the given faces to the side set.

        The parameter 'new_side_set_members' should be a list of tuples of the
        form:
        * '(element_block_id, local_element_index, element_side_index)'

        Example:
        >>> model.add_faces_to_side_set(1, [(2, 0, 3), (3, 0, 4)])

        """
        if not self.side_set_exists(side_set_id):
            self._warning('Side set does not exist.'
                          'The specified side set "%s% does not exist.  It '
                          'will be created.' % side_set_id)
            self.create_side_set(side_set_id, new_side_set_members)
            return
        members = self.get_side_set_members(side_set_id)
        fields = self._get_side_set_fields(side_set_id)
        # remove duplicate faces
        new_members = self._remove_duplicates(
            new_side_set_members,
            other_list=members,
            preserve_order=True)
        if len(new_members) != len(new_side_set_members):
            self._warning('Duplicates faces in set',
                          'The node set already contains some nodes of the '
                          'faces in the list to add.  These members will not '
                          'be duplicated.')
        # add the members
        members.extend(new_members)
        # add new values to each field
        for name, all_values in fields.items():
            value = self._get_default_field_value(name)
            for values in all_values:
                values.extend([value] * len(new_members))

    def delete_node_field(self, node_field_names):
        """
        Delete one or more node fields.

        Examples:
        >>> model.delete_node_field('temperature')
        >>> model.delete_node_field('all')
        >>> model.delete_node_field('disp_*')

        """
        node_field_names = self._format_id_list(
            node_field_names,
            self.get_node_field_names(),
            'node field')
        for node_field_name in node_field_names:
            del self.node_fields[node_field_name]

    def delete_timestep(self, timesteps):
        """
        Delete one or more timesteps.

        Because fields are defined on each timestep, this also deletes field
        values for the corresponding timestep.

        Examples:
        >>> model.delete_timestep(0.0)
        >>> model.delete_timestep('all')

        """
        timesteps = self._format_id_list(
            timesteps,
            self.get_timesteps(),
            'timestep')
        # get list of deleted indices in descending order
        deleted_indices = reversed(sorted([self.timesteps.index(x)
                                           for x in timesteps]))
        for index in deleted_indices:
            # process element fields
            for element_block_id in self.get_element_block_ids():
                fields = self._get_element_block_fields(element_block_id)
                for values in fields.values():
                    del values[index]
            # process node fields
            for values in self.node_fields.values():
                del values[index]
            # process side set fields
            for side_set_id in self.get_side_set_ids():
                fields = self._get_side_set_fields(side_set_id)
                for values in fields.values():
                    del values[index]
            # process node set fields
            for node_set_id in self.get_node_set_ids():
                fields = self._get_node_set_fields(node_set_id)
                for values in fields.values():
                    del values[index]
            # process global variables
            for values in self.global_variables.values():
                del values[index]
            # delete the timestep
            del self.timesteps[index]

    def delete_element_block(self,
                             element_block_ids,
                             delete_orphaned_nodes=True):
        """
        Delete one or more element blocks.

        This will also delete any references to elements in that block in
        side sets.

        By default, this will delete any nodes that become unused as a result
        of deleting the element blocks.  To prevent this, set
        'delete_orphaned_nodes=False'.

        Examples:
        >>> model.delete_element_block(1)
        >>> model.delete_element_block([1, 3, 4])

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        # if we're not deleting anything, skip some computations.
        if not element_block_ids:
            return
        # find unreferenced nodes
        if delete_orphaned_nodes:
            unreferenced_nodes = self._get_unreferenced_nodes()
        # delete the element blocks and associated data
        for element_block_id in element_block_ids:
            # delete the element block itself
            del self.element_blocks[element_block_id]
            # delete faces of that element block from side set
            for side_set_id in self.get_side_set_ids():
                members = self.get_side_set_members(side_set_id)
                fields = self._get_side_set_fields(side_set_id)
                # find indices to delete
                deleted_indices = []
                for index, (id_, _, _) in enumerate(members):
                    if id_ == element_block_id:
                        deleted_indices.append(index)
                # delete them from members
                for index in reversed(deleted_indices):
                    del members[index]
                # delete them from the fields
                for all_values in fields.values():
                    for values in all_values:
                        for index in reversed(deleted_indices):
                            del values[index]
        # now find the new unreferenced nodes
        if delete_orphaned_nodes:
            new_unreferenced_nodes = self._get_unreferenced_nodes()
            nodes_to_delete = sorted(set(new_unreferenced_nodes) -
                                     set(unreferenced_nodes))
            if nodes_to_delete:
                self._delete_nodes(nodes_to_delete)

    def delete_element_field(self,
                             element_field_names,
                             element_block_ids='all'):
        """
        Delete one or more element fields.

        Examples:
        >>> model.delete_element_field('eqps')
        >>> model.delete_element_field('all')

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        element_field_names = self._format_id_list(
            element_field_names,
            self.get_element_field_names(),
            'element field')
        # for each field
        for element_field_name in element_field_names:
            any_deleted = False
            # for each element block
            for element_block_id in element_block_ids:
                fields = self._get_element_block_fields(element_block_id)
                if element_field_name in fields:
                    any_deleted = True
                    del fields[element_field_name]
            if not any_deleted:
                self._warning('Element field not defined.',
                              'The element field "%s" was not defined on any '
                              'of the given element blocks.  It cannot be '
                              'deleted.' % element_field_name)

    def delete_node_set_field(self, node_set_field_names, node_set_ids='all'):
        """
        Delete one or more node set fields.

        Examples:
        >>> model.delete_node_set_field('contact_pressure', 1)

        """
        node_set_ids = self._format_node_set_id_list(node_set_ids)
        node_set_field_names = self._format_id_list(
            node_set_field_names,
            self.get_node_set_field_names(),
            'node set field')
        # for each field
        for name in node_set_field_names:
            any_deleted = False
            # for each node set
            for node_set_id in node_set_ids:
                fields = self._get_node_set_fields(node_set_id)
                if name in fields:
                    any_deleted = True
                    del fields[name]
            if not any_deleted:
                self._warning('Node set field not defined.',
                              'The node set field "%s" was not defined on any '
                              'of the given node sets.  It cannot be '
                              'deleted.' % name)

    def delete_side_set_field(self, side_set_field_names, side_set_ids='all'):
        """
        Delete one or more side set fields.

        Examples:
        >>> model.delete_side_set_field('contact_pressure', 1)

        """
        side_set_ids = self._format_side_set_id_list(side_set_ids)
        side_set_field_names = self._format_id_list(
            side_set_field_names,
            self.get_side_set_field_names(),
            'side set field')
        # for each field
        for name in side_set_field_names:
            any_deleted = False
            # for each node set
            for id_ in side_set_ids:
                fields = self._get_side_set_fields(id_)
                if name in fields:
                    any_deleted = True
                    del fields[name]
            if not any_deleted:
                self._warning('Side set field not defined.',
                              'The side set field "%s" was not defined on any '
                              'of the given side sets.  It cannot be '
                              'deleted.' % name)

    def global_variable_exists(self, global_variable_name):
        """
        Return 'True' if the given global variable exists.

        Examples:
        >>> model.global_variable_exists('timestep')

        """
        return global_variable_name in self.global_variables

    def side_set_field_exists(self, side_set_field_name, side_set_ids='all'):
        """
        Return 'True' if the side set field exists on the side sets.

        Examples:
        >>> model.side_set_field_exists('contact_pressure')

        """
        side_set_ids = self._format_side_set_id_list(side_set_ids)
        for id_ in side_set_ids:
            fields = self._get_side_set_fields(id_)
            if side_set_field_name not in fields:
                return False
        return True

    def node_set_field_exists(self, node_set_field_name, node_set_ids='all'):
        """
        Return 'True' if the node set field is defined on the node sets.

        Examples:
        >>> model.node_set_field_exists('contact_pressure')

        """
        node_set_ids = self._format_node_set_id_list(node_set_ids)
        for node_set_id in node_set_ids:
            fields = self._get_node_set_fields(node_set_id)
            if node_set_field_name not in fields:
                return False
        return True

    def element_field_exists(self,
                             element_field_name,
                             element_block_ids='all'):
        """
        Return 'True' if the element field exists on the element blocks.

        Examples:
        >>> model.element_field_exists('eqps')

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        for id_ in element_block_ids:
            fields = self.element_blocks[id_][-1]
            if element_field_name not in fields:
                return False
        return True

    def node_field_exists(self, node_field_name):
        """
        Return 'True' if the given node field exists.

        Examples:
        >>> model.node_field_exists('temperature')

        """
        return node_field_name in self.node_fields

    def element_block_exists(self, element_block_id):
        """
        Return 'True' if the given element block exists.

        Examples:
        >>> model.element_block_exists(1)
        >>> model.element_block_exists('block_name')

        """
        if type(element_block_id) is str:
            return element_block_id in self.get_all_element_block_names()
        else:
            return element_block_id in self.get_element_block_ids()

    def node_set_exists(self, node_set_id):
        """
        Return 'True' if the given node set exists.

        Examples:
        >>> model.node_set_exists(1)
        >>> model.node_set_exists('nodeset_name')

        """
        if type(node_set_id) is str:
            return node_set_id in self.get_all_node_set_names()
        else:
            return node_set_id in self.get_node_set_ids()

    def side_set_exists(self, side_set_id):
        """
        Return 'True' if the given side set exists.

        Examples:
        >>> model.side_set_exists(1)
        >>> model.node_set_exists('sideset_name')

        """
        if type(side_set_id) is str:
            return side_set_id in self.get_all_side_set_names()
        else:
            return side_set_id in self.get_side_set_ids()

    def timestep_exists(self, timestep):
        """
        Return 'True' if the given timestep exists.

        Examples:
        >>> model.timestep_exists(0.0)

        """
        return timestep in self.timesteps

    def get_node_field_values(self, node_field_name, timestep='last'):
        """
        Return the list of node field values for the given field and timestep.

        This returns the actual list of values, so any modifications to the
        list will be retained in the model.

        Examples:
        >>> model.get_node_field_values('disp_x')
        >>> model.get_node_field_values('disp_x', 0.0)

        """
        [node_field_name] = self._format_id_list(
            [node_field_name],
            self.get_node_field_names(),
            'node field',
            single=True)
        timestep_index = self._get_internal_timestep_index(timestep)
        return self.node_fields[node_field_name][timestep_index]

    def get_node_set_name(self, node_set_id):
        """Return the name of the given node set."""
        [node_set_id] = self._format_node_set_id_list(
            [node_set_id],
            single=True)
        return self.node_sets[node_set_id][0]

    def get_all_node_set_names(self):
        """Return a list of all node set names."""
        return sorted(x[0]
                      for x in self.node_sets.values()
                      if x[0])

    def get_node_set_members(self, node_set_id):
        """
        Return the list of node indices that belong to the given node set.

        Example:
        >>> model.get_node_set_members(1)

        """
        [node_set_id] = self._format_node_set_id_list(
            [node_set_id],
            single=True)
        return self.node_sets[node_set_id][1]

    def _get_face_mapping(self, element_name):
        """
        Return the mapping from face number to node number for an element.

        A list is returned with members of the form '(type, numbering)' where
        'type' is a string giving the face element type and 'numbering' is a
        list if element indices.

        """
        element_type = self._get_standard_element_type(element_name)
        if element_type not in self.FACE_MAPPING:
            self._error('Undefined face mapping.',
                        'The mapping defining which nodes belong to a given '
                        'face has not been defined for the \"%s\" element '
                        'type.  This could be because it has not been '
                        'implemented or because the element type is '
                        'invalid.  Mappings for the following element types '
                        'are defined: %s' %
                        (element_name, ', '.join(self.FACE_MAPPING.keys())))
        return self.FACE_MAPPING[element_type]

    def _get_face_mapping_from_id(self, element_block_id):
        """Return the face mapping from the element block id."""
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        element_type = self._get_element_type(element_block_id)
        element_type = self._get_standard_element_type(element_type)
        return self._get_face_mapping(element_type)

    def get_nodes_in_side_set(self, side_set_id):
        """
        Return a list of node indices which belong to the given side set.

        Example:
        >>> model.get_nodes_in_side_set(1)

        """
        [side_set_id] = self._format_side_set_id_list(
            [side_set_id],
            single=True)
        included_nodes = list()
        # store for element number in each element block
        last_element_block_id = 0
        nodes_per_element = 0
        face_mapping = None
        members = self.get_side_set_members(side_set_id)
        for block_id, element_index, face_index in members:
            if block_id != last_element_block_id or face_mapping is None:
                last_element_block_id = block_id
                nodes_per_element = self.get_nodes_per_element(block_id)
                face_mapping = self._get_face_mapping(
                    self._get_element_type(block_id))
                connectivity = self.get_connectivity(block_id)
            for local_index in face_mapping[face_index][1]:
                included_nodes.append(
                    connectivity[
                        element_index * nodes_per_element + local_index])
        # return list after deleting duplicates
        included_nodes = sorted(set(included_nodes))
        return included_nodes

    def _invert_face(self, face, face_type):
        """
        Return the connectivity of the inverted 2D face.

        Example:
        >>> model._invert_face(tuple(7, 2, 1, 3), 'quad4')

        """
        return tuple(face[x] for x in self.INVERTED_CONNECTIVITY[face_type])

    def _rotate_face(self, face, face_type):
        """
        Return a rotation of the given 2D face.

        Example:
        >>> model._invert_face(tuple(7, 2, 1, 3), 'quad4')

        """
        return tuple(face[x] for x in self.ROTATED_CONNECTIVITY[face_type])

    def _minimum_face(self, face, face_type):
        """
        Find a unique identifier for the given face.

        This finds the local face connectivity list and rotates it until the
        minimum such list is found.

        """
        original = face
        best = face
        best_count = 0
        face = self._rotate_face(face, face_type)
        rotate_count = 0
        while face != original:
            rotate_count += 1
            if face < best:
                best = face
                best_count = rotate_count
            face = self._rotate_face(face, face_type)
        return best, best_count

    def _get_internal_faces(self, element_block_ids='all', set_of_nodes=None):
        """
        Return a list of internal element faces.

        If 'element_block_ids' is specified, only faces within those elements
        will be considered.  If 'set_of_nodes' is given, only faces which
        contain one or more nodes in that set will be returned.

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        # hold internal faces
        internal_faces = dict()
        # hold faces which need matched
        faces_to_match = dict()
        # form node mask if possible
        if set_of_nodes is not None:
            node_mask = [x in set_of_nodes
                         for x in xrange(len(self.nodes))]
        # go through each element block
        for id_ in self.get_element_block_ids():
            connectivity = self.get_connectivity(id_)
            # if we have a node mask, skip this block if no nodes in the block
            # are part of the set of nodes we care about
            if set_of_nodes is not None:
                if set(connectivity).isdisjoint(set_of_nodes):
                    continue
            face_mappings = self._get_face_mapping_from_id(id_)
            element_count = self.get_element_count(id_)
            nodes_per_element = self.get_nodes_per_element(id_)
            for element_index in xrange(element_count):
                local_node = connectivity[
                    element_index * nodes_per_element:
                    (element_index + 1) * nodes_per_element]
                for face_index, (face_type, face_mapping) in enumerate(
                        face_mappings):
                    this_member = (id_, element_index, face_index)
                    this_face = tuple(local_node[x] for x in face_mapping)
                    # skip the face if no nodes are relevant
                    if set_of_nodes is not None:
                        if not [True
                                for x in this_face
                                if node_mask[x]]:
                            continue
                    # transform into lowest form
                    unique_face, unique_offset = self._minimum_face(
                        this_face, face_type)
                    # match if possible
                    match = faces_to_match.pop((face_type, unique_face), None)
                    # if no match, add it
                    if match is None:
                        # invert and transform into lowest form
                        inverted_face = self._invert_face(this_face, face_type)
                        inverted_face, inverted_offset = self._minimum_face(
                            inverted_face, face_type)
                        # add it and continue
                        faces_to_match[(face_type, inverted_face)] = (
                            this_member,
                            inverted_offset)
                        continue
                    # set first mapping
                    transform = (match[1] - unique_offset) % len(this_face)
                    internal_faces[this_member] = (match[0], transform)
                    # set opposite mapping
                    internal_faces[match[0]] = (this_member, transform)
        return internal_faces

    def _get_face_indices(self, member):
        """
        Return the indices of the given face.

        'member = (element_block_id, element_index, face_index)'

        This function is slow.

        """
        face_mapping = self._get_face_mapping(
            self._get_element_type(member[0]))
        nodes_per_element = self.get_nodes_per_element(member[0])
        connectivity = self.get_connectivity(member[0])
        return [connectivity[nodes_per_element * member[1] + x]
                for x in face_mapping[member[2]][1]]

    def convert_side_set_to_cohesive_zone(self,
                                          side_set_ids,
                                          new_element_block_id):
        """
        Convert the given side sets into a block of cohesive zone elements.

        The given side set must contain faces which are internal (i.e. shared
        by exactly two elements).  Support for the following face types is
        implemented:
        * 'quad4' creates 'hex8' cohesive elements
        * 'tri3' creates 'wedge6' cohesive elements
        * 'tri6' creates 'wedge12' cohesive elements

        Example:
        >>> model.convert_side_set_to_cohesive_zone(1, 2)

        """
        side_set_ids = self._format_side_set_id_list(side_set_ids)
        # if nothing to operate on, just return
        if not side_set_ids:
            return
        # order all faces by element block
        members_by_block = self._order_element_faces_by_block(
            itertools.chain(*[self.get_side_set_members(x)
                              for x in side_set_ids]))
        # find face shape
        face_type = set()
        element_types = set()
        for id_, members in members_by_block.items():
            element_types.add(self._get_element_type(id_))
            face_mapping = self._get_face_mapping_from_id(id_)
            numbers = set(x for _, x in members)
            face_type.update(face_mapping[x][0] for x in numbers)
        # ensure single face type
        if len(face_type) != 1:
            self._error('Multiple face types',
                        'The side set contains multiple types of faces and '
                        'cannot be converted into a cohesive zone.')
        # ensure face type is supported
        face_type = sorted(face_type)[0]
        if face_type not in self.COHESIVE_FORMULA:
            self._bug('Unsupported face type',
                      'The side set chosen has elements of face type "%s".  '
                      'At this time, only the following face types are '
                      'supported:\n* %s'
                      % (face_type,
                         '\n* '.join(sorted(self.COHESIVE_FORMULA.keys()))))
        # for each element type, find adjacent faces
        # (adjacent faces in this context are faces which share one or more
        # nodes)
        adjacent_faces = dict()
        element_types = sorted(element_types)
        for element_type in element_types:
            face_mapping = self._get_face_mapping(element_type)
            face_count = len(face_mapping)
            adjacent_face = [[] for i in xrange(face_count)]
            for face_one in xrange(face_count):
                for face_two in xrange(face_count):
                    if set(face_mapping[face_one][1]).intersection(
                            face_mapping[face_two][1]):
                        adjacent_face[face_one].append(face_two)
            adjacent_faces[element_type] = adjacent_face
        # create a set of all nodes in the side set
        nodes_on_surface = set()
        for side_set_id in side_set_ids:
            nodes_on_surface.update(self.get_nodes_in_side_set(side_set_id))
        # find all relevant internal faces
        internal_faces = self._get_internal_faces(
            list(members_by_block.keys()),
            set_of_nodes=nodes_on_surface)
        # create a set of all free faces
        free_faces = set()
        for id_, members in members_by_block.items():
            for element_index, face_index in members:
                face = (id_, element_index, face_index)
                free_faces.add(face)
                matching_face = internal_faces.get(face, None)
                if not matching_face:
                    self._error('Invalid face',
                                'One or more members of the side set are '
                                'external faces which cannot be converted '
                                'to cohesive elements.')
                free_faces.add(matching_face[0])
        # count the number of times each node is used
        node_use_count = [0] * len(self.nodes)
        for id_ in self.get_element_block_ids():
            connectivity = self.get_connectivity(id_)
            for node_index in connectivity:
                node_use_count[node_index] += 1
        # hold a list of all nodes to duplicate
        nodes_to_duplicate = []
        for node_index in nodes_on_surface:
            nodes_to_duplicate.extend([node_index] *
                                      node_use_count[node_index])
        nodes_to_duplicate.sort()
        # duplicate nodes
        first_new_node_index = len(self.nodes)
        new_indices = []
        self._duplicate_nodes(nodes_to_duplicate, new_indices)
        # for each duplicated node, find the first new node index
        next_index = [None] * first_new_node_index
        for old_index, new_index in itertools.izip(nodes_to_duplicate,
                                                   new_indices):
            if next_index[old_index] is None:
                next_index[old_index] = new_index
        # for faster lookup, create a vector to see if a given node
        # is on the cohesive zone
        node_in_zone = [x in nodes_on_surface
                        for x in xrange(first_new_node_index)]
        # for each element, unmerge nodes on the cohesive zone
        for id_ in self.get_element_block_ids():
            connectivity = self.get_connectivity(id_)
            for index, value in enumerate(connectivity):
                if node_in_zone[value]:
                    connectivity[index] = next_index[value]
                    next_index[value] += 1
        # find pairs of nodes which need merged
        node_pairs = set()
        for old_face, (new_face, transform_count) in internal_faces.items():
            # since the list is duplicated, only operate on one
            if old_face > new_face:
                continue
            # skip faces on the free surface
            if old_face in free_faces:
                continue
            # merge these nodes
            original = tuple(self._get_face_indices(old_face))
            target = tuple(self._get_face_indices(new_face))
            original = self._invert_face(original, face_type)
            for _ in xrange(transform_count):
                original = self._rotate_face(original, face_type)
            node_pairs.update((min(x, y), max(x, y))
                              for x, y in itertools.izip(original, target)
                              if x != y)
        # merge nodes
        self._merge_node_pairs(node_pairs)
        # delete nodes we created that are no longer referenced as well as the
        # original nodes which are no longer used
        nodes_to_delete = self._get_unreferenced_nodes()
        nodes_to_delete = [x
                           for x in nodes_to_delete
                           if x >= first_new_node_index]
        nodes_to_delete.extend(nodes_on_surface)
        self._delete_nodes(nodes_to_delete)
        # get the formula for converting two faces into a cohesive element
        formula = self.COHESIVE_FORMULA[face_type]
        nodes_per_element = self.NODES_PER_ELEMENT[formula[0]]
        # create cohesive zone element block connectivity by
        connectivity = []
        for face in free_faces:
            other_face, transform_count = internal_faces[face]
            if other_face < face:
                continue
            bottom = self._get_face_indices(face)
            top = self._get_face_indices(other_face)
            bottom = self._invert_face(bottom, face_type)
            for _ in xrange(transform_count):
                bottom = self._rotate_face(bottom, face_type)
            top.extend(bottom)
            # apply formula to convert two faces into an element
            this_connectivity = [top[x]
                                 for x in formula[1]]
            # add this element to the connectivity list
            connectivity.extend(this_connectivity)
        # create the actual element block
        self.create_element_block(
            new_element_block_id,
            [formula[0],
             len(connectivity) / nodes_per_element,
             nodes_per_element,
             0],
            connectivity)

    def create_averaged_element_field(self,
                                      from_element_field_names,
                                      new_element_field_name,
                                      element_block_ids='all'):
        """
        Create an element field by averaging the given field values.

        Examples:
        >>> model.create_averaged_element_field(['temp_1', 'temp_2'],
        ...                                     'temp_avg')
        >>> model.create_averaged_element_field('temp_*', 'temp_avg')

        """
        # format the arguments
        from_element_field_names = self._format_id_list(
            from_element_field_names,
            self.get_element_field_names(),
            'element field')
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        # for each timestep
        for element_block_id in element_block_ids:
            fields = self._get_element_block_fields(element_block_id)
            # see if all fields are defined on this block
            all_defined = True
            for name in from_element_field_names:
                if not self.element_field_exists(
                        name, element_block_id):
                    all_defined = False
                    break
            if not all_defined:
                self._warning('Fields not defined.',
                              'Not all of the requested element fields are '
                              'defined on element block %s.  The averaged '
                              'field will not be created.' % element_block_id)
                continue
            # create the field if it doesn't exist
            if self.element_field_exists(new_element_field_name,
                                         element_block_id):
                self._exists_on_entity_warning(
                    new_element_field_name,
                    'element field',
                    element_block_id,
                    'element block')
            else:
                self.create_element_field(new_element_field_name,
                                          element_block_id)
            # get lists to average
            field_list = []
            for name in from_element_field_names:
                field_list.append(fields[name])
            new_element_field = fields[new_element_field_name]
            # for each timestep, average the values
            field_count = len(from_element_field_names)
            for timestep_index in xrange(len(self.timesteps)):
                lists_to_average = [x[timestep_index] for x in field_list]
                averaged_field = [sum(x) / field_count
                                  for x in zip(*lists_to_average)]
                new_element_field[timestep_index] = averaged_field

    def convert_element_field_to_node_field(self,
                                            element_field_name,
                                            node_field_name='auto'):
        """
        Convert an element field to a node field by performing.

        For each node, the value of the field at that node will be the average
        of the values in every element which shares that node for elements on
        which the field is defined.

        By default, the name of the node field will be the same as the element
        field.

        Example:
        >>> model.convert_element_field_to_node_field('temperature')

        """
        [element_field_name] = self._format_id_list(
            [element_field_name],
            self.get_element_field_names(),
            'element field',
            single=True)
        # format the arguments
        if node_field_name == 'auto':
            node_field_name = element_field_name
        # create the field
        self.create_node_field(node_field_name)
        # store field indices
        node_field = self.node_fields[node_field_name]
        # store default value
        default_value = self._get_default_field_value(node_field_name)
        # process each timestep
        for timestep_index in xrange(len(self.timesteps)):
            # initialize node field
            node_field_values = [0.0] * len(self.nodes)
            node_field_elements = [0] * len(self.nodes)
            # for each element block
            for element_block_id in self.get_element_block_ids():
                # if field doesn't exist on this element, skip it
                if not self.element_field_exists(element_field_name,
                                                 element_block_id):
                    continue
                # store connectivity
                connectivity = self.get_connectivity(element_block_id)
                # for each node within each element
                fields = self._get_element_block_fields(element_block_id)
                element_field_values = fields[
                    element_field_name][timestep_index]
                nodes_per_element = self.get_nodes_per_element(
                    element_block_id)
                element_count = self.get_element_count(element_block_id)
                for element_index in xrange(element_count):
                    for node_index in xrange(nodes_per_element):
                        this_node = connectivity[
                            element_index * nodes_per_element + node_index]
                        node_field_values[this_node] += element_field_values[
                            element_index]
                        node_field_elements[this_node] += 1
            # average each value, or replace with default value
            for node_index in xrange(len(self.nodes)):
                if node_field_elements[node_index] == 0:
                    node_field_values[node_index] = default_value
                else:
                    node_field_values[node_index] /= node_field_elements[
                        node_index]
            # replace field values with generated values
            node_field[timestep_index] = node_field_values

    def convert_node_field_to_element_field(self,
                                            node_field_name,
                                            element_field_name='auto',
                                            element_block_ids='all'):
        """
        Convert a node field to an element field.

        For each element, the value of the field at that element will be the
        average of the values in each node of that element.

        Example:
        >>> model.convert_node_field_to_element_field('temperature')

        """
        [node_field_name] = self._format_id_list(
            [node_field_name],
            self.get_node_field_names(),
            'node field',
            single=True)
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        # format the arguments
        if element_field_name == 'auto':
            element_field_name = node_field_name
        # store the node field
        node_field = self.node_fields[node_field_name]
        # for each element block
        for element_block_id in element_block_ids:
            element_count = self.get_element_count(element_block_id)
            connectivity = self.get_connectivity(element_block_id)
            nodes_per_element = self.get_nodes_per_element(element_block_id)
            fields = self._get_element_block_fields(element_block_id)
            # create the field
            self.create_element_field(element_field_name, element_block_id)
            element_field = fields[element_field_name]
            # for each timestep
            for timestep_index in xrange(len(self.timesteps)):
                # initialize element field
                element_field_values = [0.0] * element_count
                # for each node within each element
                for connectivity_index, node_index in enumerate(connectivity):
                    element_index = connectivity_index / nodes_per_element
                    element_field_values[element_index] += node_field[
                        timestep_index][node_index]
                # average each value
                element_field_values = [x / float(nodes_per_element)
                                        for x in element_field_values]
                # replace field values with generated values
                element_field[timestep_index] = element_field_values

    def _find_member_indices_to_keep(self, all_members, members_to_delete):
        """
        Return a list of indices to keep from the list.

        This will handle duplicate entries correctly.  Duplicate entries will
        either all be deleted or all be kept.  Each index will appear in the
        resulting list at most once.

        """
        # find entries which are not in either list
        all_members_set = set(all_members)
        members_to_delete_set = set(members_to_delete)
        invalid_entries = members_to_delete_set - all_members_set
        if invalid_entries:
            self._warning('Invalid members',
                          'The member list contains entries that were not '
                          'members.  There were %d such occurrences.  These '
                          'will be ignored.'
                          % (len(invalid_entries)))
        # find members to delete
        indices_to_keep = [i
                           for i, x in enumerate(all_members)
                           if x not in members_to_delete_set]
        return indices_to_keep

    def _delete_side_set_members(self, side_set_id, side_set_members):
        """
        Delete the specified members from the given node set.

        Node set fields are updates.

        Members are given as (node_index).

        """
        # validate input
        [side_set_id] = self._format_side_set_id_list(
            [side_set_id],
            single=True)
        members = self.get_side_set_members(side_set_id)
        fields = self._get_side_set_fields(side_set_id)
        # find members to keep
        indices_to_keep = self._find_member_indices_to_keep(members,
                                                            side_set_members)
        # delete members from node set fields
        for all_values in fields.values():
            all_values[:] = [[values[x] for x in indices_to_keep]
                             for values in all_values]
        # delete members
        members[:] = [members[x] for x in indices_to_keep]

    def _delete_node_set_members(self, node_set_id, node_set_members):
        """
        Delete the specified members from the given node set.

        Node set fields are updates.

        Members are given as (node_index).

        """
        # validate input
        [node_set_id] = self._format_node_set_id_list(
            [node_set_id],
            single=True)
        node_set = self.node_sets[node_set_id]
        # find members to keep
        indices_to_keep = self._find_member_indices_to_keep(node_set[0],
                                                            node_set_members)
        # delete members from node set fields
        node_set[-1] = [[these_values[x] for x in indices_to_keep]
                        for these_values in node_set[1]]
        # delete members
        node_set[0] = [node_set[0][x] for x in indices_to_keep]

    def _delete_nodes(self, node_list):
        """
        Delete the given nodes.

        This will also delete all references to those nodes in node sets.  If a
        node is still being used by an element, it cannot be deleted.

        Example:
        >>> model.delete_nodes([0, 1, 2, 3])

        """
        node_list = self._remove_duplicates(node_list, preserve_order=False)
        # find node mapping
        # old node i refers to new node node_map[i]
        keep_node = [True] * len(self.nodes)
        for node_index in node_list:
            keep_node[node_index] = False
        node_map = [None] * len(self.nodes)
        reverse_node_map = []
        next_index = 0
        for index, keep in enumerate(keep_node):
            if keep:
                reverse_node_map.append(index)
                node_map[index] = next_index
                next_index += 1
        # delete nodes
        new_nodes = [self.nodes[x] for x in reverse_node_map]
        self.nodes = new_nodes
        # update connectivity in each element block
        for element_block_id in self.get_element_block_ids():
            connectivity = self.get_connectivity(element_block_id)
            new_connectivity = [node_map[x] for x in connectivity]
            if None in new_connectivity:
                self._error('Node still used.',
                            'A node in the list of nodes to delete is still '
                            'used by elements in element block %d and cannot '
                            'be deleted.' % element_block_id)
            connectivity[:] = new_connectivity
        # update node fields
        for field in self.node_fields.values():
            for timestep_index in xrange(len(self.timesteps)):
                new_values = [field[timestep_index][x]
                              for x in reverse_node_map]
                field[timestep_index] = new_values
        # delete nodes from node sets and fields
        for node_set_id in self.get_node_set_ids():
            members = self.get_node_set_members(node_set_id)
            fields = self._get_node_set_fields(node_set_id)
            # find new mapping
            value_map = []
            new_members = []
            for index, member in enumerate(members):
                if node_map[member] is not None:
                    value_map.append(index)
                    new_members.append(node_map[member])
            # update member list
            members[:] = new_members
            # delete these nodes from the field
            for field in fields.values():
                for timestep_index in xrange(len(self.timesteps)):
                    new_values = [field[timestep_index][x]
                                  for x in value_map]
                    field[timestep_index] = new_values

    def _get_unreferenced_nodes(self):
        """Return a list of node indices which are not used by any element."""
        used_node = [False] * len(self.nodes)
        for id_ in self.get_element_block_ids():
            connectivity = self.get_connectivity(id_)
            for node_index in connectivity:
                used_node[node_index] = True
        unused_nodes = [index
                        for index, used in enumerate(used_node)
                        if not used]
        return unused_nodes

    def delete_unused_nodes(self):
        """
        Delete nodes which are not used by any elements.

        Example:
        >>> model.delete_unused_nodes()

        """
        nodes_to_delete = self._get_unreferenced_nodes()
        self._delete_nodes(nodes_to_delete)

    def delete_global_variable(self, global_variable_names):
        """
        Delete one or more global variables.

        This function may also be called by its plural form
        'delete_global_variables()'.

        Examples:
        >>> model.delete_global_variable('internal_energy')
        >>> model.delete_global_variable('all')

        """
        global_variable_names = self._format_id_list(
            global_variable_names,
            self.get_global_variable_names(),
            'global variable')
        for name in global_variable_names:
            del self.global_variables[name]

    def delete_side_set(self, side_set_ids):
        """
        Delete one or more side sets.

        This function may also be called by its plural form
        'delete_side_sets()'.

        Examples:
        >>> model.delete_side_set(1)
        >>> model.delete_side_set('all')

        """
        side_set_ids = self._format_side_set_id_list(side_set_ids)
        for side_set_id in side_set_ids:
            del self.side_sets[side_set_id]

    def delete_empty_side_sets(self):
        """
        Delete all side sets with zero members.

        Example:
        >>> model.delete_empty_side_sets()

        """
        for id_ in self.get_side_set_ids():
            if not self.get_side_set_members(id_):
                self.delete_side_set(id_)

    def delete_node_set(self, node_set_ids):
        """
        Delete the given node set.

        This function may also be called by its plural form
        'delete_node_sets()'.

        Example:
        >>> model.delete_node_set(1)

        """
        node_set_ids = self._format_node_set_id_list(node_set_ids)
        for node_set_id in node_set_ids:
            del self.node_sets[node_set_id]

    def delete_empty_node_sets(self):
        """
        Delete all side sets with zero members.

        Example:
        >>> model.delete_empty_node_sets()

        """
        for id_ in self.get_node_set_ids():
            if not self.get_node_set_members(id_):
                self.delete_node_set(id_)

    def summarize(self):
        """
        Print a summary of the information in the current model.

        Example:
        >>> model.summarize()

        """
        print('The model contains the following:')
        print('  %d elements' % (self.get_element_count()))
        print('  %d timesteps' % (len(self.timesteps)))
        print('  %d global variables' % (len(self.global_variables)))
        print('  %d nodes' % (len(self.nodes)))
        print('  %d node fields' % (len(self.node_fields)))
        print('  %d element blocks' % (len(self.element_blocks)))
        print('  %d element fields' % (len(self.get_element_field_names())))
        print('  %d node sets' % (len(self.node_sets)))
        print('  %d node set fields' % (len(self.get_node_set_field_names())))
        print('  %d side sets' % (len(self.side_sets)))
        print('  %d side set fields' % (len(self.get_side_set_field_names())))

    def _get_element_block_fields(self, element_block_id):
        """Return the dictionary of element block field values."""
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        return self.element_blocks[element_block_id][-1]

    def _get_node_set_fields(self, node_set_id):
        """Return the dictionary of node set field values."""
        [node_set_id] = self._format_node_set_id_list(
            [node_set_id],
            single=True)
        return self.node_sets[node_set_id][-1]

    def get_element_count(self, element_block_ids='all'):
        """
        Return the total number of elements in the given element blocks.

        Example:
        >>> print object.get_element_count()

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        return sum(self.element_blocks[x][1][1]
                   for x in element_block_ids)

    def process_element_fields(self, element_block_ids='all'):
        """
        Process element field information to create node based fields.

        For element fields with 8 integration points, this takes the average.
        This is useful for fully integrated or selective deviatoric elements
        within the SIERRA/SM code.

        For element fields with 9 integration points, this takes the first one.
        This is useful for q1p0 elements within the SIERRA/SM code.

        This function is provided as a convenience for post processing hex8
        elements with multiple integration points.

        Example:
        >>> model.process_element_fields()

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        # for each element block
        for block_id in element_block_ids:
            element_field_names = self.get_element_field_names(block_id)
            for element_field_name in element_field_names:
                # if it's the first integration point, count the total number
                # and process accordingly
                if re.match('^.*_1$', element_field_name):
                    prefix = element_field_name[:-1]
                    points = 1
                    while self.element_field_exists(
                            prefix + str(points + 1),
                            block_id):
                        points += 1
                    if points == 9:
                        self.rename_element_field(prefix + '1',
                                                  prefix[:-1],
                                                  block_id)
                        for i in xrange(2, 10):
                            self.delete_element_field(prefix + str(i),
                                                      block_id)
                    elif points == 8:
                        self.create_averaged_element_field(
                            [prefix + str(x) for x in xrange(1, 9)],
                            prefix[:-1],
                            block_id)
                        for i in xrange(1, 9):
                            self.delete_element_field(prefix + str(i),
                                                      block_id)
        # nodal average all fields
        for element_field_name in self.get_element_field_names():
            self.convert_element_field_to_node_field(element_field_name)
        # delete element fields
        self.delete_element_field('all')

    def combine_element_blocks(self,
                               element_block_ids,
                               target_element_block_id='auto'):
        """
        Combine multiple element blocks into a single block.

        By default, the target element block id will be the smallest of the
        merged element block ids.

        The element blocks to combine must have the same element type.

        Example:
        >>> model.combine_element_blocks('all', 1)

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids,
            empty_list_okay=False)
        if target_element_block_id == 'auto':
            target_element_block_id = min(element_block_ids)
        if len(element_block_ids) == 1:
            self.rename_element_block(element_block_ids[0],
                                      target_element_block_id)
            return
        # ensure all blocks have the same number of nodes per element
        nodes_per_element = set(self.get_nodes_per_element(x)
                                for x in element_block_ids)
        if len(nodes_per_element) != 1:
            self._error('Incompatible element types.',
                        'The number of nodes per element on each element '
                        'block to be merged must by the same.  This is an '
                        'ExodusII file requirement.')
        # ensure all blocks are the same type
        element_types = set(
            self._get_standard_element_type(self._get_element_type(x))
            for x in element_block_ids)
        if len(element_types) != 1:
            self._warning('Element types vary.',
                          'The element types of the merged blocks differ.  '
                          'This may cause issues.')
        # create new connectivity
        new_connectivity = list(
            itertools.chain(*[self.get_connectivity(x)
                              for x in element_block_ids]))
        # create new info
        new_info = list(self._get_block_info(element_block_ids[0]))
        # set element count to zero
        new_info[1] = 0
        # find a temporary element block id
        temporary_element_block_id = self._new_element_block_id()
        # create translation list for elements
        # new_block[old_block_id] = (new_block_id, element_offset)
        new_block_info = dict((id_, (id_, 0))
                              for id_ in self.get_element_block_ids())
        for id_ in element_block_ids:
            new_block_info[id_] = (temporary_element_block_id, new_info[1])
            new_info[1] += self.element_blocks[id_][1][1]
        # for fields which only exist on certain block, define them on all
        element_field_names = self.get_element_field_names(element_block_ids)
        give_nonexistant_field_warning = True
        for name in element_field_names:
            for id_ in element_block_ids:
                if not self.element_field_exists(name, id_):
                    if give_nonexistant_field_warning:
                        self._warning('Inconsistent element fields.',
                                      'The element field "%s" is not defined '
                                      'on all element blocks included in the '
                                      'merging operation.  It will be '
                                      'created.  Future warnings of this type '
                                      'will be suppressed.' % name)
                        give_nonexistant_field_warning = False
                    self.create_element_field(name, id_)
        # create new field information
        new_fields = dict()
        for name in element_field_names:
            new_values = [[] for _ in self.timesteps]
            for id_ in element_block_ids:
                for index, values in enumerate(
                        self.element_blocks[id_][-1][name]):
                    new_values[index].extend(values)
            new_fields[name] = new_values
        # create new block
        self.create_element_block(temporary_element_block_id,
                                  new_info,
                                  new_connectivity)
        self.element_blocks[temporary_element_block_id][-1] = new_fields
        # translate side set members
        for side_set_id in self.get_side_set_ids():
            members = self.get_side_set_members(side_set_id)
            new_members = [
                (new_block_info[element_block_id][0],
                 element_index + new_block_info[element_block_id][1],
                 face_index)
                for element_block_id, element_index, face_index in members]
            members[:] = new_members
        # delete old blocks
        # (nodes can not be orphaned by this procedure)
        for element_block_id in element_block_ids:
            self.delete_element_block(element_block_id,
                                      delete_orphaned_nodes=False)
        # rename element block
        self.rename_element_block(temporary_element_block_id,
                                  target_element_block_id)

    def duplicate_element_block(self,
                                element_block_id,
                                new_element_block_id,
                                duplicate_nodes=True):
        """
        Create an duplicate of the given element block.

        Nodes are duplicated by default.  The new element block references
        these duplicated nodes, not the original ones.

        Example:
        >>> model.duplicate_element_block(1, 2)

        """
        [element_block_id] = self._format_element_block_id_list(
            [element_block_id],
            single=True)
        info = list(self._get_block_info(element_block_id))
        old_connectivity = self.get_connectivity(element_block_id)
        # create new nodes
        if duplicate_nodes:
            unique_node_indices = sorted(set(old_connectivity))
            new_node_indices = []
            self._duplicate_nodes(unique_node_indices, new_node_indices)
            new_node_indices = dict((x, y)
                                    for x, y in zip(unique_node_indices,
                                                    new_node_indices))
            new_connectivity = [new_node_indices[x]
                                for x in old_connectivity]
        else:
            new_connectivity = list(old_connectivity)
        # create the block
        self.create_element_block(new_element_block_id, info, new_connectivity)
        # copy fields
        fields = dict()
        for name, all_values in self._get_element_block_fields(
                element_block_id).items():
            fields[name] = [list(x) for x in all_values]
        self.element_blocks[new_element_block_id][-1] = fields
        # update side sets and side set fields
        for side_set_id in self.get_side_set_ids():
            members = self.get_side_set_members(side_set_id)
            fields = self._get_side_set_fields(side_set_id)
            new_members = []
            source_face = []
            for index, member in enumerate(members):
                if member[0] == element_block_id:
                    new_members.append((new_element_block_id,
                                        member[1],
                                        member[2]))
                    source_face.append(index)
            members.extend(new_members)
            for all_values in fields.values():
                for values in all_values:
                    new_values = [values[x] for x in source_face]
                    values.extend(new_values)

    def _get_inverted_connectivity(self, element_type):
        """Return the connectivity pemutation to invert an element."""
        element_type = self._get_standard_element_type(element_type)
        if element_type not in self.INVERTED_CONNECTIVITY:
            self._error('Unknown inversion algorithm.',
                        'The required numbering to invert an element of type '
                        '"%s" is not defined.' % element_type)
        return self.INVERTED_CONNECTIVITY[element_type]

    def _get_inverted_face_mapping(self, element_type):
        """Return the connectivity pemutation to invert an element face."""
        face_mapping = self._get_face_mapping(element_type)
        inversion_mapping = self._get_inverted_connectivity(element_type)
        original_face = [set(face) for _, face in face_mapping]
        new_face = [set(inversion_mapping[x] for x in face)
                    for _, face in face_mapping]
        try:
            new_face = [original_face.index(x) for x in new_face]
        except ValueError:
            self._bug('Unable to determine face mapping.',
                      'We were unable to determine how to map faces when '
                      'inverting an element of type "%s".' % element_type)
        return new_face

    def _invert_element_blocks(self, element_block_ids):
        """Invert all elements within one or more element blocks."""
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        for id_ in element_block_ids:
            element_count = self.get_element_count(id_)
            nodes_per_element = self.get_nodes_per_element(id_)
            # invert the connectivity
            inverted_mapping = self._get_inverted_connectivity(
                self._get_element_type(id_))
            connectivity = self.get_connectivity(id_)
            new_connectivity = [
                connectivity[element_index * nodes_per_element + x]
                for element_index in xrange(element_count)
                for x in inverted_mapping]
            connectivity[:] = new_connectivity
            # adjust side set members
            new_face_indices = self._get_inverted_face_mapping(
                self._get_element_type(id_))
            for side_set_id in self.get_side_set_ids():
                members = self.get_side_set_members(side_set_id)
                for index, member in enumerate(members):
                    if member[0] == id_:
                        members[index] = (id_,
                                          member[1],
                                          new_face_indices[member[2]])

    def create_element_block(self, element_block_id, info, connectivity=None):
        """
        Create a new element block.

        The nodes for the elements in the block must have already been defined.

        The info list should be comprised of the following information.
        * '[element_type, element_count, nodes_per_element, 0]'

        For example, the following would be valid.
        * '["hex8", elements, 8, 0]'

        The connectivity list should be a shallow list of element connectivity
        and must be of length 'element_count * nodes_per_element'.

        Element blocks are unnamed when created.  To name them, use the
        'rename_element_block' function.

        Example:
        >>> model.create_element_block(1, ['hex8', 0, 8, 0])

        """
        # make sure it doesn't exist already
        if self.element_block_exists(element_block_id):
            self._exists_error(element_block_id, 'element block')
        # set up an empty connectivity if none is given
        if not connectivity:
            connectivity = []
        # create the actual block
        self.element_blocks[element_block_id] = ['',
                                                 info,
                                                 connectivity,
                                                 {}]

    def create_nodes(self, new_nodes):
        """
        Create nodes corresponding to the given coordinate list.

        The list must contain coordinate triples '[x, y, z]' defined for each
        new node.  The new nodes are assigned a local index starting with the
        current length of the nodes list.

        Example:
        >>> model.create_nodes([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])

        """
        # add nodes
        self.nodes += new_nodes
        # adjust node_fields
        for name, all_values in self.node_fields.items():
            default_value = self._get_default_field_value(name)
            for values in all_values:
                values.extend([default_value] * len(new_nodes))

    def _exists_error(self, name, entity):
        """Warn the user something already exists and exit."""
        self._error(
            entity[0].upper() + entity[1:] + ' already exists.',
            'The specified %s "%s" already exists.  This operation cannot '
            'be completed.' % (entity, str(name)))

    def _exists_warning(self, name, entity):
        """Warn the user something already exists."""
        self._warning(
            entity[0].upper() + entity[1:] + ' already exists.',
            'The specified %s "%s" already exists.  Information may be '
            'overwritten.' % (entity, str(name)))

    def _exists_on_entity_warning(self, name, entity, base_name, base_entity):
        """Warn the user something already exists on a given entity."""
        self._warning(
            entity[0].upper() + entity[1:] + ' already exists.',
            'The specified %s "%s" already exists on %s %s.  Information may '
            'be overwritten.'
            % (entity, str(name), base_entity, str(base_name)))

    def _missing_error(self, name, entity):
        """Tell the user something does not exist and exit."""
        self._error(
            entity[0].upper() + entity[1:] + ' does not exist.',
            'The specified %s "%s" does not exist.'
            % (entity, str(name)))

    def _missing_warning(self, name, entity):
        """Warn the user something does not exist."""
        self._error(
            entity[0].upper() + entity[1:] + ' does not exist.',
            'The specified %s "%s" does not exist and will be ignored.'
            % (entity, str(name)))

    def _missing_on_entity_error(self, name, entity, base_name, base_entity):
        """Tell the user something does not exist on an entity and exit."""
        self._error(
            entity[0].upper() + entity[1:] + ' does not exist.',
            'The specified %s "%s" does not exist on %s %s.'
            % (entity, str(name), base_entity, str(base_name)))

    def create_global_variable(self, global_variable_name, value='auto'):
        """
        Create a new global variable.

        Example:
        >>> model.create_global_variable('gravitational_acceleration', 9.8)

        """
        if self.global_variable_exists(global_variable_name):
            self._exists_warning(global_variable_name, 'global variable')
        if value == 'auto':
            value = self._get_default_field_value(global_variable_name)
        values = [value] * len(self.timesteps)
        self.global_variables[global_variable_name] = values

    def calculate_element_centroids(self,
                                    element_field_name_prefix='centroid_',
                                    element_block_ids='all'):
        """
        Calculate and store the centroid of each element.

        This will approximate the element centroid as the nodal average of each
        element and will store that value in an element field.  Since a
        timestep must be defined in order for element fields to exist, one will
        be created if none exist.

        By default, the centroid will be stored in the fields 'centroid_x',
        'centroid_y', and 'centroid_z'.  Alternatively, a prefix can be given
        or a list of three strings can be given.

        Example:
        >>> model.calculate_element_centroids()

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        if not self.timesteps:
            self.create_timestep(0.0)
        if isinstance(element_field_name_prefix, str):
            centroid_field_names = [element_field_name_prefix + x
                                    for x in ['x', 'y', 'z']]
        else:
            centroid_field_names = element_field_name_prefix
        for element_block_id in element_block_ids:
            for name in centroid_field_names:
                if self.element_field_exists(name, element_block_id):
                    self._exists_on_entity_warning(name,
                                                   'element field',
                                                   element_block_id,
                                                   'element block')
            # calculate centroids
            centroid = [[], [], []]
            element_count = self.get_element_count(element_block_id)
            nodes_per_element = self.get_nodes_per_element(element_block_id)
            connectivity = self.get_connectivity(element_block_id)
            for element_index in xrange(element_count):
                this_centroid = [0.0, 0.0, 0.0]
                for connectivity_index in xrange(nodes_per_element):
                    node_index = connectivity[
                        element_index * nodes_per_element + connectivity_index]
                    for i in xrange(3):
                        this_centroid[i] += self.nodes[node_index][i]
                for i in xrange(3):
                    centroid[i].append(this_centroid[i] / nodes_per_element)
            fields = self._get_element_block_fields(element_block_id)
            for index, name in enumerate(centroid_field_names):
                values = []
                for _ in xrange(len(self.timesteps)):
                    values.append(list(centroid[index]))
                fields[name] = values

    def calculate_element_volumes(self,
                                  element_field_name='volume',
                                  element_block_ids='all'):
        """
        Calculate and store the volume of each element.

        This will approximate the element volume.  Since a timestep must be
        defined in order for element fields to exist, one will be created if
        none exist.

        For two dimensional elements, this calculates the area.  For one
        dimensional elements, this calculates the length.

        Example:
        >>> model.calculate_element_volumes()

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        if not self.timesteps:
            self.create_timestep(0.0)
        for element_block_id in element_block_ids:
            # get the element type
            element_type = self._get_standard_element_type(
                self._get_element_type(element_block_id))
            if element_type not in self.VOLUME_FORMULA:
                self._warning('Unrecognized element type',
                              'The formula for calculating the volume of '
                              'element type "%s" is not implemented or not '
                              'known.  This block will be skipped.'
                              % (element_type))
                continue
            # warn if field already exists
            if self.element_field_exists(element_field_name, element_block_id):
                self._exists_on_entity_warning(element_field_name,
                                               'element field',
                                               element_block_id,
                                               'element block')
            # get the formula
            formula = self.VOLUME_FORMULA[element_type]
            if len(formula) == 2:
                # distance between two points
                equation = ('(x[0] ** 2.0 + '
                            'x[1] ** 2.0 + '
                            'x[2] ** 2.0) ** 0.5')
            elif len(formula) == 3:
                # magnitude of cross product
                equation = ('((x[0] * x[4] - x[1] * x[3]) ** 2.0 + '
                            '(x[2] * x[3] - x[0] * x[5]) ** 2.0 + '
                            '(x[1] * x[5] - x[2] * x[4]) ** 2.0) ** 0.5')
            elif len(formula) == 4:
                # triple product
                equation = ('(x[1] * x[5] - x[2] * x[4]) * x[6] + '
                            '(x[2] * x[3] - x[0] * x[5]) * x[7] + '
                            '(x[0] * x[4] - x[1] * x[3]) * x[8]')
            else:
                self._bug('Unknown case')
            # process the formula to account for averaged nodes
            transforms = []
            for i, rule in enumerate(formula[1:]):
                for d in range(3):
                    value = '('
                    node_list = rule[1]
                    if isinstance(node_list, int):
                        value += 'n[%d]' % (node_list * 3 + d)
                    else:
                        value += '('
                        value += 'n[%d]' % (node_list[0] * 3 + d)
                        for x in node_list[1:]:
                            value += ' + n[%d]' % (x * 3 + d)
                        value += ') / %d.0' % len(node_list)
                    value += ' - '
                    node_list = rule[0]
                    if isinstance(node_list, int):
                        value += 'n[%d]' % (node_list * 3 + d)
                    else:
                        value += '('
                        value += 'n[%d]' % (node_list[0] * 3 + d)
                        for x in node_list[1:]:
                            value += ' + n[%d]' % (x * 3 + d)
                        value += ') / %d.0' % len(node_list)
                    value += ')'
                    transforms.append(('x[%d]' % (i * 3 + d), value))
            for transform in transforms:
                equation = equation.replace(transform[0], transform[1])
            # store coordinates of nodes into x, y, z lists
            function = eval('lambda n: ' + equation)
            # store some element block values
            element_count = self.get_element_count(element_block_id)
            nodes_per_element = self.get_nodes_per_element(element_block_id)
            connectivity = self.get_connectivity(element_block_id)
            volumes = []
            for element_index in xrange(element_count):
                # get local node values
                local_node = connectivity[
                    element_index * nodes_per_element:
                    (element_index + 1) * nodes_per_element]
                # get local coordinates
                n = [x for i in local_node for x in self.nodes[i]]
                # add the volume
                volumes.append(float(function(n)) * formula[0])
            # now make the field for each timestep
            values = [list(volumes) for _ in range(len(self.timesteps))]
            # assign the field
            fields = self._get_element_block_fields(element_block_id)
            fields[element_field_name] = values

    def get_element_block_volume(self,
                                 element_block_ids,
                                 element_volume_field_name=None):
        """Return the total volume of the given element blocks.

        Example:
        >>> model.get_element_block_volume('all')

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        # create a timestep if none exist
        created_timestep = False
        if not self.timesteps:
            created_timestep = True
            self.create_timestep(0.0)
        # calculate temporary field with element volumes
        created_volume_field = False
        if element_volume_field_name is None:
            created_volume_field = True
            element_volume_field_name = self._new_element_field_name()
            self.calculate_element_volumes(element_volume_field_name,
                                           element_block_ids)
        # add up the volumes
        volume = 0.0
        for id_ in element_block_ids:
            fields = self._get_element_block_fields(id_)
            volume += sum(fields[element_volume_field_name][0])
        # delete the temporary timestep
        if created_timestep:
            self.delete_timestep(0.0)
        # delete the temporary field
        if created_volume_field:
            self.delete_element_field(element_volume_field_name,
                                      element_block_ids)
        return volume

    def get_element_block_centroid(self,
                                   element_block_ids,
                                   element_volume_field_name=None,
                                   element_centroid_field_names=None):
        """Return the centroid of the given element blocks.

        Example:
        >>> model.get_element_block_centroid('all')
        [1.0, 2.0, 3.0]

        """
        # format inputs
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        # create a timestep if none exist
        created_timestep = False
        if not self.timesteps:
            created_timestep = True
            self.create_timestep(0.0)
        # calculate temporary field with element volumes
        created_volume_field = False
        if element_volume_field_name is None:
            created_volume_field = True
            element_volume_field_name = self._new_element_field_name()
            self.calculate_element_volumes(element_volume_field_name,
                                           element_block_ids)
        # calculate temporary field with element centroids
        created_centroid_fields = False
        if element_centroid_field_names is None:
            created_centroid_fields = True
            element_centroid_field_names = self._new_element_field_name(3)
            self.calculate_element_centroids(element_centroid_field_names,
                                             element_block_ids)
        # calculate the centroid
        total_volume = 0.0
        centroid = [0.0, 0.0, 0.0]
        for id_ in element_block_ids:
            fields = self._get_element_block_fields(id_)
            volumes = fields[element_volume_field_name][0]
            centroids = [fields[element_centroid_field_names[x]][0]
                         for x in range(3)]
            for x in range(3):
                centroid[x] += sum(a * b
                                   for a, b in itertools.izip(volumes,
                                                              centroids[x]))
            total_volume += sum(volumes)
        if total_volume == 0.0:
            centroid = [float('nan')] * 3
        else:
            centroid = [x / total_volume for x in centroid]
        # delete the temporary timestep
        if created_timestep:
            self.delete_timestep(0.0)
        # delete the temporary fields
        if created_volume_field:
            self.delete_element_field(element_volume_field_name,
                                      element_block_ids)
        if created_centroid_fields:
            for name in element_centroid_field_names:
                self.delete_element_field(name, element_block_ids)
        return centroid

    def create_interpolated_timestep(self, timestep, interpolation='cubic'):
        """
        Create a new timestep by interpolating neighboring steps.

        This does not extrapolate, so the given timestep to be interpolated
        must lie within the range of timesteps already defined.

        By default, 'cubic' interpolation is used.  This can be changed to
        'linear' interpolation if desired.

        Examples:
        >>> model.create_interpolated_timestep(0.5)
        >>> model.create_interpolated_timestep(0.5, interpolation='linear')

        """
        if timestep in self.timesteps:
            self._exists_warning(timestep, 'timestep')
            return
        # find neighboring timesteps and how to form the new step
        steps = self.get_timesteps()
        # make sure enough steps exist
        if len(steps) < 2:
            self._error('Invalid interpolation.',
                        'A timestep cannot be interpolated unless 2 or '
                        'more timesteps already exist.  There are %d defined '
                        'timesteps: %s' %
                        (len(self.timesteps),
                         ', '.join([str(x) for x in self.get_timesteps()])))
        # make sure value lies within time bounds
        if steps[0] > timestep or steps[-1] < timestep:
            self._error('Invalid interpolation.',
                        'The specified timestep %s does not lie within the '
                        'range of timesteps already defined: [%s, %s].'
                        % (str(timestep), str(steps[0]), str(steps[-1])))
        if interpolation == 'linear':
            # find two bounding timesteps
            nearby = min(steps, key=lambda x: abs(x - timestep))
            if nearby < timestep:
                nearby = [nearby,
                          steps[steps.index(nearby) + 1]]
            else:
                nearby = [steps[steps.index(nearby) - 1],
                          nearby]
            # find proportion of each one to use
            phi = (timestep - nearby[0]) / (nearby[1] - nearby[0])
            formula = [(nearby[0], 1.0 - phi), (nearby[1], phi)]
        elif interpolation == 'cubic':
            # find four bounding timesteps
            # if step is within first or last segment, create
            # an imaginary point to do the interpolation
            nearby = min(steps, key=lambda x: abs(x - timestep))
            index = steps.index(nearby)
            four_steps = [0.0] * 4
            if nearby > timestep:
                index -= 1
            four_steps[1] = steps[index]
            four_steps[2] = steps[index + 1]
            if index == 0:
                four_steps[0] = 2 * four_steps[1] - four_steps[2]
            else:
                four_steps[0] = steps[index - 1]
            if index + 2 == len(steps):
                four_steps[3] = 2 * four_steps[2] - four_steps[1]
            else:
                four_steps[3] = steps[index + 2]
            # find interpolation coefficients
            coefficients = self._cubic_interpolation(timestep, *four_steps)
            formula = [[steps[index], coefficients[1]],
                       [steps[index + 1], coefficients[2]]]
            if index == 0:
                formula[0][1] += 2 * coefficients[0]
                formula[1][1] -= coefficients[0]
            else:
                formula.append([steps[index - 1], coefficients[0]])
            if index + 2 == len(steps):
                formula[0][1] -= coefficients[3]
                formula[1][1] += 2 * coefficients[3]
            else:
                formula.append([steps[index + 2], coefficients[3]])
        else:
            self._error('Unknown interpolation technique',
                        'The specified interpolation technique "%s" is not '
                        'recognized.' % interpolation)
        # create the new timestep
        self.create_timestep(timestep)
        # use the given formula to create the new step
        # formula = list of (timestep_index, contribution)
        formula = [(self._get_internal_timestep_index(x[0]), x[1])
                   for x in formula]
        this_index = self._get_internal_timestep_index(timestep)
        # element fields, side set fields, node set fields
        for thing in itertools.chain(self.element_blocks.values(),
                                     self.node_sets.values(),
                                     self.side_sets.values()):
            fields = thing[-1]
            for values in fields.values():
                new_values = [0.0] * len(values[0])
                for index, phi in formula:
                    for i in xrange(len(values[0])):
                        new_values[i] += values[index][i] * phi
                values[this_index] = new_values
        # node fields
        for values in self.node_fields.values():
            new_values = [0.0] * len(values[0])
            for index, phi in formula:
                for i in xrange(len(values[0])):
                    new_values[i] += values[index][i] * phi
            values[this_index] = new_values
        # global variables
        for values in self.global_variables.values():
            new_value = sum([values[x[0]] * x[1] for x in formula])
            values[this_index] = new_value

    def _find_new_timestep_index(self, new_timestep):
        """
        Return the index at which to create the new timestep.

        This returns the index to keep the list in monotonically increasing
        order.

        """
        # find the index at which to to insert this
        # if the new step is above all the others, put it at the end
        higher_times = [x
                        for x in self.timesteps
                        if x > new_timestep]
        # if the new time is higher than all existing times, then add it to the
        # end
        if not higher_times:
            return len(self.timesteps)
        return self.timesteps.index(min(higher_times))

    def copy_timestep(self, timestep, new_timestep):
        """
        Create a copy of an existing timestep.

        This copies all information from the old timestep including values
        of global variables, node fields, element fields, node set fields and
        side set fields.

        Example:
        object.copy_timestep(0.0, 1.0)

        """
        [timestep] = self._format_id_list(
            [timestep],
            self.get_timesteps(),
            'timestep',
            single=True)
        if self.timestep_exists(new_timestep):
            self._exists_warning(new_timestep, 'timestep')
            return
        # find the index at which to insert the new timestep
        new_index = self._find_new_timestep_index(new_timestep)
        # get the index of the existing timestep prior to adding the new
        # timestep
        old_index = self._get_internal_timestep_index(timestep)
        # insert the new timestep
        self.timesteps.insert(new_index, new_timestep)
        # adjust node_fields
        for field in self.node_fields.values():
            field.insert(new_index, list(field[old_index]))
        # adjust self.global_variables
        for values in self.global_variables.values():
            values.insert(new_index, values[old_index])
        # adjust self.element_blocks, self.node_sets, self.side_sets
        for thing in itertools.chain(self.element_blocks.values(),
                                     self.node_sets.values(),
                                     self.side_sets.values()):
            fields = thing[-1]
            for field in fields.values():
                field.insert(new_index, list(field[old_index]))

    def create_timestep(self, timestep):
        """
        Create a new timestep.

        Field information at the created timestep is set to the default value.

        Example:
        >>> model.create_timestep(0)

        """
        timestep = float(timestep)
        if self.timestep_exists(timestep):
            self._exists_warning(timestep, 'timestep')
            return
        # find the index at which to to insert this
        timestep_index = self._find_new_timestep_index(timestep)
        # add the given timestep
        self.timesteps.insert(timestep_index, timestep)
        # adjust self.node_fields
        for name, field in self.node_fields.items():
            value = self._get_default_field_value(name)
            field.insert(timestep_index, [value] * len(self.nodes))
        # adjust element blocks fields
        for id_ in self.get_element_block_ids():
            element_count = self.get_element_count(id_)
            fields = self._get_element_block_fields(id_)
            for name, field_values in fields.items():
                value = self._get_default_field_value(name)
                field_values.insert(timestep_index, [value] * element_count)
        # adjust side set fields
        for id_ in self.get_side_set_ids():
            members = self.get_side_set_members(id_)
            fields = self._get_side_set_fields(id_)
            member_count = len(members)
            for name, values in fields.items():
                value = self._get_default_field_value(name)
                values.insert(timestep_index, [value] * member_count)
        # adjust node set fields
        for id_ in self.get_node_set_ids():
            members = self.get_node_set_members(id_)
            fields = self._get_node_set_fields(id_)
            member_count = len(members)
            for name, values in fields.items():
                value = self._get_default_field_value(name)
                values.insert(timestep_index, [value] * member_count)
        # adjust self.global_variables
        for name, values in self.global_variables.items():
            values.insert(timestep_index, self._get_default_field_value(name))

    def _replace_name_case(self, new_list, original_list):
        """
        Return the lowercase version of all strings in the given list.

        Example:
        >>> model._replace_name_case(['x', 'z', 'fred'], ['X', 'Fred', 'Z'])
        ['X', 'Z', 'Fred']

        """
        original_case = dict((x.lower(), x)
                             for x in original_list)
        if len(original_case) != len(original_list):
            self._warning('Ambiguous string case.',
                          'There are multiple strings in the list which have '
                          'identical lowercase representations.  One will be '
                          'chosen at random.')
        for item in new_list:
            if not item.lower() in original_case:
                self._bug('Unrecognized string.',
                          'The string "%s" appears in the new list but '
                          'not in the original list.' % item)
        return [original_case[x.lower()] for x in new_list]

    def _sort_field_names(self, original_field_names):
        """
        Return field names sorted in a SIERRA-friendly manner.

        In order for SIERRA to recognize vectors, tensors, and element fields
        with multiple integration points, fields must be sorted in a specific
        order.  This function provides that sort order.

        As fields within exomerge are stored in a set, exomerge has no internal
        or natural field order.  This routine is only necessary for writing to
        ExodusII files.

        """
        field_names = [x.lower() for x in original_field_names]
        # Look through all fields to find multi-component fields and store
        # these as tuples of the following form:
        # ('base_name', 'component', integration_points)
        multicomponent_fields = set()
        for name in field_names:
            # see if it has an integration point
            if re.match('.*_[0-9]+$', name):
                (name, integration_point) = name.rsplit('_', 1)
                integration_point = int(integration_point)
            else:
                integration_point = None
            # see if it possibly has a component
            if re.match('.*_.+$', name):
                component = name.rsplit('_', 1)[1]
                if component in self.ALL_MULTI_COMPONENT_FIELD_SUBSCRIPTS:
                    name = name.rsplit('_', 1)[0]
                    multicomponent_fields.add(
                        (name, component, integration_point))
        # now sort multi-component fields
        base_names = set(x for x, _, _ in multicomponent_fields)
        sorted_field_names = dict()
        field_names = set(field_names)
        for base_name in base_names:
            # find all components of this form
            components = set(x
                             for name, x, _ in multicomponent_fields
                             if name == base_name)
            # find max integration point value
            integration_point_count = max(
                x
                for name, _, x in multicomponent_fields
                if name == base_name)
            # see if the components match the form of something
            matching_form = None
            for form, included_components in (
                    self.MULTI_COMPONENT_FIELD_SUBSCRIPTS.items()):
                if set(included_components) == components:
                    matching_form = form
            if not matching_form:
                continue
            # see if all components and integration points are present
            mid = ['_' + x
                   for x
                   in self.MULTI_COMPONENT_FIELD_SUBSCRIPTS[matching_form]]
            if integration_point_count is None:
                last = ['']
            else:
                last = ['_' + str(x + 1)
                        for x in xrange(integration_point_count)]
            all_names = [base_name + m + l
                         for l in last
                         for m in mid]
            if set(all_names).issubset(field_names):
                sorted_field_names[all_names[0]] = all_names
                field_names = field_names - set(all_names)
        # sort field names which are not part of multicomponent fields
        field_names = sorted(field_names)
        # for each list of field names, find place to splice into list
        place_to_insert = dict()
        for name in sorted_field_names.keys():
            place = bisect.bisect_left(field_names, name)
            if place not in place_to_insert:
                place_to_insert[place] = [name]
            else:
                place_to_insert[place].append(name)
        # splice them in
        for place in sorted(place_to_insert.keys(), reverse=True):
            for name in place_to_insert[place]:
                field_names[place:place] = sorted_field_names[name]
        return self._replace_name_case(field_names, original_field_names)

    def _reorder_list(self, the_list, new_index):
        """
        Reorder a list by permuting items.

        The item at index 0 in the original list is moved to index
        'new_index[0]', and so on.

        """
        # ensure arguments are valid
        self._input_check(new_index, [list, len(the_list), int])
        # create new list
        the_list[:] = [the_list[i] for i in new_index]

    def _sort_timesteps(self):
        """Sort timesteps so that they are in ascending order."""
        index_map = [self.timesteps.index(x) for x in sorted(self.timesteps)]
        self._reorder_list(self.timesteps, index_map)
        for element_block_id in self.get_element_block_ids():
            fields = self._get_element_block_fields(element_block_id)
            for values in fields.values():
                self._reorder_list(values, index_map)
        for values in self.node_fields.values():
            self._reorder_list(values, index_map)
        for values in self.global_variables.values():
            self._reorder_list(values, index_map)
        for node_set_id in self.get_node_set_ids():
            fields = self._get_node_set_fields(node_set_id)
            for values in fields.values():
                self._reorder_list(values, index_map)
        for side_set_id in self.get_side_set_ids():
            fields = self._get_side_set_fields(side_set_id)
            for values in fields.values():
                self._reorder_list(values, index_map)

    def _apply_node_map(self, node_map):
        """
        Apply the given node map to reorder all nodes.

        Node i becomes node node_map[i].

        """
        assert sorted(node_map) == range(len(self.nodes))
        # find reverse_node_map
        reverse_node_map = self._get_reverse_index_map(node_map)
        # reorder nodes
        self.nodes = [self.nodes[x] for x in reverse_node_map]
        # reorder element connectivity
        for element_block_id in self.get_element_block_ids():
            connectivity = self.get_connectivity(element_block_id)
            connectivity[:] = [node_map[x] for x in connectivity]
        # reorder node fields
        for field in self.node_fields.values():
            for timestep_index in xrange(len(self.timesteps)):
                values = field[timestep_index]
                field[timestep_index] = [values[x] for x in reverse_node_map]
        # reorder nodes in node sets
        for node_set_id in self.get_node_set_ids():
            members = self.get_node_set_members(node_set_id)
            members[:] = [node_map[x] for x in members]

    def _order_nodes(self):
        """
        Redefine node ordering in a consistent manner.

        This was introduced because running exodiff on two different
        models with the same information but differing node numbering will
        fail.  By running this routine before export_model, the numbering
        will be consistent and will avoid this problem.

        """
        # Define new node numbering by the order in which nodes are encountered
        # in element blocks by ascending element block id.  For nodes not used
        # in any element blocks, they are sorted by x coordinate, then y
        # coordinate, then z coordinate.
        # node i on old model becomes node_map[i]
        # node i in new model was reverse_node_map[i]
        next_index = 0
        node_map = [None] * len(self.nodes)
        for element_block_id in self.get_element_block_ids():
            connectivity = self.get_connectivity(element_block_id)
            for node_index in connectivity:
                if node_map[node_index] is None:
                    node_map[node_index] = next_index
                    next_index += 1
        still_to_sort = []
        for index, map_ in enumerate(node_map):
            if map_ is None:
                still_to_sort.append((self.nodes[index], index))
        still_to_sort.sort()
        for _, node_index in still_to_sort:
            node_map[node_index] = next_index
            next_index += 1
        self._apply_node_map(node_map)

    def get_input_deck(self):
        """
        Return each line of the input deck stored in the file.

        Many SIERRA applications, when running a problem, will store the input
        deck within the results file.  This function retrieves that
        information, if it exists.  Note that due to format restriction, the
        retrieved input deck may not exactly match the original file.

        Example:
        >>> object = exomerge.import_model('results.e')
        >>> model.get_input_deck()

        """
        continuation = False
        input_deck = []
        begin_depth = 0
        first_word = None
        for line in self.info_records:
            if not continuation:
                first_word = line.strip().split(' ', 1)[0].lower()
            if not continuation and first_word == 'begin':
                begin_depth += 1
            if begin_depth > 0:
                if continuation:
                    input_deck[-1] = input_deck[-1][:-1] + line
                else:
                    input_deck.append(line)
            continuation = (begin_depth > 0 and
                            len(line) == 80 and
                            line[-1] == '\\')
            if not continuation and first_word == 'end':
                begin_depth -= 1
        return '\n'.join(input_deck)

    def get_length_scale(self):
        """
        Return the length scale of the model.

        The length scale is defined as the largest of the following:
        * absolute nodal coordinate component
        * total range in nodal coordinate component

        Example:
        >>> model.get_length_scale()

        """
        if not self.nodes:
            return 0.0
        bounds = [[self.nodes[0][d], self.nodes[0][d]]
                  for d in range(3)]
        for node in self.nodes:
            for d in xrange(3):
                if node[d] < bounds[d][0]:
                    bounds[d][0] = node[d]
                if node[d] > bounds[d][1]:
                    bounds[d][1] = node[d]
        relative = max([bound[1] - bound[0] for bound in bounds])
        absolute = max([abs(x) for x in bound for bound in bounds])
        return max(relative, absolute)

    @staticmethod
    def _get_reverse_index_map(index_map):
        """Return the reverse map for the given index mapping scheme."""
        reverse_index_map = [None] * len(index_map)
        for index, new_index in enumerate(index_map):
            reverse_index_map[new_index] = index
        return reverse_index_map

    @staticmethod
    def _dot(vector_one, vector_two):
        """Return the dot product of two things."""
        return sum(a * b for a, b in zip(vector_one, vector_two))

    @staticmethod
    def _distance_squared_between(point_one, point_two):
        """Return the distance squared between two three-dimensional points."""
        return ((point_two[0] - point_one[0]) ** 2 +
                (point_two[1] - point_one[1]) ** 2 +
                (point_two[2] - point_one[2]) ** 2)

    @staticmethod
    def _distance_between(point_one, point_two):
        """Return the distance between two three-dimensional points."""
        return math.sqrt((point_two[0] - point_one[0]) ** 2 +
                         (point_two[1] - point_one[1]) ** 2 +
                         (point_two[2] - point_one[2]) ** 2)

    @staticmethod
    def _values_match(value, list_of_values):
        """
        Return True if the list contains no items apart from the given value.

        Example:
        >>> _values_match(0, [0, 0, 0])
        True
        >>> _values_match(1, [0, 1, 0])
        False

        """
        if type(value) is not float or not math.isnan(value):
            return not any(x != value for x in list_of_values)
        return not any(not math.isnan(x) for x in list_of_values)

    def _merge_node_groups(self, node_groups, suppress_warnings=False):
        """
        Merge nodes in the given node groups.

        This updates node field and node set field values with the average of
        the values at the merged nodes.  If they differ, a warning is output in
        addition.

        Node sets and node set fields are updated accordingly.

        Example:
        >>> model._merge_node_groups({1: [2, 3], 5: [7, 8]})

        """
        # ensure slave nodes are not duplicated
        slave_nodes_set = set()
        duplicated_slave_nodes = set()
        for master, slaves in node_groups.items():
            for slave in slaves:
                if slave in slave_nodes_set:
                    duplicated_slave_nodes.add(slave)
                else:
                    slave_nodes_set.add(slave)
        if duplicated_slave_nodes:
            problem_groups = dict()
            for master, slaves in node_groups.items():
                for index in slaves:
                    if index in duplicated_slave_nodes:
                        problem_groups[master] = slaves
            groups = '\n  '.join(str(x) + ': ' + str(y)
                                 for x, y in problem_groups.items())
            self._bug('Invalid merged node groups.',
                      'Slaves nodes were found in multiple merged groups.  '
                      'Conflicting merged node groups:\n  %s' % (groups))
        # ensure validity of input
        # slave nodes are not repeated
        # master nodes are never slave nodes
        master_nodes = sorted(node_groups.keys())
        slave_nodes = sorted(list(itertools.chain(*node_groups.values())))
        slave_nodes_set = set(slave_nodes)
        # ensure master nodes are not slave nodes
        for master in master_nodes:
            if master in slave_nodes_set:
                self._bug('Invalid merged node groups.',
                          'The master node %d is also found in a slave node '
                          'group.' % (master))
        # ensure slave nodes are not in multiple groups
        if not sorted(slave_nodes_set) == slave_nodes:
            self._bug('Invalid merged node groups.',
                      'Slave nodes are duplicated in multiple groups.')
        # First, remap all nodes such that slave nodes appear at the very
        # end.
        next_master_index = 0
        first_slave_index = len(self.nodes) - len(slave_nodes)
        next_slave_index = first_slave_index
        node_map = []
        index = 0
        for slave_node in slave_nodes:
            if slave_node != index:
                count = slave_node - index
                assert count > 0
                node_map.extend(xrange(next_master_index,
                                       next_master_index + count))
                index += count
                next_master_index += count
            node_map.append(next_slave_index)
            next_slave_index += 1
            index += 1
        count = first_slave_index - next_master_index
        node_map.extend(xrange(next_master_index,
                               next_master_index + count))
        next_master_index += count
        assert next_master_index == first_slave_index
        assert next_slave_index == len(self.nodes)
        for master, slaves in node_groups.items():
            assert node_map[master] < first_slave_index
            assert min([node_map[x] for x in slaves]) >= first_slave_index
        # apply this node map
        self._apply_node_map(node_map)
        # apply the mapping to node_groups
        node_groups = dict((node_map[key], [node_map[x] for x in values])
                           for key, values in node_groups.items())
        for master, slaves in node_groups.items():
            assert master < first_slave_index
            assert min(slaves) >= first_slave_index
        # get connectivity mapping
        connectivity_map = range(len(self.nodes))
        for master, slaves in node_groups.items():
            for slave in slaves:
                connectivity_map[slave] = master
        # change connectivity in element_blocks
        for element_block_id in self.get_element_block_ids():
            connectivity = self.get_connectivity(element_block_id)
            connectivity[:] = [connectivity_map[x] for x in connectivity]
        # change self.node_fields
        node_field_value_warnings = 0
        for name, all_values in self.node_fields.items():
            for values in all_values:
                for master, slaves in node_groups.items():
                    master_value = values[master]
                    slave_values = [values[x] for x in slaves]
                    if not self._values_match(master_value,
                                              slave_values):
                        if (not node_field_value_warnings and
                                not suppress_warnings):
                            self._warning(
                                'Node field values do not match.',
                                'Nodes are being merged but values at these '
                                'nodes for node field "%s" do not match.  An '
                                'averaged value will be used.\n'
                                '\n'
                                'Future warnings of this type will be '
                                'suppressed.' % (name))
                        node_field_value_warnings += 1
                        values[master] = (values[master] + sum(slave_values)
                                          ) / float(1 + len(slaves))
                del values[first_slave_index:]
        # change self.node_sets
        node_set_member_warnings = 0
        node_set_value_warnings = 0
        for node_set_id in self.get_node_set_ids():
            members = self.get_node_set_members(node_set_id)
            fields = self._get_node_set_fields(node_set_id)
            members_set = set(members)
            member_indices_to_delete = []
            for master, slaves in node_groups.items():
                master_included = master in members_set
                slaves_included = []
                slaves_not_included = []
                for slave in slaves:
                    if slave in members_set:
                        slaves_included.append(slave)
                    else:
                        slaves_not_included.append(slave)
                if not master_included and not slaves_included:
                    continue
                # warn if not all nodes are in the set
                if not master_included or slaves_not_included:
                    if not node_set_member_warnings and not suppress_warnings:
                        self._warning(
                            'Ambiguous merge of nodes.',
                            'Node are being merged, but only some of these '
                            'nodes belong to a given node set.  The operation '
                            'is therefore ambiguous on whether or not to '
                            'include the merged node in the set.  The merged '
                            'node will be included in the set.\n\nFuture '
                            'warnings of this type will be suppressed.')
                    node_set_member_warnings += 1
                # if master node is not included, steal the position of a
                # slave node
                if not master_included:
                    members[members.index(slaves_included[0])] = master
                    del slaves_included[0]
                    master_included = True
                    if not slaves_included:
                        continue
                master_index = members.index(master)
                slave_indices = [members.index(x) for x in slaves_included]
                # mark slaves to delete
                member_indices_to_delete.extend(slave_indices)
                # average values, warn if they are not the same
                for name, all_values in fields.items():
                    for values in all_values:
                        slave_values = [values[x] for x in slave_indices]
                        if not self._values_match(values[master_index],
                                                  slave_values):
                            if (not node_set_value_warnings and
                                    not suppress_warnings):
                                self._warning(
                                    'Node set field values do not match.',
                                    'Nodes are being merged but values at '
                                    'these nodes for node set field %s do '
                                    'not match.  An averaged values will '
                                    'be used.\n\nFuture warnings of this '
                                    'type will be suppressed.' % (name))
                            node_set_value_warnings += 1
                        new_value = (values[master_index] + sum(slave_values)
                                     ) / float(1 + len(slave_indices))
                        values[master_index] = new_value
            # delete slave members
            for index in sorted(member_indices_to_delete, reverse=True):
                del members[index]
                for all_values in fields.values():
                    for values in all_values:
                        del values[index]
        # delete those nodes
        del self.nodes[first_slave_index:]

    def _merge_node_pairs(self, node_pairs):
        """
        Merge the given node pairs.

        Example:
        >>> model._merge_node_pairs([(1, 3), (2, 5), (1, 5)])

        """
        # create groups of nodes to merge
        node_group = [None] * len(self.nodes)
        merge_group = []
        for pair in node_pairs:
            master, slave = sorted(pair)
            master_group = node_group[master]
            slave_group = node_group[slave]
            if master_group is None and slave_group is None:
                node_group[master] = len(merge_group)
                node_group[slave] = len(merge_group)
                merge_group.append([master, slave])
            elif master_group is None:
                node_group[master] = slave_group
                merge_group[slave_group].append(master)
            elif slave_group is None:
                node_group[slave] = master_group
                merge_group[master_group].append(slave)
            elif master_group != slave_group:
                for index in merge_group[slave_group]:
                    node_group[index] = master_group
                merge_group[master_group].extend(merge_group[slave_group])
                merge_group[slave_group] = []
        group_to_merge = dict()
        for group in merge_group:
            if not group:
                continue
            master = min(group)
            assert master not in group_to_merge
            group_to_merge[master] = sorted(set(group) - set([master]))
        self._merge_node_groups(group_to_merge)

    def delete_duplicate_elements(self, element_block_ids='all'):
        """
        Delete duplicate elements.

        For this calculation, a duplicate element is an element which shares
        all of its nodes with another element.

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        # go through each block and delete elements
        for id_ in element_block_ids:
            nodes_per_element = self.get_nodes_per_element(id_)
            element_count = self.get_element_count(id_)
            connectivity = self.get_connectivity(id_)
            assert len(connectivity) == nodes_per_element * element_count
            # sort elements into sets which contain
            elements = set()
            duplicates = []
            for x in xrange(0, element_count):
                start = x * nodes_per_element
                element = set(connectivity[start: start + nodes_per_element])
                element = tuple(sorted(element))
                if element in elements:
                    duplicates.append(x)
                else:
                    elements.add(element)
            self._delete_elements(id_, duplicates)

    def _find_close_nodes(self, tolerance):
        """
        Return groups of nodes that are close to one another.

        This will search for all nodes which are separated by at most
        'tolerance' and return a list of node groups.

        If node A is close to node B, and node B is close to node C, then each
        A, B, and C will belong to the same node group, regardless of whether
        or not node A is close to node C.

        Results are returned as a dict of node groups with the lowest index
        node being the key and the value as the list of slave nodes.  For
        example, if nodes 0 and 1 are close, and 2 and 3 are close, it would
        return:
        >>> model._find_close_nodes()
        {0: [1], 2: [3]}

        """
        # create a sorting vector
        sorting_vector = [math.pi, 2.0, 0.5 * math.sqrt(2)]
        scale = math.sqrt(sum([x * x for x in sorting_vector]))
        sorting_vector = [x / scale for x in sorting_vector]
        # create a list of node indices sorted by distance along that vector
        sorted_nodal_distances = sorted([(self._dot(x, sorting_vector), i)
                                         for i, x in enumerate(self.nodes)])
        # go through the list to find node groups
        node_count = len(self.nodes)
        lower_index = 0
        upper_index = 1
        master_nodes = range(len(self.nodes))
        while lower_index < node_count:
            if (upper_index >= node_count or
                    sorted_nodal_distances[upper_index][0] -
                    sorted_nodal_distances[lower_index][0] > tolerance):
                lower_index += 1
                upper_index = lower_index + 1
                continue
            lower_element = sorted_nodal_distances[lower_index]
            upper_element = sorted_nodal_distances[upper_index]
            if self._distance_between(
                    self.nodes[lower_element[1]],
                    self.nodes[upper_element[1]]) <= tolerance:
                one, two = sorted([lower_element[1], upper_element[1]])
                master_nodes[two] = one
            upper_index += 1
        # create a list of close node groups
        close_node_groups = dict()
        for index, master_node in enumerate(master_nodes):
            if index != master_node:
                while master_nodes[master_node] != master_node:
                    assert master_nodes[master_node] < master_node
                    master_node = master_nodes[master_node]
                if master_node not in close_node_groups:
                    close_node_groups[master_node] = []
                close_node_groups[master_node].append(index)
        # return the result
        return close_node_groups

    def merge_nodes(self,
                    tolerance=1e-6,
                    relative_tolerance=True,
                    suppress_warnings=False):
        """
        Merge nodes that are closer than the given tolerance.

        Node fields, node sets and node set fields are updated accordingly.
        Using a tolerance of 0 will merge nodes at the exact same location.

        If 'relative_tolerace=True', this will multiply the tolerance by the
        length scale of the model as obtained from 'get_length_scale()'.

        Examples:
        >>> model.merge_nodes(0)
        >>> model.merge_nodes()

        """
        # find the absolute tolerance
        if relative_tolerance:
            tolerance *= self.get_length_scale()
        # find the node groups
        close_node_groups = self._find_close_nodes(tolerance)
        # merge these node groups
        self._merge_node_groups(close_node_groups,
                                suppress_warnings=suppress_warnings)

    def _duplicate_nodes(self, node_indices, new_node_indices):
        """
        Create duplicates of the given nodes.

        This updates node fields, node sets, and node set fields consistently.
        If a node is part of a node set, its duplicated node will be added to
        the node set.

        """
        new_node_indices[:] = xrange(len(self.nodes),
                                     len(self.nodes) + len(node_indices))
        # update self.nodes
        new_nodes = [list(self.nodes[x]) for x in node_indices]
        self.create_nodes(new_nodes)
        # update self.node_fields
        for all_values in self.node_fields.values():
            for values in all_values:
                for master, slave in itertools.izip(node_indices,
                                                    new_node_indices):
                    values[slave] = values[master]
        # update self.node_sets
        for node_set_id in self.get_node_set_ids():
            members = self.get_node_set_members(node_set_id)
            fields = self._get_node_set_fields(node_set_id)
            for master, slave in itertools.izip(node_indices,
                                                new_node_indices):
                if master in members:
                    members.append(slave)
                    master_index = members.index(master)
                    for all_values in fields.values():
                        for values in all_values:
                            values.append(values[master_index])

    def _create_averaged_nodes(self, node_indices, new_node_indices):
        """
        Create duplicates nodes which are averaged from the given nodes.

        If all nodes which are being averaged are in a given node set, the
        newly created node is also added the node set.

        For example, 'node_indices = [(0, 1, 2)]' will create a new node which
        is the average of nodes 1, 2, and 3.  Passing
        'node_indices = [((0, 0.25), (1, 0.75))]' will create a node which is
        a weighted combination of nodes 0 and 1.

        """
        first_node_index = len(self.nodes)
        new_node_indices[:] = range(len(self.nodes),
                                    len(self.nodes) + len(node_indices))
        # format node_indices such that each value is of the given form:
        # ((index, weight), (index2, weight2), ...)
        new_node_indices = []
        for index_list in node_indices:
            new_index_list = []
            if isinstance(index_list, int):
                new_index_list = [[index_list, 1.0]]
            else:
                length = float(len(index_list))
                for specifier in index_list:
                    if isinstance(specifier, int):
                        new_index_list.append([specifier, 1.0 / length])
                    else:
                        new_index_list.append(specifier)
            new_node_indices.append(new_index_list)
        node_indices = new_node_indices
        # ensure weights add up to 1
        for index_list in node_indices:
            total = sum(x[1] for x in index_list)
            if abs(total - 1.0) > 1e-14:
                self._bug('Incorrect weights.',
                          'The given node averaging weights do not add up '
                          'to 1.')
        # create the new nodes
        first_new_node_index = len(self.nodes)
        new_nodes = []
        for index_list in node_indices:
            this_new_node = [sum(weight * self.nodes[index][d]
                                 for index, weight in index_list)
                             for d in xrange(3)]
            new_nodes.append(this_new_node)
        self.create_nodes(new_nodes)
        # update self.node_fields
        for all_values in self.node_fields.values():
            for values in all_values:
                new_values = [sum(weight * values[index]
                                  for index, weight in index_list)
                              for index_list in node_indices]
                values[first_node_index:] = new_values
        # update self.node_sets
        for node_set_id in self.get_node_set_ids():
            members = self.get_node_set_members(node_set_id)
            fields = self._get_node_set_fields(node_set_id)
            members_set = set(members)
            index_from_member_index = dict(
                (node_index, index)
                for index, node_index in enumerate(members))
            for node_index, index_list in enumerate(node_indices):
                check = set(x in members_set for x, _ in index_list)
                if False in check:
                    continue
                # okay, we need to add this node to the set
                members.append(first_new_node_index + node_index)
                # and calculate its value for all fields
                for all_values in fields.values():
                    for values in all_values:
                        new_value = sum(
                            weight * values[index_from_member_index[index]]
                            for index, weight in index_list)
                        values.append(new_value)

    def unmerge_element_blocks(self, element_block_ids='all'):
        """
        Duplicate nodes to unmerge element blocks.

        For elements blocks that share one or more nodes, duplicate these nodes
        and unmerge the element blocks.  For example, if element block A and B
        share node 1, that node will be duplicated, block A will use the
        original node and block B will use the duplicate.

        Node fields, node sets, and node set fields are updated accordingly.

        Example:
        >>> model.unmerge_element_blocks()

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        if len(element_block_ids) <= 1:
            return
        # Separate each pair of blocks.  Keep nodes on block id1 the same and
        # create new nodes for block id2.  Update node fields.  Update node
        # sets such that if a shared node is part of the set, the duplicated
        # node is also part of the set.
        for id1, id2 in itertools.combinations(element_block_ids, 2):
            shared_nodes_set = set.intersection(
                set(self.get_nodes_in_element_block(id1)),
                set(self.get_nodes_in_element_block(id2)))
            shared_nodes = sorted(shared_nodes_set)
            if not shared_nodes:
                continue
            # on block id2, old node x becomes new node
            # node_offset + shared_nodes.index(x)
            # create the new nodes
            new_shared_nodes = []
            self._duplicate_nodes(shared_nodes, new_shared_nodes)
            # update element block 2 connectivity
            connectivity = self.get_connectivity(id2)
            for index, node_index in enumerate(connectivity):
                if node_index in shared_nodes_set:
                    connectivity[index] = new_shared_nodes[
                        shared_nodes.index(node_index)]
            # blocks should no longer share nodes
            assert not set.intersection(
                set(self.get_nodes_in_element_block(id1)),
                set(self.get_nodes_in_element_block(id2)))

    def import_model(self,
                     filename,
                     element_block_ids='all',
                     timesteps='all',
                     node_field_names='all',
                     element_field_names='all',
                     side_set_ids='all',
                     node_set_ids='all',
                     global_variable_names='all',
                     node_set_field_names='all',
                     side_set_field_names='all'):
        """
        Import information from an ExodusII file.

        This will add to the current model in memory, so multiple calls will
        act to merge files.  As a shortcut, one may use the
        'exomerge.import_model' to create a new model from a file.
        >>> model = exomerge.import_model('output.e')

        By default, this will import all information from the given ExodusII
        results file.  To import only part of a mesh file, or to load only a
        particular timestep, one can use the following options.
        >>> model.import_model('output.e', element_block_ids=[1, 2])
        >>> model.import_model('output.e', timesteps='last')
        >>> model.import_model('output.e', node_field_names='disp_*')

        To import only the mesh without any field information, one can use the
        following syntax.
        >>> model.import_model('output.e', timesteps='none')

        Example:
        >>> model.import_model('mesh_file.g')
        >>> model.import_model('results_file.e')

        """
        # open the file for read
        if not os.path.isfile(filename):
            self._error('File not found',
                        'The specified file "%s" was not found.' % filename)
        if SUPPRESS_EXODUS_OUTPUT:
            save_stdout = sys.stdout
            sys.stdout = DummyFile()
        exodus_file = exodus.exodus(filename, mode='r')
        if SUPPRESS_EXODUS_OUTPUT:
            sys.stdout = save_stdout
        # format timesteps to retrieve
        file_timesteps = list(exodus_file.get_times())
        if timesteps == 'last_if_any':
            if file_timesteps:
                timesteps = 'last'
            else:
                timesteps = 'none'
        timesteps = self._format_id_list(
            timesteps,
            file_timesteps,
            'timestep')
        # format block id list
        file_element_block_ids = list(exodus_file.get_elem_blk_ids())
        file_element_block_names = []
        if file_element_block_ids:
            file_element_block_names = list(exodus_file.get_elem_blk_names())
        element_block_ids = self._format_id_list(
            element_block_ids,
            sorted(file_element_block_ids),
            'element block')
        # format other id lists
        node_field_names = self._format_id_list(
            node_field_names,
            sorted(list(exodus_file.get_node_variable_names())),
            'node field')
        file_element_field_names = list(
            exodus_file.get_element_variable_names())
        element_field_names = self._format_id_list(
            element_field_names,
            sorted(file_element_field_names),
            'element field')
        file_side_set_ids = list(exodus_file.get_side_set_ids())
        file_side_set_names = []
        if file_side_set_ids:
            file_side_set_names = list(exodus_file.get_side_set_names())
        side_set_ids = self._format_id_list(
            side_set_ids,
            sorted(file_side_set_ids),
            'side set')
        file_node_set_ids = list(exodus_file.get_node_set_ids())
        file_node_set_names = []
        if file_node_set_ids:
            file_node_set_names = list(exodus_file.get_node_set_names())
        node_set_ids = self._format_id_list(
            node_set_ids,
            sorted(file_node_set_ids),
            'node set')
        global_variable_names = self._format_id_list(
            global_variable_names,
            sorted(list(exodus_file.get_global_variable_names())),
            'global variable')
        file_node_set_field_names = list(
            exodus_file.get_node_set_variable_names())
        node_set_field_names = self._format_id_list(
            node_set_field_names,
            sorted(file_node_set_field_names),
            'node set field')
        file_side_set_field_names = list(
            exodus_file.get_side_set_variable_names())
        side_set_field_names = self._format_id_list(
            side_set_field_names,
            sorted(file_side_set_field_names),
            'node set field')
        # ensure element blocks in this file do not already exist
        for element_block_id in element_block_ids:
            if self.element_block_exists(element_block_id):
                self._error(
                    'Element block already exists.',
                    'Cannot import element block \"%d\" since it already '
                    'exists in the model.' % element_block_id)
        # create new nodes used in this file
        node_offset = len(self.nodes)
        # store nodes we are going to import
        # if we're importing all element blocks, then import all nodes
        # else only import nodes in the element blocks we're choosing
        # note that if there are no element blocks, this means we import all
        # nodes
        # note: new_used_nodes is 1-based
        if element_block_ids == sorted(file_element_block_ids):
            new_used_nodes = range(1, exodus_file.num_nodes() + 1)
        else:
            new_used_nodes = [False] * exodus_file.num_nodes()
            # add nodes in blocks we are importing
            for element_block_id in element_block_ids:
                connectivity = exodus_file.get_elem_connectivity(
                    element_block_id)[0]
                for i in connectivity:
                    new_used_nodes[i - 1] = True
            # save indices for nodes we want to import
            new_used_nodes = [i + 1
                              for i, x in enumerate(new_used_nodes)
                              if x]
        # get new node coordinates
        new_nodes = []
        file_coords = [list(x) for x in exodus_file.get_coords()]
        for i in new_used_nodes:
            new_nodes.append([file_coords[x][i - 1] for x in xrange(3)])
        # create node mapping
        # node i in the file refers to node node_map[i] in the model
        node_map = [None] * (exodus_file.num_nodes() + 1)
        for i in xrange(len(new_used_nodes)):
            node_map[new_used_nodes[i]] = i + node_offset
        # create new nodes
        self.create_nodes(new_nodes)
        # find element mapping
        # element i in the file refers to
        # - element block id file_block_from_element[i]
        # - local element index file_index_from_element[i]
        # in the model
        elements_in_block = []
        for i in file_element_block_ids:
            elements_in_block.append(exodus_file.elem_blk_info(i)[1])
        file_block_from_element = [None]
        file_index_from_element = [None]
        for i in xrange(len(file_element_block_ids)):
            if file_element_block_ids[i] in element_block_ids:
                file_block_from_element.extend(
                    [file_element_block_ids[i]] * elements_in_block[i])
                file_index_from_element.extend(xrange(elements_in_block[i]))
            else:
                file_block_from_element.extend([None] * elements_in_block[i])
                file_index_from_element.extend([None] * elements_in_block[i])
        # create new element blocks
        for element_block_id in element_block_ids:
            new_connectivity = list(exodus_file.get_elem_connectivity(
                element_block_id))[0]
            # apply node map
            new_connectivity = [node_map[x] for x in new_connectivity]
            # get the element block name
            element_block_name = file_element_block_names[
                file_element_block_ids.index(element_block_id)]
            # create the element block
            self.create_element_block(
                element_block_id,
                list(exodus_file.elem_blk_info(element_block_id)),
                connectivity=new_connectivity)
            if element_block_name:
                self.rename_element_block(element_block_id, element_block_name)
        # get indices of each timestep and create new timesteps
        # list of (timestep_index_in_file, timestep_index_in_model)
        timestep_indices = []
        duplicate_timesteps = False
        for timestep in timesteps:
            if not self.timestep_exists(timestep):
                self.create_timestep(timestep)
            else:
                duplicate_timesteps = True
            timestep_indices.append(
                (file_timesteps.index(timestep) + 1,
                 self._get_internal_timestep_index(timestep)))
        # get list of timestep indices present in the model but not included
        # in the import
        included_timestep_indices = [x[1] for x in timestep_indices]
        excluded_timestep_indices = []
        for i in xrange(len(self.timesteps)):
            if i not in included_timestep_indices:
                excluded_timestep_indices.append(i)
        # create new node fields
        for node_field_name in node_field_names:
            if not self.node_field_exists(node_field_name):
                self.create_node_field(node_field_name)
        # populate node field info
        for node_field_name in node_field_names:
            for timestep_index in timestep_indices:
                file_node_field_values = list(
                    exodus_file.get_node_variable_values(
                        node_field_name,
                        timestep_index[0]))
                model_values = self.node_fields[
                    node_field_name][timestep_index[1]]
                for i in new_used_nodes:
                    model_values[node_map[i]] = file_node_field_values[i - 1]
        # populate node sets
        for node_set_id in node_set_ids:
            file_node_set_members = exodus_file.get_node_set_nodes(
                node_set_id)
            model_node_set_members = []
            # add nodes which were included
            for i in file_node_set_members:
                if node_map[i] is not None:
                    model_node_set_members.append(node_map[i])
            # create the node set or add to an existing set
            if not self.node_set_exists(node_set_id):
                self.create_node_set(node_set_id, model_node_set_members)
            else:
                self.add_nodes_to_node_set(node_set_id,
                                           model_node_set_members)
            node_set_name = file_node_set_names[
                file_node_set_ids.index(node_set_id)]
            if node_set_name:
                self.rename_node_set(node_set_id, node_set_name)
        # populate side sets
        for side_set_id in side_set_ids:
            file_side_set_members = exodus_file.get_side_set(
                side_set_id)
            model_side_set_members = []
            # add nodes which were included
            for element_index, side_number in zip(*file_side_set_members):
                if file_block_from_element[element_index] is not None:
                    model_side_set_members.append(
                        (file_block_from_element[element_index],
                         file_index_from_element[element_index],
                         side_number - 1))
            # create the side set or add to an existing set
            if not self.side_set_exists(side_set_id):
                self.create_side_set(side_set_id,
                                     model_side_set_members)
            else:
                self.add_faces_to_side_set(side_set_id,
                                           model_side_set_members)
            side_set_name = file_side_set_names[
                file_side_set_ids.index(side_set_id)]
            if side_set_name:
                self.rename_side_set(side_set_id, side_set_name)
        # store truth table for node set field info
        if node_set_field_names:
            file_node_set_truth_table = []
            for file_node_set_id in file_node_set_ids:
                file_node_set_truth_table.append(
                    exodus_file.get_node_set_variable_truth_table(
                        file_node_set_id))
        # populate node set fields
        for node_set_field_name in node_set_field_names:
            for node_set_id in node_set_ids:
                # don't process if field does not exist in file
                field_exists = file_node_set_truth_table[
                    file_node_set_ids.index(node_set_id)][
                        file_node_set_field_names.index(
                            node_set_field_name)]
                if not field_exists:
                    continue
                # process each included timestep
                model_values = [None] * len(self.timesteps)
                for timestep_index in timestep_indices:
                    file_values = list(
                        exodus_file.get_node_set_variable_values(
                            node_set_id,
                            node_set_field_name,
                            timestep_index[0]))
                    model_values[timestep_index[1]] = file_values
                # add default field value to excluded timestep
                for timestep_index in excluded_timestep_indices:
                    default_value = self._get_default_field_value(
                        node_set_field_name)
                    node_count = len(self.get_node_set_members(node_set_id))
                    field = [default_value] * node_count
                    model_values[timestep_index] = field
                # assign all values
                fields = self._get_node_set_fields(node_set_id)
                fields[node_set_field_name] = model_values
        # store truth table for side set field info
        if side_set_field_names:
            file_side_set_truth_table = []
            for file_side_set_id in file_side_set_ids:
                file_side_set_truth_table.append(
                    exodus_file.get_side_set_variable_truth_table(
                        file_side_set_id))
        # populate side set fields
        for side_set_field_name in side_set_field_names:
            for side_set_id in side_set_ids:
                # don't process if field does not exist in file
                field_exists = file_side_set_truth_table[
                    file_side_set_ids.index(side_set_id)][
                        file_side_set_field_names.index(
                            side_set_field_name)]
                if not field_exists:
                    continue
                # process each included timestep
                model_values = [None] * len(self.timesteps)
                for timestep_index in timestep_indices:
                    file_values = list(
                        exodus_file.get_side_set_variable_values(
                            side_set_id,
                            side_set_field_name,
                            timestep_index[0]))
                    model_values[timestep_index[1]] = file_values
                # add default field value to excluded timestep
                for timestep_index in excluded_timestep_indices:
                    default_value = self._get_default_field_value(
                        side_set_field_name)
                    side_count = len(self.get_side_set_members(side_set_id))
                    field = [default_value] * side_count
                    model_values[timestep_index] = field
                # assign all values
                fields = self._get_side_set_fields(side_set_id)
                fields[side_set_field_name] = model_values
        # store truth table for element field info
        if element_field_names:
            file_truth_table = []
            for file_element_block_id in file_element_block_ids:
                file_truth_table.append(
                    exodus_file.get_element_variable_truth_table(
                        file_element_block_id))
        # populate element fields
        for element_field_name in element_field_names:
            for element_block_id in element_block_ids:
                # don't process if field does not exist in file
                field_exists = file_truth_table[
                    file_element_block_ids.index(element_block_id)][
                        file_element_field_names.index(element_field_name)]
                if not field_exists:
                    continue
                # process each included timestep
                model_values = [None] * len(self.timesteps)
                for timestep_index in timestep_indices:
                    file_values = list(
                        exodus_file.get_element_variable_values(
                            element_block_id,
                            element_field_name,
                            timestep_index[0]))
                    model_values[timestep_index[1]] = file_values
                # add default field value to excluded timestep
                for timestep_index in excluded_timestep_indices:
                    default_value = self._get_default_field_value(
                        element_field_name)
                    element_count = self.get_element_count(element_block_id)
                    field = [default_value] * element_count
                    model_values[timestep_index] = field
                # assign all values
                fields = self._get_element_block_fields(element_block_id)
                fields[element_field_name] = model_values
        # get global variables
        for global_variable_name in global_variable_names:
            # create the variable
            if not self.global_variable_exists(global_variable_name):
                self.create_global_variable(global_variable_name)
            else:
                if duplicate_timesteps:
                    self._exists_warning(global_variable_name,
                                         'global variable')
            # get values
            model_values = self.global_variables[global_variable_name]
            for timestep_index in timestep_indices:
                file_value = exodus_file.get_global_variable_value(
                    global_variable_name, timestep_index[0])
                model_values[timestep_index[1]] = file_value
        # add info records
        self.info_records += exodus_file._ex_get_info_recs_quietly()
        # add qa records
        self.qa_records += exodus_file.get_qa_records()
        # add title if one does not already exist
        # else add it to an info record
        if not self.title:
            self.title = exodus_file.title()
        else:
            self.info_records.append('Discarded title from the '
                                     'following file:')
            #split filename string if filename is larger than 79 characters
            filename_wrap_list = textwrap.fill(filename, width=79).split('\n')
            #append interpreter continuation char "\\" to end of continuation line while splitting
            for i in range(len(filename_wrap_list)-1):
                filename_wrap_list[i] += "\\"
            #append multiple split list to records
            self.info_records.extend(filename_wrap_list)
            self.info_records.append(exodus_file.title())
        # run a check on the model to ensure arrays are correct sizes
        self._verify()
        # close the file
        if SUPPRESS_EXODUS_OUTPUT:
            save_stdout = sys.stdout
            sys.stdout = DummyFile()
        exodus_file.close()
        if SUPPRESS_EXODUS_OUTPUT:
            sys.stdout = save_stdout

    def export_model(self,
                     filename='output_exomerge.e',
                     element_block_ids='all',
                     timesteps='all',
                     side_set_ids='all',
                     node_set_ids='all',
                     global_variable_names='auto',
                     node_field_names='auto',
                     element_field_names='auto',
                     side_set_field_names='auto',
                     node_set_field_names='auto'):
        """
        Export the current model to an ExodusII file.

        Examples:
        >>> model.export_model('output.g')

        """
        # verify information is valid
        self._verify()
        # format subset of data to export
        element_block_ids = self._format_element_block_id_list(
            element_block_ids)
        side_set_ids = self._format_side_set_id_list(side_set_ids)
        node_set_ids = self._format_node_set_id_list(node_set_ids)
        timesteps = self._format_id_list(
            timesteps,
            self.get_timesteps(),
            'timestep')
        # If no timesteps are exported, or if no nodes are defined, then no
        # fields can be exported.
        if element_field_names == 'auto':
            if timesteps:
                element_field_names = self.get_element_field_names(
                    element_block_ids)
            else:
                element_field_names = 'none'
        if global_variable_names == 'auto':
            if timesteps:
                global_variable_names = 'all'
            else:
                global_variable_names = 'none'
        if node_field_names == 'auto':
            if timesteps and self.nodes:
                node_field_names = 'all'
            else:
                node_field_names = 'none'
        if node_set_field_names == 'auto':
            if timesteps and self.nodes:
                node_set_field_names = self.get_node_set_field_names(
                    node_set_ids)
            else:
                node_set_field_names = 'none'
        if side_set_field_names == 'auto':
            if timesteps:
                side_set_field_names = self.get_side_set_field_names(
                    side_set_ids)
            else:
                side_set_field_names = 'none'
        # format each list
        global_variable_names = self._format_id_list(
            global_variable_names,
            self.get_global_variable_names(),
            'global variable')
        element_field_names = self._format_id_list(
            element_field_names,
            self.get_element_field_names(),
            'element field')
        node_field_names = self._format_id_list(
            node_field_names,
            self.get_node_field_names(),
            'node field')
        node_set_field_names = self._format_id_list(
            node_set_field_names,
            self.get_node_set_field_names(),
            'node set field')
        side_set_field_names = self._format_id_list(
            side_set_field_names,
            self.get_side_set_field_names(),
            'side set field')
        # delete the file if it exists
        if os.path.isfile(filename):
            os.remove(filename)
        # create file
        if SUPPRESS_EXODUS_OUTPUT:
            save_stdout = sys.stdout
            sys.stdout = DummyFile()
        new_file = exodus.exodus(filename,
                                 'w',
                                 'ctype',
                                 self.title,
                                 3,
                                 len(self.nodes),
                                 self.get_element_count(element_block_ids),
                                 len(element_block_ids),
                                 len(node_set_ids),
                                 len(side_set_ids))
        if SUPPRESS_EXODUS_OUTPUT:
            sys.stdout = save_stdout
        # put timesteps into chronological order
        timestep_indices = [(self.timesteps.index(x), index + 1)
                            for index, x in enumerate(timesteps)]
        # write times
        for index, timestep in enumerate(timesteps):
            new_file.put_time(index + 1, timestep)
        # write global variables
        new_file.set_global_variable_number(len(global_variable_names))
        for index, name in enumerate(global_variable_names):
            new_file.put_global_variable_name(name, index + 1)
            values = self.global_variables[name]
            for timestep_index in timestep_indices:
                new_file.put_global_variable_value(
                    name,
                    timestep_index[1],
                    values[timestep_index[0]])
        # write nodes
        new_file.put_coords([x[0] for x in self.nodes],
                            [x[1] for x in self.nodes],
                            [x[2] for x in self.nodes])
        # write node fields
        new_file.set_node_variable_number(len(node_field_names))
        for index, name in enumerate(node_field_names):
            new_file.put_node_variable_name(name, index + 1)
            values = self.node_fields[name]
            for timestep_index in timestep_indices:
                new_file.put_node_variable_values(
                    name,
                    timestep_index[1],
                    values[timestep_index[0]])
        # write element blocks
        for id_ in element_block_ids:
            name, info, connectivity, fields = self.element_blocks[id_]
            new_file.put_elem_blk_info(id_, *info)
            # connectivity in file must be 1-based
            temp_connectivity = [x + 1 for x in connectivity]
            new_file.put_elem_connectivity(id_, temp_connectivity)
            if name:
                new_file.put_elem_blk_name(id_, name)
        # write element fields
        new_file.set_element_variable_number(len(element_field_names))
        for index, name in enumerate(element_field_names):
            new_file.put_element_variable_name(name, index + 1)
        if element_field_names:
            truth_table = self._create_element_field_truth_table(
                element_block_ids, element_field_names)
            new_file.set_element_variable_truth_table(truth_table)
        for block_id in element_block_ids:
            fields = self._get_element_block_fields(block_id)
            for name in element_field_names:
                if name not in fields:
                    continue
                field = fields[name]
                for timestep_index in timestep_indices:
                    new_file.put_element_variable_values(
                        block_id,
                        name,
                        timestep_index[1],
                        field[timestep_index[0]])
        # get first element in each block
        element_count = [self.get_element_count(id_)
                         for id_ in element_block_ids]
        first_element = [sum(element_count[:i]) + 1
                         for i in xrange(len(element_count))]
        first_element = dict(
            (element_block_ids[i], first_element[i])
            for i in xrange(len(element_block_ids)))
        # write side sets
        for id_ in side_set_ids:
            name = self.get_side_set_name(id_)
            members = self.get_side_set_members(id_)
            new_file.put_side_set_params(id_, len(members), 0)
            if members:
                elements = [first_element[block_id] + index
                            for block_id, index, _ in members]
                sides = [x[2] + 1 for x in members]
                new_file.put_side_set(id_, elements, sides)
            if name:
                new_file.put_side_set_name(id_, name)
        # write side set fields
        new_file.set_side_set_variable_number(len(side_set_field_names))
        for index, name in enumerate(side_set_field_names):
            new_file.put_side_set_variable_name(name, index + 1)
        if side_set_field_names:
            truth_table = self._create_side_set_field_truth_table(
                side_set_ids, side_set_field_names)
            new_file.set_side_set_variable_truth_table(truth_table)
        for side_set_id in side_set_ids:
            members = self.get_side_set_members(side_set_id)
            fields = self._get_side_set_fields(side_set_id)
            if not members:
                continue
            for name in side_set_field_names:
                if name not in fields:
                    continue
                field = fields[name]
                for timestep_index in timestep_indices:
                    new_file.put_side_set_variable_values(
                        side_set_id,
                        name,
                        timestep_index[1],
                        field[timestep_index[0]])
        # write node sets
        for id_ in node_set_ids:
            name = self.get_node_set_name(id_)
            members = self.get_node_set_members(id_)
            new_file.put_node_set_params(id_, len(members))
            if members:
                # convert to 1-based indexing
                temp_members = [x + 1 for x in members]
                new_file.put_node_set(id_, temp_members)
            if name:
                new_file.put_node_set_name(id_, name)
        # write node set fields
        new_file.set_node_set_variable_number(len(node_set_field_names))
        for index, name in enumerate(node_set_field_names):
            new_file.put_node_set_variable_name(name, index + 1)
        if node_set_field_names:
            truth_table = self._create_node_set_field_truth_table(
                node_set_ids, node_set_field_names)
            new_file.set_node_set_variable_truth_table(truth_table)
        for node_set_id in node_set_ids:
            members = self.get_node_set_members(node_set_id)
            fields = self._get_node_set_fields(node_set_id)
            if not members:
                continue
            for name in node_set_field_names:
                if name not in fields:
                    continue
                field = fields[name]
                for timestep_index in timestep_indices:
                    new_file.put_node_set_variable_values(
                        node_set_id,
                        name,
                        timestep_index[1],
                        field[timestep_index[0]])
        # write info records
        new_file._exodus__ex_put_info_recs(self.info_records)
        # write qa records (and append one for this program)
        self.qa_records.append(self._get_qa_record())
        new_file.put_qa_records(self.qa_records)
        del self.qa_records[-1]
        # close file
        if SUPPRESS_EXODUS_OUTPUT:
            save_stdout = sys.stdout
            sys.stdout = DummyFile()
        new_file.close()
        if SUPPRESS_EXODUS_OUTPUT:
            sys.stdout = save_stdout

    def get_side_set_area(self, side_set_ids):
        """
        Return the total area of the given side sets.

        Example:
        >>> model.get_side_set_area('all')

        """
        side_set_ids = self._format_side_set_id_list(side_set_ids)
        # we do this by first creating one or more new element block out of
        # the faces within the side sets and then calculating their "volume."
        new_blocks = self._create_element_blocks_from_side_sets(side_set_ids)
        total_area = 0.0
        for id_ in new_blocks:
            total_area += self.get_element_block_volume(id_)
            self.delete_element_block(id_,
                                      delete_orphaned_nodes=False)
        return total_area

    def _get_thickness_from_volume_and_area(self, volume, area):
        """Return the thickness of a disk given the volume and area."""
        # max phi (for a sphere) is 6^(-1/3) * pi^(-1/6)
        phi = math.pow(volume, 1 / 3.0) / math.pow(area, 1 / 2.0)
        # find a solution, if possible
        max_alpha = math.pi ** (-1.0 / 6) * 6.0 ** (-1.0 / 3)
        low = 0.0
        high = max_alpha
        a = 1.0
        b = 1.5 - 1 / (16 * math.pi * phi ** 6)
        c = 0.75
        d = 0.125
        high_value = a * high ** 3 + b * high ** 2 + c * high + d
        for _ in range(53):
            mid = (low + high) / 2
            mid_value = a * high ** 3 + b * high ** 2 + c * high + d
            if (high_value > 0) == (mid_value > 0):
                high_value = mid_value
                high = mid
            else:
                low = mid
        # in the case this fails, return NaN
        if low == 0 or high == max_alpha:
            return float('nan')
        alpha = (low + high) / 2.0
        height = (volume * 4 * alpha ** 2 / math.pi) ** (1.0 / 3)
        return height

    def _calculate_element_block_thickness(self, element_block_ids):
        """
        Return the approximate thickness of the given element blocks.

        This approximates the thickness by calculating the volume and exposed
        surface area of the element blocks, then returning the thickness of a
        short disk with the same characteristics.

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids,
            empty_list_okay=False)
        # delete any element blocks which are not of dimension 3
        new_ids = []
        for id_ in element_block_ids:
            if self.get_element_block_dimension(id_) != 3:
                self._warning('Unexpected element block dimension',
                              'Element blocks used to determine element edge '
                              'length are expected to be of dimension 3.')
            else:
                new_ids.append(id_)
        # get the volume
        volume = self.get_element_block_volume(element_block_ids)
        # find the external faces and put them into a side set
        external_faces = self._get_external_element_faces(element_block_ids)
        # find all face types
        side_set_id = self._new_side_set_id()
        self.create_side_set(side_set_id, external_faces)
        area = self.get_side_set_area(side_set_id)
        self.delete_side_set(side_set_id)
        return self._get_thickness_from_volume_and_area(volume, area)

    def get_element_block_extents(self, element_block_ids='all'):
        """
        Return the extents of the element blocks as a list.

        The results are returned in the following format:
        [[min_x, max_x], [min_y, max_y], [min_z, max_z]]

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids,
            empty_list_okay=False)
        # get a set of all nodes within the given element blocks
        all_nodes = set()
        for id_ in element_block_ids:
            connectivity = self.get_connectivity(id_)
            all_nodes.update(set(connectivity))
        # find the extents of that set
        extents = [[min(self.nodes[x][d] for x in all_nodes),
                    max(self.nodes[x][d] for x in all_nodes)]
                   for d in range(3)]
        return extents

    def _create_element_blocks_from_side_sets(self, side_set_ids):
        """
        Create new element block from the members of the given side sets.

        A list of new element block ids is returned.  Typically, only one
        ID is in the list, but since multiple face types can be present within
        the sideset, and element blocks can only contain a single element type,
        it is sometimes necessary to make more than one element block.

        """
        side_set_ids = self._format_side_set_id_list(side_set_ids)
        # create a dict for new face elements to create
        # e.g. new_elements['quad4'] = [[1, 2, 3, 4], [5, 6, 7, 8], ...]
        new_elements = {}
        for id_ in side_set_ids:
            members = self.get_side_set_members(id_)
            these_members = self._order_element_faces_by_block(members)
            for element_block_id, members in these_members.items():
                connectivity = self.get_connectivity(element_block_id)
                nodes_per_element = self.get_nodes_per_element(
                    element_block_id)
                face_mapping = self._get_face_mapping_from_id(element_block_id)
                for element_index, face_index in members:
                    face_type = face_mapping[face_index][0]
                    if face_type not in new_elements:
                        new_elements[face_type] = []
                    local_node = connectivity[
                        element_index * nodes_per_element:
                        (element_index + 1) * nodes_per_element]
                    new_elements[face_type].extend(
                        [local_node[x]
                         for x in face_mapping[face_index][1]])
        # now that we have the face elements in a local format, we create
        # the elements
        new_block_ids = []
        for element_type, connectivity in new_elements.items():
            new_id = self._new_element_block_id()
            self._assert(element_type in self.NODES_PER_ELEMENT)
            nodes_per_element = self.NODES_PER_ELEMENT[element_type]
            info = [element_type,
                    len(connectivity) / nodes_per_element,
                    nodes_per_element,
                    0]
            self.create_element_block(new_id, info, connectivity)
            new_block_ids.append(new_id)
        return new_block_ids

    def _get_element_edge_indices(self, element_type):
        """
        Return a list of element edges for the given element type.

        The returned list is a list of 2-tuples of local element node indices.

        """
        element_type = self._get_standard_element_type(element_type)
        # if not a standard type, then we don't know the edge indices
        if not self._is_standard_element_type(element_type):
            return []
        # if dimension is zero, it doesn't have edges
        if self._get_dimension(element_type) == 0:
            return []
        # create a mock element
        elements = dict()
        elements[element_type] = [range(self.NODES_PER_ELEMENT[element_type])]
        iterations = self._get_dimension(element_type) - 1
        # iterate until dimensions are zero
        for _ in range(iterations):
            new_elements = dict()
            for element_type, connectivities in elements.items():
                face_mapping = self._get_face_mapping(element_type)
                for face_type, indices in face_mapping:
                    if face_type not in new_elements:
                        new_elements[face_type] = []
                    for local_nodes in connectivities:
                        new_elements[face_type].append([local_nodes[x]
                                                        for x in indices])
            elements = new_elements
        # now find the endpoints using the volume formula
        edges = []
        for element_type, values in elements.items():
            assert element_type in self.VOLUME_FORMULA
            assert self.DIMENSION[element_type] == 1
            formula = self.VOLUME_FORMULA[element_type][-1]
            for local_nodes in values:
                edges.append([local_nodes[x] for x in formula])
        # now get the unique edges
        unique_edges = set(tuple(sorted(x)) for x in edges)
        return unique_edges

    def get_element_edge_length_info(self, element_block_ids='all'):
        """
        Return the minimum and average element edge lengths.

        Only edges within element in the specified element blocks are counted.
        All element blocks specified must be 3-dimensional.

        The result is returned as [minimum, average].

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids,
            empty_list_okay=False)
        minimum = sys.float_info.max
        total = 0.0
        edge_count = 0
        for element_block_id in element_block_ids:
            # get the edge endpoint info
            endpoints = self._get_element_edge_indices(
                self._get_element_type(element_block_id))
            # form all element edges
            element_count = self.get_element_count(element_block_id)
            connectivity = self.get_connectivity(element_block_id)
            nodes_per_element = self.get_nodes_per_element(element_block_id)
            edge_count += element_count * len(endpoints)
            for element_index in xrange(element_count):
                local_node = connectivity[
                    element_index * nodes_per_element:
                    (element_index + 1) * nodes_per_element]
                for edge in endpoints:
                    this_distance = self._distance_between(
                            self.nodes[local_node[edge[0]]],
                            self.nodes[local_node[edge[1]]])
                    total += this_distance
                    if this_distance < minimum:
                        minimum = this_distance
        if edge_count == 0:
            return [float('nan')] * 2
        return [minimum, total / edge_count]

    def _get_closest_point_distance_brute(self, points):
        """Return the distance between the two closest points."""
        point_count = len(points)
        dist = sys.float_info.max
        for i in range(1, point_count):
            for j in range(i):
                dist = min(dist, self._distance_between(points[i],
                                                        points[j]))
        return dist

    def _get_closest_point_distance(self, points, sorting_coord=0, trials=0):
        """Return the distance between the two closest points."""
        # handle cases of 3 or less points by brute force
        point_count = len(points)
        if point_count <= 3:
            return self._get_closest_point_distance_brute(points)
        # find min and max values for the sorting coord
        low = min(x[sorting_coord] for x in points)
        high = max(x[sorting_coord] for x in points)
        # if min and max values are the same, change the coord and try again
        mid = (low + high) / 2.0
        if low == mid or mid == high:
            if trials >= 2:
                return self._get_closest_point_distance_brute(points)
            else:
                return self._get_closest_point_distance(
                    points,
                    (sorting_coord + 1) % 3,
                    trials + 1)
        # sort into two lists
        low_points = [x
                      for x in points
                      if x[sorting_coord] <= mid]
        high_points = [x
                       for x in points
                       if x[sorting_coord] > mid]
        assert len(low_points) < point_count
        assert len(high_points) < point_count
        # find closest pair within each set
        dist = min(self._get_closest_point_distance(low_points),
                   self._get_closest_point_distance(high_points))
        del low_points
        del high_points
        # find points sufficiently close to centerline
        mid_low = mid - dist
        mid_high = mid - dist
        mid_points = [x
                      for x in points
                      if x[sorting_coord] >= mid_low and
                      x[sorting_coord] <= mid_high]
        if len(mid_points) == point_count:
            return self._get_closest_point_distance_brute(points)
        dist = min(dist,
                   self._get_closest_point_distance(mid_points,
                                                    (sorting_coord + 1) % 3))
        return dist

    def get_closest_node_distance(self):
        """Return the distance between the two closest nodes."""
        # create list of all nodes
        point_list = list(tuple(x) for x in self.nodes)
        return self._get_closest_point_distance(point_list)

    def _input_check_error(self, argument, format):
        """
        Report an error that the argument is not of the expected format.

        This is a helper function which gets called by _input_check in several
        places.

        """
        argument_string = '%s' % argument
        if len(argument_string) > 50:
            argument_string = argument_string[:47] + '...'
        format_string = '%s' % format
        if len(format_string) > 50:
            format_string = format_string[:47] + '...'
        self._error('Unexpected argument type',
                    'The argument to a function failed an input check.  This '
                    'is most likely due to passing an invalid argument to '
                    'the function.\n\n'
                    'Argument: %s\n'
                    'Expected type: %s'
                    % (argument_string, format_string))

    def _input_check(self, argument, format):
        """
        Verify an argument is of the expected format.

        For example, to check if an argument is a list of 3 integers, set
        format=[list, 3, int]

        For example, to check if an argument is a list of at least 3 integers,
        set format=[list, -3, int]

        To check if an argument is a list of 3 lists, each with 2 integers, set
        format=[list, 3, list, 2, int].

        Example:
        >>> model._argument_check([1, 2, 3], [list, 3, int])

        """
        self._assert(type(format) is list)
        if format[0] is list:
            self._assert(type(format[1]) is int)
            self._assert(len(format) > 2)
            if type(argument) is not list:
                self._input_check_error(argument, format)
            if format[1] > 0 and len(argument) != format[1]:
                self._input_check_error(argument, format)
            if format[1] <= 0 and len(argument) < -format[1]:
                self._input_check_error(argument, format)
            for x in argument:
                self._input_check(x, format[2:])
        else:
            self._assert(len(format) == 1)
            if format[0] is int:
                if type(argument) is float:
                    if not argument.is_integer():
                        self._input_check_error(argument, format)
                elif type(argument) is not int:
                    self._input_check_error(argument, format)
            elif format[0] is float:
                if type(argument) is not float and type(argument) is not int:
                    self._input_check_error(argument, format)
            else:
                if type(argument) is not format[0]:
                    self._input_check_error(argument, format)

    def build_hex8_cube(self,
                        element_block_id='auto',
                        extents=1.0,
                        divisions=3):
        """
        Create an element block in the shape of a cuboid.

        Extent defines the extent of the block in the form
        '[[minx, maxx], [miny, maxy], [minz, maxz]]'.  If only one value per
        dimension is given, the minimum is assumed to be zero.  If a single
        value is given, it is assumed a cube with that edge length, with min=0
        for all dimensions.

        Divisions is the number of divisions per dimension.  This can be
        defined as a single number for all dimensions or as a list of
        3 numbers.

        """
        # process shortcuts on arguments
        if element_block_id == 'auto':
            element_block_id = self._new_element_block_id()
        if type(extents) is not list:
            self._input_check(extents, [float])
            extents = [[0.0, float(extents)]] * 3
        if type(extents) is list:
            for i, x in enumerate(extents):
                if type(x) is not list:
                    self._input_check(x, [float])
                    extents[i] = [0.0, float(x)]
        if type(divisions) is int:
            divisions = [divisions] * 3
        # verify the input
        self._input_check(extents, [list, 3, list, 2, float])
        self._input_check(divisions, [list, 3, int])
        dimx, dimy, dimz = divisions
        [[minx, maxx], [miny, maxy], [minz, maxz]] = extents
        # create the nodes
        new_nodes = []
        node_offset = len(self.nodes)
        for k in range(dimz + 1):
            z = minz + (maxz - minz) * float(k) / dimz
            for j in range(dimy + 1):
                y = miny + (maxy - miny) * float(j) / dimy
                for i in range(dimx + 1):
                    x = minx + (maxx - minx) * float(i) / dimx
                    new_nodes.append([x, y, z])
        del i, j, k
        del x, y, z
        # create the connectivity matrix
        connectivity = []
        formula = [0, 1, 1 + (dimx + 1), (dimx + 1)]
        formula.extend([x + (dimx + 1) * (dimy + 1) for x in formula])
        formula = [x + node_offset for x in formula]
        for i in range(dimx):
            for j in range(dimy):
                for k in range(dimz):
                    first = (k * (dimy + 1) + j) * (dimx + 1) + i
                    connectivity.extend(first + x for x in formula)
        del i, j, k
        # now create the actual nodes
        self.create_nodes(new_nodes)
        # now create the actual block
        self.create_element_block(element_block_id,
                                  ['hex8', dimx * dimy * dimz, 8, 0],
                                  connectivity=connectivity)

    def count_degenerate_elements(self, element_block_ids='all'):
        """
        Return the number of degenerate elements in the given element blocks.

        A degenerate element is an element which contains one or more nodes
        which are a duplicate of another node within the same element.

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids,
            empty_list_okay=False)
        degenerate_element_count = 0
        for element_block_id in element_block_ids:
            _, info, connectivity, _ = self.element_blocks[element_block_id]
            nodes_per_element = self.get_nodes_per_element(element_block_id)
            element_count = len(connectivity) / nodes_per_element
            for element_index in range(element_count):
                local_node = connectivity[
                    element_index * nodes_per_element:
                    (element_index + 1) * nodes_per_element]
                if len(set(local_node)) != nodes_per_element:
                    degenerate_element_count += 1
        return degenerate_element_count

    def count_disconnected_blocks(self, element_block_ids='all'):
        """
        Return the number of disconnected blocks.

        A disconnected block is a group of elements which are connected to
        each other through one or more nodes.

        """
        element_block_ids = self._format_element_block_id_list(
            element_block_ids,
            empty_list_okay=False)
        nodes = self.get_nodes_in_element_block(element_block_ids)
        # for each node, find the lowest index node that it's connected to
        master = range(len(self.nodes))
        for element_block_id in element_block_ids:
            connectivity = self.get_connectivity(element_block_id)
            nodes_per_element = self.get_nodes_per_element(element_block_id)
            element_count = self.get_element_count(element_block_id)
            for i in xrange(element_count):
                local_node = connectivity[i * nodes_per_element:
                                          (i + 1) * nodes_per_element]
                # find lowest index master out of these
                low = min(local_node)
                for x in local_node:
                    this_low = x
                    while this_low != master[this_low]:
                        this_low = master[this_low]
                    low = min(low, this_low)
                # now set the current master to the lowest index found
                for x in local_node:
                    this_low = x
                    while this_low != master[this_low]:
                        this_low = master[this_low]
                        master[this_low] = low
                    master[this_low] = low
        # now make sure master node list is one-deep
        for i in nodes:
            master[i] = master[master[i]]
        # double check that master node list is one-deep
        for i in nodes:
            assert master[i] == master[master[i]]
        # count the number of master nodes
        block_count = sum(1 for x in nodes if master[x] == x)
        return block_count

    def _get_mating_faces(self, side_set_members_one, side_set_members_two):
        """
        Return the mating faces between the side sets.

        This will only return the faces within 'side_set_members_one'.

        """
        # order by element block
        members_one_by_block = self._order_element_faces_by_block(
            side_set_members_one)
        members_two_by_block = self._order_element_faces_by_block(
            side_set_members_two)
        # find the nodes within set two and create a set of them
        faces_to_match = set()
        for id_, members in members_two_by_block.items():
            nodes_per_element = self.get_nodes_per_element(id_)
            connectivity = self.get_connectivity(id_)
            face_mapping = self._get_face_mapping_from_id(id_)
            for element_index, face_index in members:
                local_node = connectivity[
                    element_index * nodes_per_element:
                    (element_index + 1) * nodes_per_element]
                face_nodes = [local_node[x]
                              for x in face_mapping[face_index][1]]
                faces_to_match.add(tuple(sorted(face_nodes)))
        # now look through the original list for duplicates
        mating_faces = []
        for id_, members in members_one_by_block.items():
            nodes_per_element = self.get_nodes_per_element(id_)
            connectivity = self.get_connectivity(id_)
            face_mapping = self._get_face_mapping_from_id(id_)
            for element_index, face_index in members:
                local_node = connectivity[
                    element_index * nodes_per_element:
                    (element_index + 1) * nodes_per_element]
                face_nodes = [local_node[x]
                              for x in face_mapping[face_index][1]]
                if tuple(sorted(face_nodes)) in faces_to_match:
                    mating_faces.append((id_, element_index, face_index))
        return mating_faces
        # now create a sorted list with (nodes, member)
        # now find duplicate entries

    def _detailed_summary(self):
        """Print a detailed summary of the model."""
        # if no timesteps are defined, create one so we can store field values
        # it will be deleted at the nd
        timestep_created = False
        if not self.timesteps:
            timestep_created = True
            self.create_timestep(0.0)
        # print out general info
        print('\nMODEL SUMMARY\n')
        ids = self.get_element_block_ids()
        element_block_count = len(ids)
        element_count = self.get_element_count(ids)
        print('- Model contains %d element blocks and %d elements'
              % (element_block_count, element_count))
        # if any elements exist...
        if element_count:
            # calculate element volume field
            element_volume_field_name = self._new_element_field_name(1)
            self.calculate_element_volumes(element_volume_field_name)
            # calculate element centroid field
            element_centroid_field_names = self._new_element_field_name(3)
            self.calculate_element_centroids(element_centroid_field_names)
        # print out element block info
        if element_count:
            extents = self.get_element_block_extents(ids)
            print('- Extents are:')
            for d in range(3):
                print('  - %s: %g to %g, range of %g'
                      % ('XYZ'[d],
                         extents[d][0],
                         extents[d][1],
                         extents[d][1] - extents[d][0]))
            # print center of mass
            cg = self.get_element_block_centroid(ids,
                                                 element_volume_field_name,
                                                 element_centroid_field_names)
            cg = ['%g' % x for x in cg]
            print('- Center of volume is at [%s]' % (', '.join(cg)))
            # print total volume
            volume = self.get_element_block_volume(ids,
                                                   element_volume_field_name)
            print('- Total volume is %g' % volume)
            # print total surface area
            side_set_id = self._new_side_set_id()
            self.create_side_set(side_set_id,
                                 self._get_external_element_faces(ids))
            area = self.get_side_set_area(side_set_id)
            self.delete_side_set(side_set_id)
            print('- Total surface area is %d' % area)
            # print number of disconnected blocks
            connected_blocks = self.count_disconnected_blocks('all')
            print('- Contains %d disconnected blocks' % (connected_blocks))
            # print element edge length stats
            minimum, average = self.get_element_edge_length_info(ids)
            print('- Average element edge length is %g' % (average))
            print('- Smallest element edge length is %g' % (minimum))
        if self.nodes:
            node_distance = self.get_closest_node_distance()
            print('- The closest node pair is %g apart' % (node_distance))
        print
        print('ELEMENT BLOCK INFO')
        # find external faces for each element block
        if element_count:
            external_faces = dict((id_, self._get_external_element_faces(id_))
                                  for id_ in self.get_element_block_ids())
        # print info on each element block
        for id_ in self.get_element_block_ids():
            print
            name = self.get_element_block_name(id_)
            print('Element block ID %d%s:'
                  % (id_, (' "%s"' % (name)) if name else ''))
            dim = self.get_element_block_dimension(id_)
            element_count = self.get_element_count(id_)
            element_type = self._get_element_type(id_)
            dim_name = ' %d-dimensional' % dim if dim != -1 else ''
            print('- Contains %d "%s"%s elements'
                  % (element_count,
                     element_type,
                     dim_name))
            # if no elements, skip detailed info on this block
            if not element_count:
                continue
            extents = self.get_element_block_extents(id_)
            print('- Extents are:')
            for d in range(3):
                print('  - %s: %g to %g, range of %g'
                      % ('XYZ'[d],
                         extents[d][0],
                         extents[d][1],
                         extents[d][1] - extents[d][0]))
            # print center of mass
            cg = self.get_element_block_centroid(id_,
                                                 element_volume_field_name,
                                                 element_centroid_field_names)
            cg = ['%g' % x for x in cg]
            print('- Center of volume is at [%s]' % (', '.join(cg)))
            # print total volume
            volume = self.get_element_block_volume(id_,
                                                   element_volume_field_name)
            print('- Total volume is %g' % (volume))
            # print total surface area
            side_set_id = self._new_side_set_id()
            self.create_side_set(side_set_id, external_faces[id_])
            area = self.get_side_set_area(side_set_id)
            self.delete_side_set(side_set_id)
            print('- Total surface area is %g' % (area))
            # print number of disconnected blocks
            connected_blocks = self.count_disconnected_blocks(id_)
            print('- Contains %d disconnected blocks' % (connected_blocks))
            # print element edge length stats
            if dim == 3:
                minimum, average = self.get_element_edge_length_info(id_)
                print('- Average element edge length is %g' % (average))
                print('- Smallest element edge length is %g' % (minimum))
                # print thickness
                thickness = self._get_thickness_from_volume_and_area(volume,
                                                                     area)
                print('- Approximate thickness is %g' % (thickness))
            # print surface connectivity to other blocks
            # find faces connected to other blocks
            remaining_faces = set(external_faces[id_])
            header_output = False
            for other_id in self.get_element_block_ids():
                if other_id == id_:
                    continue
                mating_faces = self._get_mating_faces(external_faces[id_],
                                                      external_faces[other_id])
                if mating_faces:
                    remaining_faces -= set(mating_faces)
                    side_set_id = self._new_side_set_id()
                    self.create_side_set(side_set_id, list(mating_faces))
                    area = self.get_side_set_area(side_set_id)
                    self.delete_side_set(side_set_id)
                    if not header_output:
                        print('- Connected to the following element blocks:')
                        header_output = True
                    print('  - To element block %d though %d faces '
                          '(area of %g)'
                          % (other_id, len(mating_faces), area))
            if header_output and remaining_faces:
                remaining_faces = list(remaining_faces)
                side_set_id = self._new_side_set_id()
                self.create_side_set(side_set_id, remaining_faces)
                area = self.get_side_set_area(side_set_id)
                self.delete_side_set(side_set_id)
                print('  - To the outside through %d faces '
                      '(area of %g)'
                      % (len(remaining_faces), area))
            if not header_output:
                print('- Not connected to any element blocks')
        # print node set info
        print('\nNODE SET INFO\n')
        ids = self.get_node_set_ids()
        print('There are %d node sets defined.' % (len(ids)))
        for id_ in ids:
            print
            name = self.get_node_set_name(id_)
            print('Node set ID %d%s:'
                  % (id_, (' "%s"' % (name)) if name else ''))
            print('- Contains %d members'
                  % (len(self.get_node_set_members(id_))))
            field_names = self.get_node_set_field_names(id_)
            if field_names:
                print('- Has %d fields defined:' % (len(field_names)))
                for name in field_names:
                    print('  - "%s"' % (name))
        # print node set info
        print('\nSIDE SET INFO\n')
        ids = self.get_side_set_ids()
        print('There are %d side sets defined.' % (len(ids)))
        for id_ in ids:
            print
            name = self.get_side_set_name(id_)
            print('Side set ID %d%s:'
                  % (id_, (' "%s"' % (name)) if name else ''))
            members = self.get_side_set_members(id_)
            member_count = len(members)
            print('- Contains %d members' % (member_count))
            parent_blocks = sorted(set(x[0] for x in members))
            parent_string = ', '.join('%s' % (x) for x in parent_blocks)
            print('- Parent element block IDs: %s' % (parent_string))
            face_types = []
            members_by_block = self._order_element_faces_by_block(members)
            for block_id, these_members in members_by_block.items():
                element_type = self._get_element_type(block_id)
                if not self._is_standard_element_type(element_type):
                    face_types.append('unknown')
                    continue
                face_mapping = self._get_face_mapping_from_id(block_id)
                face_indices = set(x[1] for x in these_members)
                face_types.extend(face_mapping[x][0] for x in face_indices)
            face_types = sorted(set(face_types))
            print('- Face types: %s' % (', '.join(face_types)))
            area = self.get_side_set_area(id_)
            print('- Total area of %g' % (area))
            field_names = self.get_side_set_field_names(id_)
            if field_names:
                print('- Has %d fields defined:' % (len(field_names)))
                for name in field_names:
                    print('  - "%s"' % (name))
        # delete temporary timestep if created
        if timestep_created:
            self.delete_timestep(0.0)
