#!/usr/bin/env python
"""
This file performs unit tests of functions within Exomerge.

Copyright 2018, 2021, 2022 National Technology and Engineering
Solutions of Sandia.  Under the terms of Contract DE-NA-0003525, there
is a non-exclusive license for use of this work by or on behalf of the
U.S. Government.  Export of this program may require a license from
the United States Government.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

Author: Tim Kostka (tdkostk@sandia.gov)
Created: May 3, 2012

This module runs unit tests in a randomized fashion.  To run:

$ python exomerge_unit_test.py 1000

Although the tests performed and randomly generated, since the random number is
seeded with a static value, the *same* sequence of randomly generated tests is
run each time.  Processor speed can dictate how many tests are run, so setting
min_tests and max_tests to the same value or calling this script with an
argument can be used to prevent this.

"""

# standard libraries
import random
import inspect
import time
import copy
import os
import sys
import re
import math

# import the exomerge module
ACCESS = os.getenv("ACCESS", "@ACCESSDIR@")
sys.path.append(os.path.join(ACCESS, "lib"))
import exomerge


def _random_element(the_list):
    """
    Return a random element from the list or None if none exist.

    """
    if not the_list:
        return None
    return the_list[random.randint(0, len(the_list) - 1)]


def _random_boolean():
    """
    Return one of True or False at random.

    """
    return [True, False][random.randint(0, 1)]


def _random_subset(the_list, count=None):
    """
    Return a random subset of elements from the list.

    """
    if count is None:
        if len(the_list) == 0:
            count = 0
        else:
            count = random.randint(1, len(the_list))
    return random.sample(the_list, count)


def _random_scalar():
    """
    Return a random scalar.

    """
    return random.uniform(0.5, 2)


def _random_vector():
    """
    Return a random vector.

    """
    vector = [_random_scalar() for _ in range(3)]
    return vector


def _new_id():
    """
    Return a new id.

    """
    return random.randint(0, 9999)


def compares_equal_with_nan(one, two):
    """
    Return True if the two objects are equal, assuming NaN == NaN.

    Objects can include lists, dictionaries, floats, integers, strings, and/or
    floats.  They can be nested.

    """
    if type(one) is not type(two):
        return False
    if isinstance(one, dict):
        if sorted(one.keys()) != sorted(two.keys()):
            return False
        return all(compares_equal_with_nan(one[x], two[x]) for x in list(one.keys()))
    elif isinstance(one, list):
        return len(one) == len(two) and all(
            compares_equal_with_nan(x, y) for x, y in zip(one, two)
        )
    else:
        # string, integer, or float
        if isinstance(one, float):
            return one == two or (math.isnan(one) and math.isnan(two))
        else:
            return one == two


class DummyFile(object):
    """Dummy class used to suppress output."""

    def write(self, x):
        """Ignore the write command."""
        pass


class OutputSuppression:
    """
    This class is used to suppress output to stdout.


    Example:
    >>> std.stdout('Hello world!\n')
    Hello world!
    >>> print('Hello world!')
    Hello world!
    >>> with OutputSuppression():
    >>>     std.stdout('This will not appear\n')
    >>>     print('Hello world!')

    """

    def __enter__(self):
        """Initialization routine."""
        self.old_stdout = sys.stdout
        sys.stdout = DummyFile()

    def __exit__(self, type_, value, tracebackself):
        """Destructor routine."""
        sys.stdout = self.old_stdout


class ExomergeUnitTester:
    """
    A class to perform unit tests of the exomerge module.

    """

    # identifiers used to create names
    RANDOM_IDENTIFIERS = (
        "abcdefghijklmnopqrstuvwxyz" "ABCDEFGHIJKLMNOPQRSTUVWXYZ" "0123456789"
    )

    def __init__(self):
        # store the exodus model
        self.model = exomerge.ExodusModel()
        # run tests for up to this amount of time in seconds
        self.time_limit = 60
        # run at least this many tests
        self.min_tests = 3000
        # run up to this many tests
        self.max_tests = 3000
        # only this many IO tests will be run
        self.remaining_io_tests = 40
        # set the maximum number of each object to prevent self.model from
        # growing too large
        self.maximum_objects = 10

    def _random_identifier(self):
        """
        Return a random name.

        """
        return "".join([_random_element(self.RANDOM_IDENTIFIERS) for _ in range(6)])

    def _random_timestep(self):
        """
        Return a random timestep.

        """
        if not self.model.timesteps:
            return None
        return random.choice(self.model.timesteps)

    def _random_field_name(self):
        """
        Return a random field name.

        """
        return "field_" + self._random_identifier()

    def _new_node_field_name(self):
        """
        Return a new node field name.

        """
        name = self._random_field_name()
        while name in self.model.get_node_field_names():
            name = self._random_field_name()
        return name

    def _new_global_variable_name(self):
        """
        Return a new global variable name.

        """
        name = self._random_field_name()
        while name in self.model.get_global_variable_names():
            name = self._random_field_name()
        return name

    def _new_element_field_name(self):
        """
        Return a new element field name.

        """
        name = self._random_field_name()
        while name in self.model.get_element_field_names():
            name = self._random_field_name()
        return name

    def _new_node_set_field_name(self):
        """
        Return a new node set field name.

        """
        name = self._random_field_name()
        while name in self.model.get_node_set_field_names():
            name = self._random_field_name()
        return name

    def _new_side_set_field_name(self):
        """
        Return a new side set field name.

        """
        name = self._random_field_name()
        while name in self.model.get_side_set_field_names():
            name = self._random_field_name()
        return name

    def _new_timestep(self):
        """
        Return a new timestep.

        """
        value = random.random()
        while value in self.model.timesteps:
            value = random.random()
        return value

    def _new_element_block_id(self):
        """
        Return a new element block id.

        """
        value = _new_id()
        while value in self.model.get_element_block_ids():
            value = _new_id()
        return value

    def _new_node_set_id(self):
        """
        Return a new node set id.

        """
        value = _new_id()
        while value in self.model.get_node_set_ids():
            value = _new_id()
        return value

    def _new_side_set_id(self):
        """
        Return a new side set id.

        """
        value = _new_id()
        while value in self.model.get_side_set_ids():
            value = _new_id()
        return value

    def _random_element_field_name(self):
        """
        Return one element field name or None if none exist in the form
        (id_, name).

        """
        ids = self.model.get_element_block_ids()
        random.shuffle(ids)
        for id_ in ids:
            fields = self.model._get_element_block_fields(id_)
            names = list(fields.keys())
            random.shuffle(names)
            if names:
                return (id_, names[0])
        return None

    def _random_node_field_name(self):
        """
        Return one node field name or None.

        """
        return _random_element(self.model.get_node_field_names())

    def _random_node_set_field_name(self):
        """
        Return one node set field name or None in the form (id_, name).

        """
        ids = self.model.get_node_set_ids()
        random.shuffle(ids)
        for id_ in ids:
            fields = self.model._get_node_set_fields(id_)
            names = list(fields.keys())
            random.shuffle(names)
            if names:
                return (id_, names[0])
        return None

    def _random_side_set_field_name(self):
        """
        Return one side set field name or None in the form (id_, name).

        """
        ids = self.model.get_side_set_ids()
        random.shuffle(ids)
        for id_ in ids:
            fields = self.model._get_side_set_fields(id_)
            names = list(fields.keys())
            random.shuffle(names)
            if names:
                return (id_, names[0])
        return None

    def _random_global_variable_name(self):
        """
        Return one global variable name or None.

        """
        return _random_element(self.model.get_global_variable_names())

    def _random_element_block_id(self):
        """Return a random element block id or name, or None if none exist."""
        ids = self.model.get_element_block_ids()
        ids.extend(self.model.get_all_element_block_names())
        if not ids:
            return None
        return _random_element(ids)

    def _random_element_block_ids(self):
        """Return a random subset of element block ids or names."""
        ids = _random_subset(self.model.get_element_block_ids())
        for i, id_ in enumerate(ids):
            if random.randint(0, 1) == 0:
                continue
            name = self.model.get_element_block_name(id_)
            if not name:
                continue
            ids[i] = name
        return ids

    def _random_side_set_id(self):
        """Return a random side set id or name, or None if none exist."""
        ids = self.model.get_side_set_ids()
        ids.extend(self.model.get_all_side_set_names())
        return _random_element(ids)

    def _random_side_set_ids(self):
        """Return a random subset of side set ids or name."""
        ids = _random_subset(self.model.get_side_set_ids())
        for i, id_ in enumerate(ids):
            if random.randint(0, 1) == 0:
                continue
            name = self.model.get_side_set_name(id_)
            if not name:
                continue
            ids[i] = name
        return ids

    def _random_node_set_id(self):
        """Return a random side set id or name, or None if none exist."""
        ids = self.model.get_node_set_ids()
        ids.extend(self.model.get_all_node_set_names())
        return _random_element(ids)

    def _random_node_set_ids(self):
        """Return a random subset of node set ids or name."""
        ids = _random_subset(self.model.get_node_set_ids())
        for i, id_ in enumerate(ids):
            if random.randint(0, 1) == 0:
                continue
            name = self.model.get_node_set_name(id_)
            if not name:
                continue
            ids[i] = name
        return ids

    def _truncate_element_blocks(self, number=None):
        """
        Delete element blocks until only the given number exist.

        """
        if number is None:
            number = self.maximum_objects
        objects = self.model.get_element_block_ids()
        count = len(objects)
        if count > number:
            self.model.delete_element_block(_random_subset(objects, count - number))

    def _truncate_timesteps(self, number=None):
        """
        Delete timesteps until only the given number exist.

        """
        if number is None:
            number = self.maximum_objects
        objects = self.model.get_timesteps()
        count = len(objects)
        if count > number:
            self.model.delete_timestep(_random_subset(objects, count - number))

    def _truncate_node_sets(self, number=None):
        """
        Delete node sets until only the given number exist.

        """
        if number is None:
            number = self.maximum_objects
        objects = self.model.get_node_set_ids()
        count = len(objects)
        if count > number:
            self.model.delete_node_set(_random_subset(objects, count - number))

    def _truncate_side_sets(self, number=None):
        """
        Delete side sets until only the given number exist.

        """
        if number is None:
            number = self.maximum_objects
        objects = self.model.get_side_set_ids()
        count = len(objects)
        if count > number:
            self.model.delete_side_set(_random_subset(objects, count - number))

    def _truncate_node_fields(self, number=None):
        """
        Delete random node fields until only the given number exist.

        """
        if number is None:
            number = self.maximum_objects
        objects = self.model.get_node_field_names()
        count = len(objects)
        if count > number:
            self.model.delete_node_field(_random_subset(objects, count - number))

    def _truncate_global_variables(self, number=None):
        """
        Delete random global variables until only the given number exist.

        """
        if number is None:
            number = self.maximum_objects
        objects = self.model.get_global_variable_names()
        count = len(objects)
        if count > number:
            self.model.delete_global_variable(_random_subset(objects, count - number))

    def _truncate_element_fields(self, number=None):
        """
        Delete random element fields until only the given number exist.

        """
        if number is None:
            number = self.maximum_objects
        while len(self.model.get_element_field_names()) > number:
            (id_, name) = self._random_element_field_name()
            self.model.delete_element_field(name, id_)

    def _truncate_node_set_fields(self, number=None):
        """
        Delete random node set fields until only the given number exist.

        """
        if number is None:
            number = self.maximum_objects
        while len(self.model.get_node_set_field_names()) > number:
            (id_, name) = self._random_node_set_field_name()
            self.model.delete_node_set_field(name, id_)

    def _truncate_side_set_fields(self, number=None):
        """
        Delete random side set fields until only the given number exist.

        """
        if number is None:
            number = self.maximum_objects
        while len(self.model.get_side_set_field_names()) > number:
            (id_, name) = self._random_side_set_field_name()
            self.model.delete_side_set_field(name, id_)

    def _create_hex8_element_block(self, origin=None):
        """
        Create a 2-element hex8 element block of dimension 2x1x1 with the
        given minimum coordinates and return the element block id.

        """
        new_nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 0.0, 1.0],
            [2.0, 1.0, 1.0],
        ]
        if origin is None:
            origin = [0, 0, 0]
        new_nodes = [[x + y for x, y in zip(node, origin)] for node in new_nodes]
        connectivity = list(range(8))
        connectivity.extend([1, 8, 9, 2, 5, 10, 11, 6])
        connectivity = [x + len(self.model.nodes) for x in connectivity]
        self.model.create_nodes(new_nodes)
        info = ["hex8", 2, 8, 0]
        element_block_id = self._new_element_block_id()
        self.model.create_element_block(element_block_id, info, connectivity)
        return element_block_id

    def _topology_test(self):
        """
        Test the topology definitions to ensure they are consistent.

        FACE_MAPPING, INVERTED_CONNECTIVITY, DIMENSION

        """
        # see if keys are defined for all places
        self.model._assert(
            (
                set(self.model.FACE_MAPPING.keys())
                == set(self.model.INVERTED_CONNECTIVITY.keys())
            )
        )
        self.model._assert(
            (set(self.model.FACE_MAPPING.keys()) == set(self.model.DIMENSION.keys()))
        )
        self.model._assert(
            (
                set(self.model.FACE_MAPPING.keys())
                == set(self.model.NODES_PER_ELEMENT.keys())
            )
        )
        # faces should be 1 dimension less than their element
        for element_type, face_list in list(self.model.FACE_MAPPING.items()):
            element_dimension = self.model._get_dimension(element_type)
            for face_type, _ in face_list:
                self.model._assert(
                    self.model._get_dimension(face_type) + 1 == element_dimension
                )
        # all 2D elements should be able to be triangulated
        for element_type in self.model.STANDARD_ELEMENT_TYPES:
            if self.model.DIMENSION[element_type] == 2:
                self.model._assert(element_type in self.model.TRIANGULATED_FACES)

    # The following functions are unit tests of public functions within
    # exomerge.  For each public_function, a unit test of the name
    # _test_public_function should exist here.  If there is not a one to
    # one mapping between the two, warning messages will appear when running
    # this module.
    #
    # Tests should return None if successful (no return statement needed)
    # Tests should return False if the test was unable to be run.
    # Tests should raise an exception or exit(1) if unsuccessful.

    def _test_calculate_element_volumes(self):
        ids = self.model._get_standard_element_block_ids()
        if not ids:
            return False
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        self.model.calculate_element_volumes(self._new_element_field_name(), ids)

    def _test_get_element_block_centroid(self):
        ids = self._random_element_block_ids()
        if not ids:
            return False
        self.model.get_element_block_centroid(ids)

    def _test_get_element_block_volume(self):
        ids = self._random_element_block_ids()
        if not ids:
            return False
        self.model.get_element_block_volume(ids)

    def _test_convert_element_blocks(self):
        id_ = _random_element(self.model._get_standard_element_block_ids())
        if not id_:
            return False
        element_type = self.model._get_element_type(id_)
        if element_type not in self.model.ELEMENT_CONVERSIONS:
            return False
        scheme = self.model.ELEMENT_CONVERSIONS[element_type]
        if not scheme:
            return False
        new_element_type = _random_element(list(scheme.keys()))
        self.model.convert_element_blocks(id_, new_element_type)

    def _test_make_elements_linear(self):
        ids = []
        for id_ in self.model._get_standard_element_block_ids():
            element_type = self.model._get_element_type(id_)
            if self.model.ELEMENT_ORDER[element_type] == 1:
                continue
            if element_type not in self.model.ELEMENT_CONVERSIONS:
                continue
            scheme = self.model.ELEMENT_CONVERSIONS[element_type]
            if 1 not in [self.model.ELEMENT_ORDER[x] for x in list(scheme.keys())]:
                continue
            ids.append(id_)
        ids = _random_subset(ids)
        if not ids:
            return False
        self.model.make_elements_linear(ids)

    def _test_make_elements_quadratic(self):
        ids = []
        for id_ in self.model._get_standard_element_block_ids():
            element_type = self.model._get_element_type(id_)
            if self.model.ELEMENT_ORDER[element_type] == 2:
                continue
            if element_type not in self.model.ELEMENT_CONVERSIONS:
                continue
            scheme = self.model.ELEMENT_CONVERSIONS[element_type]
            if 2 not in [self.model.ELEMENT_ORDER[x] for x in list(scheme.keys())]:
                continue
            ids.append(id_)
        ids = _random_subset(ids)
        if not ids:
            return False
        self.model.make_elements_quadratic(ids)

    def _test_convert_hex8_block_to_tet4_block(self):
        ids = [
            x
            for x in self.model.get_element_block_ids()
            if self.model._get_element_type(x) == "hex8"
        ]
        if not ids:
            return False
        self.model.convert_hex8_block_to_tet4_block(_random_element(ids))

    def _test_duplicate_element_block(self):
        old_id = self._random_element_block_id()
        if old_id is None:
            return False
        new_id = self._new_element_block_id()
        self.model.duplicate_element_block(old_id, new_id)

    def _test_to_lowercase(self):
        self.model.to_lowercase()

    def _test_reflect_element_blocks(self):
        id_ = self._create_hex8_element_block()
        self.model.reflect_element_blocks(id_, [0, 0, 0], [1, -1, 0])

    def _test_combine_element_blocks(self):
        id_one = self._create_hex8_element_block()
        id_two = self._create_hex8_element_block([2, 0, 0])
        self.model.combine_element_blocks([id_one, id_two])

    def _test_create_side_set_from_expression(self):
        new_nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
        connectivity = list(range(8))
        connectivity = [x + len(self.model.nodes) for x in connectivity]
        self.model.create_nodes(new_nodes)
        info = ["hex8", 1, 8, 0]
        element_block_id = self._new_element_block_id()
        self.model.create_element_block(element_block_id, info, connectivity)
        self.model.create_side_set_from_expression(
            self._new_side_set_id(),
            "Z = 0 || Z == 1 || X == 1",
            element_block_ids=element_block_id,
            timesteps="none",
        )

    def _test_convert_side_set_to_cohesive_zone(self):
        new_nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 0.0, 1.0],
            [2.0, 1.0, 1.0],
        ]
        connectivity = list(range(8))
        connectivity.extend([1, 8, 9, 2, 5, 10, 11, 6])
        connectivity = [x + len(self.model.nodes) for x in connectivity]
        self.model.create_nodes(new_nodes)
        info = ["hex8", 2, 8, 0]
        element_block_id = self._new_element_block_id()
        self.model.create_element_block(element_block_id, info, connectivity)
        side_set_id = self._new_side_set_id()
        self.model.create_side_set(side_set_id, [(element_block_id, 1, 3)])
        self.model.convert_side_set_to_cohesive_zone(
            side_set_id, self._new_element_block_id()
        )

    def _test_calculate_global_variable(self):
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        variable_names = set(["time"])
        variable_names.update(self.model.get_global_variable_names())
        expression = "%s = %s" % (
            self._new_global_variable_name(),
            " + ".join(variable_names),
        )
        self.model.calculate_global_variable(expression)

    def _test_calculate_node_field(self):
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        variable_names = set(["time"])
        variable_names.update(self.model.get_global_variable_names())
        variable_names.update(self.model.get_node_field_names())
        expression = "%s = %s" % (
            self._new_node_field_name(),
            " + ".join(variable_names),
        )
        self.model.calculate_node_field(expression)

    def _test_calculate_element_field(self):
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        info = self._random_element_field_name()
        if info is None:
            return False
        (id_, _) = info
        variable_names = set(["time"])
        variable_names.update(self.model.get_global_variable_names())
        variable_names.update(self.model.get_element_field_names(id_))
        expression = "%s = %s" % (
            self._new_element_field_name(),
            " + ".join(variable_names),
        )
        self.model.calculate_element_field(expression, id_)

    def _test_calculate_element_field_maximum(self):
        info = self._random_element_field_name()
        if info is None or len(info[1]) + 13 > 32:
            return False
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        (id_, name) = info
        block_ids = [
            x
            for x in self.model.get_element_block_ids()
            if self.model.element_field_exists(name, x)
        ]
        block_ids = _random_subset(block_ids)
        if self.model.get_element_count(block_ids) == 0:
            return False
        self.model.delete_global_variables(name + "_*")
        self.model.calculate_element_field_maximum(
            name,
            element_block_ids=block_ids,
            calculate_location=_random_boolean(),
            calculate_block_id=_random_boolean(),
        )

    def _test_calculate_element_field_minimum(self):
        info = self._random_element_field_name()
        if info is None or len(info[1]) + 13 > 32:
            return False
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        (id_, name) = info
        block_ids = [
            x
            for x in self.model.get_element_block_ids()
            if self.model.element_field_exists(name, x)
        ]
        block_ids = _random_subset(block_ids)
        if self.model.get_element_count(block_ids) == 0:
            return False
        self.model.delete_global_variables(name + "_*")
        self.model.calculate_element_field_minimum(
            name,
            element_block_ids=block_ids,
            calculate_location=_random_boolean(),
            calculate_block_id=_random_boolean(),
        )

    def _test_calculate_node_field_minimum(self):
        if not self.model.nodes:
            return False
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        name = self._random_node_field_name()
        if name is None or len(name) + 6 > 32:
            return False
        self.model.delete_global_variables(name + "_*")
        if _random_boolean():
            self.model.calculate_node_field_minimum(
                name, calculate_location=_random_boolean()
            )
        else:
            block_ids = self._random_element_block_ids()
            if not self.model.get_nodes_in_element_block(block_ids):
                return False
            self.model.calculate_node_field_minimum(
                name, element_block_ids=block_ids, calculate_location=_random_boolean()
            )

    def _test_calculate_node_field_maximum(self):
        if not self.model.nodes:
            return False
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        name = self._random_node_field_name()
        if name is None or len(name) + 6 > 32:
            return False
        self.model.delete_global_variables(name + "_*")
        self.model.calculate_node_field_maximum(
            name, calculate_location=_random_boolean()
        )

    def _test_output_global_variables(self):
        filename = "temp.csv" if _random_boolean() else None
        if self.remaining_io_tests <= 0:
            filename = None
        if filename is not None:
            self.remaining_io_tests -= 1
        with OutputSuppression():
            self.model.output_global_variables(filename)
        if filename:
            os.remove(filename)

    def _test_calculate_node_set_field(self):
        info = self._random_node_set_field_name()
        if info is None:
            return False
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        (id_, name) = info
        variable_names = set(["time"])
        variable_names.update(self.model.get_global_variable_names())
        variable_names.update(self.model.get_node_field_names())
        variable_names.update(self.model.get_node_set_field_names(id_))
        expression = "%s = %s" % (
            self._new_node_set_field_name(),
            " + ".join(variable_names),
        )
        self.model.calculate_node_set_field(expression, id_)

    def _test_calculate_side_set_field(self):
        info = self._random_side_set_field_name()
        if info is None:
            return False
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        (id_, name) = info
        variable_names = set(["time"])
        variable_names.update(self.model.get_global_variable_names())
        variable_names.update(self.model.get_side_set_field_names(id_))
        expression = "%s = %s" % (
            self._new_side_set_field_name(),
            " + ".join(variable_names),
        )
        self.model.calculate_side_set_field(expression, id_)

    def _test_merge_nodes(self):
        if not self.model.nodes:
            return False
        x = len(self.model.nodes)
        self.model._duplicate_nodes(
            _random_subset(list(range(len(self.model.nodes))), 1), []
        )
        self.model.merge_nodes(suppress_warnings=True)
        assert len(self.model.nodes) <= x

    def _test_unmerge_element_blocks(self):
        if len(self.model.element_blocks) < 2:
            return False
        (id1, id2) = _random_subset(self.model.get_element_block_ids(), count=2)
        connectivity_1 = self.model.get_connectivity(id1)
        connectivity_2 = self.model.get_connectivity(id2)
        if not connectivity_1 or not connectivity_2:
            return False
        connectivity_2[0] = connectivity_1[0]
        self.model.unmerge_element_blocks([id1, id2])
        assert not (
            set(self.model.get_connectivity(id1))
            & set(self.model.get_connectivity(id2))
        )

    def _test_get_length_scale(self):
        self.model.get_length_scale()

    def _test_get_input_deck(self):
        self.model.info_records = [
            "asdf",
            "begin",
            "  begin statement",
            "  end statement",
            "end sierra",
            "other stuff",
        ]
        assert len(self.model.get_input_deck().split("\n")) == 4

    def _test_create_interpolated_timestep(self):
        if len(self.model.timesteps) < 2:
            return False
        steps = self.model.get_timesteps()
        new_time = steps[0] + random.random() * (steps[-1] - steps[0])
        if self.model.timestep_exists(new_time):
            return False
        self.model.create_interpolated_timestep(new_time)

    def _test_rename_element_block(self):
        id_ = self._random_element_block_id()
        if id_ is None:
            return False
        self.model.rename_element_block(id_, self._new_element_block_id())

    def _test_rename_side_set(self):
        id_ = self._random_side_set_id()
        if id_ is None:
            return False
        self.model.rename_side_set(id_, self._new_side_set_id())

    def _test_rename_node_set(self):
        id_ = self._random_node_set_id()
        if id_ is None:
            return False
        self.model.rename_node_set(id_, self._new_node_set_id())

    def _test_rename_global_variable(self):
        name = self._random_global_variable_name()
        if name is None:
            return False
        self.model.rename_global_variable(name, self._new_global_variable_name())

    def _test_rename_node_field(self):
        name = self._random_node_field_name()
        if name is None:
            return False
        self.model.rename_node_field(name, self._new_node_field_name())

    def _test_rename_node_set_field(self):
        info = self._random_node_set_field_name()
        if info is None:
            return False
        (id_, name) = info
        self.model.rename_node_set_field(name, self._new_node_set_field_name(), id_)

    def _test_rename_side_set_field(self):
        info = self._random_side_set_field_name()
        if info is None:
            return False
        (id_, name) = info
        self.model.rename_side_set_field(name, self._new_side_set_field_name(), id_)

    def _test_delete_node_set_field(self):
        info = self._random_node_set_field_name()
        if info is None:
            self._test_create_node_set_field()
            info = self._random_node_set_field_name()
            if info is None:
                return False
        (id_, name) = info
        self.model.delete_node_set_field(name, id_)

    def _test_delete_side_set_field(self):
        info = self._random_side_set_field_name()
        if info is None:
            return False
        (id_, name) = info
        self.model.delete_side_set_field(name, id_)

    def _test_calculate_element_centroids(self):
        if not self.model.element_blocks:
            return False
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        for id_ in self.model.get_element_block_ids():
            for name in ["centroid_x", "centroid_y", "centroid_z"]:
                if self.model.element_field_exists(name, id_):
                    self.model.delete_element_field(name, id_)
        self.model.calculate_element_centroids()
        self._truncate_element_fields()

    def _test_convert_element_field_to_node_field(self):
        field = self._random_element_field_name()
        if field is None:
            return False
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        self.model.convert_element_field_to_node_field(
            field[1], self._new_node_field_name()
        )
        self._truncate_node_fields()

    def _test_convert_node_field_to_element_field(self):
        field = self._random_node_field_name()
        if field is None:
            return False
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        self.model.convert_node_field_to_element_field(
            field, self._new_element_field_name()
        )
        self._truncate_element_fields()

    def _test_create_averaged_element_field(self):
        self._truncate_element_fields()
        id_ = self._random_element_block_id()
        if id_ is None:
            return False
        field_names = self.model.get_element_field_names(id_)
        if not field_names:
            return False
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        field = _random_subset(field_names)
        self.model.create_averaged_element_field(
            field, self._new_element_field_name(), id_
        )

    def _test_create_displacement_field(self):
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        self.model.create_displacement_field()
        self._truncate_node_fields()

    def _test_create_element_block(self):
        for _ in range(random.randint(1, self.maximum_objects)):
            new_nodes = [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 1.0],
                [0.0, 1.0, 1.0],
                [0.0, 1.0, 1.0],
                [2.0, 0.0, 0.0],
                [2.0, 1.0, 0.0],
                [2.0, 0.0, 1.0],
                [2.0, 1.0, 1.0],
            ]
            connectivity = list(range(8))
            connectivity.extend([1, 8, 9, 2, 5, 10, 11, 6])
            connectivity = [x + len(self.model.nodes) for x in connectivity]
            self.model.create_nodes(new_nodes)
            info = ["hex8", 2, 8, 0]
            self.model.create_element_block(
                self._new_element_block_id(), info, connectivity
            )
        self._truncate_element_blocks()

    def _test_create_element_field(self):
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        for _ in range(random.randint(1, self.maximum_objects)):
            ids = self._random_element_block_ids()
            if ids:
                self.model.create_element_field(
                    self._new_element_field_name(), ids, random.random()
                )
        self._truncate_element_fields()

    def _test_create_global_variable(self):
        for _ in range(random.randint(1, self.maximum_objects)):
            self.model.create_global_variable(self._new_global_variable_name())
        self._truncate_global_variables()

    def _test_create_node_field(self):
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        for _ in range(random.randint(1, self.maximum_objects)):
            self.model.create_node_field(self._new_node_field_name(), random.random())
        self._truncate_node_fields()

    def _test_add_nodes_to_node_set(self):
        self._truncate_node_sets()
        if not self.model.nodes or not self.model.node_sets:
            return False
        id_ = self._random_node_set_id()
        current_members = self.model.get_node_set_members(id_)
        members = [random.randint(0, len(self.model.nodes) - 1) for _ in range(20)]
        members = sorted(set(members))
        new_members = []
        for member in members:
            if member not in current_members:
                new_members.append(member)
        self.model.add_nodes_to_node_set(id_, new_members)

    def _test_create_node_set(self):
        if not self.model.nodes:
            return False
        for _ in range(random.randint(1, self.maximum_objects)):
            members = [random.randint(0, len(self.model.nodes) - 1) for _ in range(20)]
            members = sorted(set(members))
            self.model.create_node_set(self._new_node_set_id(), members)
        self._truncate_node_sets()

    def _test_create_node_set_field(self):
        for _ in range(random.randint(1, self.maximum_objects)):
            id_ = self._random_node_set_id()
            if id_ is None:
                return False
            if not self.model.timesteps:
                self.model.create_timestep(0.0)
            self.model.create_node_set_field(
                self._new_node_set_field_name(), id_, random.random()
            )
        self._truncate_node_set_fields()

    def _test_create_node_set_from_side_set(self):
        id_ = self._random_side_set_id()
        if id_ is not None:
            self.model.create_node_set_from_side_set(self._new_node_set_id(), id_)
        self._truncate_node_sets()

    def _test_create_nodes(self):
        new_nodes = []
        for _ in range(7):
            new_nodes.append([random.random() for _ in range(3)])
        self.model.create_nodes(new_nodes)

    def _test_create_side_set(self):
        if not self.model.element_blocks:
            return False
        for _ in range(random.randint(1, self.maximum_objects)):
            members = []
            for _ in range(20):
                id_ = _random_element(self.model.get_element_block_ids())
                element_count = self.model.get_element_count(id_)
                if element_count == 0:
                    continue
                element_index = random.randint(0, element_count - 1)
                side_index = random.randint(0, 3)
                members.append((id_, element_index, side_index))
            members = sorted(set(members))
            self.model.create_side_set(self._new_side_set_id(), members)
        self._truncate_side_sets()

    def _test_add_faces_to_side_set(self):
        self._truncate_side_sets()
        if not self.model.element_blocks:
            return False
        if not self.model.side_sets:
            return False
        members = []
        for _ in range(20):
            id_ = _random_element(self.model.get_element_block_ids())
            element_count = self.model.get_element_count(id_)
            if element_count == 0:
                continue
            element_index = random.randint(0, element_count - 1)
            side_index = random.randint(0, 3)
            members.append((id_, element_index, side_index))
        members = sorted(set(members))
        id_ = self._random_side_set_id()
        current_members = self.model.get_side_set_members(id_)
        new_members = []
        for member in members:
            if member not in current_members:
                new_members.append(member)
        self.model.add_faces_to_side_set(id_, new_members)

    def _test_create_side_set_field(self):
        for _ in range(random.randint(1, self.maximum_objects)):
            id_ = self._random_side_set_id()
            if id_ is None:
                return False
            if not self.model.timesteps:
                self.model.create_timestep(0.0)
            self.model.create_side_set_field(
                self._new_side_set_field_name(), id_, random.random()
            )
        self._truncate_side_set_fields()

    def _test_create_timestep(self):
        for _ in range(random.randint(1, self.maximum_objects)):
            self.model.create_timestep(self._new_timestep())
        self._truncate_timesteps()

    def _test_copy_timestep(self):
        step = self._random_timestep()
        if step is None:
            return False
        self.model.copy_timestep(step, self._new_timestep())
        self._truncate_timesteps()

    def _test_delete_element_block(self):
        ids = self._random_element_block_ids()
        if not ids:
            return False
        self.model.delete_element_block(ids)

    def _test_delete_element_field(self):
        names = _random_subset(self.model.get_element_field_names())
        if not names:
            return False
        self.model.delete_element_field(names)

    def _test_delete_empty_node_sets(self):
        self.model.create_node_set(self._new_node_set_id(), [])
        self.model.delete_empty_node_sets()

    def _test_delete_empty_side_sets(self):
        self.model.create_side_set(self._new_side_set_id(), [])
        self.model.delete_empty_side_sets()

    def _test_delete_global_variable(self):
        names = _random_subset(self.model.get_global_variable_names())
        if not names:
            return False
        self.model.delete_global_variable(names)

    def _test_delete_node_field(self):
        names = _random_subset(self.model.get_node_field_names())
        if not names:
            return False
        self.model.delete_node_field(names)

    def _test_delete_node_set(self):
        ids = self._random_node_set_ids()
        if not ids:
            return False
        self.model.delete_node_set(ids)

    def _test_delete_side_set(self):
        ids = self._random_side_set_ids()
        if not ids:
            return False
        self.model.delete_side_set(ids)

    def _test_delete_timestep(self):
        times = _random_subset(self.model.get_timesteps())
        if not times:
            return False
        self.model.delete_timestep(times)

    def _test_delete_unused_nodes(self):
        self.model.create_nodes([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        node_count = len(self.model.nodes)
        self.model.delete_unused_nodes()
        assert len(self.model.nodes) <= node_count - 2

    def _test_displace_element_blocks(self):
        ids = self._random_element_block_ids()
        if not ids:
            return False
        if not self.model.timesteps:
            return False
        self.model.unmerge_element_blocks()
        self.model.displace_element_blocks(ids, _random_vector())

    def _test_displacement_field_exists(self):
        self.model.displacement_field_exists()

    def _test_element_block_exists(self):
        id_ = self._random_element_block_id()
        if id_ is None:
            assert not self.model.element_block_exists(0)
        else:
            assert self.model.element_block_exists(id_)

    def _test_element_field_exists(self):
        id_ = self._random_element_block_id()
        if not id_:
            return False
        name = _random_element(self.model.get_element_field_names(id_))
        if name is None:
            assert not self.model.element_field_exists("temp", id_)
        else:
            assert self.model.element_field_exists(name, id_)

    def _test_export_model(self):
        if self.remaining_io_tests <= 0:
            return False
        self.remaining_io_tests -= 1
        self.model.export_model("temp.e")
        # make sure the exported model is equal
        model2 = exomerge.import_model("temp.e")
        model2.qa_records = model2.qa_records[:-1]
        one = self.model
        two = model2
        assert one.nodes == two.nodes
        assert one.timesteps == two.timesteps
        # There is some odd bug here where characters get added sometimes
        # but not often.  To avoid this, we don't compare info records.  This
        # does not appear to be a bug within exomerge.py or exodus.py.  It may
        # be within ctypes, but more likely is within the exodus library
        # itself.
        # exported: ['asdf', 'begin', '  begin statement', ...]
        # imported: ['asdf\x7f\xf8', 'begin', '  begin statement', ...]
        # assert one.info_records == two.info_records
        assert one.title == two.title
        assert one.qa_records == two.qa_records
        assert one.title == two.title
        assert compares_equal_with_nan(one.node_fields, two.node_fields)
        assert compares_equal_with_nan(one.global_variables, two.global_variables)
        assert compares_equal_with_nan(one.side_sets, two.side_sets)
        assert compares_equal_with_nan(one.node_sets, two.node_sets)
        # There is a bug/feature within exodus.py where element types are
        # renamed to uppercase.  To avoid this, we don't compare the info
        # field within each block
        for id_ in list(one.element_blocks.keys()):
            assert compares_equal_with_nan(
                one.element_blocks[id_][0], two.element_blocks[id_][0]
            )
            # assert compares_equal_with_nan(one.element_blocks[id_][1],
            #                                two.element_blocks[id_][1])
            assert compares_equal_with_nan(
                one.element_blocks[id_][2], two.element_blocks[id_][2]
            )
            assert compares_equal_with_nan(
                one.element_blocks[id_][3], two.element_blocks[id_][3]
            )
        os.remove("temp.e")

    def _test_export_stl_file(self):
        if self.remaining_io_tests <= 0:
            return False
        self.remaining_io_tests -= 1
        self.model.export_stl_file("temp.stl")
        os.remove("temp.stl")

    def _test_export_wrl_model(self):
        if self.remaining_io_tests <= 0:
            return False
        self.remaining_io_tests -= 1
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        if not self.model.element_blocks:
            self._create_hex8_element_block()
        node_field_name = self._new_node_field_name()
        self.model.calculate_node_field(node_field_name + " = X + Y + Z")
        self.model.export_wrl_model("temp.wrl", node_field_name)
        os.remove("temp.wrl")
        os.remove("temp_wrl.e")

    def _test_export(self):
        if self.remaining_io_tests <= 0:
            return False
        self.remaining_io_tests -= 1
        filename = "temp." + _random_element(["e", "stl"])
        self.model.export(filename)
        os.remove(filename)

    def _test_get_element_block_connectivity(self):
        id_ = self._random_element_block_id()
        if id_ is None:
            return False
        self.model.get_element_block_connectivity(id_)

    def _test_get_nodes_per_element(self):
        id_ = self._random_element_block_id()
        if id_ is None:
            return False
        self.model.get_nodes_per_element(id_)

    def _test_get_element_block_ids(self):
        self.model.get_element_block_ids()

    def _test_get_element_count(self):
        ids = self._random_element_block_ids()
        if ids is None:
            return False
        self.model.get_element_count(ids)

    def _test_get_element_field_names(self):
        ids = self._random_element_block_ids()
        if ids is None:
            return False
        self.model.get_element_field_names()

    def _test_get_element_field_values(self):
        field = self._random_element_field_name()
        if field is None:
            return False
        if not self.model.timesteps:
            return False
        self.model.get_element_field_values(field[1], field[0])

    def _test_get_global_variable_names(self):
        self.model.get_global_variable_names()

    def _test_get_node_field_names(self):
        self.model.get_node_field_names()

    def _test_get_node_field_values(self):
        name = self._random_node_field_name()
        if name is None:
            return False
        if not self.model.timesteps:
            return False
        self.model.get_node_field_values(name)

    def _test_get_node_set_field_names(self):
        self.model.get_node_set_field_names()

    def _test_get_node_set_ids(self):
        self.model.get_node_set_ids()

    def _test_get_node_set_members(self):
        id_ = self._random_node_set_id()
        if id_ is None:
            return False
        self.model.get_node_set_members(id_)

    def _test_get_nodes_in_element_block(self):
        ids = self._random_element_block_ids()
        if ids is None:
            return False
        self.model.get_nodes_in_element_block(ids)

    def _test_get_nodes_in_side_set(self):
        id_ = self._random_side_set_id()
        if id_ is None:
            return False
        self.model.get_nodes_in_side_set(id_)

    def _test_get_side_set_ids(self):
        self.model.get_side_set_ids()

    def _test_get_side_set_field_names(self):
        self.model.get_side_set_field_names()

    def _test_node_set_field_exists(self):
        info = self._random_node_set_field_name()
        if info is None:
            return False
        id_, name = info
        assert self.model.node_set_field_exists(name, id_)

    def _test_get_node_set_field_values(self):
        if not self.model.timesteps:
            return False
        info = self._random_node_set_field_name()
        if info is None:
            return False
        id_, name = info
        self.model.get_node_set_field_values(name, id_)

    def _test_get_side_set_field_values(self):
        if not self.model.timesteps:
            return False
        info = self._random_side_set_field_name()
        if info is None:
            return False
        id_, name = info
        self.model.get_side_set_field_values(name, id_)

    def _test_side_set_field_exists(self):
        info = self._random_side_set_field_name()
        if info is None:
            return False
        id_, name = info
        assert self.model.side_set_field_exists(name, id_)

    def _test_get_timesteps(self):
        self.model.get_timesteps()

    def _test_global_variable_exists(self):
        name = self._random_global_variable_name()
        if name is None:
            assert not self.model.global_variable_exists("temp")
        else:
            assert self.model.global_variable_exists(name)

    def _test_import_model(self):
        if self.remaining_io_tests <= 0:
            return False
        self.remaining_io_tests -= 1
        new_id = random.randint(1, 5)
        if self.model.element_block_exists(new_id):
            self.model.delete_element_block(new_id)
        input_dir = os.path.dirname(__file__)
        temp_exo_path = os.path.join(input_dir, "exomerge_unit_test.e")

        self.model.import_model(
            temp_exo_path,
            element_block_ids=new_id,
            side_set_ids="none",
            node_set_ids="none",
            global_variable_names="none",
            element_field_names="none",
            node_field_names="none",
            node_set_field_names="none",
            side_set_field_names="none",
            timesteps="none",
        )
        self._truncate_element_blocks()

    def _test_node_field_exists(self):
        name = self._random_node_field_name()
        if name is None:
            assert not self.model.node_field_exists("temp")
        else:
            assert self.model.node_field_exists(name)

    def _test_node_set_exists(self):
        id_ = self._random_node_set_id()
        if id_ is None:
            assert not self.model.node_set_exists(0)
        else:
            assert self.model.node_set_exists(id_)

    def _test_process_element_fields(self):
        if not self.model.element_blocks:
            return False
        if not self.model.timesteps:
            self.model.create_timestep(0.0)
        self.model.delete_element_field("all")
        self.model.create_element_field(
            self.model._new_element_field_name(), self._random_element_block_id()
        )
        for name in ["eqps_" + str(x) for x in range(1, 9)]:
            self.model.create_element_field(name, "all")
        self.model.create_element_field(
            self.model._new_element_field_name(), self._random_element_block_id()
        )
        for name in ["theta_" + str(x) for x in range(1, 10)]:
            self.model.create_element_field(name, "all")
        self.model.create_element_field(
            self.model._new_element_field_name(), self._random_element_block_id()
        )
        if self.model.get_node_field_names():
            self.model.delete_node_field("all")
        self.model.process_element_fields()
        self._truncate_node_fields()

    def _test_rename_element_field(self):
        field = self._random_element_field_name()
        if field is None:
            return False
        self.model.rename_element_field(
            field[1], self._new_element_field_name(), field[0]
        )

    def _test_scale_element_blocks(self):
        ids = self._random_element_block_ids()
        if ids is None:
            return False
        self.model.unmerge_element_blocks()
        self.model.scale_element_blocks(ids, _random_scalar())

    def _test_rotate_element_blocks(self):
        ids = self._random_element_block_ids()
        if ids is None:
            return False
        self.model.unmerge_element_blocks()
        self.model.rotate_element_blocks(ids, [1, 0, 0], 90)

    def _test_rotate_geometry(self):
        self.model.rotate_geometry([1, 0, 0], 90)

    def _test_scale_geometry(self):
        self.model.scale_geometry(_random_scalar())

    def _test_side_set_exists(self):
        id_ = self._random_side_set_id()
        if id_ is None:
            assert not self.model.side_set_exists(0)
        else:
            assert self.model.side_set_exists(id_)

    def _test_summarize(self):
        with OutputSuppression():
            self.model.summarize()

    def _test_timestep_exists(self):
        timestep = _random_element(self.model.timesteps)
        if timestep is None:
            assert not self.model.timestep_exists(0.0)
        else:
            assert self.model.timestep_exists(timestep)

    def _test_translate_element_blocks(self):
        ids = self._random_element_block_ids()
        if not ids:
            return False
        self.model.unmerge_element_blocks()
        self.model.translate_element_blocks(ids, _random_vector())

    def _test_translate_geometry(self):
        self.model.translate_geometry(_random_vector())

    def _test_threshold_element_blocks(self):
        id_ = self._create_hex8_element_block()
        assert self.model.get_element_count(id_) == 2
        self.model.calculate_element_centroids(element_block_ids=[id_])
        self.model.threshold_element_blocks("centroid_x>1", element_block_ids=[id_])
        assert self.model.get_element_count(id_) == 1

    def _test_count_degenerate_elements(self):
        ids = self._random_element_block_ids()
        if not ids:
            return False
        self.model.count_degenerate_elements(ids)

    def _test_delete_duplicate_elements(self):
        id_ = self._random_element_block_id()
        if id_ is None:
            return False
        new_id = self._new_element_block_id()
        element_count = self.model.get_element_count(id_)
        if element_count == 0:
            return False
        self.model.duplicate_element_block(id_, new_id, duplicate_nodes=False)
        self.model.combine_element_blocks([id_, new_id], target_element_block_id=id_)
        self.model.delete_duplicate_elements(id_)
        assert self.model.get_element_count(id_) <= element_count

    def _test_get_closest_node_distance(self):
        self.model.get_closest_node_distance()

    def _test_get_connectivity(self):
        id_ = self._random_element_block_id()
        if id_ is None:
            return False
        self.model.get_connectivity(id_)

    def _test_get_element_block_dimension(self):
        id_ = self._random_element_block_id()
        if id_ is None:
            return False
        self.model.get_element_block_dimension(id_)

    def _test_get_element_block_extents(self):
        ids = self._random_element_block_ids()
        if not ids:
            return False
        self.model.get_element_block_extents(ids)

    def _test_get_element_block_name(self):
        id_ = self._random_element_block_id()
        if id_ is None:
            return False
        self.model.get_element_block_name(id_)

    def _test_get_all_element_block_names(self):
        self.model.get_all_element_block_names()

    def _test_get_element_edge_length_info(self):
        ids = self._random_element_block_ids()
        if not ids:
            return False
        self.model.get_element_edge_length_info(ids)

    def _test_get_node_set_name(self):
        id_ = self._random_node_set_id()
        if id_ is None:
            return False
        self.model.get_node_set_name(id_)

    def _test_get_all_node_set_names(self):
        self.model.get_all_node_set_names()

    def _test_get_side_set_area(self):
        id_ = self._random_side_set_id()
        if id_ is None:
            return False
        self.model.get_side_set_area(id_)

    def _test_get_side_set_members(self):
        id_ = self._random_side_set_id()
        if id_ is None:
            return False
        self.model.get_side_set_members(id_)

    def _test_get_side_set_name(self):
        id_ = self._random_side_set_id()
        if id_ is None:
            return False
        self.model.get_side_set_name(id_)

    def _test_get_all_side_set_names(self):
        self.model.get_all_side_set_names()

    def _test_build_hex8_cube(self):
        new_id = self._new_element_block_id()
        self.model.build_hex8_cube(new_id)

    def _test_count_disconnected_blocks(self):
        ids = self._random_element_block_ids()
        if not ids:
            return False
        self.model.count_disconnected_blocks(ids)

    # End of unit test functions.

    @staticmethod
    def _source_calls_target(source, target):
        """
        Return True if the source function calls target somewhere.

        """
        source_code = inspect.getsource(source)
        return bool(
            re.search("[^A-Za-z0-9_]" + target.__name__ + "[ \t\n\r]*\(", source_code)
        )

    def test(self):
        """
        Run the unit tests.

        """
        random.seed(0)
        exomerge.EXIT_ON_WARNING = True
        # get list of all unit tests
        tests = inspect.getmembers(self, inspect.ismethod)
        unit_tests = []
        for test in tests:
            if test[0].startswith("_test_"):
                unit_tests.append(test)
        # get a list of all public functions in exomerge
        public_functions = []

        for function, _ in inspect.getmembers(exomerge.ExodusModel, inspect.isfunction):
            if not function.startswith("_"):
                public_functions.append(function)
        print(
            (
                "We found %d unit tests and %d public functions."
                % (len(unit_tests), len(public_functions))
            )
        )
        # If a test exists that doesn't match a public function, issue a
        # warning message and remove that unit test.
        print("\nVerifying each unit test matches a public exomerge function.")
        unmatched = []
        matched_unit_tests = []
        for unit_test in unit_tests:
            (test, _) = unit_test
            if test[6:] not in public_functions:
                unmatched.append(test)
            else:
                matched_unit_tests.append(unit_test)
        if unmatched:
            print(
                (
                    "\nWARNING: Found %d unit test(s) without a matching "
                    "public function\nUnit tests:\n  %s"
                    % (len(unmatched), "\n  ".join(unmatched))
                )
            )
            print("")
        unit_tests = matched_unit_tests
        # If a public function exists without a unit test, issue a warning
        # message
        print("\nVerifying public function has a matching unit test.")
        unmatched = []
        unit_test_names = [x[0] for x in unit_tests]
        for name in public_functions:
            if "_test_" + name not in unit_test_names:
                unmatched.append(name)
        if unmatched:
            print(
                (
                    "\nWARNING: Found %d public functions without a matching "
                    "unit test\nPublic functions:\n  %s"
                    % (len(unmatched), "\n  ".join(unmatched))
                )
            )
            print("")
        # make sure unit tests call the appropriate function
        # i.e. _test_create_nodes better call create_nodes somewhere
        print("\nVerifying each unit test calls the appropriate function.")
        bad_unit_tests = []
        for name, function in unit_tests:
            target = eval("exomerge.ExodusModel." + name[6:])
            if not self._source_calls_target(function, target):
                bad_unit_tests.append(name)
        if bad_unit_tests:
            print(
                (
                    "\nWARNING: Found %d unit tests which do not call their "
                    "corresponding public functions:\n  %s"
                    % (len(bad_unit_tests), "\n  ".join(bad_unit_tests))
                )
            )
            print("")
        # start off the model
        input_dir = os.path.dirname(__file__)
        temp_exo_path = os.path.join(input_dir, "exomerge_unit_test.e")
        self.model = exomerge.import_model(temp_exo_path)
        # run for the given amount of walltime or the number of tests
        # is sufficient
        start_time = time.time()
        end_time = time.time() + self.time_limit
        tests = 0
        next_tests_to_run = []
        passed_tests = set()
        while tests < self.max_tests and (
            time.time() < end_time or tests < self.min_tests
        ):
            if not next_tests_to_run:
                for _ in range(5):
                    next_tests_to_run.extend(copy.copy(unit_tests))
                random.shuffle(next_tests_to_run)
            tests += 1
            test = next_tests_to_run[0]
            del next_tests_to_run[0]
            print(("[%d] %s" % (tests, test[0])))
            if test[1]() is None:
                passed_tests.add(test[0])
            if tests % 123 == 0:
                self.model._verify(allow_aliased_lists=False)
        # if some tests did not successfully run, not that
        if len(passed_tests) != len(unit_tests):
            untested = []
            for x, _ in unit_tests:
                if x not in passed_tests:
                    untested.append(x)
            print(
                (
                    "\nWARNING: Some unit tests were unable to be run:\n  "
                    + "\n  ".join(untested)
                )
            )
        print(("\nRan %d tests in %g seconds." % (tests, time.time() - start_time)))
        print("\nSuccess")


# if this module is executed (as opposed to imported), run the tests
if __name__ == "__main__":
    if len(sys.argv) > 2:
        sys.stderr.write("Invalid syntax.\n")
        exit(1)
    tester = ExomergeUnitTester()
    tester._topology_test()
    if len(sys.argv) == 2:
        tester.min_tests = int(sys.argv[1])
        tester.max_tests = tester.min_tests
    tester.test()
