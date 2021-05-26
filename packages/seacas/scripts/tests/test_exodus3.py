"""
Copyright(C) 1999-2021 National Technology & Engineering Solutions
of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
NTESS, the U.S. Government retains certain rights in this software.

See packages/seacas/LICENSE for details

"""

import unittest
import sys
import os
import tempfile

path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
sys.path.append(os.path.join(path, "lib"))
import exodus as exo


class MyTestCase(unittest.TestCase):

    def setUp(self):
        exofile = exo.exodus("test-assembly.exo", mode='r')
        self.tempdir = tempfile.TemporaryDirectory()
        self.temp_exo_path = os.path.join(self.tempdir.name, "temp-test-assembly.exo")
        self.temp_exofile = exofile.copy(self.temp_exo_path)
        exofile.close()

    def tearDown(self):
        self.tempdir.cleanup()

    def test_ex_obj_to_inq(self):
        self.assertEqual("EX_INQ_ASSEMBLY", exo.ex_obj_to_inq('EX_ASSEMBLY'))

    def test_ex_obj_to_inq_wrong_obj_type(self):
        self.assertEqual(-1, exo.ex_obj_to_inq('fake_obj'))

    def test_ex_entity_type_to_objType(self):
        self.assertEqual("assembly", exo.ex_entity_type_to_objType(exo.get_entity_type('EX_ASSEMBLY')))

    def test_ex_entity_type_to_objType_wrong_entity_type(self):
        self.assertEqual("EX_INVALID", exo.ex_entity_type_to_objType('fake_entity'))

    def test_get_name_assembly(self):
        temp_exofile = exo.exodus("test-assembly.exo", mode='r')
        assembly_ids = temp_exofile.get_ids("EX_ASSEMBLY")
        names = [temp_exofile.get_name("EX_ASSEMBLY", assembly) for assembly in assembly_ids]
        self.assertEqual(['Root', 'Child2', 'Child3', 'Child4', 'NewAssembly', 'FromPython'], names)

    def test_inquire_assembly(self):
        temp_exofile = exo.exodus("test-assembly.exo", mode='r')
        assem_count = temp_exofile.inquire("EX_INQ_ASSEMBLY")
        self.assertEqual(6, assem_count)

    def test_get_assembly(self):
        temp_exofile = exo.exodus("test-assembly.exo", mode='r')
        assembly_ids = temp_exofile.get_ids("EX_ASSEMBLY")
        assemblies = [temp_exofile.get_assembly(assembly) for assembly in assembly_ids]
        root = exo.assembly(name='Root', type='assembly', id=100)
        root.entity_list = [100, 200, 300, 400]
        self.assertEqual(str(root), str(assemblies[0]))

    def test_get_assemblies(self):
        temp_exofile = exo.exodus("test-assembly.exo", mode='r')
        assembly_ids = temp_exofile.get_ids("EX_ASSEMBLY")
        assemblies = temp_exofile.get_assemblies(assembly_ids)
        expected = [exo.assembly(name='Root', type='assembly', id=100),
                    exo.assembly(name='Child2', type='element block', id=200),
                    exo.assembly(name='Child3', type='element block', id=300),
                    exo.assembly(name='Child4', type='element block', id=400),
                    exo.assembly(name='NewAssembly', type='assembly', id=222),
                    exo.assembly(name='FromPython', type='assembly', id=333)]
        for i, x in enumerate(expected):
            entity_lists = [[100, 200, 300, 400],
                            [10, 11, 12, 13],
                            [14, 15, 16],
                            [10, 16],
                            [100, 200, 300, 400],
                            [100, 222]]
            x.entity_list = entity_lists[i]
        self.maxDiff=None
        self.assertEqual(str(expected), str(assemblies))

    def test_put_assembly(self):
        new = exo.assembly(name='Unit_test', type='EX_ASSEMBLY', id=444)
        new.entity_list = [100, 222]
        self.temp_exofile.put_assembly(new)
        self.temp_exofile.close()
        temp_exofile = exo.exodus(self.temp_exo_path)
        assembly_ids = temp_exofile.get_ids("EX_ASSEMBLY")
        assemblies = [temp_exofile.get_assembly(assembly) for assembly in assembly_ids]
        self.assertEqual(new.name, assemblies[6].name)
        self.assertEqual('assembly', assemblies[6].type)
        self.assertEqual(new.id, assemblies[6].id)
        self.assertEqual(new.entity_list, assemblies[6].entity_list)

    def test_get_reduction_variables_assembly(self):
        temp_exofile = exo.exodus("test-assembly.exo", mode='r')
        red_var = temp_exofile.get_reduction_variable_number("EX_ASSEMBLY")
        self.assertEqual(4, red_var)
        names = temp_exofile.get_reduction_variable_names("EX_ASSEMBLY")
        self.assertIn("Momentum_X", names)
        self.assertIn("Momentum_Y", names)
        self.assertIn("Momentum_Z", names)
        self.assertIn("Kinetic_Energy", names)

    def test_get_reduction_variable_values_assembly(self):
        temp_exofile = exo.exodus("test-assembly.exo", mode='r')
        assembly_ids = temp_exofile.get_ids("EX_ASSEMBLY")
        assemblies = [temp_exofile.get_assembly(assembly) for assembly in assembly_ids]
        values = temp_exofile.get_reduction_variable_values('EX_ASSEMBLY', assemblies[0].id, 1)
        self.assertListEqual([0.02, 0.03, 0.04, 0.05], list(values))

    def test_get_reduction_variable_values_assembly_no_values(self):
        temp_exofile = exo.exodus("test-assembly.exo", mode='r')
        times = temp_exofile.get_times()
        assembly_ids = temp_exofile.get_ids("EX_ASSEMBLY")
        assemblies = [temp_exofile.get_assembly(assembly) for assembly in assembly_ids]
        values = temp_exofile.get_reduction_variable_values('EX_ASSEMBLY', assemblies[5].id, 1)
        self.assertListEqual([0.00, 0.00, 0.00, 0.00], list(values))


    def test_put_assemblies(self):
        new = exo.assembly(name='Unit_test', type='EX_ASSEMBLY', id=444)
        new.entity_list = [100, 222]
        self.temp_exofile.put_assemblies([new])
        self.temp_exofile.close()
        temp_exofile = exo.exodus(self.temp_exo_path)
        assembly_ids = temp_exofile.get_ids("EX_ASSEMBLY")
        assemblies = [temp_exofile.get_assembly(assembly) for assembly in assembly_ids]
        print(assemblies)
        self.assertEqual(new.name, assemblies[6].name)
        self.assertEqual('assembly', assemblies[6].type)
        self.assertEqual(new.id, assemblies[6].id)
        self.assertEqual(new.entity_list, assemblies[6].entity_list)


    def test_put_assemblies_multiple_assemblies(self):
        new = exo.assembly(name='Unit_test1', type='EX_ASSEMBLY', id=444)
        new.entity_list = [100, 222]
        new2 = exo.assembly(name='Unit_test2', type='EX_ASSEMBLY', id=448)
        new2.entity_list = [102, 224]
        self.temp_exofile.put_assemblies([new, new2])
        self.temp_exofile.close()
        temp_exofile = exo.exodus(self.temp_exo_path)
        assembly_ids = temp_exofile.get_ids("EX_ASSEMBLY")
        assemblies = [temp_exofile.get_assembly(assembly) for assembly in assembly_ids]
        print(assemblies)
        self.assertEqual(new.name, assemblies[6].name)
        self.assertEqual('assembly', assemblies[6].type)
        self.assertEqual(new.id, assemblies[6].id)
        self.assertEqual(new.entity_list, assemblies[6].entity_list)

        self.assertEqual(new2.name, assemblies[7].name)
        self.assertEqual('assembly', assemblies[7].type)
        self.assertEqual(new2.id, assemblies[7].id)
        self.assertEqual(new2.entity_list, assemblies[7].entity_list)


if __name__ == '__main__':
    unittest.main()
