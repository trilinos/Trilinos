import unittest
from PyROL.PyROL import Teuchos

class TestParameterList(unittest.TestCase):
    def test_all(self):
        a = Teuchos.ParameterList()
        b = Teuchos.ParameterList()
        print(a)
        print(a.numParams())
        a.set('relaxation: sweeps', 5)
        a.set('relaxation: damping factor', 0.95)
        print(dir(Teuchos.ParameterList))
        print(a.numParams())
        b.set('Ifpack2 Settings', a)
        print(b)
        print('just before')
        print(b.sublist('Ifpack2 Settings'))
        #b.sublist('Ifpack2 Settings').set('relaxation: sweeps', 6)
        print(b)
        #print(a['relaxation: damping factor'])
        a['relaxation: damping factor'] = .65
        #print(a['relaxation: damping factor'])
        print(a.get('relaxation: sweeps'))
        print(a)
        print(a['relaxation: sweeps'])
        print(a['relaxation: damping factor'])
        print(b.get('Ifpack2 Settings'))

        print(b['Ifpack2 Settings']['relaxation: damping factor'])
        b['Ifpack2 Settings']['relaxation: damping factor'] = 0.5
        b['Ifpack2 Settings']['damped'] = True
        b['Ifpack2 Settings']['not damped'] = False
        print(b['Ifpack2 Settings']['relaxation: damping factor'])
        self.assertEqual(b['Ifpack2 Settings']['relaxation: damping factor'], 0.5)
        self.assertEqual(b['Ifpack2 Settings']['not damped'], False)

if __name__ == '__main__':
    unittest.main()