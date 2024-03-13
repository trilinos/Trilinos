# @HEADER
# ***********************************************************************
#
#          PyTrilinos2: Automatic Python Interfaces to Trilinos Packages
#                 Copyright (2022) Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia
# Corporation, the U.S. Government retains certain rights in this
# software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Kim Liegeois (knliege@sandia.gov)
#
# ***********************************************************************
# @HEADER

import unittest
from mpi4py import MPI
from PyTrilinos2.PyTrilinos2 import Teuchos

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