#!/usr/bin/env python

import os
import sys

impl_path = os.path.join('..', 'trilinos-test2')
sys.path.append(os.path.join(os.path.abspath(os.path.dirname(__file__)), impl_path))
import nightly_create_tarball

packages = [
    ('Teuchos', True),
    ('Epetra', True),
]

options, args = nightly_create_tarball.get_options(os.getcwd())
options.enableBoost = False
options.enableNetcdf = False
nightly_create_tarball.main(packages, options)
