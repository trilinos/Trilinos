#!/usr/bin/env python
import os
import sys
import argparse
import logging
import subprocess


THIS_DIR = os.path.dirname(os.path.realpath(__file__))


def find_trilinos_makefile(trilinos_install):
    """Find the Trilinos makefile.  The makefile is looked for in order of
    preference"""

    makefile = os.path.join(trilinos_install,
                            'Makefile.export.Trilinos_install')
    if os.path.isfile(makefile):
        return makefile

    makefile = os.path.join(trilinos_install,
                            'include/Makefile.export.Trilinos')
    if os.path.isfile(makefile):
        return makefile

    makefile = os.path.join(trilinos_install,
                            'Makefile.export.Trilinos')
    if os.path.isfile(makefile):
        return makefile

    return None


def make_example(trilinos_install, example):
    assert os.path.isdir(trilinos_install)

    # Find the Trilinos makefile
    makefile = find_trilinos_makefile(trilinos_install)
    if makefile is None:
        raise Exception('Could not find Trilinos makefile in {0!r}, '
                        'is Trilinos built?'.format(trilinos_install))

    # The build environment
    env = dict(os.environ)
    env['TRILINOS_INSTALL'] = trilinos_install
    env['TRILINOS_MAKEFILE'] = makefile

    command = ['make', '-f', 'Examples.mk', example]

    logging.info('Making example {0!r}'.format(example))
    p = subprocess.Popen(command, env=env, cwd=THIS_DIR)
    p.wait()
    return p.returncode


def main():
    p = argparse.ArgumentParser()
    p.add_argument('example', help='Name of example file to build')
    p.add_argument('--trilinos-install', default=os.getenv('TRILINOS_INSTALL'),
                   help='Trilinos installation [default: %(default)s]')
    args = p.parse_args()

    if not args.trilinos_install:
        logging.error("The '--trilinos-install' option must be passed to this "
                      "script or the 'TRILINOS_INSTALL' environment variable must "
                      "exist in your environment to run this script.")
        return 1

    trilinos_install = os.path.realpath(args.trilinos_install)
    if not os.path.isdir(trilinos_install):
        logging.error("{0!r} is not a directory".format(trilinos_install))
        return 1

    root, ext = os.path.splitext(os.path.basename(args.example))
    if not os.path.isfile(os.path.join(THIS_DIR, root + '.cpp')):
        logging.error('{0!r} is not a valid example to '
                      'build'.format(args.example))
        return 1

    return make_example(trilinos_install, root)


if __name__ == '__main__':
    sys.exit(main())
