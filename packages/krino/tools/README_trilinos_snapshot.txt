This script, trilinos_snapshot.sh, is used to sync sierra krino to Trilinos.

It is executed using:
SIERRA=/path/sierra/code TRILINOS=/path/trilinos trilinos_snapshot.sh

For my current setup, 
SIERRA=/fgs/drnoble/projects/code
TRILINOS=/fgs/drnoble/projects/Trilinos

I think these are the defaults, so you can just run:
cd /fgs/drnoble/projects/code
krino/tools/trilinos_snapshot.sh

The Trilinos repository in /fgs/drnoble/projects/Trilinos was created off
of a clone of Trilinos using the commands:

cd /fgs/drnoble/projects
git clone git@github.com:drnobleabq/Trilinos
cd Trilinos
git remote add upstream git@github.com:trilinos/Trilinos

Where drnobleabq/Trilinos is a github fork of trilinos/Trilinos.

After this script is run, a PR can be created.  Be sure it is relative
to the develop branch of trilinos/Trilinos.
