These scripts are used to sync sierra krino to Trilinos.

They are executed using:
SIERRA=/path/sierra TRILINOS=/path/trilinos cmd.sh

For my current setup, 
SIERRA=/fgs/drnoble/projects/code_votd
TRILINOS=/fgs/drnoble/projects/Trilinos

These are the defaults, so you can just run:
cd
krino/tools/trilinos_prepare_snapshot_branch.sh
krino/tools/trilinos_create_snapshot_commit.sh

If more iterations are needed, make changes in sierra, merge them and rerun
krino/tools/trilinos_create_snapshot_commit.sh

The Trilinos repository in /fgs/drnoble/projects/Trilinos was created off
of a clone of Trilinos using the commands:

cd /fgs/drnoble/projects
git clone git@github.com:drnobleabq/Trilinos
cd Trilinos
git remote add upstream git@github.com:trilinos/Trilinos

Where drnobleabq/Trilinos is a github fork of trilinos/Trilinos.

After the initial commit, a link is printed where a PR can be created.  Be sure it is relative
to the develop branch of trilinos/Trilinos.
