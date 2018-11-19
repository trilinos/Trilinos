This atdm/common/ directory is used to put files that are shared by more than
one env.  For example, the 'chama' and 'serrano' envs share files in the
atdm/common/toss3/ subdir.  This atdm/common/ directory is always installed on
every platform.  But since these files don't represent a complete env and
can't be directly loaded from the base atdm/load-env.sh script, the existence
of the files in this subdir will not allow a user to try to load an env that
is incompatible with the installed configuration of Trilinos.  But yet we
don't need to be constantly updating the install rules in the file
atdm/ATDMDevEnvSettings.cmake to accommodate shared files like these.
