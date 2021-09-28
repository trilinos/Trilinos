include(TribitsAddInstallGroupAndPermsFixups)

tribits_add_install_group_and_perms_fixups()

# NOTE: This is a workaround for a problem with an issue with CMake where you
# must call install(SCRIPT ...) in a subdirectory added with
# add_subdirectory() in order for CMake to run that install script after all
# of the other subdirectories (which are TriBITS packages).  As of CMake 3.14,
# you can set the policy set_policy(CMP0082 NEW) which would have allowed us
# to put this in the base tribits_project_impl() macro.  But since TriBITS is
# not allowed to require CMake 3.14 yet, we must use this workaround.  To make
# super sure that this install(SCRIPTS ...) script will get called last, the
# policy CMP0082 is set to NEW in the TribitsCMakePolicies.cmake file.  In
# automated testing of TriBITS, this seems to show that this custom script
# runs after all of the other files get installed (even with CMake 3.11) and
# passes the automated tests that ensures that the last file installed is
# given the correct permissions!  At least that is what happened with TriBITS.
