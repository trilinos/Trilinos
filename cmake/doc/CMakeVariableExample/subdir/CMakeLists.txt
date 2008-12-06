
MESSAGE("Set the value of SOMEVAR in the parent to newName1")
SET(SOMEVAR newName1 PARENT_SCOPE)
MESSAGE("Value of SOMEVAR here is '${SOMEVAR}'")
# Above will print the orginal value 'baseName' that was set in the
# base scope and copied into this scope.  This is because SOMEVAR in
# this scope is separate and distict from SOMEVAR in the base scope.

MESSAGE("Now set the value of SOMEVAR in this scope to newName2")
SET(SOMEVAR newName2)
MESSAGE("Value of SOMEVAR here is '${SOMEVAR}'")
# Above will print 'newName2' and will not affect the value SOMEVAR in
# the base scope.

MESSAGE("Set the value of SOMEGLOBALVAR in global scope to newName3")
SET(SOMEGLOBALVAR newName3 CACHE INTERNAL "")
MESSAGE("Value of SOMEGLOBALVAR here is '${SOMEGLOBALVAR}'")
# Above will set the value 'newName3' that will be set in all scopes,
# including this one.

MESSAGE("Create a local version SOMEGLOBALVAR is this scope set to newName4")
SET(SOMEGLOBALVAR newName4)
MESSAGE("Value of now local SOMEGLOBALVAR here is '${SOMEGLOBALVAR}'")

# Above will set the value 'newName4' and will not affect the separate
# and district global varaible SOMEGLOBALVAR.
