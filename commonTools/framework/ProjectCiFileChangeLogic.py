# Specialized logic for what file changes should trigger a global build in CI
# testing where testing should only involve packages impacted by the change.
#
class ProjectCiFileChangeLogic:

  def isGlobalBuildFileRequiringGlobalRebuild(self, modifiedFileFullPath):
    modifiedFileFullPathArray = modifiedFileFullPath.split('/')
    lenPathArray = len(modifiedFileFullPathArray)

    if lenPathArray==1:  # Base project directory
      # Files directly under <projectDir>/
      if modifiedFileFullPathArray[0] == "CMakeLists.txt":
        return True
      elif modifiedFileFullPathArray[0].rfind(".cmake") != -1:
        # Any <projectDir>/*.cmake files should trigger testing all packages
        return True
      else:
        # Other files under <projectDir>/ likely don't require a global
        # rebuild
        return False
    elif modifiedFileFullPathArray[0] == 'cmake':
      # cmake/
      if lenPathArray==2:
        # Files directly under cmake/
        if modifiedFileFullPathArray[1]=='ExtraRepositoriesList.cmake':
          return False
        elif modifiedFileFullPathArray[1].rfind(".cmake") != -1:
          # All CMake files directly under cmake/ should trigger
          # testing of all packages.
          return True
        else:
          # Any other files directly under cmake/ don't require a global
          # rebuild
          return False
      elif modifiedFileFullPathArray[1] == 'std':
        # cmake/std/
        if lenPathArray==3:
          # Files directly under Trilinos/cmake/std/.  These may be used for
          # CI or PR testing so we need to test all packages if any of these
          # files change.
          return True
        elif modifiedFileFullPathArray[2] == "atdm":
          # Files under Trilinos/cmake/std/atdm/.  Currently, none of these
          # files are being used for CI or PR testing so we don't need to
          # enable all packages when these change.
          return False
        elif modifiedFileFullPathArray[2] == "sems":
          # Files under Trilinos/cmake/std/sems/.  Currently, none of these
          # files are being used for CI or PR testing so we don't need to
          # enable all packages when these change.
          return False
      elif modifiedFileFullPathArray[1] == 'TPLs':
        # cmake/TPLs/.  Any TPL file modules that change really needs to
        # trigger a global build to be safe
        return True
      elif modifiedFileFullPathArray[1] == 'ctest':
        # cmake/ctest/
        if lenPathArray==3:
          # Files directly under <project>/cmake/ctest/.  These could impact
          # CI and PR testing if they use the TriBITS CTest driver so to be
          # safe, we should test all packages if these changes.
          return True
        elif modifiedFileFullPathArray[2] == 'drivers':
          # cmake/ctest/drivers/
          if lenPathArray==4:
            # Files directly under <project>/cmake/ctest/drivers/.
            # These files could be used by CI or PR testing, and they don't
            # change very often, so changes to these files should trigger a
            # global rebuild.
            return True
          else:
            # This is a file under
            # cmake/ctest/<something>/[...something...] so this
            # is likely build files for a specific system that CI and PR
            # testing does not depend on.  Therefore, we should not need to
            # test all packages.
            return False
      elif modifiedFileFullPath.rfind(".cmake") != -1:
        # All other *.cmake files under any subdir of cmake/ should trigger
        # a global rebuild to be safe.
       return True
    elif lenPathArray >= 2 and modifiedFileFullPathArray[0] == 'packages' and \
      modifiedFileFullPathArray[1] == 'framework' \
      :
      # Changes under packages/framework/ likely impact the GenConfig PR build
      # configurations and therefore to be safe, everything needs to be tested.
      return True
    elif lenPathArray >= 2 and modifiedFileFullPathArray[0] == '.github' and \
      modifiedFileFullPathArray[1] == 'workflows' \
      :
      # Changes under .github/workflows/ impact CI-type runs on GitHub Actions
      # and therefore to be safe, everything needs to be tested.
      return True
    # Any other files not already covered above should *not* trigger a global
    # build
    return False
