#
# Specialized logic for what file changes should trigger a global build in CI
# testing where testing should only occur package impacted by the change.
#

class ProjectCiFileChangeLogic:

  def isGlobalBuildFileRequiringGlobalRebuild(self, modifiedFileFullPath):
    modifiedFileFullPathArray = modifiedFileFullPath.split('/')
    lenPathArray = len(modifiedFileFullPathArray)
    if lenPathArray==1:
      # Files sitting directly under <projectDir>/
      if modifiedFileFullPathArray[0] == "CMakeLists.txt":
        return True
      if modifiedFileFullPathArray[0].rfind(".cmake") != -1:
        return True
    elif modifiedFileFullPathArray[0] == 'cmake':
      # Files under <projectDir>/cmake/
      if modifiedFileFullPathArray[1]=='ExtraRepositoriesList.cmake':
        return False
      elif modifiedFileFullPathArray[1] == 'ctest' and lenPathArray >= 3:
        if lenPathArray > 3:
          # This is a file
          # <projectDir>/cmake/ctest/<something>/[...something...]  so this is
          # for a specific machine and should not trigger a global build.
          return False
        else:
          # Any other file directly under cmake/ctest/ should trigger a global
          # build.
          return True
      else:
        # All other files under cmake/
        if modifiedFileFullPath.rfind(".cmake") != -1:
          # All other *.cmake files under cmake/ trigger a global build.
          return True
    # Any other files should not trigger a global build
    return False
