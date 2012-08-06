doxygen_ver = "1.4.7"
def checkDoxygenVersion(root):
    if root.get("version") != doxygen_ver:
      msgStr = 'XML files must be generated with Doxygen ' + doxygen_ver + '.'
      raise Exception(msgStr)
