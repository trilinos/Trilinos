from FindTribitsCiSupportDir import *
import GeneralScriptSupport as GSS


# Standard set of build stats fields we want to read in
#
def getStdBuildStatsColsAndTypesList():
  return [
    ColNameAndType('max_resident_size_Kb', 'float'),
    ColNameAndType('elapsed_real_time_sec', 'float'),
    ColNameAndType('FileName', 'string'),
    ColNameAndType('FileSize', 'float'),
    ]
# NOTE: Above, we use type 'float' instead of 'int' for fields that are ints
# because we want to allow a very large size.


def getColNameTypeIdxListGivenColNameAndTypeList(csvFileName, columnHeadersList,
    colNameAndTypesToGetList,
  ):
  colNameTypeIdxList = []
  for colNameAndTypeToGet in colNameAndTypesToGetList:
    colIdx = GSS.findInSequence(columnHeadersList, colNameAndTypeToGet.colName())
    if colIdx != -1:
      colNameTypeIdxList.append(ColNameTypeIdx(colNameAndTypeToGet, colIdx))
    else:
      raise Exception(
        "Error, the CSV file column header '"+colNameAndTypeToGet.colName()+"'"+\
        " does not exist in the list of column headers "+str(columnHeadersList)+\
        " from the CSV file '"+csvFileName+"'!")
  return colNameTypeIdxList


class ColNameAndType(object):
  def __init__(self, colName, colType):
    self.__colName = colName
    self.__colType = colType
    self.assertType()
  def colName(self):
    return self.__colName
  def colType(self):
    return self.__colType
  def __repr__(self):
    myStr = "ColNameAndType{"+self.__colName+","+str(self.__colType)+"}"
    return myStr
  def convertFromStr(self, strIn):
    if self.__colType == "string":
      return strIn
    elif self.__colType == "int":
      return int(strIn)
    elif self.__colType == "float":
      return float(strIn)
  def assertType(self):
    supportedTypes = [ "string", "int", "float" ]
    if -1 == GSS.findInSequence(supportedTypes, self.__colType):
      raise Exception(
        "Error, type '"+str(self.__colType)+"' is not supported!  Supported types"+\
        " include "+str(supportedTypes)+"!")
  def __eq__(self, other):
    return((self.__colName,self.__colType)==(other.__colName,other.__colType))
  def __ne__(self, other):
    return((self.__colName,self.__colType)!=(other.__colName,other.__colType))


class ColNameTypeIdx(object):
  def __init__(self, colNameAndType, colIdx):
    self.__colNameAndType = colNameAndType
    self.__colIdx = colIdx
  def colName(self):
    return self.__colNameAndType.colName()
  def getIdx(self):
    return self.__colIdx
  def convertFromStr(self, strIn):
    return self.__colNameAndType.convertFromStr(strIn)
  def __repr__(self):
    myStr = "ColNameTypeIdx{"+str(self.__colNameAndType)+","+str(self.__colIdx)+"}"
    return myStr
  def __eq__(self, other):
    return ((self.__colNameAndType,self.__colIdx)==(other.__colNameAndType,other.__colIdx))
  def __ne__(self, other):
    return ((self.__colNameAndType,self.__colIdx)!=(other.__colNameAndType,other.__colIdx))
