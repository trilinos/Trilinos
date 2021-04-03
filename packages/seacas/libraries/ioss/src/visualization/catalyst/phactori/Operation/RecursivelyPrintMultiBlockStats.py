from phactori import *

def PrintMultiBlockStatsOneBlock(inCsdata, inMetaData, ioFlatIndexCounter):
  icsdClassname = inCsdata.GetClassName()
  myDebugPrint3("leaf block " + str(icsdClassname) + " index " + str(ioFlatIndexCounter[0]) + "\n")
  if icsdClassname == "vtkStructuredGrid":
    wholeExtent = inMetaData.Get(vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT())
    myDebugPrint3("this leaf is structured grid: " + str(ioFlatIndexCounter[0]) + "\n"
      "vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(): " + str(wholeExtent) + "\n"
      "B inCsdata.GetExtent(): " + str(inCsdata.GetExtent()) + "\n")

def PrintMultiBlockStatsRecurse1(inCsdata, inMetaData, recursionLevel, ioFlatIndexCounter):
  myDebugPrint3("PrintMultiBlockStatsRecurse1 entered recursionLevel " + str(recursionLevel) + "\n")

  icsdClassname = inCsdata.GetClassName()
  myDebugPrint3("inCsdata.GetClassName() " + str(icsdClassname) + " ioFlatIndexCounter[0] " + str(ioFlatIndexCounter[0]) + "\n")
  myDebugPrint3("numpoints: " + str(inCsdata.GetNumberOfPoints()) + " numcells: " + str(inCsdata.GetNumberOfCells()) + "\n")
  vtkInfo1 = inCsdata.GetInformation()
  myDebugPrint3("vtkInfo1 (inCsdata.Getinformation()): " + str(vtkInfo1) + "\n")


  ioFlatIndexCounter[0] += 1

  if icsdClassname == "vtkMultiBlockDataSet":
    thisIsLeaf = False
    nonLeafType = 0
  elif icsdClassname == "vtkExodusIIMultiBlockDataSet":
    thisIsLeaf = False
    nonLeafType = 1
  elif icsdClassname == "vtkMultiPieceDataSet":
    thisIsLeaf = False
    nonLeafType = 2
  elif icsdClassname == "vtkPartitionedDataSet":
    thisIsLeaf = False
    nonLeafType = 3
  else:
    thisIsLeaf = True
    nonLeafType = -1

  myDebugPrint3("thisIsLeaf " + str(thisIsLeaf) + " nonLeafType " + str(nonLeafType) + "\n")
  if not thisIsLeaf:
    if nonLeafType == 2:
      numBlocks = inCsdata.GetNumberOfPieces()
    elif nonLeafType == 3:
      numBlocks = inCsdata.GetNumberOfPartitions()
    else:
      numBlocks = inCsdata.GetNumberOfBlocks()

    myDebugPrint3("recursing, flat index: " + \
        str(ioFlatIndexCounter[0]) + \
        "    num blocks: " + \
        str(numBlocks) + "\n")
    for ii in range(0, numBlocks):
      myDebugPrint3("doing sub block " + str(ii) + " of " + str(numBlocks) + "\n")
      if nonLeafType == 2:
        oneBlock = inCsdata.GetPiece(ii)
      elif nonLeafType == 3:
        oneBlock = inCsdata.GetPartition(ii)
      else:
        oneBlock = inCsdata.GetBlock(ii)
      oneBlockMetaData = inCsdata.GetMetaData(ii)
      myDebugPrint3("oneBlockMetaData: " + str(oneBlockMetaData) + "\n")
      #myDebugPrint3("name: " + \
      #    oneBlockMetaData.Get(vtk.vtkCompositeDataSet.NAME()) + "\n")
      if oneBlock != None:
        PrintMultiBlockStatsRecurse1(oneBlock, oneBlockMetaData, recursionLevel + 1, ioFlatIndexCounter)
      else:
        myDebugPrint3("this sub block is None, index " + str(ioFlatIndexCounter[0]) + "\n")
        ioFlatIndexCounter[0] += 1
  else:
    PrintMultiBlockStatsOneBlock(inCsdata, inMetaData, ioFlatIndexCounter)

  myDebugPrint3("PrintMultiBlockStatsRecurse1 returning recursionLevel " + str(recursionLevel) + "\n")

def RecursivelyPrintMultiBlockStats(inFilter):
  if PhactoriDbg(100) == False:
    return

  myDebugPrint3("RecursivelyPrintMultiBlockStats entered\n")
  csdata = inFilter.GetClientSideObject().GetOutputDataObject(0)

  vtkInfo1 = csdata.GetInformation()
  myDebugPrint3("top level vtkInfo1 (inCsdata.Getinformation()): " + str(vtkInfo1) + "\n")
  #csdata.CopyInformationToPipeline(vtkInfo1)
  #csdata.CopyInformationFromPipeline(vtkInfo1)

  flatIndexCounter = [0]
  recursionLevel = 0
  PrintMultiBlockStatsRecurse1(csdata, None, recursionLevel, flatIndexCounter)

  myDebugPrint3("RecursivelyPrintMultiBlockStats returning\n")

