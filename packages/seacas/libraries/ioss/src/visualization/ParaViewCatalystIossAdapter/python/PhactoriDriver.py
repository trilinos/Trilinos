# Copyright(C) 1999-2020 National Technology & Engineering Solutions
# of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
# NTESS, the U.S. Government retains certain rights in this software.
#
# See packages/seacas/LICENSE for details

try: paraview.simple
except: from paraview.simple import *

from paraview import coprocessing

# ----------------------- Pipeline definition -----------------------

import phactori

#do sierra/catalyst logging (from process id 0 only)
#note: logging should be disabled for testing, but enabled for production
if phactori.SmartGetLocalProcessId() == 0:
  import os
  from datetime import datetime

  loggingIsEnabled = True
  if "SNL_CATALYST_SIERRA_USAGE_LOG_FLAG" in os.environ:
    #print "SNL_CATALYST_SIERRA_USAGE_LOG_FLAG environment variable: " + \
            #os.environ["SNL_CATALYST_SIERRA_USAGE_LOG_FLAG"]
    if os.environ["SNL_CATALYST_SIERRA_USAGE_LOG_FLAG"] == "disable":
      loggingIsEnabled = False

  if loggingIsEnabled:
    #print "I am process 0 doing logging!"
    if "HOSTNAME" in os.environ:
      os.environ["LOG_PLATFORM"] = os.environ["HOSTNAME"]
    else:
      os.environ["LOG_PLATFORM"] = "HOSTNAMEnotset"
    os.environ["LOG_PRODUCT"] = "sierra-catalyst"
    os.environ["CATALYST_VERSION"] = "p4.1.0-s4.31.6-p"
    nowdate = datetime.now()
    #os.environ["DATE"] = str(nowdate.month) + "/" + str(nowdate.day) + "/" + str(nowdate.year)
    #currently, DATE environment variable is overwritten in logParserScr script
    #print "date is :  " + os.environ["DATE"]
    import subprocess
    try:
      subprocess.call("/projects/viz/catalyst/utilities/logParserScr")
    except:
      print "usage logging warning: logParserScr call failed\n"

    #print "platform: " + os.environ["LOG_PLATFORM"]
    #print "product:  " + os.environ["LOG_PRODUCT"]
    #print "version:  " + os.environ["CATALYST_VERSION"]
    #print "date:     ->" + os.environ["DATE"] + "<-"
    #print "I am process 0 done logging!"
  #else:
    #print "I am process 0, logging disabled in environment variable"

#def UseViewMapB(inViewMapB):
#  PhactoriScript.UseViewMapB_ps(inViewMapB)

#def UseDataSetupMapB(inDataSetupMapB):
#  PhactoriScript.UseDataSetupMapB_ps(inDataSetupMapB)

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      dummyItem = "blah"
      #global cp_views, cp_writers

      #print "CreatePipeline startup"
      #print "cp_views len is " + str(len(cp_views))

      #rigid_body_impact_6_ff_e = coprocessor.CreateProducer( datadescription, "input" )

      #phactori.SetUpCoProcessor(coprocessor)
      #phactori.SetCpViewsAndCpWriters(cp_views, cp_writers)

#QQQQQ->we want to do this and CreateProducer once for each output results block,
#reuse same coprocesor instance (?), then have LocalWriteImages3 only do the images for a given output results bloxk:q

      #PhactoriScript.CreatePipeline(datadescription)

    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

    def LocalWriteImages3(self, datadescription, rescale_lookuptable=False):
      phactori.WriteImagesForCurrentPipeAndViewsState(datadescription)

    def LocalExportOperationsData3(self, datadescription, rescale_lookuptable=False):
      phactori.ExportOperationsDataForCurrentPipeAndViewsState(datadescription)



  coprocessor = CoProcessor()
  freqs = {'input': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor


#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
#coprocessor = CreateCoProcessor()
coprocessor = None

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
#coprocessor.EnableLiveVisualization(False)

global gDoNewScatterPlotsC
gDoNewScatterPlotsC = True
#gDoNewScatterPlotsC = False

def PerFrameUpdate(datadescription):
  phactori.PerRendersetInitialization(datadescription)
  phactori.UpdateAllOperationsWhichMayChangeWithData()
  phactori.UpdateAllImagesetViewsWhichMayChangeWithData()

  #we are now updating color range stuff right before WriteImage
  #phactori.UpdateDataRangesForColorValues()

  phactori.myDebugPrint3('PerFrameUpdate entered\n')
  #GetAndReduceViewControl()
  global tearListPersistent
  global deathListPersistent

  global gDoNewScatterPlotsC
  if gDoNewScatterPlotsC:
    phactori.myDebugPrint3('doing phactori.UpdateAllScatterPlots (C)\n')
    phactori.UpdateAllScatterPlots()
    phactori.myDebugPrint3('doing phactori.UpdateAllPlotsOverTime\n')
    phactori.UpdateAllPlotsOverTime()
    phactori.myDebugPrint3('did plot updates\n')

  #currentFrameTearList = phactori.CollectCells1('TEAR_DOUBLE', 0.01, 1)
  #tearListPersistent = mergeCurrentIntoPersistent(tearListPersistent,
  #                       currentFrameTearList)
  #currentFrameDeathList = phactori.CollectCells1('STATUS', 0.8, -1)
  #deathListPersistent = mergeCurrentIntoPersistent(deathListPersistent,
  #                        currentFrameDeathList)
  #compareTearDeath(tearListPersistent, deathListPersistent)

  phactori.myDebugPrint3('PerFrameUpdate exiting\n')


#begin tear/death persistence; not used now but may be useful later
global tearListPersistent
tearListPersistent = []
global deathListPersistent
deathListPersistent = []

def mergeCipCompare(item1, item2):
  if item1[1] > item2[1]:
   return 1
  if item1[1] < item2[1]:
   return -1
  if item1[2] > item2[2]:
   return 1
  if item1[2] < item2[2]:
   return -1
  return 0

def mergeCurrentIntoPersistent(pList, cList):
  phactori.myDebugPrint2('mergeCurrentIntoPersistent entered\n')
  pListLen = len(pList)
  cListLen = len(cList)
  phactori.myDebugPrint2('  pList has ' + str(pListLen) + ' elements, cList has ' + str(cListLen) + ' elements\n')
  pIndex = 0
  cIndex = 0
  mergeList = []
  while (pIndex < pListLen) or (cIndex < cListLen):
    if(pIndex >= pListLen):
      mergeList.append(cList[cIndex])
      #phactori.myDebugPrint2('  pIndex ' + str(pIndex) + ' cIndex ' + str(cIndex) + '\n')
      #phactori.myDebugPrint2('  (p ended) [cIndex] ' + str(cList[cIndex]) + '\n')
      cIndex += 1
      continue
    if(cIndex >= cListLen):
      mergeList.append(pList[pIndex])
      #phactori.myDebugPrint2('  pIndex ' + str(pIndex) + ' cIndex ' + str(cIndex) + '\n')
      #phactori.myDebugPrint2('  (c ended) [pIndex] ' + str(cList[pIndex]) + '\n')
      pIndex += 1
      continue
    compareResult = mergeCipCompare(cList[cIndex], pList[pIndex])
    #phactori.myDebugPrint2('  pIndex ' + str(pIndex) + ' cIndex ' + str(cIndex) + '\n')
    #phactori.myDebugPrint2('  [pIndex] ' + str(pList[pIndex]) + ' [cIndex] ' + str(cList[cIndex]) + '\n')
    if compareResult == 1:
      #phactori.myDebugPrint2('  item from pList is less\n')
      mergeList.append(pList[pIndex])
      pIndex += 1
    if compareResult == -1:
      #phactori.myDebugPrint2('  item from cList is less\n')
      mergeList.append(cList[cIndex])
      cIndex += 1
    if compareResult == 0:
      #phactori.myDebugPrint2('  items from cList and pList are same element\n')
      mergeList.append(pList[pIndex])
      pIndex += 1
      cIndex += 1
  phactori.myDebugPrint2('  merged list has ' + str(len(mergeList)) + ' elements\n')
  phactori.myDebugPrint2('mergeCurrentIntoPersistent returning\n')
  return mergeList

def compareTearDeath(tList, dList):
  phactori.myDebugPrint2('compareTearDeath entered\n')
  phactori.myDebugPrint2('compareTearDeath returning\n')
#end tear/death persistence; not used now but may be useful later


# ---------------------- Data Selection method ----------------------

global gFirstTimeInDoCoProcessing
gFirstTimeInDoCoProcessing = True
global gSkipCountdown
#gSkipCountdown = 3
gSkipCountdown = 0

global gCatchAllExceptionsAndPassUpFlag
gCatchAllExceptionsAndPassUpFlag = True
#gCatchAllExceptionsAndPassUpFlag = False

def RequestDataDescription(datadescription):
  phactori.myDebugPrint3("PhactoriDriver.RequestDataDescription entered: " + str(gDoCoProcessingCount)+ "\n");

  phactori.TestUserDataForBypassScript(datadescription)

  if phactori.GetBypassUserDataFlag() == False:
    fd = datadescription.GetUserData()

    if fd == None:
      phactori.myDebugPrint2("no user data, returning {}\n")
      returnViewMapC = {}
      return returnViewMapC

  global gCatchAllExceptionsAndPassUpFlag
  if gCatchAllExceptionsAndPassUpFlag:
    try:
      return RequestDataDescriptionSub(datadescription)
    except:
      import traceback
      tb = traceback.format_exc()
      phactori.IssueErrorOrWarningThroughSierraIO(datadescription, tb, True)
  else:
    return RequestDataDescriptionSub(datadescription)

def RequestDataDescriptionSub(datadescription):
    phactori.myDebugPrint("PhactoriDriver.RequestDataDescriptionSub entered\n");
    "Callback to populate the request for current timestep"

    phactori.TestUserDataForBypassScript(datadescription)

    if phactori.GetBypassUserDataFlag() == False:
      fd = datadescription.GetUserData()

      if fd == None:
        phactori.myDebugPrint2("no user data, returning {}\n")
        returnViewMapC = {}
        return returnViewMapC

    global coprocessor

    global gFirstTimeInDoCoProcessing
    global gSkipCountdown

    if gFirstTimeInDoCoProcessing == True:
      phactori.myDebugPrint2("RequestDataDescription doing gFirstTimeInDoCoProcessing\n")
      phactori.myDebugPrint2(" skip countdown is " + str(gSkipCountdown) + "\n")
      if gSkipCountdown > 0:
        gSkipCountdown = gSkipCountdown - 1
        return 0
      coprocessor = CreateCoProcessor()
      coprocessor.EnableLiveVisualization(False)
      gFirstTimeInDoCoProcessing = False

    #import pdb
    #pdb.set_trace()

    phactori.InitializePerPipeRoot(datadescription, coprocessor)

    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return 1

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)
    return 1




# ------------------------ Processing method ------------------------

global tripleBufferCount
tripleBufferCount = 0

global gDoCoProcessingCount
gDoCoProcessingCount = 0

def DoCoProcessing(datadescription):
  phactori.myDebugPrint3("PhactoriDriver.DoCoProcessing entered: " + str(gDoCoProcessingCount)+ "\n");

  fd = datadescription.GetUserData()

  if phactori.GetBypassUserDataFlag() == False:
    if fd == None:
      phactori.myDebugPrint2("no user data, returning {}\n")
      returnViewMapC = {}
      return returnViewMapC

  global gCatchAllExceptionsAndPassUpFlag
  if gCatchAllExceptionsAndPassUpFlag:
    try:
      DoCoProcessingSub(datadescription)
    except:
      import traceback
      tb = traceback.format_exc()
      phactori.IssueErrorOrWarningThroughSierraIO(datadescription, tb, True)
  else:
    DoCoProcessingSub(datadescription)

def DoCoProcessingSub(datadescription):
    "Callback to do co-processing for current timestep"

    global gDoCoProcessingCount
    phactori.myDebugPrint3("PhactoriDriver.DoCoProcessingSub entered: " + str(gDoCoProcessingCount)+ "\n");

    gDoCoProcessingCount += 1


    fd = datadescription.GetUserData()

    if phactori.GetBypassUserDataFlag() == False:
      if fd == None:
        phactori.myDebugPrint2("no user data, returning {}\n")
        returnViewMapC = {}
        return returnViewMapC

    global coprocessor
    global gFirstTimeInDoCoProcessing
    global gSkipCountdown

    if gFirstTimeInDoCoProcessing == True:
      phactori.myDebugPrint2("DoCoProcessing doing gFirstTimeInDoCoProcessing\n")
      phactori.myDebugPrint2(" skip countdown is " + str(gSkipCountdown) + "\n")
      if gSkipCountdown > 0:
        return
      coprocessor = CreateCoProcessor()
      coprocessor.EnableLiveVisualization(False)
      gFirstTimeInDoCoProcessing = False

    #import pdb
    #pdb.set_trace()

    phactori.InitializePerPipeRoot(datadescription, coprocessor)

    "Callback to do co-processing for current timestep"
    timestep = datadescription.GetTimeStep()

    phactori.myDebugPrint("timestep is: " + str(timestep) + "\n");

    phactori.SmartGetLocalProcessId()

    # Load the Pipeline if not created yet
    #if not pipeline:
    #   phactori.myDebugPrint("PhactoriDriver.DoCoProcessing creating pipeline\n");
    #   pipeline = CreatePipeline(datadescription)
    #else:
    #   phactori.myDebugPrint("PhactoriDriver.DoCoProcessing updating pipeline\n");
    #   # update to the new input and time
    #   UpdateProducers(datadescription)
    #   PerFrameUpdate(datadescription)


    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    PerFrameUpdate(datadescription)

    # check for simulation-data-based i/o filtering--skip image creation
    # and writing if criteria has been set up to determine whether to
    # create images, such as 'maximum of variable X above 80.0'
    result = phactori.WriteOutImagesTest(datadescription, coprocessor)
    if result == False:
      #don't write images
      return

    # Write output data, if appropriate.
    #coprocessor.WriteData(datadescription);

    # Write output data
    #WriteAllData(datadescription, cp_writers, timestep);

    # Write image capture (Last arg: rescale lookup table)
    #phactori.myDebugPrint("PhactoriDriver.DoCoProcessing writing images\n");
    #LocalWriteAllImages(datadescription, cp_views, timestep, False)
    #WriteAllImages(datadescription, cp_views, timestep, False)

    # Live Visualization
    #if (len(cp_views) == 0) and live_visu_active:
    #   DoLiveInsitu(timestep, pv_host, pv_port)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription)

    coprocessor.LocalExportOperationsData3(datadescription)

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.LocalWriteImages3(datadescription,
        rescale_lookuptable=False)

    #test and allow for looping when doing user vis interaction while
    #pausing simulation
    continueWriteAndInteractionCheckLoop = True
    while continueWriteAndInteractionCheckLoop:
      interactionTestResult = \
          phactori.DoUserInteractionWithSimulationPausedIfEnabled()
      if interactionTestResult == 0:
        #looping interaction is not or is no longer on; allow simulation to
        #continue
        continueWriteAndInteractionCheckLoop = False
      elif interactionTestResult == 1:
        #looping interaction is on, but there were no changes to the vis
        #(i.e. no trigger was given to update vis).  Therefore do not write
        #images, but continue looping and waiting for vis change trigger
        continueWriteAndInteractionCheckLoop = True
      elif interactionTestResult == 2:
        #there was a vis change triggered; update the images for the new
        #vis, write out the images, and continue looping
        continueWriteAndInteractionCheckLoop = True
        imagesNeedWriting = True
        if imagesNeedWriting:
          phactori.UpdateAllImagesetViewsWhichMayChangeWithData()
          coprocessor.LocalExportOperationsData3(datadescription)
          coprocessor.LocalWriteImages3(datadescription,
          rescale_lookuptable=False)

    #coprocessor.WriteImages(datadescription, rescale_lookuptable=True)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)

