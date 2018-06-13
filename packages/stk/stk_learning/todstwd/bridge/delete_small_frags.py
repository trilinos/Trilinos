#!/usr/local/epd/1.7.4/appdata/canopy-1.7.4.3348.rh5-x86_64/User/bin/python
import sys
import os

var = 'SIERRA_SNTOOLS_PATH'
try:
  snToolsPath = os.environ[var];
except:
  raise RuntimeError('Please set environment variable: ' + var)
##snToolsPath = "/scratch/sierra/toolset"

toolset_dir = snToolsPath + "/contrib/testTools/adagio/"
module_dir = toolset_dir + "modules"
sys.path.append(module_dir)

from exodus import exodus
import geometry

def delete_small_frags():
    demo_unit = os.getcwd() + '/../../../../' + 'bin/linux-gcc-4.9.3-ip-openmpi-1.10.2/release/demo_disconnected_components'
    if not os.path.isfile(demo_unit):
      print (demo_unit)
      print( 'Did you build disconnected component finder?')
      sys.exit(1)
    os.system("module purge;module load sierra-devel;" + demo_unit + " -mesh_basename my_bridge -universal 1")
    e = exodus('my_bridge.e','r')
    coordinates = e.get_coords()
    elem_frag_ids = []
    blkIds = e.get_elem_blk_ids()
    step = e.num_times()
    fragId_numElems = dict()
    for blkId in blkIds:
      (elem_block_connectivity,num_elem_this_blk,num_nodes_per_elem) = e.get_elem_connectivity(blkId)
      eVals = e.get_element_variable_values(blkId,"elem_fragment_id",step)
      node_index = 0
      for iE in range(num_elem_this_blk):
        nodeCoords = []
        for iN in range(num_nodes_per_elem):
          node_id = elem_block_connectivity[node_index] - 1 ## index from 0
          node_index += 1
          nodeCoords.append( (coordinates[0][node_id],coordinates[1][node_id],coordinates[2][node_id]) )
        elem_center = geometry.averageCoord(nodeCoords)
        distL = geometry.SQdistance((-100,0,0),elem_center)
        distR = geometry.SQdistance(( 100,0,0),elem_center)
        fid = eVals[iE]
        if fid in fragId_numElems:
          fragId_numElems[fid][0] = fragId_numElems[fid][0] + 1
          fragId_numElems[fid][1] = geometry.plus(fragId_numElems[fid][1],elem_center)
          fragId_numElems[fid][2] = min(fragId_numElems[fid][2], distL)
          fragId_numElems[fid][3] = max(fragId_numElems[fid][3], distR)
        else:
          fragId_numElems[fid] = [ 1 , elem_center, distL, distR ]
      elem_frag_ids.extend(eVals)
    fragIdsToDelete = []
    car1 = (-24.056451612903224, 24.274193548387096, 0.0)
    car2 = (-58.01239669421488, 24.305785123966942, 0.0)
    car3 = (-89.01239669421489, 24.305785123966942, 0.0)
    print(fragId_numElems)
    for fid in fragId_numElems.keys():
      fragCenter = geometry.scale(1.0/fragId_numElems[fid][0],fragId_numElems[fid][1])
      dist1 = geometry.SQdistance(car1,fragCenter)
      dist2 = geometry.SQdistance(car2,fragCenter)
      dist3 = geometry.SQdistance(car3,fragCenter)
      if dist1 > 0.1 and dist2 > 0.1 and dist3 > 0.1 and fragId_numElems[fid][2] > 1.0 and fragId_numElems[fid][3] > 1.0:
          fragIdsToDelete.append(fid)
    nFrags = len(fragIdsToDelete)
    print('Delete ',nFrags,' frags')
    if nFrags > 0:
      elemsToDelete = []
      for i in range(len(elem_frag_ids)):
        if elem_frag_ids[i] in fragIdsToDelete:
            elemsToDelete.append(str(i+1))
      cb = open('stripElems.jou','w')
      cb.write('import mesh "my_bridge.e" no_geom\n')
      cb.write('delete hex ' + ' '.join(elemsToDelete) + '\n')
      cb.write('export mesh "my_bridge_clean.g" overwrite\n')
      cb.write('exit\n')
      cb.close()
      os.system('cubit -batch -nographics stripElems.jou')
      print (elemsToDelete)
    else:
      os.system('cp my_bridge.e my_bridge_clean.g')
    e.close()

if __name__ == "__main__":
  delete_small_frags()