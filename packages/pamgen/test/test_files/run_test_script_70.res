PARSE ERROR
    1:    mesh
    2:    brick
    3:     numz 1
    4:       zblock 1 10.0 interval 6
    5:     numx 3 initial radius 10.
    6:       xblock 1 2. interval 12
    7:       xblock 2 5. interval 6
    8:       xblock 3 5. interval 12
    9:     numy 2
   10:       yblock 1 18. interval 18
   11:       yblock 2 18. interval 18
   12:     end
   13:     set assign
   14:       block nodeset, jlo, 1 1
   15:     end
   16:     topology modification
   17:     suppress block, 1
   18:     suppress block, 2
   19:     suppress block, 6
   20:     suppress block, 600
   21:     end
   22:     decomposition strategy
   23:       sequential
   24:     end
   25: 
   26:   end
   27: 
   28: SETUP ERROR IN CHECK_BLOCKS Terminating from Inline_Mesh_Desc::Check_Blocks block 600 may not be suppressed as it does not exist.
