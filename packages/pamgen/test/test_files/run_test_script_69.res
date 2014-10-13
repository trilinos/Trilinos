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
   13: $    set assign
   14: $      block nodeset, jlo, 1 1
   15: $      block nodeset, jlo, 2 2
   16: $      block sideset, jlo, 1 1
   17: $      block sideset, jlo, 2 2
   18: $    end
   19:     topology modification
   20:     suppress block, 1
   21:     suppress block, 2
   22:     suppress block, 6
   23:     suppress block, 600
   24:     end
   25:     decomposition strategy
   26:       sequential
   27:     end
   28: 
   29:   end
   30: 
   31: SETUP ERROR IN CHECK_BLOCKS Terminating from Inline_Mesh_Desc::Check_Blocks block 600 may not be suppressed as it does not exist.
