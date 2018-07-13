begin sierra TOSDTWD Bridge

  define direction down with vector 0.0 -1.0 0.0

  begin definition for function DUMMY
    type is piecewise linear
    begin values
      0.0 0.0
      1.0 1.0
    end values
  end definition for function DUMMY

  begin material a36_steel
      density = 490
      begin parameters for model elastic
          youngs modulus = 4176000000.0003
          poissons ratio = 0.26
      end
  end

  begin material car_material
    density = {490 / 1.0 }
    begin parameters for model elastic
      youngs modulus = { 4176000000.0003 / 50 }
      poissons ratio = 0.26
    end parameters for model elastic
  end

  begin material soft
    density = {490 / 10}
    begin parameters for model elastic
      youngs modulus = {4176000000.0003 / 10}
      poissons ratio = 0.26
    end parameters for model elastic
  end

  begin finite element model bridgeFalls
    Database name = pic.g
    Database type = exodusII

    begin parameters for block block_1
      material = a36_steel
      model = elastic
    end parameters for block block_1
    begin parameters for block block_2 # red pixels
      material = a36_steel
      model = elastic
    end parameters for block block_2 
    begin parameters for block block_3 # green pixels
      material = soft
      model = elastic
    end parameters for block block_3
    begin parameters for block block_4 # blue pixels
      material = car_material
      model = elastic
    end parameters for block block_4
  end

  begin solid mechanics procedure adagio_procedure_1

    begin time control
      begin time stepping block p0
        start time = 0.0
        begin parameters for adagio region adagio_1
          number of time steps = 10
        end parameters for adagio region adagio_1
      end
      termination time = 1.0
    end time control

    begin solid mechanics region adagio_1
      use finite element model bridgeFalls
  
      begin gravity
        include all blocks
        direction = down
        function = dummy
        gravitational constant = {32 * 1.0}
      end
  
      begin fixed displacement
        include all blocks
        components = z
      end fixed displacement
  
      begin fixed displacement
        block = block_2 block_4
        components = x y 
      end fixed displacement
  

      begin results output output_adagio
        database name = bridgeLoads.exo
        database type = exodusII
        at step 0 interval = 1
        nodal variables = displacement as displ
        element variables = von_mises as stress
      end results output output_adagio
  
      begin contact definition
        skin all blocks = on
        begin interaction defaults
          general contact = on
          self contact = off
        end
      end
  
      begin solver
        begin cg
          maximum iterations = 500
          begin full tangent preconditioner
            tangent diagonal scale = 1e-6
          end
          target relative residual = 1e-8
        end
        begin control contact
          maximum iterations = 10
          acceptable relative residual = 1e-2
          target relative residual = 1e-4
        end
      end
    end
  end
  begin solid mechanics procedure presto_procedure

    begin procedural transfer migration1
      include all blocks
    end

    begin time control
      begin time stepping block p0
        start time = 0.0
        begin parameters for presto region presto1
              
        end 
      end
      termination time = 3.0
    end time control

    begin solid mechanics region presto1
      use finite element model bridgeFalls
  
      begin gravity
        include all blocks
        direction = down
        gravitational constant = {32 * 1}
      end
  
      begin fixed displacement
        include all blocks
        components = z
      end fixed displacement
  
      begin fixed displacement
        block = block_2 
        components = x y 
      end fixed displacement
  
      begin initial condition
        block = block_4
        initialize variable name = velocity
        variable type = node
        magnitude = 50 0 0 
      end

      begin results output output_adagio
        database name = bridgeFalls.exo
        database type = exodusII
        at time 0.0 interval = 1.0e-2
        nodal variables = displacement as displ
        element variables = von_mises as stress
        element variables = death_status as death
      end results output output_adagio
  
      begin contact definition
        contact surface block_1 contains block_1
        contact surface block_2 contains block_2
        contact surface block_3 contains block_3
        contact surface block_4 contains block_4

        update all surfaces for element death = on
#        begin interaction defaults
#          general contact = on
#          self contact = on
#        end

        begin interaction
            surfaces = block_1 block_2
        end
        begin interaction
            surfaces = block_1 block_3
        end
        begin interaction
            surfaces = block_1 block_4
        end
        begin interaction
            surfaces = block_2 block_3
        end
        begin interaction
            surfaces = block_2 block_4
        end
      end
 
      begin element death death1
        block = block_1 block_3 block_4
        criterion is element value of von_mises > 2e7
      end
    end
  end
end sierra
