begin sierra TOSDTWD Bridge

  define direction down with vector 0.0 -1.0 0.0

  begin definition for function DUMMY
    type is piecewise linear
    begin values
      0.0 0.0
      1.0 1.0
    end values
  end definition for function DUMMY

  begin material hard
    density = 1.0
    begin parameters for model elastic
      youngs modulus = 1000000.
      poissons ratio = 0.3
    end parameters for model elastic
  end

  begin material soft
    density = {1.0 / 10}
    begin parameters for model elastic
      youngs modulus = {1000000. / 10}
      poissons ratio = 0.3
    end parameters for model elastic
  end

  begin finite element model bridgeFalls
    Database name = pic.g
    Database type = exodusII

    begin parameters for block block_1
      material = hard
      model = elastic
    end parameters for block block_1
    begin parameters for block block_2 # red pixels
      material = hard
      model = elastic
    end parameters for block block_2 
    begin parameters for block block_3 # green pixels
      material = soft
      model = elastic
    end parameters for block block_3
    begin parameters for block block_4 # blue pixels
      material = hard
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
        gravitational constant = {9.8 * .5}
      end
  
      begin fixed displacement
        include all blocks
        components = z
      end fixed displacement
  
      begin fixed displacement
        block = block_2 
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
      termination time = 2.0
    end time control

    begin solid mechanics region presto1
      use finite element model bridgeFalls
  
      begin gravity
        include all blocks
        direction = down
        gravitational constant = {9.8 * .5}
      end
  
      begin fixed displacement
        include all blocks
        components = z
      end fixed displacement
  
      begin fixed displacement
        block = block_2 
        components = x y 
      end fixed displacement
  
      begin results output output_adagio
        database name = bridgeFalls.exo
        database type = exodusII
        at time 0.0 interval = 1.0e-2
        nodal variables = displacement as displ
        element variables = von_mises as stress
        element variables = death_status as death
      end results output output_adagio
  
      begin contact definition
        skin all blocks = on
        begin interaction defaults
          general contact = on
          self contact = off
        end
      end
 
      begin element death death1
        block = block_1 block_3 block_4
        criterion is element value of von_mises > 3e3
      end
    end
  end
end sierra
