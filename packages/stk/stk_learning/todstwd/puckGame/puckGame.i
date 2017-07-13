begin sierra TOSDTWD Puck Game

  define direction down with vector 0.0 -1.0 0.0

  begin definition for function DUMMY
    type is constant
    begin values
      1.0
    end values
  end definition for function DUMMY

  begin material hard
    density = 1.0
    begin parameters for model elastic
      youngs modulus = 1000000.
      poissons ratio = 0.3
    end parameters for model elastic
  end

  begin rigid body puck
  end

  begin solid section puck
      rigid body = puck
  end

  begin finite element model rigidBlocks
    Database name = pic.g
    Database type = exodusII

    begin parameters for block block_1
      material = hard
      model = elastic
    end parameters for block block_1
    begin parameters for block block_2
      material = hard
      model = elastic
    end parameters for block block_2
    begin parameters for block block_3
      material = hard
      model = elastic
    end parameters for block block_3
    begin parameters for block block_4
      material = hard
      model = elastic
      section = puck
    end parameters for block block_4
  end

  begin solid mechanics procedure adagio_procedure_1

    begin time control
      begin time stepping block p0
        start time = 0.0
        begin parameters for presto region adagio_1
          user time step = 1e-4
        end parameters for presto region adagio_1
      end
      termination time = 5.0
    end time control

    begin solid mechanics region adagio_1

    use finite element model rigidBlocks

    begin gravity
      block = block_4
      direction = down
      function = dummy
      gravitational constant = {9.8 * 100}
    end

    begin fixed displacement
      include all blocks
      components = z
    end fixed displacement

    begin fixed displacement
      block = block_1 block_2 block_3
      components = x y 
    end fixed displacement

    begin fixed rotation
      rigid body = puck
      components = x y 
    end

    begin results output output_adagio
      database name = game_results.exo
      database type = exodusII
      at step 0 interval = 50
      nodal variables = displacement as displ
    end results output output_adagio

    begin contact definition
      skin all blocks = on
      begin interaction defaults
        general contact = on
        self contact = off
      end
    end
    
    begin user output
      block = block_4
      compute at every step
      compute global min_x as min of nodal coordinates(x)
      compute global max_x as max of nodal coordinates(x)
      compute global min_y as min of nodal coordinates(y)
      compute global max_y as max of nodal coordinates(y)
    end

    begin user output
      block = block_2
      compute at every step
      compute global laserDeath as max of nodal force_contact
    end

    begin solution termination
      terminate global kinetic_Energy < 1e-5
      skip times = 0.0 to 0.5
      terminate global min_x < -50
      terminate global max_x > 250
      terminate global min_y < -50
      terminate global max_y > 250
      terminate global laserDeath > 0
    end
  end
end
end sierra
