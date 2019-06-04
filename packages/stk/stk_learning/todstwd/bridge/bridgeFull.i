begin sierra TOSDTWD Bridge

  define direction down with vector 0.0 -1.0 0.0

  begin definition for function DUMMY
    type is piecewise linear
    begin values
      0.0 0.0
      1.0 1.0
    end values
  end definition for function DUMMY

  begin definition for function pointsFunc
    type = analytic
    expression variable: x = nodal coordinates(x)
    expression variable: y = nodal coordinates(y)
    evaluate expression = " x > 100.0 ? ( y > 18.5 ? 1.0 : 0.0 ) : 0.0"
  end

  # { gc = 32.0 }
  # { mock_steel_mod = 4.176e9 } 
  # { mock_steel_dens = 490 }
  # { mock_steel_peak_stress = 2e7 }
  # { steel_dens = 8.26 }
  # { aluminum_dens = 2.69 }
  # { wood_dens = 0.56 }
  # { steel_stiff = 235.0 }
  # { aluminum_stiff = 68.0 }
  # { wood_stiff = 8.0 }

  begin material steel
    density = {mock_steel_dens}
    begin parameters for model elastic
      youngs modulus = {mock_steel_mod}
      poissons ratio = 0.26
    end
  end

  begin material aluminum
    density = { mock_steel_dens * aluminum_dens / steel_dens }
    begin parameters for model elastic
      youngs modulus = { mock_steel_mod * aluminum_stiff / steel_stiff }
      poissons ratio = 0.26
    end parameters for model elastic
  end

  begin material wood
    density = { mock_steel_dens * wood_dens / steel_dens }
    begin parameters for model elastic
      youngs modulus = { mock_steel_mod * wood_stiff / steel_stiff }
      poissons ratio = 0.26
    end parameters for model elastic
  end

  begin material car_mat
    density = { ( mock_steel_dens * aluminum_dens / steel_dens ) * 0.1 }
    begin parameters for model elastic
      youngs modulus = { (mock_steel_mod * aluminum_stiff / steel_stiff) * 0.1 }
      poissons ratio = 0.26
    end parameters for model elastic
  end

  begin finite element model bridgeFalls
    Database name = my_bridge_clean.g
    Database type = exodusII

    begin parameters for block block_1 # black pixels
      material = steel
      model = elastic
    end parameters for block block_1
    begin parameters for block block_2 # red pixels
      material = steel
      model = elastic
    end parameters for block block_2 
    begin parameters for block block_3 # green pixels
      material = aluminum
      model = elastic
    end parameters for block block_3
    begin parameters for block block_4 # blue pixels
      material = car_mat
      model = elastic
    end parameters for block block_4
    begin parameters for block block_5 # grey pixels
      material = steel
      model = elastic
    end parameters for block block_5
    begin parameters for block block_6 # brown pixels
      material = wood
      model = elastic
    end parameters for block block_6
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
        gravitational constant = {gc}
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
        database name = bridgeLoads.e
        database type = exodusII
        at step 0 interval = 1
        nodal variables = displacement as displ
        element variables = von_mises as stress
      end results output output_adagio

      begin contact definition
        skin all blocks = on
        begin interaction defaults
          general contact = on
          self contact = on
          constraint formulation = node_face
        end
      end

      begin solver
        begin cg
          maximum iterations = 500
          begin full tangent preconditioner
            tangent diagonal scale = 1e-6
          end
          target relative residual = 1e-8
          acceptable relative residual = 1.0
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

      begin node based time step parameters
        step interval = 100
      end

      begin gravity
        include all blocks
        direction = down
        gravitational constant = {gc}
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
        magnitude = {50*1.50} 0 0 
      end

      begin results output output_adagio
        database name = bridgeFalls.e
        database type = exodusII
        at time 0.0 interval = 1.0e-2
        nodal variables = displacement as displ
        nodal variables = points
        global variables = total_points
        element variables = von_mises as stress
        element variables = death_status as death
      end results output output_adagio

      begin contact definition
        skin all blocks = on
        begin interaction defaults
          general contact = on
          self contact = on
          constraint formulation = node_face
        end
        ##begin enforcement options
        ##  momentum balance iterations = 1
        ##end
      end

      ##begin mass scaling
      ##  include all blocks
      ##  target time step = 3.0e-4
      ##end

      # { check_death = 10 }
      begin element death steel
        block = block_1 block_5
        criterion is element value of von_mises > {mock_steel_peak_stress}
        check step interval = {check_death}
      end
      begin element death aluminum
        block = block_3
        criterion is element value of von_mises > {mock_steel_peak_stress * aluminum_stiff/steel_stiff}
        check step interval = {check_death}
      end
      begin element death wood
        block = block_6
        criterion is element value of von_mises > {mock_steel_peak_stress * wood_stiff/steel_stiff}
        check step interval = {check_death}
      end
      begin element death points
        block = block_4
        criterion is max nodal value of points >=  0.5
      end
      begin element death bounds
        include all blocks
        check step interval = {check_death}
        criterion is avg nodal value of coordinates(x) >=  101
        criterion is avg nodal value of coordinates(x) <= -101
        #criterion is avg nodal value of coordinates(y) >=  40
        criterion is avg nodal value of coordinates(y) <= -40
      end
      begin element death inversion
        include all blocks
        death on inversion = on
      end

      begin user output
        block = block_4
        compute nodal points as function pointsFunc
        compute at every step
      end

    end
  end
end sierra
