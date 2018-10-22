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
end sierra
