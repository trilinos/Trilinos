      integer         m_max   , n_max
            parameter(m_max=512, n_max=2)
      integer         p_max
            parameter(p_max = m_max)
      integer         q_max
            parameter(q_max = 2)
      integer         buff
            parameter(buff = 36)
      integer         NSTEP_max
            parameter(NSTEP_max = 36)

      integer         m_max_act, n_max_act
            parameter(m_max_act  = (m_max/3)*2 + 2)
            parameter(n_max_act  = (n_max/3)*2 + 2)

      integer         m2max
            parameter(m2max = 2*m_max)

