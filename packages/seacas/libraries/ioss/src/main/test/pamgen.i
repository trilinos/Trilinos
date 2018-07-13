  mesh
   radial trisection
    trisection blocks, 4
    transition radius, 0.02
    numz 1
      zblock 1 1.0 interval 20
    numr 3
      rblock 1 0.04 first size 0.004  last size 0.006
      rblock 2 0.06 first size 0.006  last size 0.0075
      rblock 3 0.04 first size 0.0075 last size 0.01
    numa 1
      ablock 1 360.0  interval 24
  end
  set assign
     block sideset, ihi, 1, 3
     block sideset, klo, 10, 1
     block sideset, khi, 11, 1
     block sideset, klo, 20, 2
     block sideset, khi, 21, 2
     block sideset, klo, 30, 3
     block sideset, khi, 31, 3
  end
end
