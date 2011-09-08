      subroutine fixmap(numel, map)
      integer map(numel)
      
      do 10 i=1, numel
        if (map(i) .eq. 0) map(i) = i
 10   continue
      return
      end
