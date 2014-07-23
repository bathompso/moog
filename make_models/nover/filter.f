	program filter
c
c a routine to filter out unwanted stuff from the Kurucz model files
c
	character*80 line
	do 100 i = 1,100000
	   read(5,'(a80)',end=999)line
	   if (line(:4).eq.'TEFF') then
	      read(line(5:12),*)teff
              if (teff.gt.8750.) then
                 stop
	      else
                 write(6,'(a80)')line
              endif
           elseif (line(:4).eq.'READ') then
              write(6,'(a80)')line
	      do 200 j = 1,72
	         read(5,'(a80)',end=999)line
                 write(6,'(a44)')line(:44)
 200	      continue
	   endif
 100	continue
 999	continue
	stop
	end
