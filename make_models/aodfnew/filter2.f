	program filter
c
c a routine to filter out unwanted stuff from the Kurucz model files
c
	character*80 line
	real logg
	do 100 i = 1,100000
	   read(5,'(a80)',end=999)line
	   if (line(:4).eq.'TEFF') then
	      read(line(5:12),*)teff
	      read(line(23:29),*)logg
              write(6,'(5HMODEL)')
              write(6,'(f10.0,5x,f6.2)')teff,logg
	      read(5,'(a80)',end=999)line
           else
              write(6,'(a44)')line(:44)
	   endif
 100	continue
 999	continue
	stop
	end
