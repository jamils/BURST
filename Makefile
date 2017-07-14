#CF90    = mpif90
CF90    = gfortran
#SCFLAGS = -O0 -cpp -fbounds-check # serial
SCFLAGS = -O0 -cpp # serial
#PCFLAGS = -O0 -cpp -Dprllel # parallel w/o gprof profiling data
PCFLAGS = -O0 -cpp -Dprllel # parallel w/o gprof profiling data
#PCFLAGS = -pg -O0 -Dprllel # parallel w/ gprof profiling data
#PCFLAGS = -O0 -cpp -Dprllel -fbounds-check # parallel w/o gprof profiling data
CFLAGS  = $(SCFLAGS) # set to $(PCFLAGS) or $(SCFLAGS)
CMD     = ./burst

.SUFFIXES: .f .f90 .F .f95

%.o : %.f
	$(CF90) $(CFLAGS) -c $<

%.o : %.f90
	$(CF90) $(CFLAGS) -c $<

%.o : %.F
	$(CF90) $(CFLAGS) -c $<

%.o : %.f95
	$(CF90) $(CFLAGS) -c $<

%.mod : %.f90
	$(CF90) $(CFLAGS) -c $<

#%.mod : %.o

include .depend

clean:
	/bin/rm -f $(CMD) core *.mod *.o .depend
depend .depend:
	makedepf90 -W -m %f.mod \
		-l"$(CF90) $(CFLAGS) -o $(CMD) $(FOBJ)" -o $(CMD) ./*.f90 > .depend

