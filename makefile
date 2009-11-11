FC=ifort
FFLAGS= -c -132 -O0 -g -traceback -CB
LIBS=-L/opt/intel/mkl/10.0.2.018/lib/em64t -lmkl_lapack -lmkl 
LDFLAGS= -parallel 


OBJ= m_eulerangle.o  m_search_lma.o  m_wancenterdef.o cross.o mate.o rotdlmm.o linpackdummy.o writemesh_l.o 
main1: main1.o $(OBJ)
	$(FC) -o $@  main1.o $(OBJ) $(LIBS) $(LDFLAGS) 


main1.o : $(OBJ) 

clean:
	rm -f *.o 
