FC=gfortran

gwsrc=/home/kino/kit/GW/ecal3ad7c59d8abb1fe0553f3f0f415ad37d7e36b674/fpgw/gwsrc/

FFLAGS= -c -132 -O0 -g -traceback -CB -module /home/kino/kit/GW/ecal3ad7c59d8abb1fe0553f3f0f415ad37d7e36b674/fpgw/exec 

FFLAGS= -ffixed-line-length-132 -c -g -O0 \
 -frange-check \
 -I/home/kino/kit/GW/ecal3ad7c59d8abb1fe0553f3f0f415ad37d7e36b674/fpgw/exec
# -frange-check \# -ffpe-trap=invalid,overflow,zero,underflow,precision,denormal \

LD=ifort
LIBS=-L/opt/intel/mkl/10.0.2.018/lib/em64t -lmkl_lapack -lmkl  -lgfortran
LDFLAGS=  -parallel 



OBJ0= m_eulerangle.o  m_search_lma.o  m_wancenterdef.o writemesh_l.o  m_lmlabel.o 

OBJ= $(OBJ0) cross.o mate.o rotdlmm.o linpackdummy.o

para_g=.o 

 GW0 = \
 $(gwsrc)readpomat.o \
 $(gwsrc)keyvalue.o \
 $(gwsrc)rppovl.o \
 $(gwsrc)nocctotg.o \
 $(gwsrc)ppbafp.fal$(para_g) \
 $(gwsrc)psi2b_v2$(para_g) \
 $(gwsrc)psi2b_v3$(para_g) \
 $(gwsrc)wfacx.o \
 $(gwsrc)sortea.o \
 $(gwsrc)rydberg.o \
 $(gwsrc)polinta.o \
 $(gwsrc)efsimplef.o \
 $(gwsrc)extension.o \
 $(gwsrc)rangedq.o \
 $(gwsrc)nword.o \
 $(gwsrc)scg.o \
 $(gwsrc)matm.o \
 $(gwsrc)rdpp.o \
 $(gwsrc)mptauof.o \
 $(gwsrc)genallcf_mod.o \
 $(gwsrc)rgwinf_mod.o \
 $(gwsrc)rotdlmm.o \
 $(gwsrc)iopen.o \
 $(gwsrc)cputid.o \
 $(gwsrc)rw.o \
 $(gwsrc)ext.o \
 $(gwsrc)ext2.o \
 $(gwsrc)cross.o \
 $(gwsrc)mate.o \
 $(gwsrc)mate1.o \
 $(gwsrc)icopy.o \
 $(gwsrc)bib1.o \
 $(gwsrc)index.o \
 $(gwsrc)idxk.o \
 $(gwsrc)maxnn.o \
 $(gwsrc)reindx.o \
 $(gwsrc)pointops.o \
 $(gwsrc)iolib.o \
 $(gwsrc)iprint.o \
 $(gwsrc)bz.o \
 $(gwsrc)bzmesh.o \
 $(gwsrc)genqbz.o \
 $(gwsrc)linpackdummy.o \
 $(gwsrc)switches.o \
 $(gwsrc)rwbzdata.o \
 $(gwsrc)llnew.o  \
 $(gwsrc)readeigen.o \
 $(gwsrc)readqg.o \
 $(gwsrc)iqindx.o 



all: main2

main1: main1.o $(OBJ)
	$(LD) -o $@  main1.o $(OBJ) $(LIBS) $(LDFLAGS) 

main2: main2.o $(OBJ0)
	$(LD) -o $@ main2.o $(OBJ0) $(GW0)  $(LIBS) $(LDFLAGS)

main1.o : $(OBJ) 

main2.o: $(OBJ) 

m_wancenterdef.o: m_lmlabel.o
m_search_lma.o: m_lmlabel.o 

clean:
	rm -f *.o 
