PETSC_ARCH=intel-opt
PETSC_DIR=${HOME}/project-nufgroup/petsc/src/petsc-3.11.1

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

FFILES = parameters myfft comm grid user variables fsinterface operators_struct operators_fluid timestepper operators_pressure main
MFILES = getdata makeplot make_control001 plotcontrol

IBFILES = ${FFILES:=.o}

#FC = mpif90 
#FLAGS1 = -I/usr/local/fftw/3.3.8/intel-18.0/mvapich2-2.3/include 
#FLAGS2 = -L/usr/local/fftw/3.3.8/intel-18.0/mvapich2-2.3/lib -lfftw3_mpi -lfftw3

#FLAGS1 = -I/usr/local/fftw/3.3.8/intel-18.0/openmpi-3.1.1/include#-I/usr/local/fftw/3.3.8/intel-18.0/mvapich2-2.3/include 
#FLAGS2 = -L/usr/local/fftw/lib -lfftw3_mpi -lfftw3#-L/usr/local/fftw/3.3.8/intel-18.0/mvapich2-2.3/lib -lfftw3_mpi -lfftw3


BASE = /home/nirmal/Studies/Thesis_PhD/Codes/IB_Parallel
SRC = src
OBJ = obj
BIN = bin
OUT = output
IN  = input

DIR = ${SRC} ${OBJ} ${BIN} ${OUT} ${IN}
DIRS = ${DIR:=.dir}

IBFILES_SRC := $(foreach p, $(IBFILES), $(SRC)/$(p))
IBFILES_OBJ := $(foreach p, $(IBFILES), $(OBJ)/$(p))

ib :  ${IBFILES_OBJ}
	-${FLINKER} -o $(BIN)/ib $(IBFILES_OBJ) ${PETSC_LIB} -no-wrap-margin -heap-arrays 0

setup: ${DIRS}
	@for i in ${DIR} ; do \
               mkdir -p $${i} ;\
	done
	@for i in ${FFILES} ; do \
               echo "cp  -p ${BASE}/${SRC}/$${i}.f90 ${SRC}/$${i}.f90 " ;\
               cp  -p ${BASE}/${SRC}/$${i}.f90 ${SRC}/$${i}.f90 ;\
	done
	@for i in ${MFILES} ; do \
               echo "cp  -p ${BASE}/${OUT}/$${i}.m ${OUT}/$${i}.m " ;\
               cp  -p ${BASE}/${OUT}/$${i}.m ${OUT}/$${i}.m ;\
	done

realclean: ${DIRS} # 
	@for i in ${DIR} ; do \
               rm -rf $${i} ;\
	done

clean :: 
	rm -f $(BIN)/ib ${OBJ}/*.o ${OBJ}/*.mod *.mod

$(OBJ)/%.o: $(SRC)/%.F90
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} -o $@ $< $(FC_MODULE_OUTPUT_FLAG)$(OBJ) -no-wrap-margin -heap-arrays 0

%.dir : 	
	$< 
