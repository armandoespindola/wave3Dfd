#Directories
S=src
O=obj
E=bin


_OBJS1=geometry3D.o model.o mpi_trans.o pml.o sdm.o show.o source.o parameters.o\
	source.o receptor.o kernels.o fdt.o utilities.o wave3Dfd.o

OBJS1=$(patsubst %,$(O)/%,$(_OBJS1))


_OBJS2=geometry3D.o model.o mpi_trans.o pml.o sdm.o show.o source.o parameters.o\
	source.o receptor.o kernels.o fdt.o utilities.o MergeSGT.o

OBJS2=$(patsubst %,$(O)/%,$(_OBJS2))

CFLAGS = -std=c++11 -w -O2 -I   
LIBS = -lm

detected_OS := $(shell sh -c 'uname -s 2>/dev/null || echo not')

ifeq ($(detected_OS),Darwin)  # Mac OS X
    CC = g++-7
endif
ifeq ($(detected_OS),Linux)
    CC = mpic++
endif


all: wave3Dfd.out MergeSGT.out

wave3Dfd.out: $(OBJS1) 
	${CC} -o ${E}/wave3Dfd.out $^ $(CFLAGS) $(LIBS)

MergeSGT.out: $(OBJS2)
	${CC} -o ${E}/MergeSGT.out $^ $(CFLAGS) $(LIBS)

$O/%.o: $S/%.cpp $S/%.hpp
	${CC} -c -o $@ $< ${CFLAGS}

clean:
	rm ${O}/*.o ${E}/*.out

