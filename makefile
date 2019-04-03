#Directories
S=src
O=obj
E=bin


_OBJS=geometry3D.o model.o mpi_trans.o pml.o sdm.o show.o source.o parameters.o\
	source.o receptor.o kernels.o fdt.o kernelsw.o utilities.o wave3Dfd.o

OBJS=$(patsubst %,$(O)/%,$(_OBJS))

CFLAGS = -std=c++11 -w -O2 -I   
LIBS = -lm
TARGET = wave3Dfd.out

detected_OS := $(shell sh -c 'uname -s 2>/dev/null || echo not')

ifeq ($(detected_OS),Darwin)  # Mac OS X
    CC = g++-7
endif
ifeq ($(detected_OS),Linux)
    CC = mpic++
endif


all: ${TARGET}

${TARGET}: $(OBJS)
	${CC} -o ${E}/${TARGET} $^ $(CFLAGS) $(LIBS)

$O/%.o: $S/%.cpp $S/%.hpp
	${CC} -c -o $@ $< ${CFLAGS}

clean:
	rm ${O}/*.o ${E}/*.out

