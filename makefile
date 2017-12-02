#Directories
S = src
O = obj
E = bin

SRC = ${S}/*.cpp
OBJS = geometry3D.o  model.o  mpi_trans.o  pml.o  sdm.o  show.o  source.o\
  wave3Dfd.o

detected_OS := $(shell sh -c 'uname -s 2>/dev/null || echo not')

ifeq ($(detected_OS),Darwin)  # Mac OS X
    CC = g++-7
endif
ifeq ($(detected_OS),Linux)
    CC = mpic++
endif

CFLAGS =  -fopenmp -std=c++11 -w -O2   
LFLAGS = -w -fopenmp -O2 
TARGET = wave3Dfd.out

all: ${TARGET}

${TARGET}: $O/geometry3D.o  $O/model.o  $O/mpi_trans.o  $O/pml.o  $O/sdm.o \
$O/show.o  $O/source.o $O/wave3Dfd.o
	${CC} $(LFLAGS) -o ${E}/${TARGET} $O/*.o

$O/%.o: $S/%.cpp
	${CC} ${CFLAGS} -c -o $@ $<

clean:
	rm ${O}/*.o ${E}/*.out

