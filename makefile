DIR = src
SRC = ${DIR}/*.cpp
OBJS = *.o

detected_OS := $(shell sh -c 'uname -s 2>/dev/null || echo not')

ifeq ($(detected_OS),Darwin)  # Mac OS X
    CC = g++-7
endif
ifeq ($(detected_OS),Linux)
    CC = g++
endif

CFLAGS = -c -std=c++11 
LFLAGS =  
TARGET = wave3Dfd.out

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(TARGET)

$(OBJS): $(SRC)
	$(CC) $(CFLAGS) $(SRC)

clean:
	\rm *.o  $(TARGET)
 
run:
	./$(TARGET)

