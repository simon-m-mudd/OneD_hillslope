# make with make -f test_hillslope.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=Test_hillslope.cpp OneDImplicitHillslope.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Test_hillslope.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
