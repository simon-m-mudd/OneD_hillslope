# make with make -f CRN_oneD_hs.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=CRN_oneD_hs.cpp OneDImplicitHillslope.cpp LSDParticle.cpp LSDCRNParameters.cpp CRN_funcs.cpp LSDStatsTools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=CRN_oneD_hs.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
