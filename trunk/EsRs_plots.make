# make with make -f EsRs_plots.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=EsRs_plot.cpp OneDImplicitHillslope.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=EsRs_plot.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
