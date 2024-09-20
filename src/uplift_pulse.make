# upllift_pulse.make
# makes the upllift_pulse program.
# make with: make -f upllift_pulse.make

CC = g++
CFLAGS= -c -Wall -O3 -std=c++11 -fPIC
OFLAGS = -Wall -O3 -std=c++11 -fPIC
LDFLAGS= -Wall -fPIC
SOURCES = uplift_pulse.cpp \
		  OneDImplicitHillslope.cpp \
          LSDStatsTools.cpp \
          LSDParameterParser.cpp
OBJ = $(SOURCES:.cpp=.o)
#LIBS = -g -O0 -D_GLIBCXX_DEBUG
LIBS = -Wwrite-strings
EXEC = uplift_puse.out

all: $(SOURCES) $(SCRIPTS) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(OFLAGS) $(OBJ) $(LIBS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@