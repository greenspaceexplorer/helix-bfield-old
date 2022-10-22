IDIR=include
SDIR=src
DEPS_H=token_parser.h helix_magnet.h coil.h
SOURCES_C=helix_magnet.cpp coil.cpp
MAIN_C=field_calc.cpp

DEPS=$(addprefix $(IDIR)/,$(DEPS_H))
SOURCESC=$(addprefix $(SDIR)/,$(SOURCES_C) $(MAIN_C))
CC=g++
EXECUTABLE=BFIELD
OBJECTSC=$(SOURCESC:.cpp=.o)
CFLAGS=-I$(IDIR) -I/opt/homebrew/include `gsl-config --cflags` -std=c++11
LFLAGS=`gsl-config --libs` -L/opt/homebrew/lib

all: $(EXECUTABLE) $(SOURCESC) $(DEPS)

$(EXECUTABLE): $(OBJECTSC)
	$(CC) -o $@ $(OBJECTSC) $(CFLAGS) $(LFLAGS)

.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -f $(SDIR)/*.o  $(EXECUTABLE)
