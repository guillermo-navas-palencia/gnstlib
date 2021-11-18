SHELL=/bin/sh

# GNSTLIB project
PROJECT = gnstlib
PROJECTROOT = ../$(PROJECT)

# debugmode
ifndef
	debugmode = 1
endif

# default compiler
ifndef
	CXX = g++
endif

$(info COMPILER = $(CXX))
$(info DEBUGMODE = $(debugmode))

# OPERATING SYSTEM INFORMATION
# ----------------------------
uname_os = $(shell uname -s)

ifeq ($(findstring CYGWIN,$(uname_os)),CYGWIN)
  platos = cygwin
else ifeq ($(uname_os),Linux)
  platos = linux
else
  platos = $(uname_os)
endif

$(info INFO OS $(platos))

# COMPILER OPTIONS
# ----------------
# compiler flags (g++)
CXXFLAGSTD = -std=gnu++11 -fpic 
CXXFLAGWARNING = -Wall -Wextra -Wunused-parameter -Wno-unused \
                 -Wno-uninitialized -Wconversion -Winline

# optimize
ifneq ($(debugmode),0)
	CXXOPT = -g -ggdb -fopenmp -lquadmath
else
	CXXOPT = -g -O3 -fopenmp -lquadmath 
endif

CFLAGS = $(CXXOPT) $(CXXFLAGSTD) $(CXXFLAGWARNING)
CFLAGSCOVERAGE = $(CXXOPT) $(CXXFLAGSTD)

# coverage
ifeq ($(coverage),1)
	CFLAGS = $(CFLAGSCOVERAGE) -fprofile-arcs -ftest-coverage
	LFLAGS = -lgcov --coverage
else
	LFLAGS = -fopenmp -lquadmath
endif

SOURCEDIR  = $(PROJECTROOT)/src
SOURCETESTDIR = $(PROJECTROOT)/examples
INCLUDEDIR = $(PROJECTROOT)/include

SOURCES = $(wildcard $(SOURCEDIR)/*.cpp) $(wildcard $(SOURCETESTDIR)/*.cpp)

OBJECTS = $(SOURCES:.cpp=.o)

INCL_DIR = -I$(INCLUDEDIR)
LINKER = $(CXX) -c

all: $(SOURCES) $(PROJECT)

$(PROJECT): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) $(INCL_DIR) -c $< -o $@

clean:
	rm -f $(PROJECT) $(OBJECTS)
