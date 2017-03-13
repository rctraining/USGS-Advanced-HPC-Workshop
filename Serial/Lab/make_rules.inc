OBJS = $(SRCS:.f90=.o)

CC = gcc
CXX = g++
FC = gfortran
LD = g++

# for intel using icc and ifort
#CC = icc
#FC = ifort

# Note:
# The opt-report levels for Intel compilers are:
# 0: No vectorization report
# 1: Reports when vectorization has occurred.
# 2: Adds diagnostics why vectorization did not occur.
# 3: Adds vectorization loop summary diagnostics.
# 4: Additional detail, e.g. on data alignment
# 5: Adds detailed data dependency information
#
# Compiler flags
#
FCFLAGS = -Wall -Werror -Wextra
FCFLAGS_DBG = -g -O0 -DDEBUG
FCFLAGS_REL = -Ofast -ftree-vectorize -fopt-info-vec-all

CXX_FLAGS =  -Wall -Werror -Wextra
CXXFLAGS_DBG = -g -O0 -DDEBUG
CXXFLAGS_REL = -Ofast -DNDEBUG -ftree-vectorize -fopt-info-vec-all

ifeq (.f90 ,$(findstring .f90, $(SRCS)))
	LD = $(FC)
endif
#
# Debug build settings
#
DBG_DIR = Debug
DBG_EXE = $(DBG_DIR)/$(EXE)
DBG_OBJS = $(addprefix $(DBG_DIR)/, $(OBJS))

#
# Release build settings
#
REL_DIR = Release
REL_EXE = $(REL_DIR)/$(EXE)
REL_OBJS = $(addprefix $(REL_DIR)/, $(OBJS))

.SUFFIXES: .c .cpp .f90 .o
.PHONEY: all clean debug prep release remake



# Default build
all: prep release

debug: $(DBG_EXE)

$(DBG_EXE): $(DBG_OBJS)
	$(LD) $(LDFLAGS) -o $(DBG_EXE) $^

$(DBG_DIR)/%.o: %.f90
	$(FC) -c $(FCFLAGS) $(FCFLAGS_DBG) -o $@ $<

$(DBG_DIR)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(CXXFLAGS_DBG) -o $@ $<


#
# Release rules
#
release: $(REL_EXE)

$(REL_EXE): $(REL_OBJS)
	$(LD) $(LDFLAGS) -o $(REL_EXE) $^

$(REL_DIR)/%.o: %.f90
	$(FC) -c $(FCFLAGS) $(FCFLAGS_REL) -o $@ $<

$(REL_DIR)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(CXXFLAGS_REL) -o $@ $<

#
# Other rules
#
prep:
	@mkdir -p $(DBG_DIR) $(REL_DIR)

remake: clean all

clean:
	rm -f $(REL_EXE) $(REL_OBJS) $(DBG_EXE) $(DBG_OBJS)

