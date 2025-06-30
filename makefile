## Makefile for desktop

CC := g++
SRCDIR := src
BUILDDIR := build
BUILDDIR_DEBUG := build_debug
TARGET := bin/model
TARGETDB := bin/model_debug
INCDIR := include
EIGEN := /Users/user/Documents/Sheffield/ActinModelling/eigen
GEO := GeometricTools/GTEngine/include
# BOOST := /opt/homebrew/opt/boost
OMP := /opt/homebrew
LIBOMP := /opt/homebrew/opt/libomp

RM := rm -f


CFLAGS := -O3 -std=c++14 -DEIGEN_NO_DEBUG -Wall -arch arm64 \
  -Xpreprocessor -fopenmp -I$(LIBOMP)/include -DEIGEN_DISABLE_FP16
CFLAGS_DEBUG := -Wall -std=c++14 -DBOOST_TEST_DYN_LINK -g -arch arm64 \
  -Xpreprocessor -fopenmp -I$(LIBOMP)/include -DEIGEN_DISABLE_FP16

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS_DEBUG := $(patsubst $(SRCDIR)/%,$(BUILDDIR_DEBUG)/%,$(SOURCES:.$(SRCEXT)=.o))

INC := -I$(INCDIR) -I$(EIGEN) -I$(GEO) -I$(BOOST)/include
LIB := -L$(BOOST)/lib -lboost_program_options -lutil \
  -L$(LIBOMP)/lib -lomp

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(BUILDDIR_DEBUG)/%.o: $(SRCDIR)/%.$(SRCEXT)
	g++ -O3 -std=c++14 -Wall -arch arm64 -I/opt/homebrew/include -Iinclude -c $< -o $@
	@echo "check 4"
	@mkdir -p $(BUILDDIR_DEBUG)
	@echo " $(CC) $(CFLAGS_DEBUG) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS_DEBUG) $(INC) -c -o $@ $<
	@echo "check 5"
	
debug: $(OBJECTS_DEBUG)
	$(CC) $(CFLAGS_DEBUG) $(OBJECTS_DEBUG) $(INC) $(LIB) -o $(TARGETDB)

.PHONY: clean
clean:
	@echo " Cleaning..."; \
	$(RM) -r $(BUILDDIR) $(TARGET)

clean_debug:
	@echo " Cleaning..."; \
	$(RM) -r $(BUILDDIR_DEBUG) $(TARGETDB)
