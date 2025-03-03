# Copyright 2025 Bianca Maria Laudenzi, Caterina Dalmaso, Lucas Omar Muller
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



# define compiler options
# GNU compiler family:
# Compiler
CPP=g++
## Optimization flags
CPPFLAGS= -O3 #-ffast-math #-march=native -freorder-blocks-and-partition -floop-interchange -floop-strip-mine -floop-block  -ftree-loop-if-convert -ftree-loop-if-convert-stores -ftree-loop-distribution -ftree-vectorize -fno-exceptions -fno-rtti 

# macros
RM = rm -f

# You have to modify this dsirectory!!!
DIRSOURCE=../source_files
# libraries
LIBDEPS = $(DIRSOURCE)/libcardiorespsolver.a
LDFLAGS = $(DIRSOURCE)/libcardiorespsolver.a -lm 
# paths to include files (out of current directory)
INCLUDES = -I$(DIRSOURCE)
# output
EXEC=cardiorespiratory.x

# input

INPUT=main.cc


DIRS := ./
FILES:= $(INPUT) $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.h))
SRCS := $(notdir $(INPUT))
# rule to create object files
OBJS := $(SRCS:.cc=.o)

 
.PHONY: all clean

all: $(EXEC) $(TARGET) 

$(EXEC): $(OBJS) $(LIBDEPS)
	$(CPP) $(CPPFLAGS) $(INCLUDES) $(OBJS) $(LDFLAGS) -o $@

%.o : %.cc
	$(CPP) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

clean: 	
	rm -f *.o cardiorespiratory.x 
