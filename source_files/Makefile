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
# Optimization flags
<<<<<<< HEAD
CPPFLAGS= -O3 
=======
CPPFLAGS= -Ofast #-ffast-math #-march=native -freorder-blocks-and-partition -floop-interchange -floop-strip-mine -floop-block  -ftree-loop-if-convert -ftree-loop-if-convert-stores -ftree-loop-distribution -ftree-vectorize -fno-exceptions -fno-rtti
>>>>>>> origin/main

# GIT
GIT_VERSION := $(shell git describe --dirty --always --tags)
CPPFLAGS+= -DVERSION=\"$(GIT_VERSION)\"

# macros
RM = rm -f
# libraries
LDFLAGS = -lm 

# paths to include files (out of current directory)
INCLUDES = -I/usr/include

HASRANLIB = t

# output
TARGET:= ./libcardiorespsolver.a
TARGETLIBDIR = ./


DIRS := ./
TEMP := $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.cc))
FILES:= $(TEMP) $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.h))
SRCS := $(notdir $(TEMP))
# rule to create object files
OBJS := $(SRCS:.cc=.o)
OBJSD:= $(addprefix ./$(TARGETLIBDIR),$(OBJS))

.PHONY: all clean

all: $(TARGET) 

./$(TARGETLIBDIR)%.o : %.cc
	$(CPP) $(CPPFLAGS) $(MUDFLAP_OPTIONS) $(CPPFLAGSDEF) $(INCLUDES) -c $< -o $@

$(TARGET): $(OBJSD) 
	@ ar rv $(TARGET) $?
	@ case x$(HASRANLIB) in xt ) echo ranlib; ranlib $(TARGET);; esac

%.o : %.cc
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

clean: 	
	rm -f *.o libcardiorespsolver.a *.opari.inc *.mod.cc
