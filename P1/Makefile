NAME:=GNUMakefile
LIBS= -lm
CXX= g++
#Use version c++11 to compile, for functions like printf etc.
CXXFLAGS= -std=c++11

#every directory with one project in each, ignore .directories
DIRS := $(shell find ./ -maxdepth 1 -type d -name "[a-z]*" | tr "./" " " )
#Find every object file in the subdirectories for every subtask
find_objects = $(patsubst %.cpp, %.o, $(wildcard $(dir)/*.cpp))
OBJECTS = $(foreach dir, $(DIRS), $(find_objects))

all: zeta0/main mach0/main

zeta0/main: $(patsubst %.cpp, %.o, $(wildcard zeta0/*.cpp))
	$(CXX) $(CXXFLAGS) $< -o $@

zeta0/main.cpp: $(patsubst %.cpp, %.o, $(wildcard zeta0/*.cpp))
	$(CXX) $(CXXFLAGS) $< -o $@

mach0/main: $(patsubst %.cpp, %.o, $(wildcard mach0/*.cpp))
	$(CXX) $(CXXFLAGS) $< -o $@

veryclean:
	$(RM) zeta0/main mach0/main

clean: #clean all folders of object files.
	rm -rf $(OBJECTS)

#For debugging makefile, simply 'make print-VARIABLE' to check content of VARIABLE
print-%: ; @echo $*=$($*)
