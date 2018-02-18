NAME = Gnu
LIBS= -lm
CXX= g++
#Use version c++11 to compile, for functions like printf etc.
CXXFLAGS= -std=c++11 -Wall

#every directory with one project in each
DIRS = $(shell find ./ -maxdepth 1 -type d -name "[a-z]*" | tr "./" " " )
#Find every object file in the subdirectories for every subtask
find_objects = $(patsubst %.cpp, %.o, $(wildcard $(dir)/*.cpp))
OBJECTS_T = $(foreach dir, $(DIRS), $(find_objects))

#find every target to clean
find_target = $(shell find -maxdepth 1 -name "main" | tr "./" " " )

OBJECTS_ZETA0 = $(patsubst %.cpp, %.o, $(wildcard zeta0/*.cpp))

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

zeta0/main.cpp: $(OBJECTS_ZETA0)
	$(CXX) $(CXXFLAGS) -c $< -o $@
	$< 100

clean: #clean all folders of object files.
	rm -rf $(OBJECTS_T)

#For debugging makefile, simply 'make print-VARIABLE' to check content of VARIABLE
print-%: ; @echo $*=$($*)

