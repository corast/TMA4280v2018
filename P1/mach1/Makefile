NAME=mach1/Makefile #For debugging purposes
CEXE:=.cpp
MPI=mpirun

SRC:=$(wildcard *$(CEXE)) #Collect all .cpp files.

OBJS:=$(patsubst %$(CEXE),%.o, $(SRC)) #collect all object files

MAINSRC:=main.cpp

#our target is a file called main, remove .cpp extention to MAINSRC
TARGET:=main
#$(patsubst %$(CEXE),%,$(MAINSRC))
LIB=libutils.a #library file containing all the compiled libraries dependensies.

CXX:=mpicxx
#Use mpic++ compiler, version c++11 to compile, for functionality like printf etc.
CXXFLAGS:= -std=c++11

#LDFLAGS = -lGL -lglut -lpng -lz -lm
MAINOBJ=$(patsubst %$(CEXE),%.o,$(MAINSRC))# target object file.
LIBOBJ=$(patsubst $(MAINOBJ),,$(OBJS))# all object file dependant on the library file

.PHONY: all clean
#rule: Every object file is dependant on corresponding .cpp file
%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<

#default function running make.
all: $(LIB) $(MAINOBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(MAINOBJ) $(LIB)
#Copmile the test programs. (don't know a better way to do this)
	#$(CXX) $(CXXFLAGS) -o utestzeta0 utestzeta0.o $(LIB)
	#$(CXX) $(CXXFLAGS) -o vtestzeta0 vtestzeta0.o $(LIB)

$(LIB): $(LIBOBJ)
	@ar r $@ $(LIBOBJ)

#Try running the compiled target file with an parameter
test: #we run with default 4 processes
	$(MPI) -np 4 $(TARGET) 100

atest: #to run make NP="4" N="100" atest
	$(MPI) -np $(NP) $(TARGET) $(N)

testm: #to run make NP="4" N="100" atest
	@$(MPI) -np $(NP) $(TARGET) 8
	@$(MPI) -np $(NP) $(TARGET) 16
	@$(MPI) -np $(NP) $(TARGET) 128
	@$(MPI) -np $(NP) $(TARGET) 256

clean:
	$(RM) $(TARGET) $(OBJS)
	#rm -f *.o 

#For debugging makefile, simply 'make print-VARIABLE' to check content of VARIABLE
print-%: ; @echo $*=$($*)