lsNAME=reduc/Makefile #For debugging purposes
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
CXXFLAGS:= -std=c++11 -fopenmp -Wall
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

$(LIB): $(LIBOBJ)
	@ar r $@ $(LIBOBJ) 

#Try running the compiled target file with an parameter
test: #we run with 4 process, since i only have 4 cores.
	$(MPI) -np 4 $(TARGET) 100 1

atest: #run with: make atest N: num of tasks, T: num of threads.
	$(MPI) -np $(NP) $(TARGET) $(N) $(T)

testm:
	@$(MPI) -np $(NP) $(TARGET) $(N) 1
	@$(MPI) -np $(NP) $(TARGET) $(N) 2
	@$(MPI) -np $(NP) $(TARGET) $(N) 3
	@$(MPI) -np $(NP) $(TARGET) $(N) 4
	@$(MPI) -np $(NP) $(TARGET) $(N) 8
	@$(MPI) -np $(NP) $(TARGET) $(N) 16
	@$(MPI) -np $(NP) $(TARGET) $(N) 32

clean:
	$(RM) $(TARGET) $(OBJS)
	#rm -f *.o 

#For debugging makefile, simply 'make print-VARIABLE' to check content of VARIABLE
print-%: ; @echo $*=$($*)