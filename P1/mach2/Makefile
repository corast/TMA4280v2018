NAME=reduc/Makefile #For debugging purposes
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

$(LIB): $(LIBOBJ)
	@ar r $@ $(LIBOBJ)

#Try running the compiled target file with an parameter
test: #we run with 4 process, since i only have 4 cores.
	$(MPI) -np 4 $(TARGET) 100 1 0

atest: #run with: make atest NP="" N="" M="" T="" # M: "0"=zeta, "1"=mach ; T:"0"=Allreduce, "1"=recursive-doubling.
	$(MPI) -np $(NP) $(TARGET) $(N) $(M) $(T)

clean:
	$(RM) $(TARGET) $(OBJS) utestzeta0 vtestzeta0
	#rm -f *.o 

#For debugging makefile, simply 'make print-VARIABLE' to check content of VARIABLE
print-%: ; @echo $*=$($*)