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

CXX:=g++
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
	./$(TARGET) 100 4

atest: #run with: make atest N: num of tasks, T: num of threads.
	 ./$(TARGET) $(N) $(T)

mtest: 
	@./$(TARGET) 1000000 1
	@./$(TARGET) 1000000 2
	@./$(TARGET) 1000000 3
	@./$(TARGET) 1000000 4
	@./$(TARGET) 1000000 5
	@./$(TARGET) 1000000 6
	@./$(TARGET) 1000000 7
	@./$(TARGET) 1000000 8
	@./$(TARGET) 1000000 16
	@./$(TARGET) 1000000 32
	@./$(TARGET) 1000000 64
clean:
	$(RM) $(TARGET) $(OBJS)
	#rm -f *.o 

#For debugging makefile, simply 'make print-VARIABLE' to check content of VARIABLE
print-%: ; @echo $*=$($*)