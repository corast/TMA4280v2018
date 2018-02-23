NAME=mach0/Makefile #For debugging purposes
CEXE:=.cpp
SRC:=$(wildcard *$(CEXE)) #Collect all .cpp files.

OBJS:=$(patsubst %$(CEXE),%.o, $(SRC))

MAINSRC:=main.cpp

#our target is a file called main, remove .cpp extention to MAINSRC
TARGET:=main
#$(patsubst %$(CEXE),%,$(MAINSRC))
LIB=libutils.a

CXX:=g++
#Use version c++11 to compile, for functionality like printf etc.
CXXFLAGS:= -std=c++11

#LDFLAGS = -lGL -lglut -lpng -lz -lm
MAINOBJ=$(patsubst %$(CEXE),%.o,$(MAINSRC))
LIBOBJ=$(patsubst $(MAINOBJ),,$(OBJS))

.PHONY: all clean

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<

all: $(LIB) $(MAINOBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(MAINOBJ) $(LIB)
#Copmile the utestzeta0 program. (don't know a better way to do this)
	$(CXX) $(CXXFLAGS) -o utestmach0 utestmach0.o $(LIB)
	$(CXX) $(CXXFLAGS) -o vtestmach0 vtestmach0.o $(LIB)

$(LIB): $(LIBOBJ)
	@ar r $@ $(LIBOBJ)

#Try running the compiled target file with an parameter
run:
	./$(TARGET) 100

utest:
	./utestmach0

vtest:
	./vtestmach0

clean:
	$(RM) $(TARGET) $(OBJS) utestmach0 vtestmach0

#For debugging makefile, simply 'make print-VARIABLE' to check content of VARIABLE
print-%: ; @echo $*=$($*)