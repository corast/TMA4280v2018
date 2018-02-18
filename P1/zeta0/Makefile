NAME=zeta0/Makefile
SRC:=$(wildcard *.cpp)
OBJ:=$(patsubst %.cpp,%.o, $(SRC))
MAINSRC=main.cpp
#our target is a file called main
TARGET:=$(patsubst %.cpp,%,$(MAINSRC))
CXX=g++
#Use version c++11 to compile, for functions like printf etc.
CXXFLAGS= -std=c++11 -Wall

#LDFLAGS = -lGL -lglut -lpng -lz -lm
MAINOBJ=$(patsubst %.c,%.o,$(MAINSRC))
LIBOBJ=$(patsubst $(MAINOBJ),,$(OBJ))

‰.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

all: $(LIB) $(MAINOBJ)
	$(CXX) -o $(CXXFLAGS) -c $<

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -c $< -o $@
#include $(OBJ:.o=.d)

#$(LIB): $(LIBOBJ)
#	@ar r $@ $(LIBOBJ)

.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJ)

#For debugging makefile, simply 'make print-VARIABLE' to check content of VARIABLE
print-%: ; @echo $*=$($*)