TARGET = main
#MSK=$(MAKE) -C
MODULES = mach0 zeta0
LIBS = -lm
CXX=g++
CXXFLAGS= -std=c+11 -g -Wall
#Use version c++11 to compile, for functions like printf etc.

.PHONY: default all clean
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
#store the object files from corresponding .c files.
HEADERS = $(wildcard *.h)
#Store the header files.

%.o: %.c $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -Wall $(LIBS) -o $@

#Clean rules
clean:
	rm -f *.o
	rm -f $(TARGET)