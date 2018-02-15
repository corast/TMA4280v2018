MSK=$(MAKE) -C

CXX=g++
CXXFLAGS= -std=c+11
BIN=zeta0

zeta0:
	$(CXX) $@ -c

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.o
	rm $(BIN)