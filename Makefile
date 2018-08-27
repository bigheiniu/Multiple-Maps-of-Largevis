CXX = /usr/local/Cellar/gcc/7.3.0_1/bin/g++-7
CXXFLAGS = -lm -pthread -lgsl -lgslcblas -Ofast -march=native -ffast-math
DFLAG += -DUSEOMP
CXXFLAGS += -fopenmp
all: mmtsne
mmtsne: LargeVis.cpp  mmtsne.o
	$(CXX) $(CXXFLAGS) -o $@ $^
mmtsne.o: LargeVis.cpp	main.cpp	LargeVis.h
	$(CXX) $(CXXFLAGS) $(DFLAG) -c -o $@ $<
clean:
	rm -f mmtsne mmtsne.o
