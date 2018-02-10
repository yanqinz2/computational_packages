EXEC = $(wildcard ./include/*.hpp)
EXEC := $(EXEC:./include/%.hpp=%)

CXX = g++-7
CXXFLAGS = -c -Wall -std=c++11 -O3 -I./include -DBOOST -DOMP -fopenmp -I/usr/local/include
MAINFLAGS = -DOMP -fopenmp -I/usr/local/include

.PHONY: all
all: MakeBin $(EXEC)

$(EXEC): %: ./src/%_main.o ./src/%.o
	$(CXX) $(MAINFLAGS) $^ -o ./bin/$@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

MakeBin:
	@ mkdir -p ./bin

.PHONY: clean clean-o
clean:
	rm -rf ./src/*.o ./bin/*
clean-o:
	rm -rf ./src/*.o
