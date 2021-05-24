.PHONY: greedy fill_gaps random_subset examine build_FM_index tests

INCLUDEPATHS=-I sdsl-lite/build/include/ -I sdsl-lite/build/external/libdivsufsort/include
LIBS=-lsdsl -ldivsufsort -ldivsufsort64 -fopenmp
LIBPATHS=-L sdsl-lite/build/lib/ -L sdsl-lite/build/external/libdivsufsort/lib/
BINARY_DIR=./bin

all:

greedy:
	g++ greedy.cpp -O3 -o $(BINARY_DIR)/greedy -g -Wall -Wno-sign-compare $(LIBPATHS) $(INCLUDEPATHS) $(LIBS)

fill_gaps:
	g++ fill_gaps.cpp -O3 -o $(BINARY_DIR)/fill_gaps -g -Wall -Wno-sign-compare $(LIBPATHS) $(INCLUDEPATHS) $(LIBS)

random_subset:
	g++ random_subset.cpp -O3 -o $(BINARY_DIR)/random_subset

examine:
	g++ examine.cpp -O3 -o $(BINARY_DIR)/examine -g -Wall -Wno-sign-compare $(LIBPATHS) $(INCLUDEPATHS) $(LIBS)

build_FM_index:
	g++ build_FM_index.cpp -O3 -o $(BINARY_DIR)/build_FM_index -g -Wall -Wno-sign-compare $(LIBPATHS) $(INCLUDEPATHS) $(LIBS)

tests:
	g++ test_FM_index.cpp -lgtest -lgtest_main -lpthread -o $(BINARY_DIR)/test_FM_index $(LIBPATHS) $(INCLUDEPATHS) $(LIBS)
	g++ test_greedy.cpp -lgtest -lgtest_main -lpthread -o $(BINARY_DIR)/test_greedy $(LIBPATHS) $(INCLUDEPATHS) $(LIBS)