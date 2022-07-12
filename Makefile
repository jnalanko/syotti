.PHONY: greedy fill_gaps random_subset examine build_FM_index tests

INCLUDEPATHS=-I sdsl-lite/build/include/ -I sdsl-lite/build/external/libdivsufsort/include
LIBS=-lsdsl -ldivsufsort -ldivsufsort64 -fopenmp
LIBPATHS=-L sdsl-lite/build/lib/ -L sdsl-lite/build/external/libdivsufsort/lib/
BINARY_DIR=./bin
SRC_DIR=./src

syotti:
	g++ $(SRC_DIR)/CLI.cpp $(SRC_DIR)/examine.cpp $(SRC_DIR)/fill_gaps.cpp $(SRC_DIR)/greedy.cpp $(SRC_DIR)/random_subset.cpp $(SRC_DIR)/FM_index.cpp $(SRC_DIR)/build_FM_index.cpp $(SRC_DIR)/input_reading.cpp $(SRC_DIR)/util.cpp -O3 -o $(BINARY_DIR)/syotti -g -Wall -Wno-sign-compare $(LIBPATHS) $(INCLUDEPATHS) $(LIBS)

tests:
	g++ $(SRC_DIR)/test_FM_index.cpp -L googletest/build/lib/ -lgtest -lgtest_main -lpthread -o $(BINARY_DIR)/test_FM_index $(LIBPATHS) $(INCLUDEPATHS) $(LIBS)
	g++ $(SRC_DIR)/test_greedy.cpp -L googletest/build/lib/ -lgtest -lgtest_main -lpthread -o $(BINARY_DIR)/test_greedy $(LIBPATHS) $(INCLUDEPATHS) $(LIBS)