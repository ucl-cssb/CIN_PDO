# CCC = icc		# Using intel complier
CCC = g++
# CCC = /usr/local/opt/llvm/bin/clang
BOOST = /usr/local
FLAG = -O3
# FLAG =

all: simpdo


simpdo:
	$(CCC) $(FLAG) -std=gnu++11 sim_cin.cpp -o simpdo -L$(BOOST)/lib/ -lboost_program_options -lboost_system -lboost_filesystem -lgsl -lgslcblas -lm -lz -I$(BOOST)/include

clean:
	rm simpdo
