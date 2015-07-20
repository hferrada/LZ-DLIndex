CPP=g++ -std=c++11
CPPFLAGS=-O9 -Wall -DVERBOSE
INCLUDES=-I/home/hferrada/include/ -I/home/hferrada/drf/dir64/DRF_Utils64/includes/
LIB=/home/hferrada/lib/libsdsl.a /home/hferrada/lib/libdivsufsort.a /home/hferrada/lib/libdivsufsort64.a /home/hferrada/drf/dir64/DRF_Utils64/drflib64.a
OBJECTS=LZ78Tries64.o LZDocList64.o
BINS=build_index

%.o: %.cpp
	@echo " [C++] Compiling $<"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

all: clean stats $(OBJECTS) $(BINS)

stats:
	@echo
	@echo " COMPILING appTopkLZhq.a"
	@echo " ###################"
	@echo "  * Compiler flags: $(CPPFLAGS)"
	@echo "  * Include dirs: $(INCLUDES)"
	@echo "  * Lib dirs: $(LIB)"
	@echo

clean:
	@echo " [CLN] Removing object files"
	@rm -f $(OBJECTS) $(BINS)

build_index: 
	@echo " [BLD] Building build_index"
	ar -rvcs DL_LZ.a $(OBJECTS) $(LIB) 

build_binary:
	@echo " [BLD] Building binary buildDL_LZ"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o buildDL_LZ buildDL_LZ.cpp $(OBJECTS) $(LIB)

load_binary:
	@echo " [BLD] Building binary loadDL_LZ"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o loadDL_LZ loadDL_LZ.cpp $(OBJECTS) $(LIB)
