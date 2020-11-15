
WORKDIR = `pwd` 
 
CC = gcc 
CXX = g++ 
AR = ar 
LD = g++ 
WINDRES = windres 
 
INC = -I./
CFLAGS = -Wall -g -pthread  
RESINC = 
LIBDIR = -L./ -pthread -g 
LIB = 
LDFLAGS = 
 
INC_RELEASE = $(INC) -I./
CFLAGS_RELEASE = $(CFLAGS) -march=corei7-avx -O -std=c++11
 
RESINC_RELEASE = $(RESINC) -I./
RCFLAGS_RELEASE = $(RCFLAGS) 
LIBDIR_RELEASE = $(LIBDIR) -L./
LIB_RELEASE = $(LIB)
LDFLAGS_RELEASE = $(LDFLAGS) -s -lz
OBJDIR_RELEASE = obj/Release
DEP_RELEASE = 
OUT_RELEASE = bin/Release/quantization

 
OBJ_RELEASE = $(OBJDIR_RELEASE)/util/vstring.o $(OBJDIR_RELEASE)/cnpy/cnpy.o $(OBJDIR_RELEASE)/nn/geo.o $(OBJDIR_RELEASE)/nn/he.o $(OBJDIR_RELEASE)/info.o $(OBJDIR_RELEASE)/util/timer.o $(OBJDIR_RELEASE)/loadSift.o $(OBJDIR_RELEASE)/main.o $(OBJDIR_RELEASE)/missionagent.o $(OBJDIR_RELEASE)/nn/bruteforcenn.o $(OBJDIR_RELEASE)/nn/geohenn.o $(OBJDIR_RELEASE)/clustering/nngkmeans.o $(OBJDIR_RELEASE)/util/cleaner.o $(OBJDIR_RELEASE)/clustering/evaluator.o $(OBJDIR_RELEASE)/util/ioagent.o $(OBJDIR_RELEASE)/clustering/nng.o $(OBJDIR_RELEASE)/clustering/abstractkmeans.o $(OBJDIR_RELEASE)/util/pqmath.o $(OBJDIR_RELEASE)/util/print2scrn.o $(OBJDIR_RELEASE)/clustering/randompartition.o $(OBJDIR_RELEASE)/util/scriptparser.o 
 

all: release

 
clean: clean_release
 

before_release: 
	test -d bin/Release || mkdir -p bin/Release
	test -d $(OBJDIR_RELEASE)/clustering || mkdir -p $(OBJDIR_RELEASE)/clustering
	test -d $(OBJDIR_RELEASE)/cnpy || mkdir -p $(OBJDIR_RELEASE)/cnpy
	test -d $(OBJDIR_RELEASE)/util || mkdir -p $(OBJDIR_RELEASE)/util 
	test -d $(OBJDIR_RELEASE) || mkdir -p $(OBJDIR_RELEASE)
	test -d $(OBJDIR_RELEASE)/nn || mkdir -p $(OBJDIR_RELEASE)/nn
 

after_release: 

 
release: before_release out_release after_release 

 
out_release: before_release $(OBJ_RELEASE) $(DEP_RELEASE) 
	$(LD) $(LIBDIR_RELEASE) -o $(OUT_RELEASE) $(OBJ_RELEASE)  $(LDFLAGS_RELEASE) $(LIB_RELEASE) 

 
$(OBJDIR_RELEASE)/util/vstring.o: util/vstring.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c util/vstring.cpp -o $(OBJDIR_RELEASE)/util/vstring.o 
 
$(OBJDIR_RELEASE)/cnpy/cnpy.o: cnpy/cnpy.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c cnpy/cnpy.cpp -o $(OBJDIR_RELEASE)/cnpy/cnpy.o 
 
 
$(OBJDIR_RELEASE)/nn/geo.o: nn/geo.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c nn/geo.cpp -o $(OBJDIR_RELEASE)/nn/geo.o 
 
$(OBJDIR_RELEASE)/nn/he.o: nn/he.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c nn/he.cpp -o $(OBJDIR_RELEASE)/nn/he.o 
 
$(OBJDIR_RELEASE)/info.o: info.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c info.cpp -o $(OBJDIR_RELEASE)/info.o 
 
$(OBJDIR_RELEASE)/util/timer.o: util/timer.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c util/timer.cpp -o $(OBJDIR_RELEASE)/util/timer.o 
 
$(OBJDIR_RELEASE)/loadSift.o: loadSift.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c loadSift.cpp -o $(OBJDIR_RELEASE)/loadSift.o 
 
$(OBJDIR_RELEASE)/main.o: main.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c main.cpp -o $(OBJDIR_RELEASE)/main.o 
 
$(OBJDIR_RELEASE)/missionagent.o: missionagent.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c missionagent.cpp -o $(OBJDIR_RELEASE)/missionagent.o 
 
$(OBJDIR_RELEASE)/nn/bruteforcenn.o: nn/bruteforcenn.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c nn/bruteforcenn.cpp -o $(OBJDIR_RELEASE)/nn/bruteforcenn.o 
 
$(OBJDIR_RELEASE)/nn/geohenn.o: nn/geohenn.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c nn/geohenn.cpp -o $(OBJDIR_RELEASE)/nn/geohenn.o 
 
$(OBJDIR_RELEASE)/clustering/nngkmeans.o: clustering/nngkmeans.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c clustering/nngkmeans.cpp -o $(OBJDIR_RELEASE)/clustering/nngkmeans.o 
 
$(OBJDIR_RELEASE)/util/cleaner.o: util/cleaner.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c util/cleaner.cpp -o $(OBJDIR_RELEASE)/util/cleaner.o 
 
$(OBJDIR_RELEASE)/clustering/evaluator.o: clustering/evaluator.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c clustering/evaluator.cpp -o $(OBJDIR_RELEASE)/clustering/evaluator.o 
 
$(OBJDIR_RELEASE)/util/ioagent.o: util/ioagent.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c util/ioagent.cpp -o $(OBJDIR_RELEASE)/util/ioagent.o 
 
$(OBJDIR_RELEASE)/clustering/nng.o: clustering/nng.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c clustering/nng.cpp -o $(OBJDIR_RELEASE)/clustering/nng.o 
 
$(OBJDIR_RELEASE)/clustering/abstractkmeans.o: clustering/abstractkmeans.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c clustering/abstractkmeans.cpp -o $(OBJDIR_RELEASE)/clustering/abstractkmeans.o 
 
$(OBJDIR_RELEASE)/util/pqmath.o: util/pqmath.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c util/pqmath.cpp -o $(OBJDIR_RELEASE)/util/pqmath.o 
 
$(OBJDIR_RELEASE)/util/print2scrn.o: util/print2scrn.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c util/print2scrn.cpp -o $(OBJDIR_RELEASE)/util/print2scrn.o 
 
$(OBJDIR_RELEASE)/clustering/randompartition.o: clustering/randompartition.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c clustering/randompartition.cpp -o $(OBJDIR_RELEASE)/clustering/randompartition.o 
 
$(OBJDIR_RELEASE)/util/scriptparser.o: util/scriptparser.cpp 
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c util/scriptparser.cpp -o $(OBJDIR_RELEASE)/util/scriptparser.o 

 
clean_release: 
	rm -f $(OBJ_RELEASE) $(OUT_RELEASE) 
	rm -rf bin/Release 
	rm -rf $(OBJDIR_RELEASE)/clustering 
	rm -rf $(OBJDIR_RELEASE)/cnpy 
	rm -rf $(OBJDIR_RELEASE)/util 
	rm -rf $(OBJDIR_RELEASE) 
	rm -rf $(OBJDIR_RELEASE)/nn 
 

.PHONY: before_release after_release clean_release
