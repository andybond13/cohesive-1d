CPP       = mpic++
CPP_FLAGS = -std=c++17 -Wall -fPIC -g     #debug
#CPP_FLAGS = -std=c++17 -Wall -fPIC -O2    #production
SOFTWARE  = ./utilities
BOOST_ROOT = /usr/local/Cellar/boost/1.73.0
LIBS = -lboost_system -lboost_filesystem

INCS      = -I./inc
SRC       = ./src/
OBJ       = ./bin/
OBJMAIN   = $(OBJ)
OBJSERV   = $(OBJ)/service/
OBJCLIE   = $(OBJ)/client/
TEST      = ./test/

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	LIBS += -L $(BOOST_ROOT)stage/lib/
	INCS += -I $(BOOST_ROOT) -I$(SOFTWARE)/utilities 
endif
ifeq ($(UNAME), Darwin)
#	# do something OSX-y
    LIBS += -L $(BOOST_ROOT)/lib/
    INCS += -I $(BOOST_ROOT)/include/
endif

#check 32 or 64 bit
LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
   # do 64 bit stuff here
   CPP_FLAGS += -m64
else
   # do 32 bit stuff here
   CPP_FLAGS += -m32
endif


#all    : main.exe service client.exe
#service: $(OBJSERV)connectCR.o $(OBJSERV)connectMPG.o $(OBJSERV)libCartRing.so $(OBJSERV)libMatPropGen.so $(OBJSERV)CartRing.exe $(OBJSERV)MatPropGen.exe


# Classical compilation of the main
main.exe: $(OBJMAIN)main.o $(OBJMAIN)CartRing.o $(OBJMAIN)MatPropGen.o $(OBJMAIN)ParallelCombiner.o
	$(CPP) $(CPP_FLAGS) -o $(OBJMAIN)main.exe $(OBJMAIN)main.o $(OBJMAIN)CartRing.o $(OBJMAIN)MatPropGen.o $(OBJMAIN)ParallelCombiner.o $(LIBS)

$(OBJMAIN)main.o: $(TEST)main.cpp
	$(CPP) $(CPP_FLAGS) -c $(INCS) $(TEST)main.cpp -o $(OBJMAIN)main.o
	
$(OBJMAIN)CartRing.o: $(SRC)cartRing.C
	$(CPP) $(CPP_FLAGS) -c $(INCS) $(SRC)cartRing.C -o $(OBJMAIN)CartRing.o

$(OBJMAIN)MatPropGen.o: $(SRC)matPropGen.C
	$(CPP) $(CPP_FLAGS) -c $(INCS) $(SRC)matPropGen.C -o $(OBJMAIN)MatPropGen.o

$(OBJMAIN)ParallelCombiner.o: $(SRC)parallelCombiner.C
	$(CPP) $(CPP_FLAGS) -c $(INCS) $(SRC)parallelCombiner.C -o $(OBJMAIN)ParallelCombiner.o


clean:
	rm $(OBJMAIN)*.o $(OBJMAIN)*.exe
