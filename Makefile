CPP       = mpic++
CPP_FLAGS = -O2 -Wall -fPIC #-g -Wall
SOFTWARE  = ./utilities
BOOST_ROOT = $(SOFTWARE)/boost/
LIBS = -lboost_system -lboost_filesystem

#CTL       = $(SOFTWARE)/utils/ctl
#CTLINC    = $(CTL)/include
#CTLLIB    = $(CTL)/lib/libctl_g++Linuxx86_64/libctl.so
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
	INCS += -I$(BOOST_ROOT) -I./ci -I$(CTLINC) -I$(SOFTWARE)/utilities 
endif
ifeq ($(UNAME), Darwin)
	# do something OSX-y
	LIBS += -L /sw/opt/boost-1_55/lib/
	INCS += -I /sw/opt/boost-1_55/include/ -I /sw/include
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


# Compilation of the componant as a shared library and a remote executables
#$(OBJSERV)connectCR.o: $(SRC)connectCR.cpp
#	$(CPP) $(CPP_FLAGS) -c $(INCS) $(SRC)connectCR.cpp -o $(OBJSERV)connectCR.o

#$(OBJSERV)connectMPG.o: $(SRC)connectMPG.cpp
#	$(CPP) $(CPP_FLAGS) -c $(INCS) $(SRC)connectMPG.cpp -o $(OBJSERV)connectMPG.o

#$(OBJSERV)libCartRing.so: $(OBJSERV)connectCR.o $(OBJMAIN)CartRing.o
#	$(CPP) $(CPP_FLAGS) -shared -o $(OBJSERV)libCartRing.so $(OBJSERV)connectCR.o $(OBJMAIN)CartRing.o

#$(OBJSERV)libMatPropGen.so: $(OBJSERV)connectMPG.o $(OBJMAIN)MatPropGen.o
#	$(CPP) $(CPP_FLAGS) -shared -o $(OBJSERV)libMatPropGen.so $(OBJSERV)connectMPG.o $(OBJMAIN)MatPropGen.o

#$(OBJSERV)CartRing.exe: $(OBJSERV)connectCR.o $(OBJMAIN)CartRing.o
#	$(CPP) $(CPP_FLAGS) -o $(OBJSERV)CartRing.exe $(OBJSERV)connectCR.o $(OBJMAIN)CartRing.o $(CTLLIB)

#$(OBJSERV)MatPropGen.exe: $(OBJSERV)connectMPG.o $(OBJMAIN)MatPropGen.o
#	$(CPP) $(CPP_FLAGS) -o $(OBJSERV)MatPropGen.exe $(OBJSERV)connectMPG.o $(OBJMAIN)MatPropGen.o $(CTLLIB)


# compilation d'un client faisant appel au service
#client.exe: $(OBJCLIE)client.o
#	$(CPP) $(CPP_FLAGS) -o $(OBJCLIE)client.exe $(OBJCLIE)client.o $(CTLLIB)

#$(OBJCLIE)client.o: $(TEST)client.cpp
#	$(CPP) $(CPP_FLAGS) -c $(INCS) $(TEST)client.cpp -o $(OBJCLIE)client.o


clean:
	rm $(OBJMAIN)*.o $(OBJMAIN)*.exe
	#rm $(OBJSERV)*.o $(OBJSERV)*.so $(OBJSERV)*.exe
	#rm $(OBJCLIE)*.o $(OBJCLIE)*.exe
