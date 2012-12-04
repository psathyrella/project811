CC=g++ -O2  -g
CCFLAGS = `root-config --cflags` -I /afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms4/include
LDFLAGS = `root-config --libs` -L /afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms4/lib  -l RooFit -l RooFitCore

EXE = test.exe
OBJ = MitStyleRemix.o MyPdfV2.o
#test.cc -o bin/test.exe
all : obj exe
exe : $(EXE)
obj : $(OBJ)

MitStyleRemix.o : MitStyleRemix.cc
	$(CC) $(CCFLAGS) -c MitStyleRemix.cc -o MitStyleRemix.o $(LDFLAGS)
MyPdfV2.o : MyPdfV2.cc
	$(CC) $(CCFLAGS) -c MyPdfV2.cc -o MyPdfV2.o $(LDFLAGS)

#%.o : %.cc
#	$(CC) $(CCFLAGS) -c $< -o $@ $(LDFLAGS)

%.exe : %.cc
#	$(CC) $(CCFLAGS) $< -o $@ $(LDFLAGS) MitStyleRemix.o MyPdfV2_cc.so
	$(CC) $(CCFLAGS) -Wl,-rpath=$(PWD) $< -o $@ $(LDFLAGS) MitStyleRemix.o MyPdfV2_cc.so

clean :
	rm -f test.exe MitStyleRemix.o MyPdfV2.o
