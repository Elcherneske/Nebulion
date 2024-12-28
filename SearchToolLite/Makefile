include ../Makefile.common

# Need to handle whether $(MSTOOLKIT) is defined as absolute or relative path
# from parent directory and set MSTPATH accordingly
IS_ROOT    := $(if $(patsubst /%,,$(MSTOOLKIT)),,yes)
IS_HOME    := $(if $(patsubst ~%,,$(MSTOOLKIT)),,yes)
IS_NETWORK := $(if $(patsubst \\\\%,,$(MSTOOLKIT)),,yes)
IS_DRIVE   := $(foreach d,A B C D E F G H I J K L M N O P Q R S T U V Q X Y Z,$(if $(patsubst $(d):/%,,$(MSTOOLKIT)),,yes))
ifeq ($(strip $(IS_ROOT)$(IS_HOME)$(IS_NETWORK)$(IS_DRIVE)),yes)
MSTPATH = $(MSTOOLKIT)
else
MSTPATH = ../$(MSTOOLKIT)
endif


UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
   override CXXFLAGS += -g -static -std=c++14            -fpermissive -Wall -Wextra -Wno-write-strings -DGITHUBSHA='"$(GITHUB_SHA)"' -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -D_NOSQLITE -D__int64=off64_t -I. -I$(MSTPATH)/include -I$(MSTPATH)/src/expat-2.2.9/lib -I$(MSTPATH)/src/zlib-1.2.11
else
   override CXXFLAGS += -g -static -std=c++14 -fconcepts -fpermissive -Wall -Wextra -Wno-write-strings -DGITHUBSHA='"$(GITHUB_SHA)"' -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -D_NOSQLITE -D__int64=off64_t -I. -I$(MSTPATH)/include -I$(MSTPATH)/src/expat-2.2.9/lib -I$(MSTPATH)/src/zlib-1.2.11


endif

# COMETSEARCH = Threading.o CometInterfaces.o CometSearch.o CometPreprocess.o CometPostAnalysis.o CometMassSpecUtils.o CometWriteOut.o\
# 				CometWriteSqt.o CometWritePepXML.o CometWriteMzIdentML.o CometWritePercolator.o CometWriteTxt.o CometSearchManager.o\
# 				CombinatoricsUtils.o ModificationsPermuter.o CometFragmentIndex.o

COMETSEARCH = Threading.o CometInterfaces.o CometMassUtils.o CometSearchManager.o CometFragmentIndex.o CometCalculateUtils.o CometParam.o\
				CometSpectrumPreprocess.o CometSearch.o CometModificationUtils.o CombinatoricsUtils.o CometDecoyUtils.o CometWriteOut.o


all:  $(COMETSEARCH)
	ar rcs libcometsearch.a $(COMETSEARCH)

clean:
	rm -f *.o libcometsearch.a

oclean:
	rm -f *.o

CometSpectrumPreprocess.o: 		CometSpectrumPreprocess.cpp CometSpectrumPreprocess.h Common.h CometData.h CometParam.h $(MSTPATH)
								${CXX} ${CXXFLAGS} CometSpectrumPreprocess.cpp -c
CometMassUtils.o:				CometMassUtils.cpp CometMassUtils.h CometParam.h CometData.h Common.h
								${CXX} ${CXXFLAGS} CometMassUtils.cpp -c
Threading.o:          			Threading.cpp Threading.h
								${CXX} ${CXXFLAGS} Threading.cpp -c
CometSearch.o:        			CometSearch.cpp CometSearch.h Common.h CometData.h CometParam.h CometCalculateUtils.h CometFragmentIndex.h 
								${CXX} ${CXXFLAGS} CometSearch.cpp -c
CometCalculateUtils.o:			CometCalculateUtils.cpp CometCalculateUtils.h Common.h CometData.h CometParam.h
								${CXX} ${CXXFLAGS} CometCalculateUtils.cpp -c
CometParam.o: 					CometParam.cpp CometParam.h Common.h CometData.h Threading.h
								${CXX} ${CXXFLAGS} CometParam.cpp -c
CometFragmentIndex.o: 			CometFragmentIndex.cpp CometFragmentIndex.h Common.h CometData.h CometParam.h CometSearch.h ThreadPool.h CometStatus.h CometMassUtils.h
								${CXX} ${CXXFLAGS} CometFragmentIndex.cpp -c
CometSearchManager.o:     		CometSearchManager.cpp CometSearchManager.h Common.h CometData.h CometParam.h CometInterfaces.h CometSpectrumPreprocess.h CometStatus.h CometMassUtils.h CometFragmentIndex.h CometDecoyUtils.h CometWriteOut.h
								${CXX} ${CXXFLAGS} CometSearchManager.cpp -c
CometInterfaces.o:      		CometInterfaces.cpp CometInterfaces.h Common.h CometData.h CometSearchManager.h
								${CXX} ${CXXFLAGS} CometInterfaces.cpp -c
CometModificationUtils.o:		CometModificationUtils.cpp CometModificationUtils.h Common.h CometData.h CometParam.h CombinatoricsUtils.h
								${CXX} ${CXXFLAGS} CometModificationUtils.cpp -c
CombinatoricsUtils.o:			CombinatoricsUtils.cpp CombinatoricsUtils.h Common.h
								${CXX} ${CXXFLAGS} CombinatoricsUtils.cpp -c
CometDecoyUtils.o:				CometDecoyUtils.cpp CometDecoyUtils.h Common.h CometData.h CometParam.h
								${CXX} ${CXXFLAGS} CometDecoyUtils.cpp -c
CometWriteOut.o:				CometWriteOut.cpp CometWriteOut.h Common.h CometData.h CometParam.h
								${CXX} ${CXXFLAGS} CometWriteOut.cpp -c
