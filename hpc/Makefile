VPATH := ../../../core ../../../hpc .
#GPP := g++ -std=c++11 -O3
GPP := g++ -std=c++11 -I../../.. -g -pg -O3 -Wall -Wno-unknown-pragmas -fdiagnostics-color
all: # default target

#==========================================
# Specific targets
#==========================================

Haplo: Haplo.o \
       Analysis.o \
       AncestralChromosomes.o \
       AncestralRecombinationGraph.o \
       Blocks.o \
       BlocksRange.o \
       FileCheckpointer.o \
       FileMemorySaver.o \
       FileOutputs.o \
       FileReader.o \
       FileSNPsWriter.o \
       FileWriter.o \
       GraphDataBinned.o \
       GraphDataBinned2D.o \
       HaploBlockBoundaries.o \
       HaploBlockVariants.o \
       MutationGraph.o \
       ParentPairUrn.o \
       Parser.o \
       PolyaUrn.o \
       PopulationStructureHistory.o \

EXECS := Haplo

#==========================================
# Generic targets
#==========================================

-include *.d

%.o: %.cpp
	@ echo making $*.o 
	@ $(GPP) -MMD -c $< -o $@ 

$(EXECS): 
	@ echo making $@
	@ $(GPP) $^ -o $@ 
	@ chmod +x $@
	@ echo


#==========================================
# Declare 'phony' (non-file) targets
#==========================================
.PHONY: clean all

clean: 
	@ rm -f *.o *.d *.exe*

all: $(EXECS)


#==========================================

