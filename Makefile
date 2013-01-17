##############################################################################
# $Id: Makefile hohejd8 $
##############################################################################
#///
#/// \dir /AppKillSpin

CODE_HOME = ..

SUBDIRECTORIES = Executables Tests

OBJECT_FILES_INTO_LIB    = \
             AdmIntegralS2.o \
             AddStrahlkorperSurfaceAndMesh.o \
             ComputeStrahlkorperWithMesh.o \
             FlatspaceCKV.o \
             AKVsolver.o \
             ComputeAKV.o \
             ComputeS2ConformalFactor.o


##############################################################################
# Nothing below here needs to be changed
##############################################################################
include $(CODE_HOME)/MakefileRules/this_machine.def
include $(CODE_HOME)/MakefileRules/Rules
