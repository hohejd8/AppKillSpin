##############################################################################
# $Id: Makefile hohejd8 $
##############################################################################
#///
#/// \dir /AppKillSpin

CODE_HOME = ..

SUBDIRECTORIES = Executables Tests

OBJECT_FILES_INTO_LIB    = \
             AdmIntegralS2.o \
             AddStrahlkorperSurface.o \
             ComputeSKWM.o \
             FlatspaceCKV.o \
             AKVsolver.o \
             ComputeAKV.o


##############################################################################
# Nothing below here needs to be changed
##############################################################################
include $(CODE_HOME)/MakefileRules/this_machine.def
include $(CODE_HOME)/MakefileRules/Rules
