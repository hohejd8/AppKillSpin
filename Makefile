##############################################################################
# $Id: Makefile hohejd8 $
##############################################################################
#///
#/// \dir /hohejd8
#/// Contains code written by hohejd8

CODE_HOME = ..

SUBDIRECTORIES = Executables #SpherePackhohejd8

OBJECT_FILES_INTO_LIB    = \
             AdmIntegralS2.o \
             AddStrahlkorperSurface.o \
             ComputeSKWM.o \
             FlatspaceCKV.o \
             AKVsolver.o \
             ComputeAKV.o
             #SurfaceBasisExt.o \


##############################################################################
# Nothing below here needs to be changed
##############################################################################
include $(CODE_HOME)/MakefileRules/this_machine.def
include $(CODE_HOME)/MakefileRules/Rules
