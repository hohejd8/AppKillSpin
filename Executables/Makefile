##############################################################################
# $Id: Makefile hohejd8 $
##############################################################################
#
#/// \dir hohejd8
#/// General operators for initial value problem for black holes
#/// \ingroup Elliptic

CODE_HOME = ../..

SUBDIRECTORIES = 

LINKED_LIBS = \
	InitialValueProblem \
	EllipticSolver Preconditioning EllipticResidual \
        VolumeTerms BoundaryConditions \
	LinearizedGravity \
        ComputeItems GrItems BackgroundMetrics Observers\
        EinsteinSolutions TensorFunctions CoordinateTransformations\
	StrahlkorperItems StrahlkorperResidual StrahlkorperDataSupplier\
        StrahlkorperFinderAlgorithm Strahlkorper \
        SurfaceFinderItems SurfaceResidual SurfaceDataSupplier Surface \
        MinimizationAlgorithm  \
	TensorFunctions \
	SpatialCoordMaps \
        SpectralFilters SpectralFunctionals BasisFunctions\
        TensorSphericalHarmonics DfftPack SpherePack \
        Mover Interpolators Functionals \
        Domain Slice Subdomain CoordMap Topology DataMesh IndexMap \
        IO DataBox Vars StringParsing Tensor Math \
        MiscUtils MpiWrappers MyContainers LowLevelUtils ErrorHandling \
        hohejd8 #SpherePackhohejd8

EXECUTABLES = IVP_hohejd8 #TestSurfaceBasisExt #TestVectorHA  TestGslMultiroots

##############################################################################
# Nothing below here needs to be changed
##############################################################################
include $(CODE_HOME)/MakefileRules/this_machine.def
include $(CODE_HOME)/MakefileRules/Rules
