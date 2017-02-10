################################################################################
# Package: Rich1MDCS
################################################################################
gaudi_subdir(Rich1MDCS v1r0)

gaudi_depends_on_subdirs(Rich/RichRecBase
				AIDA
				ROOT
				COOL
				GaudiUtils)

find_package(Boost)
find_package( ROOT COMPONENTS Hist Graf Gpad RooFit RooFitCore Foam Minuit2 Tree)
find_package(AIDA)
find_package(COOL)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

gaudi_add_module(Rich1MDCS
                   src/*.cpp
                   INCLUDE_DIRS Boost ROOT RichRecBase
                   LINK_LIBRARIES Boost ROOT RichRecBase)

gaudi_add_library(libRich1MDCS
                   src/*.cpp
                   NO_PUBLIC_HEADERS
                   INCLUDE_DIRS Boost ROOT RichRecBase
                   LINK_LIBRARIES Boost ROOT RichRecBase)

gaudi_add_executable(PatternCleaning scripts/CleanRun.C
                        LINK_LIBRARIES libRich1MDCS)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable -Wno-sign-compare" ) 
