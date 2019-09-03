cmake_minimum_required(VERSION 3.13.4)

project(SlicerLineMarkerRegistration)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/SlicerLineMarkerRegistration")
set(EXTENSION_CATEGORY "IGT")
set(EXTENSION_CONTRIBUTORS "Junichi Tokuda (BWH), Eleni Siampli (Children's National), Steve Pieper (Isomics, Inc.), ")
set(EXTENSION_DESCRIPTION "Includes Line Marker Registration algorithm for registration of IGT instruments imaged in intraprocedural scans to pre-procedure scans.")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/SlicerLineMarkerRegistration.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/SlicerLineMarkerRegistration/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(LineMarkerRegistration)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
