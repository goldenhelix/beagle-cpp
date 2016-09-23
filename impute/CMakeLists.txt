# ---< I M P U T E >---

set(IMPUTE_HEADERS
  haplotypepair.h
  samples.h
  markers.h
  vcfemission.h
  iointerface.h
)

set(IMPUTE_SOURCES
  haplotypepair.cpp
  samples.cpp
  markers.cpp
  vcfemission.cpp
  iointerface.cpp
)

find_package(Qt5Core REQUIRED)
find_package(Qt5Xml REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)

set(IMPUTE_ALL
  ${IMPUTE_HEADERS} ${IMPUTE_SOURCES}
)

include_directories(${CMAKE_SOURCE_DIR}) #source root
include_directories(${CMAKE_BINARY_DIR}) #build root
include_directories(${CMAKE_CURRENT_SOURCE_DIR}) #source path
include_directories(${CMAKE_CURRENT_BINARY_DIR}) #build path

include_directories(${SQLITE_INCLUDE_DIR})
include_directories(${3RDPARTY_DIR}/jansson/)
if(WIN32)
  include_directories(${3RDPARTY_DIR}/win32)
endif()

link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
link_directories(${PRECOMPILED_LIB_DIR})

# TODO: needs to be switched to SHARED before we ship, if we directly link. If we only use a driver, than can be kept static.
add_library(impute STATIC ${IMPUTE_ALL})

set_target_properties(impute PROPERTIES COMPILE_FLAGS "${G_CPP_FLAGS}")
set_target_properties(impute PROPERTIES LINK_FLAGS "${G_LINK_FLAGS}")

set_target_properties(impute PROPERTIES AUTOMOC TRUE)

target_link_libraries(impute
  Qt5::Core
  Qt5::Xml
  ghicore
  sqlite
  jansson
)

if(APPLE)
  target_link_libraries(impute "objc" "-framework AppKit")
  target_link_libraries(impute "objc" "-framework CoreFoundation")
  target_link_libraries(impute "objc" "-framework IOKit")
endif(APPLE)

source_group("Headers" FILES ${IMPUTE_HEADERS})
source_group("Sources" FILES ${IMPUTE_SOURCES})

add_subdirectory(tests)
