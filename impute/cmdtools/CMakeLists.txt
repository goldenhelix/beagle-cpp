# ---< I M P U T E _ C M D T O O L S >---
find_package(Qt5Core REQUIRED)
find_package(Qt5Network REQUIRED)
find_package(Qt5Xml REQUIRED) #required by ghicore
set(QT_DONT_USE_QTGUI TRUE)

include_directories(${CMAKE_SOURCE_DIR}) #source root
include_directories(${CMAKE_BINARY_DIR}) #build root
include_directories(${CMAKE_CURRENT_SOURCE_DIR}) #source path
include_directories(${CMAKE_CURRENT_BINARY_DIR}) #build path

link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
link_directories(${PRECOMPILED_LIB_DIR})

set(TOOLS_LIBS Qt5::Core impute )

set(TOOLS_GENERATED )

GHI_ADD_CMD(imputec.cpp
            GENERATED ${TOOLS_GENERATED}
            LIBS ${TOOLS_LIBS})

set(TOOLS_LIBS Qt5::Network Qt5::Core impute )

GHI_ADD_CMD(svs_impute_client.cpp
            GENERATED ${TOOLS_GENERATED}
            LIBS ${TOOLS_LIBS})
