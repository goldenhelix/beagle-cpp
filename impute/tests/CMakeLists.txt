# ---< I M P U T E _ T E S T S >---

find_package(Qt5Core REQUIRED)
find_package(Qt5Test REQUIRED)
find_package(Qt5Xml)
set(QT_DONT_USE_QTGUI TRUE)

add_definitions(-DBASE_PATH="${CMAKE_CURRENT_SOURCE_DIR}")
add_definitions(-DSCRATCH_PATH="${CMAKE_BINARY_DIR}/test_scratch")

include_directories(${CMAKE_SOURCE_DIR}) #source root
include_directories(${CMAKE_BINARY_DIR}) #build root
include_directories(${CMAKE_CURRENT_SOURCE_DIR}) #source path
include_directories(${CMAKE_CURRENT_BINARY_DIR}) #build path

link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
link_directories(${PRECOMPILED_LIB_DIR})

set(TEST_LIBS Qt5::Core Qt5::Test impute)

set(EMPTY_TESTS_GENERATED
)

GHI_ADD_TEST(test_impute_datastructures.cpp
             GENERATED ${EMPTY_TESTS_GENERATED}
             LIBS ${TEST_LIBS})


