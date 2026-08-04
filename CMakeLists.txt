cmake_minimum_required(VERSION 3.24)

project(bam_seek VERSION 0.1.0 LANGUAGES CXX)

include(CTest)

option(BAM_SEEK_BUILD_GUI "Build the Qt desktop application" ON)
set(IGVCPP_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../igv-cpp" CACHE PATH
    "Path to the igv-cpp source checkout")

if(NOT TARGET igv::igv)
    if(NOT EXISTS "${IGVCPP_SOURCE_DIR}/CMakeLists.txt")
        message(FATAL_ERROR "igv-cpp was not found. Set IGVCPP_SOURCE_DIR to its source directory.")
    endif()
    set(IGVCPP_BUILD_TESTS OFF CACHE BOOL "Build igv-cpp tests" FORCE)
    add_subdirectory("${IGVCPP_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/igv-cpp")
endif()

add_library(bam_seek_core
    src/evidence_engine.cpp
    src/query_parser.cpp
)
add_library(bamseek::core ALIAS bam_seek_core)
target_compile_features(bam_seek_core PUBLIC cxx_std_20)
target_include_directories(bam_seek_core PUBLIC include)
target_link_libraries(bam_seek_core PUBLIC igv::igv)
target_compile_options(bam_seek_core PRIVATE
    $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall;-Wextra;-Wpedantic;-Wconversion;-Wsign-conversion>
)

if(BAM_SEEK_BUILD_GUI)
    find_package(Qt6 6.5 REQUIRED COMPONENTS Concurrent Widgets)
    qt_standard_project_setup()
    qt_add_executable(bam-seek
        src/main.cpp
        src/main_window.cpp
        include/bamseek/main_window.hpp
    )
    target_link_libraries(bam-seek PRIVATE bamseek::core Qt6::Concurrent Qt6::Widgets)
    target_compile_definitions(bam-seek PRIVATE BAM_SEEK_VERSION="${PROJECT_VERSION}")
    set_target_properties(bam-seek PROPERTIES MACOSX_BUNDLE TRUE WIN32_EXECUTABLE TRUE)
endif()

if(BUILD_TESTING)
    add_executable(bam_seek_tests tests/test_query_parser.cpp)
    target_link_libraries(bam_seek_tests PRIVATE bamseek::core)
    add_test(NAME bam_seek_tests COMMAND bam_seek_tests)

    find_program(SAMTOOLS_EXECUTABLE samtools)
    if(SAMTOOLS_EXECUTABLE)
        set(BAM_SEEK_TEST_DATA "${CMAKE_CURRENT_BINARY_DIR}/generated")
        add_custom_command(
            OUTPUT "${BAM_SEEK_TEST_DATA}/reads.bam" "${BAM_SEEK_TEST_DATA}/reads.bam.bai"
            COMMAND "${CMAKE_COMMAND}" -E make_directory "${BAM_SEEK_TEST_DATA}"
            COMMAND "${SAMTOOLS_EXECUTABLE}" view -b -o "${BAM_SEEK_TEST_DATA}/reads.bam" "${IGVCPP_SOURCE_DIR}/tests/data/reads.sam"
            COMMAND "${SAMTOOLS_EXECUTABLE}" index -o "${BAM_SEEK_TEST_DATA}/reads.bam.bai" "${BAM_SEEK_TEST_DATA}/reads.bam"
            DEPENDS "${IGVCPP_SOURCE_DIR}/tests/data/reads.sam"
        )
        add_custom_target(bam_seek_test_bam DEPENDS "${BAM_SEEK_TEST_DATA}/reads.bam" "${BAM_SEEK_TEST_DATA}/reads.bam.bai")
        add_executable(bam_seek_evidence_tests tests/test_evidence_engine.cpp)
        add_dependencies(bam_seek_evidence_tests bam_seek_test_bam)
        target_compile_definitions(bam_seek_evidence_tests PRIVATE
            BAM_SEEK_TEST_BAM="${BAM_SEEK_TEST_DATA}/reads.bam"
            BAM_SEEK_TEST_REFERENCE="${IGVCPP_SOURCE_DIR}/tests/data/reference.fa"
        )
        target_link_libraries(bam_seek_evidence_tests PRIVATE bamseek::core)
        add_test(NAME bam_seek_evidence_tests COMMAND bam_seek_evidence_tests)
    endif()
endif()
