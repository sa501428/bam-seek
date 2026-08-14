cmake_minimum_required(VERSION 3.24)

project(bam_seek VERSION 0.1.0 LANGUAGES CXX)

include(CTest)

option(BAM_SEEK_BUILD_GUI "Build the Qt desktop application" ON)
option(BAM_SEEK_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
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
    src/comparison.cpp
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
if(BAM_SEEK_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
    target_compile_options(bam_seek_core PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(bam_seek_core PUBLIC -fsanitize=address,undefined)
endif()

if(BAM_SEEK_BUILD_GUI)
    find_package(Qt6 6.5 REQUIRED COMPONENTS Concurrent Network Widgets)
    qt_standard_project_setup()
    qt_add_executable(bam-seek
        src/main.cpp
        src/igv_command_receiver.cpp
        src/main_window.cpp
        include/bamseek/igv_command_receiver.hpp
        include/bamseek/main_window.hpp
    )
    target_link_libraries(bam-seek PRIVATE bamseek::core Qt6::Concurrent Qt6::Network Qt6::Widgets)
    target_compile_definitions(bam-seek PRIVATE BAM_SEEK_VERSION="${PROJECT_VERSION}")
    set_target_properties(bam-seek PROPERTIES MACOSX_BUNDLE TRUE WIN32_EXECUTABLE TRUE)

    if(BUILD_TESTING)
        add_executable(bam_seek_receiver_tests
            tests/test_igv_command_receiver.cpp
            src/igv_command_receiver.cpp
            include/bamseek/igv_command_receiver.hpp
        )
        target_include_directories(bam_seek_receiver_tests PRIVATE include)
        target_link_libraries(bam_seek_receiver_tests PRIVATE Qt6::Core Qt6::Network)
        add_test(NAME bam_seek_receiver_tests COMMAND bam_seek_receiver_tests)

        add_executable(bam_seek_pileup_tests
            tests/test_pileup_view.cpp
            src/pileup_view.cpp
            include/bamseek/pileup_view.hpp
        )
        target_compile_definitions(bam_seek_pileup_tests PRIVATE
            BAM_SEEK_JAK2_PILEUP_BAM="${CMAKE_CURRENT_SOURCE_DIR}/tests/data/jak2/jak2_v617f_25pct_150x.bam")
        target_link_libraries(bam_seek_pileup_tests PRIVATE bamseek::core Qt6::Widgets)
        add_test(NAME bam_seek_pileup_tests COMMAND bam_seek_pileup_tests)
        set_tests_properties(bam_seek_pileup_tests PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen")
    endif()
endif()

if(BUILD_TESTING)
    add_executable(bam_seek_comparison_tests tests/test_comparison.cpp)
    target_link_libraries(bam_seek_comparison_tests PRIVATE bamseek::core)
    add_test(NAME bam_seek_comparison_tests COMMAND bam_seek_comparison_tests)

    add_executable(bam_seek_tests tests/test_query_parser.cpp)
    target_link_libraries(bam_seek_tests PRIVATE bamseek::core)
    target_compile_definitions(bam_seek_tests PRIVATE
        BAM_SEEK_CLINICAL_MAPPING="${CMAKE_CURRENT_SOURCE_DIR}/tests/data/clinical_mappings.tsv")
    add_test(NAME bam_seek_tests COMMAND bam_seek_tests)

    find_program(SAMTOOLS_EXECUTABLE samtools)
    if(SAMTOOLS_EXECUTABLE)
        set(BAM_SEEK_TEST_DATA "${CMAKE_CURRENT_BINARY_DIR}/generated")
        add_custom_command(
            OUTPUT "${BAM_SEEK_TEST_DATA}/reads.bam" "${BAM_SEEK_TEST_DATA}/reads.bam.bai"
                   "${BAM_SEEK_TEST_DATA}/evidence.bam" "${BAM_SEEK_TEST_DATA}/evidence.bam.bai"
                   "${BAM_SEEK_TEST_DATA}/evidence_reference.fa" "${BAM_SEEK_TEST_DATA}/evidence_reference.fa.fai"
            COMMAND "${CMAKE_COMMAND}" -E make_directory "${BAM_SEEK_TEST_DATA}"
            COMMAND "${SAMTOOLS_EXECUTABLE}" view -b -o "${BAM_SEEK_TEST_DATA}/reads.bam" "${IGVCPP_SOURCE_DIR}/tests/data/reads.sam"
            COMMAND "${SAMTOOLS_EXECUTABLE}" index -o "${BAM_SEEK_TEST_DATA}/reads.bam.bai" "${BAM_SEEK_TEST_DATA}/reads.bam"
            COMMAND "${SAMTOOLS_EXECUTABLE}" view -b -o "${BAM_SEEK_TEST_DATA}/evidence.bam" "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/evidence.sam"
            COMMAND "${SAMTOOLS_EXECUTABLE}" index -o "${BAM_SEEK_TEST_DATA}/evidence.bam.bai" "${BAM_SEEK_TEST_DATA}/evidence.bam"
            COMMAND "${CMAKE_COMMAND}" -E copy "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/evidence_reference.fa" "${BAM_SEEK_TEST_DATA}/evidence_reference.fa"
            COMMAND "${SAMTOOLS_EXECUTABLE}" faidx "${BAM_SEEK_TEST_DATA}/evidence_reference.fa"
            DEPENDS "${IGVCPP_SOURCE_DIR}/tests/data/reads.sam" tests/data/evidence.sam tests/data/evidence_reference.fa
        )
        add_custom_target(bam_seek_test_bam DEPENDS
            "${BAM_SEEK_TEST_DATA}/reads.bam" "${BAM_SEEK_TEST_DATA}/reads.bam.bai"
            "${BAM_SEEK_TEST_DATA}/evidence.bam" "${BAM_SEEK_TEST_DATA}/evidence.bam.bai"
            "${BAM_SEEK_TEST_DATA}/evidence_reference.fa" "${BAM_SEEK_TEST_DATA}/evidence_reference.fa.fai")
        add_executable(bam_seek_evidence_tests tests/test_evidence_engine.cpp)
        add_dependencies(bam_seek_evidence_tests bam_seek_test_bam)
        target_compile_definitions(bam_seek_evidence_tests PRIVATE
            BAM_SEEK_TEST_BAM="${BAM_SEEK_TEST_DATA}/reads.bam"
            BAM_SEEK_TEST_REFERENCE="${IGVCPP_SOURCE_DIR}/tests/data/reference.fa"
            BAM_SEEK_EVIDENCE_BAM="${BAM_SEEK_TEST_DATA}/evidence.bam"
            BAM_SEEK_EVIDENCE_REFERENCE="${BAM_SEEK_TEST_DATA}/evidence_reference.fa"
            BAM_SEEK_JAK2_NEGATIVE_BAM="${CMAKE_CURRENT_SOURCE_DIR}/tests/data/jak2/jak2_negative_120x.bam"
            BAM_SEEK_JAK2_25PCT_BAM="${CMAKE_CURRENT_SOURCE_DIR}/tests/data/jak2/jak2_v617f_25pct_150x.bam"
            BAM_SEEK_JAK2_28PCT_BAM="${CMAKE_CURRENT_SOURCE_DIR}/tests/data/jak2/jak2_v617f_28pct_180x.bam"
        )
        target_link_libraries(bam_seek_evidence_tests PRIVATE bamseek::core)
        add_test(NAME bam_seek_evidence_tests COMMAND bam_seek_evidence_tests)
    endif()
endif()
