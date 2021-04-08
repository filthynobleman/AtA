cmake_minimum_required(VERSION 3.16.3 FATAL_ERROR)
project(AtA_Shared VERSION 1.0 LANGUAGES C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

add_compile_options(-O3)


option(USE_DOUBLES "Use double precision for floating point values." ON)
option(PROCESS_SYNCHRONIZATION "Synchronize processes after the computation." OFF)
option(COMPUTE_ONLY "Stop the execution after the computation at leaf level." OFF)


if(USE_DOUBLES)
    set(DOUBLE_PRECISION "ON")
    message("Compiling with double precision floating point.")
else()
    message("Compiling with single precision floating point.")
endif()

if(PROCESS_SYNCHRONIZATION)
    set(SYNCHRONIZE_PROCESSES "ON")
    message("Processes will synchronize after the leaf level computation.")
endif()

if(COMPUTE_ONLY)
    set(ONLY_COMPUTE "ON")
    message("Algorithm will sop after the leaf level computation.")
endif()


add_executable(TestAtADistributed test_ata.c ATA.c matrices.c simulexec.c queue.c)

find_package(MPI REQUIRED)
if(MPI_FOUND)
    include_directories(${MPI_INCLUDE_DIRS})
    link_directories(${MPI_LIBRARIES})
else()
    message(FATAL_ERROR "MPI not found.")
endif()

find_package(MKL REQUIRED)
if (MKL_FOUND)
    set(USE_MKL "ON")
    include_directories(${MKL_INCLUDE_DIRS})
    link_directories(${MKL_LIBRARIES})
    find_package(TBB REQUIRED)
    if (TBB_FOUND)
        message("Found TBB")
        include_directories(${TBB_INCLUDE_DIRS})
        link_directories(${TBB_LIBRARIES})
        target_link_libraries(TestAtADistributed mkl_intel_ilp64 mkl_core mkl_tbb_thread pthread m ${TBB_IMPORTED_TARGETS} ${MPI_LIBRARIES})
    else()
        message(FATAL_ERROR "TBB not found.")
    endif()
else()
    message(FATAL_ERROR "MKL not found")
endif()

configure_file(common.h.in ${CMAKE_SOURCE_DIR}/common.h)