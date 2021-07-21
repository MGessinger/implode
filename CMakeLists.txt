cmake_minimum_required(VERSION 3.0)
project(PadicODE VERSION 1.0)

# Include tests, which have their own CMakeLists in a subdirectory
#enable_testing()
#add_subdirectory(tests)

set(padicSrc
	src/padic_diffeq.c
	src/solver_fuchs.c
	src/solver_frobenius.c
)

add_library(padicODE SHARED ${padicSrc})

target_link_libraries(padicODE flint)

target_compile_options(padicODE PRIVATE -Wall -Wextra -pedantic)

add_executable(padicMain src/main.c)
target_link_libraries(padicMain padicODE)
