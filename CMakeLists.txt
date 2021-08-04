cmake_minimum_required(VERSION 3.0)
project(PadicODE VERSION 0.1 LANGUAGES C)

set(DEPS gmp flint)

foreach (LIB ${DEPS})
	find_library(${LIB}_LIBRARY NAMES ${LIB})
	if (NOT ${LIB}_LIBRARY)
		message(FATAL_ERROR "~~ Library ${LIB} was not found.")
	endif()
	add_library(${LIB} UNKNOWN IMPORTED)
	set_property(TARGET ${LIB} PROPERTY IMPORTED_LOCATION ${${LIB}_LIBRARY})
	message("-- Library ${LIB} was found in ${${LIB}_LIBRARY}")
endforeach ()

# Include tests, which have their own CMakeLists in a subdirectory
enable_testing()
add_subdirectory(tests)

set(padicSrc
	src/padic_ode.c
	src/solver_fuchs.c
	src/solver_frobenius.c
)

set(padicHeaders
	src/padic_ode.h
	src/solver_fuchs.h
	src/solver_frobenius.h
)

add_library(padicODE SHARED ${padicSrc})
target_link_libraries(padicODE ${DEPS})
target_compile_options(padicODE PRIVATE -Wall -Wextra -pedantic)

set_target_properties(padicODE PROPERTIES PUBLIC_HEADER "${CascadeHeaders}")
set_target_properties(padicODE PROPERTIES VERSION ${PROJECT_VERSION})

install(TARGETS padicODE LIBRARY PUBLIC_HEADER)
