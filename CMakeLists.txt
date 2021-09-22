cmake_minimum_required(VERSION 3.0)
project(Implode VERSION 0.1 LANGUAGES C)

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

set(ExplodeSrc
	src/padic_ode.c
	src/padic_ode_solution.c
	src/fuchs_solver.c
	src/frobenius_solver.c
)

set(ExplodeHeaders
	src/padic_ode.h
	src/implode.h
)

add_library(implode SHARED ${ExplodeSrc})
target_link_libraries(implode ${DEPS})

if (MSVC)
	add_compile_options(/W4 /WX)
else()
	add_compile_options(-Wall -Wextra -Wpedantic)
endif()

set_target_properties(implode PROPERTIES PUBLIC_HEADER "${ExplodeHeaders}")
set_target_properties(implode PROPERTIES VERSION ${PROJECT_VERSION})

install(TARGETS implode LIBRARY PUBLIC_HEADER)
