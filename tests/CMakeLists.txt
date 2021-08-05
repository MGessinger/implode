set(Tests
	indicial_polynomial
	indicial_polynomial_eval
	power_series
	solution_extend
)

include_directories(${CMAKE_SOURCE_DIR}/src)

foreach (test_file IN LISTS Tests)
	add_executable(${test_file}_target ${test_file}.c)
	add_dependencies(${test_file}_target padicODE)
	target_link_libraries(${test_file}_target padicODE flint)
	add_test(NAME ${test_file}_test COMMAND ${test_file}_target)
endforeach()
