set(Tests
	indicial_polynomial
	solution_extend
	solution_update
	singleton_frobenius
	frobenius
)

include_directories(${CMAKE_SOURCE_DIR}/src)

foreach (test_file IN LISTS Tests)
	add_executable(${test_file}_target ${test_file}.c)
	add_dependencies(${test_file}_target implode)
	target_link_libraries(${test_file}_target implode flint gmp)
	add_test(NAME ${test_file}_test COMMAND ${test_file}_target)
endforeach()
