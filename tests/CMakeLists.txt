add_executable(exsqs_tests
  test_geometry.cpp
  test_correlation.cpp
  test_rng.cpp
  test_symmetry.cpp
  test_displacements.cpp
  test_config.cpp
  test_evolution.cpp
  test_e2e.cpp
  test_parallel.cpp
  test_checkpoint.cpp)
target_link_libraries(exsqs_tests PRIVATE exsqs_core Catch2::Catch2WithMain)
target_compile_definitions(exsqs_tests PRIVATE
  EXSQS_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data")
target_compile_options(exsqs_tests PRIVATE -Wall -Wextra)
add_test(NAME exsqs_tests COMMAND exsqs_tests ~[e2e])
add_test(NAME exsqs_e2e COMMAND exsqs_tests [e2e])
