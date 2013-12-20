#To keep the file list clean
set(hdrs_dir ${${PROJECT_NAME}_include_dir})
set(tests_dir ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir})
set(srcs_dir ${CMAKE_CURRENT_SOURCE_DIR}/${source_dir})

#Project header files
set(hdrs ${hdrs_dir}/tnp.hpp)

#Project source files
set(srcs ${srcs_dir}/tnp.cpp
         ${srcs_dir}/npnumber.cpp
	 ${srcs_dir}/multiplication.cpp
	 ${srcs_dir}/composition.cpp
	 ${srcs_dir}/polynomial.cpp
  )

#Project tests
set(test_sources ${tests_dir}/simpleTest.cpp ${tests_dir}/unaryAlgebraic.cpp)

