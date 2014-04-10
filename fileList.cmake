#To keep the file list clean
set(hdrs_dir ${${PROJECT_NAME}_include_dir})
set(tests_dir ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir})
set(srcs_dir ${CMAKE_CURRENT_SOURCE_DIR}/${source_dir})


#Project source files
set(${PROJECT_NAME}_srcs ${srcs_dir}/tnp.cpp
                        ${srcs_dir}/npnumber.cpp
			${srcs_dir}/multiplication.cpp
			${srcs_dir}/composition.cpp
			${srcs_dir}/polynomial.cpp
			${srcs_dir}/ops.cpp
  )

#Project tests
set(${PROJECT_NAME}_test_sources ${tests_dir}/simpleTest.cpp ${tests_dir}/unaryAlgebraic.cpp)

#Public API
set(${PROJECT_NAME}_headers ${hdrs_dir}/tnp.h
			    ${hdrs_dir}/tnp.hpp
                            ${hdrs_dir}/tnp/ops.h
			    ${hdrs_dir}/tnp/npnumber.hpp
			    ${hdrs_dir}/tnp/polynomial.hpp
			    ${hdrs_dir}/tnp/ops/multiplication.hpp
			    ${hdrs_dir}/tnp/ops/composition.hpp
			    )