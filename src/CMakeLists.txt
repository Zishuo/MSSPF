cmake_minimum_required(VERSION 2.8)
project(MSSPF)

set(CMAKE_CXX_COMPILER /etc/alternatives/g++)
set(executable_output_path ${project_binary_dir}/bin)
set(library_output_path ${project_binary_dir}/lib)

set(CMAKE_CXX_FLAGS "-std=c++0x")
find_package(Boost COMPONENTS graph)

include_directories(./)
aux_source_directory(. src_list)
add_executable(msspf.out ${src_list})

if(UNIX)
	target_link_libraries(msspf.out
		${Boost_LIBRARIES}
		)
elseif(WIN32)
elseif(APPLE)
endif(UNIX)

message(STATUS "c++ flag : "${CMAKE_CXX_FLAGS}) 
message(STATUS "this is source dir "${project_source_dir})
message(STATUS "this is binary dir "${project_binary_dir})
