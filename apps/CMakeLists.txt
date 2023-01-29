
file( GLOB swmp_fortran_app_sources *.f90)

foreach ( filepath ${swmp_fortran_app_sources})
    get_filename_component(app_source ${filepath} NAME)
    string(REPLACE ".f90" "" app_name "${app_source}")
    add_executable(${app_name} ${app_source})
    target_link_libraries(${app_name} swmp)
    target_include_directories(${app_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../modules>)
endforeach()

