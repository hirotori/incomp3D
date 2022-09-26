add_subdirectory(common)

add_library(incomp3d STATIC
 boundary_condition.f90
 case_common.f90
 fluid_field.f90
 fractional_step.f90
 fractional_step_implicit.f90
 mesh.f90
 writeout.f90
 setting_parameter.f90
 simulator.f90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    target_compile_options(incomp3d PRIVATE
        $<$<CONFIG:Debug>:-O0 -warn all -check all -traceback>
        $<$<CONFIG:Release>:-O2 -fPIC>
    )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    target_compile_options(incomp3d PRIVATE
        $<$<CONFIG:Debug>: -O0 -Wall -fcheck=all -fbacktrace>
        $<$<CONFIG:Release>: -O2>
    )
endif()

target_link_libraries(incomp3d incomp3d_util)