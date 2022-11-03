add_subdirectory(common)

add_library(incomp3d STATIC
 boundary_condition.f90
 case_common.f90
 fluid_field.f90
 fractional_step.f90
 fractional_step_implicit.f90
 mesh_reader.f90
 mesh.f90
 rectilinear_mesh.f90
 writeout.f90
 setting_parameter.f90
 simulator.f90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    target_compile_options(incomp3d PRIVATE
        $<$<CONFIG:Debug>:-O0 -warn all -check all -traceback -fpp>
        $<$<CONFIG:Release>:-O3 -fPIC -finline -fpp>
        $<$<CONFIG:RelWithDevInfo>:-O2 -fpp -p>
    )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    target_compile_options(incomp3d PRIVATE
        $<$<CONFIG:Debug>: -O0 -Wall -fcheck=all -fbacktrace -cpp>
        $<$<CONFIG:Release>: -O3 -cpp>
        $<$<CONFIG:RelWithDevInfo>:-O2 -cpp -pg>
    )
endif()

if(OpenMP_Fortran_FOUND)
    target_compile_options(incomp3d PRIVATE 
        "${OpenMP_Fortran_FLAGS}"
    )

    target_link_libraries(incomp3d 
        incomp3d_util
        OpenMP::OpenMP_Fortran 
    )
endif()