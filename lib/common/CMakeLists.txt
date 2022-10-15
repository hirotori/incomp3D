add_library(incomp3d_util STATIC
 command_line_args.f90
 floating_point_parameter.f90
 IO_operator.f90
 linear_algebra.f90
 system_operator.f90
 time_monitor.f90
 abstract_legacy_vtk_writer.f90
 vtk_field_data.f90
 vtk_structured_points_writer.f90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    target_compile_options(incomp3d_util PRIVATE
        $<$<CONFIG:Debug>:-O0 -warn all -check all -traceback>
        $<$<CONFIG:Release>:-O2 -fPIC>
    )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    target_compile_options(incomp3d_util PRIVATE
        $<$<CONFIG:Debug>: -O0 -Wall -fcheck=all -fbacktrace>
        $<$<CONFIG:Release>: -O2>
    )
endif()