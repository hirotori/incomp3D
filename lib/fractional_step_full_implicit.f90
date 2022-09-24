module fractional_step_full_implicit_m
    use floating_point_parameter_m
    use fluid_field_m, only : fluid_field_t
    use mesh_m, only : rectilinear_mesh_t
    use setting_parameter_m, only : slv_setting, case_setting
    use fractional_step_m, only : mat_a, calc_convective_and_diffusive_flux, slvr_init_common
    use fractional_step_implicit_m, only : solver_fs_imp_t, calc_pseudo_velocity_common_core

    implicit none
    private

    type,extends(solver_fs_imp_t) :: solver_fs_full_imp_t
        !!ver 1.0.0
        !!移流項と拡散項をクランク･ニコルソン法で離散化する場合.
        
        contains
        procedure :: init => init_solver_3

    end type
    
contains
subroutine init_solver_3(this, fld, grd, settings_slv, setting_case)
    !!ソルバを初期化する.
    class(solver_fs_full_imp_t),intent(inout) :: this
    type(fluid_field_t),intent(in) :: fld
    type(rectilinear_mesh_t),intent(in) :: grd
    type(slv_setting),intent(in) :: settings_slv
    type(case_setting),intent(in) :: setting_case

    integer(ip) mx(3)

    mx(:) = grd%get_extents()

    call slvr_init_common(this, fld, grd, settings_slv, setting_case)
    
    call set_matrix_(this%matrix_v, grd%ds, grd%dv, grd%dx, setting_case%reynolds_number, setting_case%dt)
    
    allocate(this%rhs_(1:3,2:mx(1),2:mx(2),2:mx(3)))

end subroutine

subroutine set_matrix_(mat)
    type(mat_a),intent(inout) :: mat

end subroutine
end module fractional_step_full_implicit_m