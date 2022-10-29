module mesh_reader_m
    use floating_point_parameter_m
    use IO_operator_m
    implicit none
    character(*),parameter,public :: mesh_form_equil = "equil-spaced"
    character(*),parameter,public :: mesh_form_rectil = "rectilinear"
contains

subroutine get_mesh_from_file(path, imx, jmx, kmx, length, r)
    character(*),intent(in) :: path
        !!ファイルパス.
    integer(ip),intent(out) :: imx, jmx, kmx
        !!格子点数.
    real(dp),intent(inout) :: length(3)
        !!寸法.
    real(dp),allocatable,intent(inout) :: r(:,:,:,:)
        !!節点座標.

    if ( file_exists(path) ) then

        select case(get_extention(path))

        case(".msh")
            call original_format_reader(path, imx, jmx, kmx, length, r)
        
        case(".vtk")
            error stop "vtk format is unsupported."
        case default
            error stop 'unknown format ""'//path//'""'
        end select
    
    else
        error stop "Mesh file does not exist in current path."
    end if

end subroutine

subroutine original_format_reader(path, imx, jmx, kmx, lengths, r)
    !!独自規格のファイルフォーマットリーダー. 拡張子はmsh, プレーンテキストのみ.
    character(*),intent(in) :: path
        !!ファイルパス    
    integer,intent(inout) :: imx, jmx, kmx
    real(dp),intent(inout) :: lengths(3)
    real(dp),allocatable,intent(inout) :: r(:,:,:,:)

    character(128) format_, comment_
    integer unit, iostat
    character(128) iomsg

    integer(ip) i, j, k

    open(newunit=unit, file=path, status="old", form="formatted")
        print "('Loading Mesh file ""', A, '""...')", path
        read(unit,"(A)",iostat=iostat,iomsg=iomsg) format_
        read(unit,"(A)",iostat=iostat,iomsg=iomsg) comment_
        read(unit,*,iostat=iostat,iomsg=iomsg) imx, jmx, kmx
        print "('Extents :: ', 3(i0,1x,:,', '))", imx, jmx, kmx
        !メッシュフォーマットによって場合分けする. 
        select case(trim(adjustl(format_)))

            case(mesh_form_equil)
                !各方向の寸法だけ必要とする.
                read(unit,*,iostat=iostat,iomsg=iomsg) lengths(1:3)

            case(mesh_form_rectil)
                !座標点.
                allocate(r(3,imx,jmx,kmx))
                do k = 1, kmx
                do j = 1, jmx
                do i = 1, imx
                    read(unit,*,iostat=iostat,iomsg=iomsg) r(1:3,i,j,k)
                end do
                end do    
                end do

        case default
            print "(A)","Wrong File format. format should be """//mesh_form_equil//""" or """//mesh_form_rectil//""". "
            print "('For further description, please see README.')"
            error stop "Error in mesh_reader.f90"
        end select
        
end subroutine

function get_extention(fname) result(extention)
    character(*),intent(in) :: fname
    character(:),allocatable :: extention

    integer loc_

    loc_ = index(fname, ".", back=.true.)
    if ( loc_ == 0 ) then
        print "('Warning :: passed argument ""',A,'"" might not have extention.')", fname
        extention = ""
    else 
        extention = fname(loc_:)
    end if
end function

end module mesh_reader_m