module vtk_field_data_m
    !/////////////////////////////////////////////////////////////////
    !! version: 1.0.0
    !! author: T.Ikeda
    !! summary:
    !! - Legacy format におけるattributes dataの書き出しに関するモジュール
    !!
    !! @note
    !! 
    !! - これ単体では動作しない. legacy formatのフィールドデータ書き出しのために用いられる.
    !! - 現状はデータを動的割り付け配列により記憶しているため, アプリケーションによってはメモリが足りなくなる
    !!   可能性がある. 
    !! 
    !! @endnote
    !/////////////////////////////////////////////////////////////////
    use floating_point_parameter_m
    implicit none
    private
    type attrib_data_holder_t
        !!ファイルに登録するフィールドデータの構造体.
        character(:),allocatable :: name
        real(dp),allocatable :: val_r1(:) 
        real(dp),allocatable :: val_r2(:,:) 
        ! real(dp),pointer :: val_r1_ptr(:) => null() 
        ! real(dp),pointer :: val_r2_ptr(:,:) => null()
        contains
        procedure,private :: register_r2, convert_r4_to_r2_and_register_
        generic,public :: register_vector => register_r2, convert_r4_to_r2_and_register_
        procedure,private :: register_r1, convert_r3_to_r1_and_register_
        generic,public :: register_scalar => register_r1, convert_r3_to_r1_and_register_
        procedure,public :: write_data => write_data_ascii
        procedure,public :: write_data_binary
    
    end type
    
    character(*),parameter,private :: LOOKUP_TABLE_NAME = "default"
    
    public attrib_data_holder_t

contains
subroutine register_r2(holder, v_r2, lbounds, ubounds)
    class(attrib_data_holder_t),intent(inout) :: holder
    real(dp),intent(in) :: v_r2(:,:)
    integer(ip),intent(in) :: lbounds
        !!データを格納する配列のインデックスの下限.
    integer(ip),intent(in) :: ubounds
        !!データを格納する配列のインデックスの上限.

    integer(ip) ist, imx, icmx
    integer(ip) astat
    character(256) errmsg

    ist = lbounds
    imx = ubounds
    icmx = imx - ist + 1

    allocate(holder%val_r2(3,icmx), source = v_r2(1:3,ist:imx))

end subroutine

subroutine convert_r4_to_r2_and_register_(holder, v_r4, lbounds, ubounds)
    class(attrib_data_holder_t),intent(inout) :: holder
    real(dp),intent(in) :: v_r4(:,:,:,:)
        !!ベクトルデータ. rank1: coordinates(x,y,z), rank2~4:index(i,j,k)
        !!@note 配列全体を指定する. 有効範囲はlbounds, uboundsで指定する.
    integer(ip),intent(in) :: lbounds(3)
        !!データを格納する配列のインデックスの下限.
    integer(ip),intent(in) :: ubounds(3)
        !!データを格納する配列のインデックスの上限.

    integer(ip) ist, jst, kst, imx, jmx, kmx, icmx
    integer(ip) astat
    character(256) errmsg

    ist = lbounds(1)
    jst = lbounds(2)
    kst = lbounds(3)
    imx = ubounds(1)
    jmx = ubounds(2)
    kmx = ubounds(3)
    icmx = product(ubounds(:) - lbounds(:) + 1)

    allocate(holder%val_r2(3,icmx), source = reshape(v_r4(1:3,ist:imx,jst:jmx,kst:kmx), shape=[3,icmx]), &
            stat = astat, errmsg = errmsg)

    if ( astat /= 0 ) error stop trim(errmsg) 

end subroutine

subroutine register_r1(holder, v_r1, lbounds, ubounds)
    class(attrib_data_holder_t),intent(inout) :: holder
    real(dp),intent(in) :: v_r1(:)
    integer(ip),intent(in) :: lbounds
        !!データを格納する配列のインデックスの下限.
    integer(ip),intent(in) :: ubounds
        !!データを格納する配列のインデックスの上限.

    integer(ip) ist, imx, icmx
    integer(ip) astat
    character(256) errmsg

    ist = lbounds
    imx = ubounds
    icmx = imx - ist + 1

    allocate(holder%val_r1(icmx), source = v_r1(ist:imx))

end subroutine

subroutine convert_r3_to_r1_and_register_(holder, v_r3, lbounds, ubounds)
    class(attrib_data_holder_t),intent(inout) :: holder
    real(dp),intent(in) :: v_r3(:,:,:)
        !!スカラーデータ配列. 
        !!@note 配列全体を指定する. 有効範囲はlbounds, uboundsで指定する.
    integer(ip),intent(in) :: lbounds(3)
        !!データを格納する配列のインデックスの下限.
    integer(ip),intent(in) :: ubounds(3)
        !!データを格納する配列のインデックスの上限.

    integer(ip) ist, jst, kst, imx, jmx, kmx, icmx
    integer(ip) astat
    character(256) errmsg

    ist = lbounds(1)
    jst = lbounds(2)
    kst = lbounds(3)
    imx = ubounds(1)
    jmx = ubounds(2)
    kmx = ubounds(3)
    icmx = product(ubounds(:) - lbounds(:) + 1)

    allocate(holder%val_r1(icmx), source = reshape(v_r3(ist:imx,jst:jmx,kst:kmx), shape=[icmx]), &
            stat = astat, errmsg = errmsg)

    if ( astat /= 0 ) error stop trim(errmsg) 
end subroutine

subroutine write_data_ascii(holder, unit, data_size)
    !!登録したデータをASCII形式で指定したunitに書き込む.
    class(attrib_data_holder_t),intent(in) :: holder
    integer(ip),intent(in) :: unit
        !!論理装置番号.
    integer(ip),intent(in) :: data_size
        !!ライターが所持するデータ点の個数. コレと一致しない場合書き込みは中止される.

    integer(ip) ic_

    if ( allocated(holder%val_r1)) then
        write(unit,"(A)") "SCALARS "//holder%name//" double"
        write(unit, "(A)") "LOOKUP_TABLE "//LOOKUP_TABLE_NAME
        
        if ( size(holder%val_r1) /= data_size ) error stop "ERROR (vtk_field_data) :: invalid data size." 
        
        do ic_ = lbound(holder%val_r1, dim=1), ubound(holder%val_r1, dim=1)
            write(unit, "(g0)") holder%val_r1(ic_)
        end do
    
    else if ( allocated(holder%val_r2)) then
        write(unit,"(A)") "VECTORS "//holder%name//" double"

        if ( size(holder%val_r2, dim=2) /= data_size ) error stop "ERROR (vtk_field_data) :: invalid data size."

        do ic_ = lbound(holder%val_r2, dim=2), ubound(holder%val_r2, dim=2)
            write(unit, "(*(g0,1x))") holder%val_r2(1:, ic_)
        end do
    end if
end subroutine

subroutine write_data_binary(holder, unit, data_size)
    !!登録したデータをBINARY形式で指定したunitに書き込む.
    !!@note ファイルはstream形式.
    class(attrib_data_holder_t),intent(in) :: holder
    integer(ip),intent(in) :: unit
        !!論理装置番号.
    integer(ip),intent(in) :: data_size
        !!ライターが所持するデータ点の個数. コレと一致しない場合書き込みは中止される.

    integer(ip) ic_
    character(*),parameter :: LF_ = new_line("")

    if ( allocated(holder%val_r1)) then
        write(unit) "SCALARS "//holder%name//" double"//LF_
        write(unit) "LOOKUP_TABLE "//LOOKUP_TABLE_NAME//LF_
        
        if ( size(holder%val_r1) /= data_size ) error stop "ERROR (vtk_field_data) :: invalid data size." 
        
        write(unit) holder%val_r1(:)
        
    else if ( allocated(holder%val_r2)) then
        write(unit) "VECTORS "//holder%name//" double"//LF_

        if ( size(holder%val_r2, dim=2) /= data_size ) error stop "ERROR (vtk_field_data) :: invalid data size."

        do ic_ = lbound(holder%val_r2, dim=2), ubound(holder%val_r2, dim=2)
            write(unit) holder%val_r2(1:, ic_)
        end do
    end if

    write(unit) LF_

end subroutine
    
end module vtk_field_data_m