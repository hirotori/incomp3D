module lin_alg_m
    !!3×3行列に対する演算. 
    use floating_point_parameter_m, only : ip, wp_ => dp
    implicit none
    private

    public Gaussian_elimination_3, Clamer_rule_3, Determinant_3, cross_3

    
contains
    subroutine Gaussian_elimination_3(A,b,x)
        !!部分ピボット付きGaussの消去法により3元連立方程式を解く. 
        !!@warning
        !!rank(A)<3に対する動作が未定義. 
        !!@endwarning
        implicit none
        real(wp_),intent(in) :: A(3,3)
        real(wp_),intent(in) :: b(3)
        real(wp_),intent(inout) :: x(3)

        real(wp_):: Mat_(4,4)
        real(wp_),parameter :: eps_ = epsilon(1.0_wp_)

        integer i,k
        real(wp_) m_ii_, m_ki_

        !拡大行列を作成する. 
        Mat_(1:3,1:3) = A(1:3,1:3)
        Mat_(1:3,4) = b(1:3)
        Mat_(4,1:4) = [1,2,3,4] !元々の列番号を保持する. 完全ピボット用. 

        !forward elimination
        do I = 1, 2
            if ( abs(Mat_(I,I)) <= eps_  ) call partial_pivotting(Mat_, I) 
            m_ii_ = Mat_(I,I)
            do K = I+1, 3
                m_ki_ = Mat_(K,I)
                Mat_(K,I+1:4) = Mat_(K,I+1:4) - m_ki_/m_ii_*Mat_(I,I+1:4)
                Mat_(K,I) = 0.0
            end do
            
        end do


        !backward substitution
        x(3) = Mat_(3,4)/Mat_(3,3)
        x(2) = (Mat_(2,4) - Mat_(2,3)*x(3))/Mat_(2,2)
        x(1) = (Mat_(1,4) - Mat_(1,2)*x(2) - Mat_(1,3)*x(3))/Mat_(1,1)


    end subroutine

    subroutine partial_pivotting(mat4, col)
        !!col行に対して部分ピボット選択を実施する. 
        implicit none
        real(wp_),intent(inout) ::mat4(4,4)
        integer,intent(in) :: col
        real(wp_) :: tmp_(4)

        integer max_col

        !col列目において絶対値の最大値を持つ行と入れ替える. 
        max_col = maxloc(abs(mat4(col:3,col)),dim=1) + col - 1
        tmp_(1:4) = mat4(col,1:4)
        mat4(col,1:4) = mat4(max_col,1:4)
        mat4(max_col,1:4) = tmp_(1:4)

        ! print*, "---pivot---"
        ! print*, mat4(1,1:4)
        ! print*, mat4(2,1:4)
        ! print*, mat4(3,1:4)

    end subroutine


    subroutine Clamer_rule_3(A,b,x,status)
        !!クラーメルの公式により連立方程式の解を求める. 
        implicit none
        real(wp_),intent(in) :: A(3,3)
            !!3×3行列
        real(wp_),intent(in) :: b(3)
            !!定数ベクトル
        real(wp_),intent(inout) :: x(3)
            !!解ベクトル
        integer,intent(out) :: status
            !!エラーのとき-1を返す. 
        real(wp_) det3_

        det3_ = Determinant_3(A)
        status = 0

        if(abs(det3_) <= epsilon(1.0_wp_)) then
            ! print"('Clamer_rule::Error: Matrix is not regular.')"
            status = -1 !他のメソッド切り替えように残しておく. 
            return
        endif

        x(1) = dot_product(b(1:3), cross_3(A(1:3,2), A(1:3,3)))/det3_
        x(2) = dot_product(A(1:3,1), cross_3(b(1:3), A(1:3,3)))/det3_
        x(3) = dot_product(A(1:3,1), cross_3(A(1:3,2), b(1:3)))/det3_

    end subroutine
    
    real(wp_) function Determinant_3(mat) result(det_)
        implicit none
        real(wp_),intent(in) :: mat(3,3)

        det_ = dot_product(mat(1:3,1), cross_3(mat(1:3,2), mat(1:3,3)))

    end function

    function cross_3(u,v) result(v_)
        implicit none
        real(wp_),intent(in) ::u(3),v(3)
        real(wp_) v_(3)

        v_(1) = u(2)*v(3) - u(3)*v(2)
        v_(2) = u(3)*v(1) - u(1)*v(3)
        v_(3) = u(1)*v(2) - u(2)*v(1)
    end

end module lin_alg_m

! program test
!     use lin_alg_m
!     implicit none

!     real(8) :: A(3,3) = reshape([1,2,3,4,5,6,7,8,9],shape=[3,3])
!     real(8) :: B(3,3)
!     real(8) :: C(3,3) = reshape([2,1,1,4,2,3,-2,1,2],shape=[3,3])
!     real(8) :: x(3), y(3) = [8.0,6.0,9.0] 
!     integer i, j, stat

!     print*, Determinant_3(A) !=0

!     !Hilbert matrix
!     do j = 1, 3
!     do i = 1, 3
!         B(i,j) = 1.0d0/(i+j-1)
!     end do
!     end do

!     print*, Determinant_3(B) !=1/2160 ~ 0.00046296296296...

!     print*, "-------"
!     print*, C(1,1:3)
!     print*, C(2,1:3)
!     print*, C(3,1:3)
!     call Clamer_rule_3(C, y, x, stat) ; print*, stat

!     call Clamer_rule_3(C, y, x, stat) ; print*, x
!     call Gaussian_elimination_3(C, y, x) ; print*, x
!     print*, matmul(C,x)-y
! end program test