module interpolation_m
    !!セル中心データから面への補間を取り扱うモジュール.
    !!@note 内部と境界とで扱う補間方法は一致しなければならない.
    use floating_point_parameter_m
    implicit none
    
contains

real(dp) elemental pure function arithmetic_average_interpolate(val1, val2) result(res)
    !!算術平均(二点のデータの平均)で補間する.
    real(dp),intent(in) :: val1, val2

    res = 0.5_dp*(val1 + val2)

end function

real(dp) elemental pure function arithmetic_average_extrapolate(val1, mean) result(res)
    !!算術平均(二点のデータの平均)に基づいて外挿する.
    real(dp),intent(in) :: val1, mean

    res = 2.0_dp*mean - val1

end function
    
real(dp) pure elemental function weighted_average_interpolate(val1, val2, w1, w2) result(res)
    !!重み付き平均(=線形補間)により内挿する.
    real(dp),intent(in) :: val1, val2
    real(dp),intent(in) :: w1, w2
    real(dp) alpha_

    alpha_ = 1.0_dp/(w1 + w2)
    res = (w1*val2 + w2*val1)*alpha_

end function

real(dp) pure elemental function weighted_average_extrapolate(val1, avg, w1, w2) result(res)
    !!重み付き平均(=線形補間)により補間値から外挿する.
    real(dp),intent(in) :: val1, avg
    real(dp),intent(in) :: w1, w2
    real(dp) alpha_

    alpha_ = w2/w1
    res = (1.0_dp + alpha_)*avg - alpha_*val1

end function

real(dp) pure elemental function harmonic_like_interpolate(val1, val2, w1, w2) result(res)
    !!逆数の重み付け平均により内挿する. w1=w2のとき調和平均となる.
    !!@note 得られる値の次元はval1(val2)の逆数と同じである.
    real(dp),intent(in) :: val1, val2
    real(dp),intent(in) :: w1, w2
    real(dp) alpha_

    alpha_ = 1.0_dp/(w1 + w2)
    !res = 1/(w1 + w2)*(w1/val1 + w2/val2)
    res = alpha_*(w1/val1 + w2/val2)

end function

real(dp) pure elemental function harmonic_like_extrapolate(val1, avg_inverse, w1, w2) result(res)
    !!逆数の重み付け平均により補間値から外挿する. w1=w2のとき調和平均となる.
    !!@note 得られる値の次元はavg_inverseと同じである.
    real(dp),intent(in) :: val1, avg_inverse
    real(dp),intent(in) :: w1, w2
    real(dp) alpha_
    
    alpha_ = 1.0_dp/w2
    res = alpha_*((w1 + w2)*avg_inverse - w1/val1)

end function

end module interpolation_m