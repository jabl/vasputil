!****h* /mathconst
! PURPOSE
!  Contains mathematical constants.
! NOTES
! by J.E. Akin: http://www.owlnet.rice.edu/~mech517/Books/oop2.pdf.
! Also see Table of Constants to 50 Places at:
! http://numbers.computation.free.fr/Constants/Miscellaneous/digits.html
! Special thanks to Howard W. Ludwig.
! Modified by Janne Blomqvist.
!****
Module mathconst

  use conf

  Implicit none
  
  real(kind=wp), public, parameter::              &
       Deg_Per_Rad  = 57.2957795130823208767981548_wp, &
       Rad_Per_Deg  = 0.01745329251994329576923691_wp, &
       
  e_Value      =  2.7182818284590452353602287_wp, &
       e_Recip      =  0.3678794411714423215955238_wp, &
       e_Squared    =  7.3890560989306502272304275_wp, &
       Log10_of_e   =  0.4342944819032518276511289_wp, &
       
  Euler        =  0.5772156649015328606065121_wp, &
       Euler_Log    = -0.5495393129816448223376617_wp, &
       Gamma        =  0.5772156649015328606065121_wp, &
       Gamma_Log    = -0.5495393129816448223376617_wp, &
       Golden_Ratio =  1.6180339887498948482045868_wp, &
       
  Ln_2         =  0.6931471805599453094172321_wp, &
       Ln_10        =  2.3025850929940456840179915_wp, &
       Log10_of_2   =  0.3010299956639811952137389_wp, &
       
  pi_Value     =  3.1415926535897932384626434_wp, &
       pi_Ln        =  1.1447298858494001741434273_wp, &
       pi_Log10     =  0.4971498726941338543512683_wp, &
       pi_Over_2    =  1.5707963267948966192313217_wp, &
       pi_Over_3    =  1.0471975511965977461542145_wp, &
       pi_Over_4    =  0.7853981633974483096156608_wp, &
       pi_Recip     =  0.3183098861837906715377675_wp, &
       pi_Squared   =  9.8696044010893586188344910_wp, &
       pi_Sq_Root   =  1.7724538509055160272981675_wp, &
       
  Sq_Root_of_2 =  1.4142135623730950488016887_wp, &
       Sq_Root_of_3 =  1.7320508075688772935274463_wp

End Module Mathconst
