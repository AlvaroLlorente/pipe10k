!***************************************************************************
! parameters
!***************************************************************************
!
!
!  i_N   Number of radial points n in [1,N]
!  i_K   Maximum k (axial), k in (−K,K)
!  i_M   Maximum m (azimuthal), m in [0,M)
!  i_Mp  Azimuthal periodicity, i.e. m=0,Mp,2Mp,…,(M−1)Mp (set =1 for no symmetry assumption)
!  d_Re    Reynolds number Re or Rmm
!  d_alpha Axial wavenumber α=2π/Lz. \alpha must be 0.25 (8pi)
!  b_const_flux    Enforce constant flux Ub=1/2.
! 
! --------------------------------------------------------  
!
!  i_save_rate1    Save frequency for snapshot data files (timesteps between saves)
!  i_save_rate2    Save frequency for time-series data
!  i_maxtstep      Maximum number of timesteps (no limit if set =-1)
!  d_cpuhours      Maximum number of cpu hours
!  d_time          Start time (taken from state.cdf.in if set =-1d0)
!  d_timestep      Fixed timestep (typically =0.01d0 or dynamically controlled if set =-1d0)
!  d_dterr         Maximum corrector norm, ∥fcorr∥ (typically =1d-5 or set =1d1 to avoid extra corrector iterations)
!  d_courant       Courant number C (unlikely to need changing)
!  d_implicit      Implicitness c (unlikely to need changing)!
!
!***************************************************************************
#include "../parallel.h"
 module parameters
!***************************************************************************
   implicit none
   save

   integer,          parameter :: i_N           = 192 ! n in [1,N] radial r_n
   integer,          parameter :: i_K           = 96 ! ! Este tiene que ser potencia de 2 y 3 y grande!!
   integer,          parameter :: i_M           = 64 ! azimutal 
   integer,          parameter :: i_Mp          = 1  ! Siempre 1. 

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   double precision            :: d_Re          = 3500d0
   double precision            :: d_alpha       = 0.4d0
   logical,          parameter :: b_const_flux  = .true.
   logical,          parameter :: b_mirrorsym   = .false.
   logical,          parameter :: b_shiftrefl   = .false.
   logical,          parameter :: b_shiftrott   = .false.
   double precision, parameter :: d_minE3d      = 1d-5
   
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   integer,          parameter :: i_save_rate1  = 1000
   integer,          parameter :: i_save_rate2  = 10
   integer,          parameter :: i_maxtstep    = 2001
   integer,          parameter :: s_step        = 5 ! interval to take statistics. 
   integer,          parameter :: n_sta         = i_save_rate1/s_step ! número de estadísticas 
   integer,          parameter :: f_step        = 1d8 ! Steps
   !double precision, parameter :: f_step        = 1d99 ! Steps, modification
   double precision, parameter :: d_maxt        = -1d0
   double precision, parameter :: d_cpuhours    = 1d99 !90d0
   double precision, parameter :: d_time        = 0d0 !-1d0 !0d0
   double precision, parameter :: d_timestep    = -1d0 !0.01d0 !-1d0
   double precision, parameter :: d_maxdt       = 1d99
   double precision, parameter :: d_dterr       = 1d-5 !1d99
   double precision, parameter :: d_courant     = 0.5d0
   double precision, parameter :: d_implicit    = 0.51d0

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   integer,          parameter :: i_KL  = 4   
   integer,          parameter :: i_K1  = i_K-1
   integer,          parameter :: i_M1  = i_M-1
   integer,          parameter :: i_Z   = 3*i_K  
   integer,          parameter :: i_Th  = 3*i_M
   integer,          parameter :: i_H1  = (2*i_K1+1)*i_M-i_K1-1
   integer,          parameter :: i_pN  = (_Nr+i_N-1)/_Nr
   integer,          parameter :: i_pH1 = (_Nr+_Hs1)/_Nr-1
   integer,          parameter :: i_pZ  = i_Z/_Ns
   double precision, parameter :: d_PI  = 3.1415926535897931d0

!---------------------------------------------------------------------------
   integer, parameter          :: i_3K = 3*i_K
   integer, parameter          :: i_3M = 3*i_M

 contains

!---------------------------------------------------------------------------
!  check params are ok
!---------------------------------------------------------------------------
   subroutine par_precompute()
      if(modulo(i_Z,_Ns)/=0) stop '_Ns must divide i_Z'
      if(modulo(i_M,_Ns)/=0) stop '_Ns must divide i_M'
   end subroutine par_precompute
 

!***************************************************************************
 end module parameters
!***************************************************************************



