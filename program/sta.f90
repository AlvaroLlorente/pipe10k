
! 1) Otf, max = 54.4
! 2) max =  52.1
! 3) max =  51.0
! 4) max =  49
! 5) max =  48.5
! 6) max =  48.4
! 7) max =  
! 8) max =  
! 9) max =  


!**************************************************************************
!  IN/OUT 
!
!**************************************************************************
#include "../parallel.h"
 module sta
!**************************************************************************
   use wksp
   use variables
   use velocity
   use h5lt
   use hdf5
   use transform


   implicit none
   save


   double precision, private :: uclm, utaum,ucl, utau, Ub
   double precision, private :: aa, bb, cc
   

! ------------------------- stats  -------------------------------

   double precision :: mean_ur(i_N,n_sta), stdv_ur(i_N,n_sta)
   double precision :: mean_ut(i_N,n_sta), stdv_ut(i_N,n_sta)
   double precision :: mean_uz(i_N,n_sta), stdv_uz(i_N,n_sta)
   double precision :: stdv_rz(i_N,n_sta), stdv_rt(i_N,n_sta), stdv_tz(i_N,n_sta)
   ! double precision :: mom_ur(i_N,10)
   ! double precision :: mom_ut(i_N,10)
   ! double precision :: mom_uz(i_N,10)
   double precision :: mean_p(i_N,n_sta), stdv_p(i_N,n_sta)

   !double precision :: piz(i_N), pit(i_N), pir(i_N)
   double precision :: durdr(i_N,n_sta), durdt(i_N,n_sta), durdz(i_N,n_sta)
   double precision :: dutdr(i_N,n_sta), dutdt(i_N,n_sta), dutdz(i_N,n_sta)
   double precision :: duzdr(i_N,n_sta), duzdt(i_N,n_sta), duzdz(i_N,n_sta)
   !double precision :: dissr(i_N,3),disst(i_N,3),dissz(i_N,3), diss(i_N,3) !, dzduzsq(i_N), dzduzcub(i_N)
   double precision :: PDT2(i_N),TDT2(i_N),DT1(i_N,n_sta), DT4(i_N,n_sta), DT5(i_N,n_sta) , DT6(i_N,n_sta), DT7(i_N,n_sta)
   double precision :: factor
   double precision :: time_sta(n_sta), utauv(n_sta)

   !double precision, private :: d(i_N) ,dd(i_N,n_sta) ! auxiliary mem
   integer :: csta
   
! ------------------------- HDF5 -------------------------------
!desde el nuevo portatil
   integer:: info,ierr
   integer(hid_t):: fid,pid, dset_id,dspace_id
   integer:: h5err
   integer(hsize_t),dimension(3):: dims
   integer(hsize_t),dimension(1):: hdims,maxdims
   integer(hsize_t),dimension(2):: hdims2
   
   integer(hid_t) :: header_id,sta_id, budget_id ! Group identifiers
   integer(size_t)::siever
   parameter (siever = 4*1024*1024)


! --------------------- Program -------------------------------------

 contains

! Compute statistics
   subroutine compute_sta()
   implicit none
   integer:: n, n_, ii, jj, kk 

   ! Compute friction velocity

   if (mpi_rnk ==0 ) then
      time_sta(csta) = tim_t
      ucl = 1d0 + dot_product(vel_uz%Re(1:1+i_KL,0),mes_D%dr0(:,0))
      utau = dot_product(vel_uz%Re(i_N-i_KL:i_N,0),mes_D%dr1(:,1))
      utau = dsqrt(dabs((Utau-2d0)/d_Re))

      ! add statistics
      uclm = uclm + ucl
      utaum = utaum + utau 
      utauv(csta) = utau 
   endif

      
      call vel_sta()

      call pressure(c1,c2,c3,p1,p2)
      !call var_coll_dissp(c1,c2,c3,c4)
       




      

   do n = 1, mes_D%pN
      n_ = mes_D%pNi + n - 1
      mean_ur(n_,csta) = mean_ur(n_,csta) + sum(vel_r%Re(:,:,n))
      stdv_ur(n_,csta) = stdv_ur(n_,csta) + sum(vel_r%Re(:,:,n)**2)
      mean_ut(n_,csta) = mean_ut(n_,csta) + sum(vel_t%Re(:,:,n))
      stdv_ut(n_,csta) = stdv_ut(n_,csta) + sum(vel_t%Re(:,:,n)**2)
      mean_uz(n_,csta) = mean_uz(n_,csta) + sum(vel_z%Re(:,:,n))
      stdv_uz(n_,csta) = stdv_uz(n_,csta) + sum(vel_z%Re(:,:,n)**2)
      stdv_rz(n_,csta) = stdv_rz(n_,csta) + sum(vel_z%Re(:,:,n)*vel_r%Re(:,:,n))
      stdv_rt(n_,csta) = stdv_rz(n_,csta) + sum(vel_r%Re(:,:,n)*vel_t%Re(:,:,n))
      stdv_tz(n_,csta) = stdv_rz(n_,csta) + sum(vel_t%Re(:,:,n)*vel_z%Re(:,:,n))
   enddo

   call compute_turb_budget()

   csta = csta + 1
  
   call mpi_barrier(mpi_comm_world, mpi_er)

end subroutine compute_sta

!------------------------------------------------------------------------
!  Initialize pressure calculation
!------------------------------------------------------------------------
   subroutine p2m_lumesh_inits(PM,BC,c11,c21, A)
      use timestep
      integer,          intent(in)  :: PM,BC
      double precision, intent(in)  :: c11,c21
      type (lumesh),    intent(inout) :: A(0:i_pH1)

      integer :: info, n,j, S
      _loop_km_vars

      _loop_km_begin
         d = -mes_D%r(:,-2)*i_Mp*m*i_Mp*m - d_alpha*k*d_alpha*k
         if(PM/=0) d = d - mes_D%r(:,-2) - 2d0*PM*i_Mp*m*mes_D%r(:,-2)
         A(nh)%M(i_KL+1:, :) = c21 * mes_D%radLap%M(:,1:)
         A(nh)%M(2*i_KL+1,:) = A(nh)%M(2*i_KL+1,:) + c21*d + c11

         ! assume symmetry on axis: S==-1 mode odd,  S==1 mode even
         S = modulo(m*i_Mp+abs(PM),2)
         S = 1 - 2*S
         do j = 1, i_KL
            do n = 1, i_KL+1-j
               A(nh)%M(2*i_KL+1+n-j, j) = A(nh)%M(2*i_KL+1+n-j, j)  &
                  + c21 * S * mes_D%radLap%M(i_KL+1+n-(1-j), (1-j))
            end do
         end do
                                        ! boundary condition
         do j = i_N-i_KL, i_N
            A(nh)%M(2*i_KL+1+i_N-j,j) = mes_D%dr1(i_KL-i_N+j+1,1)
         end do

        if(BC==1 .and. k==0 .and. m==0) then
         do j = i_N-i_KL, i_N
            A(nh)%M(2*i_KL+1+i_N-j,j) = mes_D%dr1(i_KL-i_N+j+1,0)
         end do
         end if
         call dgbtrf(i_N,i_N,i_KL,i_KL,A(nh)%M,3*i_KL+1,A(nh)%ipiv,info)
         if(info /= 0) stop 'tim_lumesh_init'
      _loop_km_end

   end subroutine p2m_lumesh_inits

! !------------------------------------------------------------------------
! !  Pressure field:
! !      Results: 
! !      
! !------------------------------------------------------------------------
subroutine pressure(c1,c2,c3,p1,p2)
implicit none
_loop_km_vars
integer:: n, n_
double precision :: BCR(0:i_pH1), BCI(0:i_pH1)
type (coll), intent(inout)    :: c1,c2,c3
type (phys), intent(inout)    :: p1,p2

! Necesitamos: 3 colls, 2 phys
! Primer cambio



         call var_coll_curl(vel_ur,vel_ut,vel_uz, c1,c2,c3)
         call var_coll_curl(c1,c2,c3, c1,c2,c3)
         
         _loop_km_begin
           BCR(nh) = - c1%Re(i_N,nh)/d_Re
           BCI(nh) = - c1%Im(i_N,nh)/d_Re
         _loop_km_end

         ! r equation 
         !durdr ----> Alvaro: he añadido dutdr y duzdr aqui porque en compute_turb_budget
         !por algun motivo ocurre segmentation fault
         call var_coll_meshmult(1,mes_D%dr(1),vel_ur,c1) !durdr
         call tra_coll2phys1d(c1,p1)
         call var_coll_meshmult(1,mes_D%dr(1),vel_ut,c3) !dutdr
         call tra_coll2phys1d(c3,p3)
         call var_coll_meshmult(1,mes_D%dr(1),vel_uz,c4) !duzdr
         call tra_coll2phys1d(c4,p4)
         do n = 1, mes_D%pN
            n_ = mes_D%pNi + n - 1
            durdr(n_,csta) = durdr(n_,csta) + sum(p1%Re(:,:,n))
            dutdr(n_,csta) = dutdr(n_,csta) + sum(p3%Re(:,:,n))
            duzdr(n_,csta) = duzdr(n_,csta) + sum(p4%Re(:,:,n))
         enddo 

         p2%Re=vel_r%Re*p1%Re
      
         ! 1/r(durdt-ut)
         _loop_km_begin
               c1%Re(:,nh) = mes_D%r(:,-1)*(-vel_ur%Im(:,nh)*m*i_Mp-vel_ut%Re(:,nh))
               c1%Im(:,nh) = mes_D%r(:,-1)*( vel_ur%Re(:,nh)*m*i_Mp-vel_ut%Im(:,nh))   
         _loop_km_end
         call tra_coll2phys1d(c1,p1)
         p2%Re=p2%Re+vel_t%Re*p1%Re
      
         ! durdz
         _loop_km_begin
               c1%Re(:,nh) = -vel_ur%Im(:,nh)*d_alpha*k
               c1%Im(:,nh) =  vel_ur%Re(:,nh)*d_alpha*k       
         _loop_km_end
         call tra_coll2phys1d(c1,p1)
         p2%Re=p2%Re+vel_z%Re*p1%Re
         do n = 1, mes_D%pN
            n_ = mes_D%pNi + n - 1
   	    p2%Re(:,:,n)=p2%Re(:,:,n)+p1%Re(:,:,n)*vel_U(n_)	
         end do

         call tra_phys2coll1d(p2,c1)

         ! theta equation
         ! dutdr
         call var_coll_meshmult(1,mes_D%dr(1),vel_ut,c2)
         call tra_coll2phys1d(c2,p1)
         p2%Re=vel_r%Re*p1%Re
      
         ! 1/r(dutdt+ur)
         _loop_km_begin   
               c2%Re(:,nh) = mes_D%r(:,-1)*(-vel_ut%Im(:,nh)*m*i_Mp+vel_ur%Re(:,nh))
               c2%Im(:,nh) = mes_D%r(:,-1)*( vel_ut%Re(:,nh)*m*i_Mp+vel_ur%Im(:,nh))
         _loop_km_end
         call tra_coll2phys1d(c2,p1)
         p2%Re=p2%Re+vel_t%Re*p1%Re
      
         ! dutdz
         _loop_km_begin          
               c2%Re(:,nh) = -vel_ut%Im(:,nh)*d_alpha*k
               c2%Im(:,nh) =  vel_ut%Re(:,nh)*d_alpha*k          
         _loop_km_end
         call tra_coll2phys1d(c2,p1)
         p2%Re=p2%Re+vel_z%Re*p1%Re
         do n = 1, mes_D%pN
            n_ = mes_D%pNi + n - 1	
            p2%Re(:,:,n)=p2%Re(:,:,n)+p1%Re(:,:,n)*vel_U(n_)
         end do
         call tra_phys2coll1d(p2,c2)

         ! z equation
         ! duzdr 
         call var_coll_meshmult(0,mes_D%dr(1),vel_uz,c3)
         call tra_coll2phys1d(c3,p1)
         p2%Re=vel_r%Re*p1%Re
      
         ! 1/r(duzdt)
         _loop_km_begin      
            c3%Re(:,nh) = mes_D%r(:,-1)*(-vel_uz%Im(:,nh)*m*i_Mp)
            c3%Im(:,nh) = mes_D%r(:,-1)*( vel_uz%Re(:,nh)*m*i_Mp)
         _loop_km_end
         call tra_coll2phys1d(c3,p1)
         p2%Re=p2%Re+vel_t%Re*p1%Re
      
         ! duzdz
         _loop_km_begin
            c3%Re(:,nh) = -vel_uz%Im(:,nh)*d_alpha*k
            c3%Im(:,nh) =  vel_uz%Re(:,nh)*d_alpha*k
        _loop_km_end
        
         call tra_coll2phys1d(c3,p1)
         p2%Re=p2%Re+vel_z%Re*p1%Re
         do n = 1, mes_D%pN
            n_ = mes_D%pNi + n - 1	
            p2%Re(:,:,n)=p2%Re(:,:,n)+p1%Re(:,:,n)*vel_U(n_)+vel_r%Re(:,:,n)*vel_Up(n_)	
         end do
         call tra_phys2coll1d(p2,c3)


         _loop_km_begin
            BCR(nh) = BCR(nh) - c1%Re(i_N,nh)
            BCI(nh) = BCI(nh) - c1%Im(i_N,nh)
         _loop_km_end
         call var_coll_div(c1,c2,c3, c1)
         c1%Re=-c1%Re;
         c1%Im=-c1%Im;
         _loop_km_begin
            c1%Re(i_N,nh) = BCR(nh)
            c1%Im(i_N,nh) = BCI(nh)
         _loop_km_end

         ! call tim_zerobc(c1)
         call p2m_lumesh_inits( 0,1,0d0,1d0, LNp)
         call tim_lumesh_invert(0,LNp, c1)
         call tra_coll2phys1d(c1,p2) !pressure field


   do n = 1, mes_D%pN
      n_ = mes_D%pNi + n - 1
      mean_p(n_,csta)  = mean_p(n_,csta)  + sum(p2%Re(:,:,n))
      stdv_p(n_,csta)  = stdv_p(n_,csta)  + sum(p2%Re(:,:,n)**2) 
   end do



end subroutine pressure

! !------------------------------------------------------------------------
! !  Derivatives & Turbulent budget:
! !     vel r,vel t, vel z     
! !------------------------------------------------------------------------
subroutine compute_turb_budget()
   
   implicit none
   integer :: n,n_
   _loop_km_vars

   
   !Estoy reservando p2 para el coll del campo de presiones
   

!!--------Derivatives-------!!

!!   vel_r
   

_loop_km_begin

 c3%Im(:,nh) = -vel_ur%Im(:,nh)*ad_k1a1(k)
 c3%Re(:,nh) =  vel_ur%Re(:,nh)*ad_k1a1(k)

 c4%Im(:,nh) = -vel_ur%Im(:,nh)*ad_m1r1(:,m)
 c4%Re(:,nh) =  vel_ur%Re(:,nh)*ad_m1r1(:,m)

_loop_km_end


call tra_coll2phys1d(c3,p3) !durdt


call tra_coll2phys1d(c4,p4) !durdz

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1

   durdt(n_,csta) = durdt(n_,csta) + sum(p3%Re(:,:,n)) 
   durdz(n_,csta) = durdz(n_,csta) + sum(p4%Re(:,:,n)) 

end do

!   vel_t

 _loop_km_begin

 c3%Im(:,nh) = -vel_ut%Im(:,nh)*ad_k1a1(k)
 c3%Re(:,nh) =  vel_ut%Re(:,nh)*ad_k1a1(k)

 c4%Im(:,nh) = -vel_ut%Im(:,nh)*ad_m1r1(:,m)
 c4%Re(:,nh) =  vel_ut%Re(:,nh)*ad_m1r1(:,m)

_loop_km_end

call tra_coll2phys1d(c1,p3) !dutdt
call tra_coll2phys1d(c3,p4) !dutdz

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
    
   dutdt(n_,csta) = dutdt(n_,csta) + sum(p3%Re(:,:,n)) 
   dutdz(n_,csta) = dutdz(n_,csta) + sum(p4%Re(:,:,n)) 

end do

!   vel_z

_loop_km_begin

 c3%Im(:,nh) = -vel_uz%Im(:,nh)*ad_k1a1(k)
 c3%Re(:,nh) =  vel_uz%Re(:,nh)*ad_k1a1(k)

 c4%Im(:,nh) = -vel_uz%Im(:,nh)*ad_m1r1(:,m)
 c4%Re(:,nh) =  vel_uz%Re(:,nh)*ad_m1r1(:,m)

_loop_km_end

call tra_coll2phys1d(c2,p3) !duzdt
call tra_coll2phys1d(c2,p4) !duzdz

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
   duzdt(n_,csta) = duzdt(n_,csta) + sum(p3%Re(:,:,n)) 
   duzdz(n_,csta) = duzdz(n_,csta) + sum(p4%Re(:,:,n)) 

end do
!--------Turbulent kinetic energy budget-------!!
!  Pressure difussion term
p1%Re = p2%Re*vel_r%Re  !presion · vel radial fisico

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
   PDT2(n_)  = PDT2(n_)  + sum(p1%Re(:,:,n)) ! saco la distribucion radial
   

end do
!En matlab derivar en r y dividir por r

!!  Turbulent difussion term 

p1%Re=vel_r%Re*(vel_r%Re**2+vel_t%Re**2+vel_z%Re**2)

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
   TDT2(n_)  = TDT2(n_)  + sum(p1%Re(:,:,n)) ! saco la distribucion radial
   TDT2(n_)  = TDT2(n_)*mes_D%r(n_,1) ! multiplico por r
end do
!En matlab derivar en r y dividir por r

!!  Viscous difussion term 

!Todos los terminos en matlab

!!  Production term

!Todos los terminos en matlab

!!   Dissipation term

!---Termino 1

   DT1(:,csta)=DT1(:,csta)+sum(duzdz(:,csta)**2)

!---Termino 4

   DT4(:,csta)=DT4(:,csta)+(dutdt(:,csta)*mes_D%r(:,-1))**2

!---Terno 5

   DT5(:,csta)=DT5(:,csta)+(durdz(:,csta)+duzdr(:,csta))**2

!---Termino 6

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
   p1%Re(:,:,n)=vel_t%Re(:,:,n)*mes_D%r(n_,-1)
enddo

 call tra_phys2coll1d(p1,c1) !pasar a coll

 call var_coll_meshmult(0,mes_D%dr(1),c1, c2) !Derivar con respecto a r

 call tra_coll2phys1d(c2,p1) !pasar a fisico

 do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
   dd(n_,csta)=dd(n_,csta)+sum(p1%Re(:,:,n))
   enddo

   DT6(:,csta)=DT6(:,csta)+(durdt(:,csta)*mes_D%r(:,-1))+(dd(:,csta)*mes_D%r(:,1))

!---Termino 7

   DT7(:,csta)=DT7(:,csta)+(dutdz(:,csta)+mes_D%r(:,-1)*duzdt(:,csta))


end subroutine compute_turb_budget


! !------------------------------------------------------------------------
! !  Turbulent budgets:
! !      Dissipation, % optimization pending
! !      Rest of the derivatives: budgets in postproc
! !------------------------------------------------------------------------
! subroutine var_coll_dissp(c1,c2,c3,c4)
!    implicit none
!   type(coll), intent(inout)  :: c1,c2,c3,c4
!
!       double precision :: factor
!       integer :: n, n_
!       _loop_km_vars
!
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       call var_coll_meshmult(0,mes_D%dr(1),vel_uz, c2) 
!
!       _loop_km_begin
!
!       c3%Im(:,nh) = -vel_uz%Im(:,nh)*ad_k1a1(k)
!       c3%Re(:,nh) =  vel_uz%Re(:,nh)*ad_k1a1(k)
!
!       c1%Im(:,nh) = -vel_uz%Im(:,nh)*ad_m1r1(:,m)
!       c1%Re(:,nh) =  vel_uz%Re(:,nh)*ad_m1r1(:,m)
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       factor = 2d0
!       if (m==0) factor = 1d0
!
!          dissr(:,3) = dissr(:,3) + factor*(c2%Re(:,nh)**2+c2%Im(:,nh)**2)
!          dissz(:,3) = dissz(:,3) + factor*(c3%Re(:,nh)**2+c3%Im(:,nh)**2)
!          disst(:,3) = disst(:,3) + factor*(c1%Re(:,nh)**2+c1%Im(:,nh)**2)
!     
!       _loop_km_end
!  
!
!       ! Lo mismo para theta
!
!      call var_coll_meshmult(1,mes_D%dr(1),vel_ut, c2) 
!
!      _loop_km_begin
!
!      c4%Re(:,nh) = -vel_ut%Im(:,nh)*ad_k1a1(k)
!      c4%Im(:,nh) =  vel_ut%Re(:,nh)*ad_k1a1(k)
!
!      c1%Im(:,nh) = -vel_ut%Im(:,nh)*ad_m1r1(:,m) + vel_ur%Re(:,nh)*mes_D%r(:,-1) !+
!      c1%Re(:,nh) =  vel_ut%Re(:,nh)*ad_m1r1(:,m) + vel_ur%Im(:,nh)*mes_D%r(:,-1) !+
!
!      factor = 2d0
!      if (m==0) factor = 1d0
!      
!         dissr(:,2) = dissr(:,2) + factor*(c2%Re(:,nh)**2+c2%Im(:,nh)**2)
!         dissz(:,2) = dissz(:,2) + factor*(c4%Re(:,nh)**2+c4%Im(:,nh)**2)
!         disst(:,2) = disst(:,2) + factor*(c1%Re(:,nh)**2+c1%Im(:,nh)**2)
!
!      _loop_km_end
!
!      
!      ! Lo mismo para r
!
!      call var_coll_meshmult(1,mes_D%dr(1),vel_ur, c2)
!
!      _loop_km_begin
!
!      c4%Im(:,nh) = -vel_ur%Im(:,nh)*ad_k1a1(k)
!      c4%Re(:,nh) =  vel_ur%Re(:,nh)*ad_k1a1(k)
!
!      c1%Im(:,nh) = -vel_ur%Im(:,nh)*ad_m1r1(:,m) - vel_ut%Re(:,nh)*mes_D%r(:,-1) !-
!      c1%Re(:,nh) =  vel_ur%Re(:,nh)*ad_m1r1(:,m) - vel_ut%Im(:,nh)*mes_D%r(:,-1) !-
!
!
!      factor = 2d0
!      if (m==0) factor = 1d0
!
!         dissr(:,1) = dissr(:,1) + factor*(c2%Re(:,nh)**2+c2%Im(:,nh)**2)
!         dissz(:,1) = dissz(:,1) + factor*(c4%Re(:,nh)**2+c4%Im(:,nh)**2)
!         disst(:,1) = disst(:,1) + factor*(c1%Re(:,nh)**2+c1%Im(:,nh)**2)
!
!   _loop_km_end
!
!                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                        !!   Rest of Turbulent budgets !! 
!                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! Pressure strain correlation term
!
!   call var_coll_meshmult(1,mes_D%dr(1),vel_ur, c1) ! r simple derivative
!   call tra_phys2coll1d(p2,c4)
!   ! call tra_phys2coll1d(vel_t,vel_ut)
!   ! call tra_phys2coll1d(vel_z,vel_uz)
!
!         _loop_km_begin
!
!                  c2%Im(:,nh) =  vel_ut%Im(:,nh) *ad_m1r1(:,m)*c4%Im(:,nh) +  vel_ut%Re(:,nh)  *ad_m1r1(:,m)*c4%Re(:,nh)
!                  ! c2%Re(:,nh) =  vel_ut%Re(:,nh)  *ad_m1r1(:,m)    
!
!                  c3%Im(:,nh) =  vel_uz%Im(:,nh)*ad_k1a1(k)*c4%Im(:,nh)  +  vel_uz%Re(:,nh)*ad_k1a1(k)*c4%Re(:,nh)
!                  ! c3%Re(:,nh) =  vel_uz%Re(:,nh)*ad_k1a1(k)
!
!                  c1%Re(:,nh) = c1%Re(:,nh)*c4%Re(:,nh) + c1%Im(:,nh)*c4%Im(:,nh)
!
!             factor = 2d0
!            if (m==0) factor = 1d0
!               pir(:) = pir(:) + factor*(c1%Re(:,nh))
!               piz(:) = piz(:) + factor*(c2%Im(:,nh))
!               pit(:) = pit(:) + factor*(c3%Im(:,nh))
!          _loop_km_end
!
!
!
!      ! call tra_coll2phys1d(c3,p3) !z
!      ! call tra_coll2phys1d(c2,p1) !t
!      ! ! ! call tra_coll2phys1d(c1,p4) !r
!
!      ! p3%Re = p3%Re * p2%Re !z
!      ! p1%Re = p1%Re * p2%Re !t
!      ! ! p4%Re = p4%Re * p2%Re !r
!
!
!      ! do n = 1, mes_D%pN
!      ! n_ = mes_D%pNi + n - 1
!      
!      ! ! p1%Re(:,:,n) =  2d0 * p1%Re(:,:,n)  ! multiplico por 2 y divido entre r
!      ! ! p3%Re(:,:,n) =  2d0 * p3%Re(:,:,n)
!      ! ! pir(n_)  = pir(n_)  + sum(p4%Re(:,:,n)) ! saco la distribucion radial
!      ! pit(n_)  = pit(n_)  + sum(p1%Re(:,:,n))
!      ! piz(n_)  = piz(n_)  + sum(p3%Re(:,:,n))
!      ! end do
!
!
!
!
!
!      ! call tra_phys2coll1d(p3,c2) !z
!      ! call tra_phys2coll1d(p1,c3) !t
!      ! call tra_phys2coll1d(p4,c1) !r
!
!
!   ! _loop_km_begin
!   !          ! z component  
!   !          c2%Re(:,nh) =  vel_uz%Re(:,nh)*ad_k1a1(k)*c4%Re(:,nh) + vel_uz%Im(:,nh)*ad_k1a1(k)*c4%Im(:,nh)
!   !          ! theta component
!   !          c3%Re(:,nh) =  vel_ut%Re(:,nh)*ad_m1r1(:,m)*c4%Re(:,nh) + vel_ut%Im(:,nh)*ad_m1r1(:,m)*c4%Im(:,nh)
!   !          ! r component
!   !          ! c1%Re(:,nh) =  c1%Re(:,nh)*c4%Re(:,nh) + c1%Im(:,nh)*c4%Im(:,nh)
!
!   !    _loop_km_end
!
!   ! _loop_km_begin
!   !    factor = 2d0
!   !    if (m==0) factor = 1d0
!
!   !       piz(:) = piz(:) + factor*(c2%Re(:,nh) + c2%Im(:,nh))!+c2%Im(:,nh)**2)
!   !       pit(:) = pit(:) + factor*(c3%Re(:,nh) + c3%Im(:,nh))!+c3%Im(:,nh)**2)
!   !       ! pir(:) = pir(:) + factor*(c1%Re(:,nh))!+c1%Im(:,nh)**2)
!   ! _loop_km_end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!           SEMI FUNCIONA
!! Pressure strain correlation term
!
!   ! call var_coll_meshmult(1,mes_D%dr(1),vel_ur, c1) ! r simple derivative
!   ! call tra_phys2coll1d(p2,c4)
!
!   !       _loop_km_begin
!
!   !                c2%Im(:,nh) =  -vel_ut%Im(:,nh)*m*i_Mp    
!   !                c2%Re(:,nh) =  vel_ut%Re(:,nh)*m*i_Mp    
!
!   !                c3%Im(:,nh) =  -vel_uz%Im(:,nh)*d_alpha*k
!   !                c3%Re(:,nh) =  vel_uz%Re(:,nh)*d_alpha*k
!
!   !                c1%Re(:,nh) = c1%Re(:,nh)*c4%Re(:,nh) + c1%Im(:,nh)*c4%Im(:,nh)
!
!   !           factor = 2d0
!   !          if (m==0) factor = 1d0
!   !             pir(:) = pir(:) + factor*(c1%Re(:,nh))
!   !    _loop_km_end
!
!
!
!   !    call tra_coll2phys1d(c3,p3) !z
!   !    call tra_coll2phys1d(c2,p1) !t
!   !    ! ! call tra_coll2phys1d(c1,p4) !r
!
!   !    p3%Re = p3%Re * p2%Re !z
!   !    p1%Re = p1%Re * p2%Re !t
!   !    ! ! p4%Re = p4%Re * p2%Re !r
!
!
!   !    ! do n = 1, mes_D%pN
!   !    ! n_ = mes_D%pNi + n - 1
!   !    ! p1%Re(:,:,n) =  2d0 * p1%Re(:,:,n) * mes_D%r(n_,-1) ! multiplico por 2 y divido entre r
!   !    ! p3%Re(:,:,n) =  2d0 * p3%Re(:,:,n)
!   !    ! ! pir(n_)  = pir(n_)  + sum(p4%Re(:,:,n)) ! saco la distribucion radial
!   !    ! pit(n_)  = pit(n_)  + sum(p1%Re(:,:,n))
!   !    ! piz(n_)  = piz(n_)  + sum(p3%Re(:,:,n))
!   !    ! end do
!
!
!
!
!
!   !    call tra_phys2coll1d(p3,c2) !z
!   !    call tra_phys2coll1d(p1,c3) !t
!   !    ! call tra_phys2coll1d(p4,c1) !r
!
!
!   ! ! _loop_km_begin
!   ! !          ! z component  
!   ! !          c2%Re(:,nh) =  vel_uz%Re(:,nh)*ad_k1a1(k)*c4%Re(:,nh) + vel_uz%Im(:,nh)*ad_k1a1(k)*c4%Im(:,nh)
!   ! !          ! theta component
!   ! !          c3%Re(:,nh) =  vel_ut%Re(:,nh)*ad_m1r1(:,m)*c4%Re(:,nh) + vel_ut%Im(:,nh)*ad_m1r1(:,m)*c4%Im(:,nh)
!   ! !          ! r component
!   ! !          ! c1%Re(:,nh) =  c1%Re(:,nh)*c4%Re(:,nh) + c1%Im(:,nh)*c4%Im(:,nh)
!
!   ! !    _loop_km_end
!
!   ! _loop_km_begin
!   !    factor = 2d0
!   !    if (m==0) factor = 1d0
!
!   !       piz(:) = piz(:) + factor*(c2%Re(:,nh) + c2%Im(:,nh))!+c2%Im(:,nh)**2)
!   !       pit(:) = pit(:) + factor*(c3%Re(:,nh) + c3%Im(:,nh))!+c3%Im(:,nh)**2)
!   !       ! pir(:) = pir(:) + factor*(c1%Re(:,nh))!+c1%Im(:,nh)**2)
!   ! _loop_km_end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         !!!!!!!!!!!!!!!!!!!     DERIVATIVES       !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! UZSQUR
!      p1%Re = vel_z%Re * vel_z%Re * vel_z%Re
!      ! do n = 1, mes_D%pN
!      ! n_ = mes_D%pNi + n - 1
!      ! ! p1%Re(:,:,n) = -2d0 * mes_D%r(n_,-1) * p1%Re(:,:,n) ! multiplico por 2 y divido entre r
!      ! ! pur(n_)  = pur(n_)  + sum(p1%Re(:,:,n)) ! saco la distribucion radial
!      ! p1%Re(:,:,n) = ( vel_z%Re(:,:,n) - vel_U(n_)) ** 3
!      ! end do
!
!      
!! UTSQUR
!      p3%re = vel_t%Re * vel_t%Re * vel_r%Re
!! URCUB
!      p4%re = vel_r%Re * vel_r%Re * vel_r%Re
!
!      do n = 1, mes_D%pN
!            n_ = mes_D%pNi + n - 1
!               
!               uzsqur(n_)  = uzsqur(n_)  + sum(p1%Re(:,:,n)) ! saco la distribucion radial
!               utsqur(n_)  = utsqur(n_)  + sum(p3%Re(:,:,n)) ! saco la distribucion radial
!               urcub(n_)   = urcub(n_)   + sum(p4%Re(:,:,n)) ! saco la distribucion radial
!      end do
!
!
!
!
!
!
!
!
!
!
!
!
!! Convection term
!   
!      ! p1%Re = vel_r%Re*vel_r%Re
!
!
!
!      ! do n = 1, mes_D%pN
!      ! n_ = mes_D%pNi + n - 1
!      ! p1%Re(:,:,n)  = p1%Re(:,:,n) * mes_D%r(n_,1)
!      ! end do
!
!      ! call tra_phys2coll1d(p1,c1)
!      ! call var_coll_meshmult(1,mes_D%dr(1),c1, c1) ! r simple derivative
!      ! call tra_coll2phys1d(c1,p1)
!
!      ! do n = 1, mes_D%pN
!      ! n_ = mes_D%pNi + n - 1
!      ! p1%Re(:,:,n) = -2d0 * mes_D%r(n_,-1) * p1%Re(:,:,n)
!      ! pur(n_)  = pur(n_)  + sum(p1%Re(:,:,n))
!      ! end do
!
!
!! Viscous diffusion term
!
!   ! p1%Re = vel_z%Re * vel_z%Re
!   ! p3%Re = vel_r%Re * vel_r%Re
!   ! p4%Re = vel_t%Re * vel_t%Re
!
!   ! call tra_phys2coll1d(p1,c1) !z
!   ! call tra_phys2coll1d(p3,c2) !r
!   ! call tra_phys2coll1d(p4,c3) !t
!
!
!   !    _loop_km_begin
!
!   !    c1%Im(:,nh) = -c1%Im(:,nh)*ad_k1a1(k)*ad_k1a1(k) !z
!   !    c1%Re(:,nh) =  c1%Re(:,nh)*ad_k1a1(k)*ad_k1a1(k)
!
!   !    c2%Im(:,nh) = -c2%Im(:,nh)*ad_k1a1(k)*ad_k1a1(k) !r
!   !    c2%Re(:,nh) =  c2%Re(:,nh)*ad_k1a1(k)*ad_k1a1(k)
!
!   !    c3%Im(:,nh) = -c3%Im(:,nh)*ad_k1a1(k)*ad_k1a1(k) !t
!   !    c3%Re(:,nh) =  c3%Re(:,nh)*ad_k1a1(k)*ad_k1a1(k)
!
!   !    _loop_km_end
!
!   !    call tra_coll2phys1d(c1,p1) !z
!   !    call tra_coll2phys1d(c2,p3) !r
!   !    call tra_coll2phys1d(c3,p4) !t
!
!   ! do n = 1, mes_D%pN
!   !    n_ = mes_D%pNi + n - 1
!   !    duzsqdz2(n_)  = duzsqdz2(n_)  + sum(p1%Re(:,:,n))
!   !    dutsqdz2(n_)  = dutsqdz2(n_)  + sum(p4%Re(:,:,n))
!   !    dursqdz2(n_)  = dursqdz2(n_)  + sum(p3%Re(:,:,n))
!       
!   ! end do
!
!      ! _loop_km_begin
!
!      !    factor = 2d0
!      !    if (m==0) factor = 1d0
!
!      !       pur(:) = pur(:) + factor*(c1%Re(:,nh))
!
!      ! _loop_km_end
!
!
!
!      ! do n = 1, mes_D%pN
!      !    n_ = mes_D%pNi + n - 1
!      
!      ! piz(n_) = piz(n_) + 2 * sqrt(stdv_p(n_)) *  sqrt(duzdz(n_))
!      ! pit(n_) = pit(n_) + 2 * mes_D%r(n_,-1) * sqrt(stdv_p(n_)) *  sqrt(dutdt(n_))
!      ! pir(n_) = pir(n_) + 2 * sqrt(stdv_p(n_)) *  sqrt(durdr(n_))
!      ! enddo
!
!
!
!! ! Z derivatives of squared velocity
!
!!          !dzduzsq
!
!!          _loop_km_begin
!!                tmpr1 = (vel_uz%Re(:,nh)**2 - vel_uz%Im(:,nh)**2) *ad_k1a1(k)
!!                tmpr2 = (2*vel_uz%Re(:,nh)*vel_uz%Im(:,nh)) *ad_k1a1(k)
!
!!                ! Operate
!!                factor = 2d0
!!                if (m==0) factor = 1d0
!!                dzduzsq(:) = dzduzsq(:) + factor*(tmpr1(:)**2+tmpr2(:)**2)
!
!!          _loop_km_end
!
!!          !dzdutsq
!
!!          _loop_km_begin
!!                tmpr1 = (vel_ut%Re(:,nh)**2 - vel_ut%Im(:,nh)**2) *ad_k1a1(k)
!!                tmpr2 = (2*vel_ut%Re(:,nh)*vel_ut%Im(:,nh)) *ad_k1a1(k)
!
!!                ! Operate
!!                factor = 2d0
!!                if (m==0) factor = 1d0
!!                dzdutsq(:) = dzdutsq(:) + factor*(tmpr1(:)**2+tmpr2(:)**2)
!
!!          _loop_km_end
!
!
!!          !dzdursq
!
!!          _loop_km_begin
!!                tmpr1 = (vel_ur%Re(:,nh)**2 - vel_ur%Im(:,nh)**2) *ad_k1a1(k)
!!                tmpr2 = (2*vel_ur%Re(:,nh)*vel_ur%Im(:,nh)) *ad_k1a1(k)
!
!!                ! Operate
!!                factor = 2d0
!!                if (m==0) factor = 1d0
!!                dzdursq(:) = dzdursq(:) + factor*(tmpr1(:)**2+tmpr2(:)**2)
!
!!          _loop_km_end
!
!
!! ! Theta derivatives of squared velocity
!
!
!!          !dtduzsq
!
!!          _loop_km_begin
!!                tmpr1 = (vel_uz%Re(:,nh)**2 - vel_uz%Im(:,nh)**2) *ad_m1r1(k)
!!                tmpr2 = (2*vel_uz%Re(:,nh)*vel_uz%Im(:,nh)) *ad_m1r1(k)
!
!!                ! Operate
!!                factor = 2d0
!!                if (m==0) factor = 1d0
!!                dtduzsq(:) = dtduzsq(:) + factor*(tmpr1(:)**2+tmpr2(:)**2)
!
!!          _loop_km_end
!
!!          !dtdutsq
!
!!          _loop_km_begin
!!                tmpr1 = (vel_ut%Re(:,nh)**2 - vel_ut%Im(:,nh)**2) *ad_m1r1(k)
!!                tmpr2 = (2*vel_ut%Re(:,nh)*vel_ut%Im(:,nh)) *ad_m1r1(k)
!
!!                ! Operate
!!                factor = 2d0
!!                if (m==0) factor = 1d0
!!                dtdutsq(:) = dtdutsq(:) + factor*(tmpr1(:)**2+tmpr2(:)**2)
!
!!          _loop_km_end
!
!
!!          !dtdursq
!
!!          _loop_km_begin
!!                tmpr1 = (vel_ur%Re(:,nh)**2 - vel_ur%Im(:,nh)**2) *ad_m1r1(k)
!!                tmpr2 = (2*vel_ur%Re(:,nh)*vel_ur%Im(:,nh)) *ad_m1r1(k)
!
!!                ! Operate
!!                factor = 2d0
!!                if (m==0) factor = 1d0
!!                dtdursq(:) = dtdursq(:) + factor*(tmpr1(:)**2+tmpr2(:)**2)
!
!!          _loop_km_end
!   
!! r derivatives of squared velocity -> Matlab
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!
!   
!!      !dzduzcub
!
!!    call var_coll_copy(vel_uz,c1)
!
!!    _loop_km_begin
!!    tmpr1 = (c1%Re(:,nh)**3 - 3*c1%Im(:,nh)**2*c1%Re(:,nh)) *ad_k1a1(k)
!!    tmpr2 = (3*c1%Im(:,nh)*c1%Re(:,nh)**2 - c1%Im(:,nh)**3) *ad_k1a1(k)
!
!!          factor = 2d0
!!       if (m==0) factor = 1d0
!
!!       dzduzcub(:) = dzduzcub(:) + factor*(tmpr1(:)**2+tmpr2(:)**2)
!
!
!!    _loop_km_end 
!
!
! end subroutine var_coll_dissp








subroutine initialiseSTD()

implicit none
   csta = 1
   utaum   = 0d0
   uclm    = 0d0
   mean_ur = 0d0
   stdv_ur = 0d0
   mean_ut = 0d0
   stdv_ut = 0d0
   mean_uz = 0d0
   stdv_uz = 0d0
   stdv_rz = 0d0
   stdv_rt = 0d0
   stdv_tz = 0d0

   ! mom_ur = 0d0
   ! mom_uz = 0d0
   ! mom_ut = 0d0

   !dissr = 0d0
   !disst = 0d0
   !dissz = 0d0

   !urf = 0d0
   mean_p = 0d0
   stdv_p = 0d0

   !diss   = 0d0

   durdr = 0d0
   durdt = 0d0
   durdz = 0d0
   dutdr = 0d0
   dutdt = 0d0
   dutdz = 0d0
   duzdr = 0d0   
   duzdt = 0d0
   duzdz = 0d0

!
   !pir  = 0d0
   !pit  = 0d0
   !piz  = 0d0

   PDT2 = 0d0
   TDT2 = 0d0
   DT1 = 0d0
   DT4 = 0d0
   DT5 = 0d0
   DT6 = 0d0
   DT7 = 0d0
   



end subroutine initialiseSTD


subroutine saveStats(fnameima)
implicit none

    character(len = 256):: fnameima
    call mpi_reduce(mean_ur, dd, i_N*n_sta, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    mean_ur = dd
    call mpi_reduce(stdv_ur, dd, i_N*n_sta, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    stdv_ur = dd
    call mpi_reduce(mean_ut, dd, i_N*n_sta, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    mean_ut = dd
    call mpi_reduce(stdv_ut, dd, i_N*n_sta, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    stdv_ut = dd
    call mpi_reduce(mean_uz, dd, i_N*n_sta, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    mean_uz = dd
    call mpi_reduce(stdv_uz, dd, i_N*n_sta, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    stdv_uz = dd
    call mpi_reduce(stdv_rz, dd, i_N*n_sta, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    stdv_rz = dd
    call mpi_reduce(stdv_rt, dd, i_N*n_sta, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    stdv_rt = dd
    call mpi_reduce(stdv_tz, dd, i_N*n_sta, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    stdv_tz = dd

    call mpi_reduce(mean_p, dd, i_N*n_sta, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    mean_p = dd
    call mpi_reduce(stdv_p, dd, i_N*n_sta, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    stdv_p = dd

   !! Dissipation
   !    diss(:,1) = dissr(:,1) + disst(:,1) + dissz(:,1)
   !    diss(:,2) = dissr(:,2) + disst(:,2) + dissz(:,2)
   !    diss(:,3) = dissr(:,3) + disst(:,3) + dissz(:,3)

   ! call mpi_reduce(diss, dissr, 3*i_N, mpi_double_precision,  & ! Usammos dissr que ya está usada y el no. de elementos a escribir es 3*i_N
   !     mpi_sum, 0, mpi_comm_world, mpi_er)
   !  diss = dissr

   !call mpi_reduce(diss(1,1), d, i_N, mpi_double_precision,  &
   !    mpi_sum, 0, mpi_comm_world, mpi_er)
   ! diss(:,1) = d

   ! call mpi_reduce(diss(1,2), d, i_N, mpi_double_precision,  &
   ! mpi_sum, 0, mpi_comm_world, mpi_er)
   ! diss(:,2) = d
   ! call mpi_reduce(diss(1,3), d, i_N, mpi_double_precision,  &
   ! mpi_sum, 0, mpi_comm_world, mpi_er)
   ! diss(:,3) = d

   !  call mpi_reduce(piz, d, i_N, mpi_double_precision,  &
   !    mpi_sum, 0, mpi_comm_world, mpi_er)
   ! piz = d
   !   call mpi_reduce(pit, d, i_N, mpi_double_precision,  &
   !    mpi_sum, 0, mpi_comm_world, mpi_er)
   ! pit = d
   !   call mpi_reduce(pir, d, i_N, mpi_double_precision,  &
   !    mpi_sum, 0, mpi_comm_world, mpi_er)
   ! pir = d

      call mpi_reduce(PDT2, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    PDT2 = d
      call mpi_reduce(TDT2, d, i_N, mpi_double_precision,  &
       mpi_sum, 0, mpi_comm_world, mpi_er)
    TDT2 = d
      call mpi_reduce(DT1, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
    DT1 = dd
      call mpi_reduce(DT4, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
    DT4 = dd
      call mpi_reduce(DT5, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
    DT5 = dd
      call mpi_reduce(DT6, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
    DT6 = dd
      call mpi_reduce(DT7, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
    DT7 = dd
    
      call mpi_reduce(durdr, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
   durdr = dd
      call mpi_reduce(durdt, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
   durdt = dd
      call mpi_reduce(durdz, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
   durdz = dd
      call mpi_reduce(dutdr, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
   dutdr = dd
      call mpi_reduce(dutdt, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
   dutdt = dd
      call mpi_reduce(dutdz, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
   dutdz = dd
      call mpi_reduce(duzdr, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
   duzdr = dd
      call mpi_reduce(duzdt, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
   duzdt = dd
      call mpi_reduce(duzdz, dd, i_N*n_sta, mpi_double_precision,  &
      mpi_sum, 0, mpi_comm_world, mpi_er)
   duzdz = dd

    

   !   call mpi_reduce(uzsqur, d, i_N, mpi_double_precision,  &
   !    mpi_sum, 0, mpi_comm_world, mpi_er)
   !   uzsqur = d
   !   call mpi_reduce(utsqur, d, i_N, mpi_double_precision,  &
   !    mpi_sum, 0, mpi_comm_world, mpi_er)
   !   utsqur = d
   !   call mpi_reduce(urcub, d, i_N, mpi_double_precision,  &
   !    mpi_sum, 0, mpi_comm_world, mpi_er)
   !   urcub = d






   ! call mpi_reduce(urf, d, i_N, mpi_double_precision,  &
   !     mpi_sum, 0, mpi_comm_world, mpi_er)
   !  urf = d

   ! or, ot, oz: vort's
   ! call mpi_reduce(or, d, i_N, mpi_double_precision,  &
   !     mpi_sum, 0, mpi_comm_world, mpi_er)
   !  or = d
   ! call mpi_reduce(ot, d, i_N, mpi_double_precision,  &
   !     mpi_sum, 0, mpi_comm_world, mpi_er)
   !  ot = d
   ! call mpi_reduce(oz, d, i_N, mpi_double_precision,  &
   !     mpi_sum, 0, mpi_comm_world, mpi_er)
   !  oz = d

   ! call mpi_reduce(dzduzsq, d, i_N, mpi_double_precision,  &
   !     mpi_sum, 0, mpi_comm_world, mpi_er)
   !  dzduzsq = d
   !     call mpi_reduce(dzduzcub, d, i_N, mpi_double_precision,  &
   !     mpi_sum, 0, mpi_comm_world, mpi_er)
   !  dzduzcub = d

   !  call mpi_reduce(mom_ur, dd, i_N*10, mpi_double_precision,  &
   !     mpi_sum, 0, mpi_comm_world, mpi_er)
   !  mom_ur = dd
   !  call mpi_reduce(mom_uz, dd, i_N*10, mpi_double_precision,  &
   !     mpi_sum, 0, mpi_comm_world, mpi_er)
   !  mom_uz = dd
   !  call mpi_reduce(mom_ut, dd, i_N*10, mpi_double_precision,  &
   !     mpi_sum, 0, mpi_comm_world, mpi_er)
   !  mom_ut = dd

    if(mpi_rnk==0)  then

       call h5fcreate_f(trim(fnameima),H5F_ACC_TRUNC_F,fid,h5err)
       call h5gcreate_f(fid, '/header', header_id, h5err)
       call h5gcreate_f(fid, '/sta'   , sta_id   , h5err)
       call h5gcreate_f(fid, '/budget', budget_id ,h5err)

       hdims = (/1/)
       call h5ltmake_dataset_double_f(header_id,"time",1,hdims,(/tim_t/),h5err)
       call h5ltmake_dataset_double_f(header_id,"Re",1,hdims,(/d_Re/),h5err)
       call h5ltmake_dataset_double_f(header_id,"alpha",1,hdims,(/d_alpha/),h5err)

       call h5ltmake_dataset_int_f(header_id,"N" ,1,hdims,(/i_N/),h5err)
       call h5ltmake_dataset_int_f(header_id,"num" ,1,hdims,(/csta/),h5err)
       call h5ltmake_dataset_int_f(header_id,"M" ,1,hdims,(/i_M/),h5err)
       call h5ltmake_dataset_int_f(header_id,"K" ,1,hdims,(/i_K/),h5err)
       call h5ltmake_dataset_int_f(header_id,"Mp",1,hdims,(/i_Mp/),h5err)

       call h5ltmake_dataset_double_f(sta_id,"utau",1,hdims,utaum,h5err)
       call h5ltmake_dataset_double_f(sta_id,"ucl",1,hdims,uclm,h5err)

       call h5ltmake_dataset_double_f(header_id,"dt"   ,1,hdims,(/tim_dt/),h5err)
       hdims = (/i_N/)
       call h5ltmake_dataset_double_f(header_id,"r"   ,1,hdims,mes_D%r(1:i_N,1),h5err)
      
       hdims2 = (/i_N,n_sta/)
       call h5ltmake_dataset_double_f(sta_id,"mean_ur",2,hdims2,mean_ur,h5err)
       call h5ltmake_dataset_double_f(sta_id,"mean_uz",2,hdims2,mean_uz,h5err)
       call h5ltmake_dataset_double_f(sta_id,"mean_ut",2,hdims2,mean_ut,h5err)
2
       call h5ltmake_dataset_double_f(sta_id,"stdv_ur",2,hdims2,stdv_ur,h5err)
       call h5ltmake_dataset_double_f(sta_id,"stdv_ut",2,hdims2,stdv_ut,h5err)
       call h5ltmake_dataset_double_f(sta_id,"stdv_uz",2,hdims2,stdv_uz,h5err)
       call h5ltmake_dataset_double_f(sta_id,"stdv_rz",2,hdims2,stdv_rz,h5err)
       call h5ltmake_dataset_double_f(sta_id,"stdv_rt",2,hdims2,stdv_rt,h5err)
       call h5ltmake_dataset_double_f(sta_id,"stdv_tz",2,hdims2,stdv_tz,h5err)

       call h5ltmake_dataset_double_f(sta_id,"stdv_p",2,hdims2,stdv_p,h5err)
       call h5ltmake_dataset_double_f(sta_id,"mean_p",2,hdims2,mean_p,h5err)

       
      ! call h5ltmake_dataset_double_f(sta_id,"disstt",1,hdims,diss(:,2),h5err)
      ! call h5ltmake_dataset_double_f(sta_id,"dissrr",1,hdims,diss(:,1),h5err)
      ! call h5ltmake_dataset_double_f(sta_id,"disszz",1,hdims,diss(:,3),h5err) 

      ! call h5ltmake_dataset_double_f(sta_id,"pir",1,hdims,pir,h5err)
      ! call h5ltmake_dataset_double_f(sta_id,"piz",1,hdims,piz,h5err)
      ! call h5ltmake_dataset_double_f(sta_id,"pit",1,hdims,pit,h5err) 

       call h5ltmake_dataset_double_f(budget_id,"pdt2",2,hdims2,PDT2,h5err)
       call h5ltmake_dataset_double_f(budget_id,"tdt2",2,hdims2,TDT2,h5err)
       call h5ltmake_dataset_double_f(budget_id,"dt1",2,hdims2,DT1,h5err)
       call h5ltmake_dataset_double_f(budget_id,"dt4",2,hdims2,DT4,h5err)
       call h5ltmake_dataset_double_f(budget_id,"dt5",2,hdims2,DT5,h5err)
       call h5ltmake_dataset_double_f(budget_id,"dt6",2,hdims2,DT6,h5err)
       call h5ltmake_dataset_double_f(budget_id,"dt7",2,hdims2,DT7,h5err)
      
      call h5ltmake_dataset_double_f(budget_id,"dt7",2,hdims2,durdr,h5err)
      call h5ltmake_dataset_double_f(budget_id,"dt7",2,hdims2,durdt,h5err)
      call h5ltmake_dataset_double_f(budget_id,"dt7",2,hdims2,durdz,h5err)
      call h5ltmake_dataset_double_f(budget_id,"dt7",2,hdims2,dutdr,h5err)
      call h5ltmake_dataset_double_f(budget_id,"dt7",2,hdims2,dutdt,h5err)
      call h5ltmake_dataset_double_f(budget_id,"dt7",2,hdims2,dutdz,h5err)
      call h5ltmake_dataset_double_f(budget_id,"dt7",2,hdims2,duzdr,h5err)
      call h5ltmake_dataset_double_f(budget_id,"dt7",2,hdims2,duzdt,h5err)
      call h5ltmake_dataset_double_f(budget_id,"dt7",2,hdims2,duzdz,h5err)



       !call h5ltmake_dataset_double_f(sta_id,"uzsqur",1,hdims,uzsqur,h5err)
       !call h5ltmake_dataset_double_f(sta_id,"utsqur",1,hdims,utsqur,h5err)
       !call h5ltmake_dataset_double_f(sta_id,"urcub",1,hdims,urcub,h5err)      
       
         

       

      !  call h5ltmake_dataset_double_f(sta_id,"mom_ur",2,hdims2,mom_ur,h5err)
      !  call h5ltmake_dataset_double_f(sta_id,"mom_uz",2,hdims2,mom_uz,h5err)
      !  call h5ltmake_dataset_double_f(sta_id,"mom_ut",2,hdims2,mom_ut,h5err)

       call h5gclose_f(header_id,h5err)
       call h5gclose_f(sta_id,h5err)
       call h5gclose_f(budget_id,h5err)

       call h5fclose_f(fid,h5err)
   endif
   ! Once saved, initialise e verything to 0 again. 
   call initialiseSTD()

end subroutine saveStats

end module sta