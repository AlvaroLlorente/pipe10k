
! 1) Otf, max = 54.4
! 2) max =  52.1
! 3) max =  51.0
! 4) max =  49
! 5) max =  48.5
! 6) max =  48.4
! 7) max =  
! 8) max =  
! 9) max =  
!Codigo optimizado

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


   double precision, private :: ucl, utau, Ub
 
 

! ------------------------- stats  -------------------------------

   double precision :: mean_ur(i_pN,n_sta), stdv_ur(i_pN,n_sta)
   double precision :: mean_ut(i_pN,n_sta), stdv_ut(i_pN,n_sta)
   double precision :: mean_uz(i_pN,n_sta), stdv_uz(i_pN,n_sta)
   double precision :: stdv_rz(i_pN,n_sta), stdv_rt(i_pN,n_sta), stdv_tz(i_pN,n_sta)
   double precision :: stdv_zzr(i_pN,n_sta),stdv_rtt(i_pN,n_sta),stdv_rrr(i_pN,n_sta)
  
   double precision :: mean_p(i_N,n_sta), stdv_p(i_N,n_sta)


   double precision :: durdrsq(i_N,n_sta), durdzsq(i_N,n_sta)
   double precision :: dutdrsq(i_N,n_sta),dutdzsq(i_N,n_sta)
   double precision :: duzdrsq(i_N,n_sta),duzdzsq(i_N,n_sta)
   double precision :: uuPST1(i_N,n_sta), uuDT3(i_N,n_sta)
   double precision :: ttPST1(i_N,n_sta), ttDT31(i_N,n_sta)
   double precision :: rrDT31(i_N,n_sta), rrPST1(i_N,n_sta), rrPDT1(i_N,n_sta)
   double precision :: kVDT1(i_N,n_sta),kTDT2(i_N,n_sta),kDT4(i_N,n_sta),kDT5(i_N,n_sta)
   double precision :: kDT61(i_N,n_sta), kDT62(i_N,n_sta), kDT63(i_N,n_sta),kDT73(i_N,n_sta)
   double precision :: time
   double precision :: time_sta(n_sta), utauv(n_sta), uclv(n_sta), dt(n_sta)

   integer :: csta
   
! ------------------------- HDF5 -------------------------------

   integer:: info,ierr
   integer(hid_t):: fid,pid, dset_id,dspace_id
   integer:: h5err
   integer(hsize_t),dimension(3):: dims
   integer(hsize_t),dimension(1):: hdims,maxdims
   integer(hsize_t),dimension(2):: hdims2
   
   !integer(hid_t) :: header_id,sta_id, budget_id, derivative_id ! Group identifiers
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
      time= tim_t + tim_dt
      time_sta(csta) = time 
      ucl = 1d0 + dot_product(vel_uz%Re(1:1+i_KL,0),mes_D%dr0(:,0))
      utau = dot_product(vel_uz%Re(i_N-i_KL:i_N,0),mes_D%dr1(:,1))
      utau = dsqrt(dabs((Utau-2d0)/d_Re))

      ! add statistics

      dt(csta)=tim_dt
      uclv(csta)=ucl
      utauv(csta) = utau 
   endif

      
      call vel_sta()

!      call pressure(c1,c2,c3,p1,p2)
         


   do n = 1, mes_D%pN
      !n_ = mes_D%pNi + n - 1
      mean_ur(n,csta) = mean_ur(n,csta) + sum(vel_r%Re(:,:,n))
      stdv_ur(n,csta) = stdv_ur(n,csta) + sum(vel_r%Re(:,:,n)**2)
      mean_ut(n,csta) = mean_ut(n,csta) + sum(vel_t%Re(:,:,n))
      stdv_ut(n,csta) = stdv_ut(n,csta) + sum(vel_t%Re(:,:,n)**2)
      mean_uz(n,csta) = mean_uz(n,csta) + sum(vel_z%Re(:,:,n))
      stdv_uz(n,csta) = stdv_uz(n,csta) + sum(vel_z%Re(:,:,n)**2)
      stdv_rz(n,csta) = stdv_rz(n,csta) + sum(vel_z%Re(:,:,n)*vel_r%Re(:,:,n))
      stdv_rt(n,csta) = stdv_rt(n,csta) + sum(vel_r%Re(:,:,n)*vel_t%Re(:,:,n))
      stdv_tz(n,csta) = stdv_tz(n,csta) + sum(vel_t%Re(:,:,n)*vel_z%Re(:,:,n))

      stdv_zzr(n,csta) = stdv_zzr(n,csta) + sum(vel_z%Re(:,:,n)*vel_z%Re(:,:,n)*vel_r%Re(:,:,n))
      stdv_rtt(n,csta) = stdv_rtt(n,csta) + sum(vel_r%Re(:,:,n)*vel_t%Re(:,:,n)*vel_t%Re(:,:,n))
      stdv_rrr(n,csta) = stdv_rrr(n,csta) + sum(vel_r%Re(:,:,n)**3)
   enddo

   !call compute_turb_budget()

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
         !durdr 
         call var_coll_meshmult(1,mes_D%dr(1),vel_ur,c1) !durdr
         call tra_coll2phys1d(c1,p1)


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
   
call var_coll_meshmult(1,mes_D%dr(1),vel_ur,c1) !durdr de todo el campo
call tra_coll2phys1d(c1,p1)

_loop_km_begin

 c4%Re(:,nh) = -vel_ur%Im(:,nh)*d_alpha*k
 c4%Im(:,nh) =  vel_ur%Re(:,nh)*d_alpha*k


_loop_km_end


call tra_coll2phys1d(c4,p4) !durdz

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
   durdrsq(n_,csta) = durdrsq(n_,csta) + sum(p1%Re(:,:,n)**2)
   durdzsq(n_,csta) = durdzsq(n_,csta) + sum(p4%Re(:,:,n)**2)  
end do


!   vel_t

call var_coll_meshmult(1,mes_D%dr(1),vel_ut,c1) !dutdr
call tra_coll2phys1d(c1,p1)
 _loop_km_begin

 c4%Re(:,nh) = -vel_ut%Im(:,nh)*d_alpha*k
 c4%Im(:,nh) =  vel_ut%Re(:,nh)*d_alpha*k


_loop_km_end


call tra_coll2phys1d(c4,p4) !dutdz

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
   dutdrsq(n_,csta) = dutdrsq(n_,csta) + sum(p1%Re(:,:,n)**2)    
   dutdzsq(n_,csta) = dutdzsq(n_,csta) + sum(p4%Re(:,:,n)**2) 
end do

!   vel_z


call var_coll_meshmult(0,mes_D%dr(1),vel_uz,c1) !duzdr
call tra_coll2phys1d(c1,p1)
_loop_km_begin

 c4%Re(:,nh) = -vel_uz%Im(:,nh)*d_alpha*k
 c4%Im(:,nh) =  vel_uz%Re(:,nh)*d_alpha*k

_loop_km_end

call tra_coll2phys1d(c4,p4) !duzdz

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
   duzdrsq(n_,csta) = duzdrsq(n_,csta) + sum(p1%Re(:,:,n)**2)
   duzdzsq(n_,csta) = duzdzsq(n_,csta) + sum(p4%Re(:,:,n)**2)


end do







!--------Turbulent kinetic energy budget-------!!




!  Pressure difussion term

!Es cero

!  Turbulent difussion term 



p1%Re=vel_r%Re*(vel_r%Re**2+vel_t%Re**2+vel_z%Re**2)

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
   kTDT2(n_,csta)  = kTDT2(n_,csta)  + sum(p1%Re(:,:,n)) ! saco la distribucion radial
end do

! Viscous diffusion term

_loop_km_begin
c2%Re(:,nh)=vel_ur%Re(:,nh)*vel_uz%Re(:,nh)*mes_D%r(:,1)
c2%Im(:,nh)=vel_ur%Re(:,nh)*vel_uz%Im(:,nh)*mes_D%r(:,1)
_loop_km_end

call var_coll_meshmult(0,mes_D%dr(1),c2,c1)

_loop_km_begin
c2%Re(:,nh)=(c1%Re(:,nh)*mes_D%r(:,-1))*d_alpha*k
c2%Im(:,nh)=(c1%Im(:,nh)*mes_D%r(:,-1))*d_alpha*k
_loop_km_end

call tra_coll2phys1d(c2,p1)

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
kVDT1(n_,csta)=kVDT1(n_,csta)+sum(p1%Re(:,:,n))
enddo

!  Dissipation term 

call var_coll_meshmult(0,mes_D%dr(1),vel_uz,c3)

_loop_km_begin
c1%Re(:,nh) = (-vel_ut%Im(:,nh)*m*i_Mp)*mes_D%r(:,-1)
c1%Im(:,nh) =  (vel_ut%Re(:,nh)*m*i_Mp)*mes_D%r(:,-1)

c4%Re(:,nh) = (-vel_ur%Im(:,nh)*d_alpha*k)*c3%Re(:,nh)
c4%Im(:,nh) =  (vel_ur%Re(:,nh)*d_alpha*k)*c3%Im(:,nh)
_loop_km_end

call tra_coll2phys1d(c1,p1) !1/r(dutdt), termino 4
call tra_coll2phys1d(c4,p4) !Termino 5

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
kDT4(n_,csta)=kDT4(n_,csta)+sum(p1%Re(:,:,n)**2)
kDT5(n_,csta)=kDT5(n_,csta)+2*sum(p4%Re(:,:,n))
enddo
!------

_loop_km_begin
c3%Re(:,nh) = vel_ut%Re(:,nh)*mes_D%r(:,-1)
c3%Im(:,nh) = vel_ut%Im(:,nh)*mes_D%r(:,-1)
_loop_km_end

call var_coll_meshmult(1,mes_D%dr(1),c3,c1)

_loop_km_begin
c3%Re(:,nh) = (-vel_ur%Im(:,nh)*m*i_Mp)
c3%Im(:,nh) = (vel_ur%Re(:,nh)*m*i_Mp)
_loop_km_end
call tra_coll2phys1d(c1,p1) !provisional para multiplicar por p3
call tra_coll2phys1d(c3,p3) !Termino 6_3

p3%Re=p3%Re*p1%Re

_loop_km_begin
c1%Re(:,nh) = (-vel_ur%Im(:,nh)*m*i_Mp)*mes_D%r(:,-1)
c1%Im(:,nh) = (vel_ur%Re(:,nh)*m*i_Mp)*mes_D%r(:,-1)
_loop_km_end

call tra_coll2phys1d(c1,p1) !Termino 6_1

_loop_km_begin
c1%Re(:,nh) = vel_ut%Re(:,nh)*mes_D%r(:,-1)
c1%Im(:,nh) = vel_ut%Im(:,nh)*mes_D%r(:,-1)
_loop_km_end

call var_coll_meshmult(1,mes_D%dr(1),c1,c2)

_loop_km_begin
c4%Re(:,nh) = c2%Re(:,nh)*mes_D%r(:,1)
c4%Im(:,nh) = c2%Im(:,nh)*mes_D%r(:,1)
_loop_km_end

call tra_coll2phys1d(c4,p4) !Termino 6_2


do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
 kDT61(n_,csta)=kDT61(n_,csta)+sum(p1%Re(:,:,n)**2) 
 kDT62(n_,csta)=kDT62(n_,csta)+sum(p4%Re(:,:,n)**2) 
 kDT63(n_,csta)=kDT63(n_,csta)+2*sum(p3%Re(:,:,n))
enddo

!kDT72 = uuDT3, no hace falta sacarlo otra vez
_loop_km_begin
c4%Re(:,nh) = (-vel_ut%Im(:,nh)*d_alpha*k)*mes_D%r(:,-1)
c4%Im(:,nh) = (vel_ut%Re(:,nh)*d_alpha*k)*mes_D%r(:,-1)

c1%Re(:,nh) = (-vel_uz%Im(:,nh)*m*i_Mp)
c1%Im(:,nh) = (vel_uz%Re(:,nh)*m*i_Mp)
_loop_km_end

call tra_coll2phys1d(c1,p1)
call tra_coll2phys1d(c4,p4) !Termino 7_3

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
 kDT73(n_,csta)=kDT73(n_,csta)+2*sum(p4%Re(:,:,n)*p1%Re(:,:,n))
enddo

!--------Reynolds normal stresses budget-------!!



!  Pressure strain term 

call var_coll_meshmult(1,mes_D%dr(1),vel_ur,c1)

_loop_km_begin
 c4%Re(:,nh) = -vel_uz%Im(:,nh)*d_alpha*k
 c4%Im(:,nh) =  vel_uz%Re(:,nh)*d_alpha*k

 c3%Re(:,nh) = -vel_ut%Im(:,nh)*m*i_Mp
 c3%Im(:,nh) =  vel_ut%Re(:,nh)*m*i_Mp
_loop_km_end

call tra_coll2phys1d(c4,p4) !duzdz uu
call tra_coll2phys1d(c3,p3) !dutdt tt
call tra_coll2phys1d(c1,p1) !durdr rr

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
uuPST1(n_,csta)=uuPST1(n_,csta)+sum(p4%Re(:,:,n)*p2%Re(:,:,n))

ttPST1(n_,csta)=ttPST1(n_,csta)+sum(p3%Re(:,:,n)*p2%Re(:,:,n))

rrPST1(n_,csta)=rrPST1(n_,csta)+sum(p1%Re(:,:,n)*p2%Re(:,:,n))

enddo

!  Pressure diffusion term 



do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
rrPDT1(n_,csta)=rrPDT1(n_,csta)+sum((p2%Re(:,:,n)*vel_r%Re(:,:,n))) ! P * Vr

enddo





!  Turbulent diffusion

!La parte de fortran es cero, derivadas en z

!  Dissipation term uzuz

_loop_km_begin
c3%Re(:,nh) = (-vel_uz%Im(:,nh)*m*i_Mp)*mes_D%r(:,-1)
c3%Im(:,nh) =  (vel_uz%Re(:,nh)*m*i_Mp)*mes_D%r(:,-1)
_loop_km_end

call tra_coll2phys1d(c3,p3) !1/r(duzdt)


do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1
uuDT3(n_,csta)=uuDT3(n_,csta)+sum(p3%Re(:,:,n)**2)

enddo

!  Dissipation term utut


_loop_km_begin   
c1%Re(:,nh) = mes_D%r(:,-1)*(-vel_ut%Im(:,nh)*m*i_Mp+vel_ur%Re(:,nh))
c1%Im(:,nh) = mes_D%r(:,-1)*( vel_ut%Re(:,nh)*m*i_Mp+vel_ur%Im(:,nh))
_loop_km_end


call tra_coll2phys1d(c1,p1)
call tra_coll2phys1d(c3,p3)
call tra_coll2phys1d(c4,p4)

do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1

ttDT31(n_,csta)=ttDT31(n_,csta)+sum(p1%Re(:,:,n)**2)   

enddo

!  Dissipation term urur


_loop_km_begin

c1%Re(:,nh) = mes_D%r(:,-1)*(-vel_ur%Im(:,nh)*m*i_Mp-vel_ut%Re(:,nh))
c1%Im(:,nh) = mes_D%r(:,-1)*( vel_ur%Re(:,nh)*m*i_Mp-vel_ut%Im(:,nh))

_loop_km_end

call tra_coll2phys1d(c1,p1)


do n = 1, mes_D%pN
   n_ = mes_D%pNi + n - 1

rrDT31(n_,csta)=rrDT31(n_,csta)+sum(p1%Re(:,:,n)**2)


enddo

end subroutine compute_turb_budget











subroutine initialiseSTD()

implicit none
   csta = 1
   time_sta=0d0
   time=0d0
   dt=0d0
   uclv=0d0
   utauv=0d0
   
   mean_ur = 0d0
   stdv_ur = 0d0
   mean_ut = 0d0
   stdv_ut = 0d0
   mean_uz = 0d0
   stdv_uz = 0d0
   stdv_rz = 0d0
   stdv_rt = 0d0
   stdv_tz = 0d0
   stdv_zzr = 0d0
   stdv_rtt = 0d0
   stdv_rrr = 0d0
   
   mean_p = 0d0
   stdv_p = 0d0




   duzdrsq= 0d0
   dutdrsq= 0d0
   durdrsq= 0d0


   duzdzsq=0d0
   dutdzsq=0d0
   durdzsq= 0d0





   uuPST1=0d0
   uuDT3=0d0


   ttPST1=0d0
   ttDT31=0d0


   rrPST1=0d0
   rrPDT1=0d0
   rrDT31=0d0




   kTDT2=0d0
   kVDT1=0d0
   kDT4=0d0
   kDT5=0d0
   kDT61=0d0
   kDT62=0d0
   kDT63=0d0
   kDT73=0d0
  


end subroutine initialiseSTD


subroutine saveStats(fnameima)
implicit none
integer:: strow
character(len = 256):: fnameima
integer :: info
integer(hid_t) :: G1, G2, G3

strow = 1

info = MPI_INFO_NULL

call MPI_BARRIER(MPI_COMM_WORLD,ierr)

call h5pcreate_f(H5P_FILE_ACCESS_F,pid,h5err)
call h5pset_fapl_mpio_f(pid,MPI_COMM_WORLD,info,h5err)
call h5pset_sieve_buf_size_f(pid, siever, h5err)
call h5fcreate_f(trim(fnameima),H5F_ACC_TRUNC_F,fid,h5err,H5P_DEFAULT_F,pid)
call h5gcreate_f(fid, '/sta', G2, h5err)

call h5pclose_f(pid,h5err)

call MPI_BARRIER(MPI_COMM_WORLD,ierr)

hdims2 = (/i_pN,n_sta/)


call h5dump_parallel2D(G2,"mean_ur",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,mean_ur,h5err)
call h5dump_parallel2D(G2,"mean_uz",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,mean_uz,h5err)
call h5dump_parallel2D(G2,"mean_ut",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,mean_ut,h5err)

call h5dump_parallel2D(G2,"stdv_ur",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,stdv_ur,h5err)
call h5dump_parallel2D(G2,"stdv_ut",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,stdv_ut,h5err)
call h5dump_parallel2D(G2,"stdv_uz",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,stdv_uz,h5err)
call h5dump_parallel2D(G2,"stdv_rz",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,stdv_rz,h5err)
call h5dump_parallel2D(G2,"stdv_rt",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,stdv_rt,h5err)
call h5dump_parallel2D(G2,"stdv_tz",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,stdv_tz,h5err)
call h5dump_parallel2D(G2,"stdv_zzr",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,stdv_zzr,h5err)
call h5dump_parallel2D(G2,"stdv_rtt",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,stdv_rtt,h5err)
call h5dump_parallel2D(G2,"stdv_rrr",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,stdv_rrr,h5err)

call h5gclose_f(G2,h5err)
call h5fclose_f(fid,h5err)


    if(mpi_rnk==0)  then
 
      call h5fopen_f(trim(fnameima),H5F_ACC_RDWR_F,fid,h5err)
      call h5gopen_f(fid, '/sta', G2, h5err)
      call h5gcreate_f(fid, '/header', G1, h5err)
       
       hdims = (/1/)
       call h5ltmake_dataset_double_f(G1,"time",1,hdims,(/tim_t/),h5err)
       call h5ltmake_dataset_double_f(G1,"Re",1,hdims,(/d_Re/),h5err)
       call h5ltmake_dataset_double_f(G1,"alpha",1,hdims,(/d_alpha/),h5err)

       call h5ltmake_dataset_int_f(G1,"N" ,1,hdims,(/i_N/),h5err)
       call h5ltmake_dataset_int_f(G1,"num" ,1,hdims,(/csta/),h5err)
       call h5ltmake_dataset_int_f(G1,"M" ,1,hdims,(/i_M/),h5err)
       call h5ltmake_dataset_int_f(G1,"K" ,1,hdims,(/i_K/),h5err)
       call h5ltmake_dataset_int_f(G1,"Mp",1,hdims,(/i_Mp/),h5err)


       hdims = (/i_N/)
       call h5ltmake_dataset_double_f(G1,"r"   ,1,hdims,mes_D%r(1:i_N,1),h5err)
       hdims = (/n_sta/)
       call h5ltmake_dataset_double_f(G1,"timev",1,hdims,time_sta,h5err)
       call h5ltmake_dataset_double_f(G2,"utauv",1,hdims,utauv,h5err)
       call h5ltmake_dataset_double_f(G2,"uclv",1,hdims,uclv,h5err)
       call h5ltmake_dataset_double_f(G1,"dt",1,hdims,dt,h5err)

       call h5gclose_f(G1,h5err)
       call h5gclose_f(G2,h5err)
       call h5fclose_f(fid,h5err)
   
   endif

      !call h5gcreate_f(fid, '/derivatives', derivative_id ,h5err)
!
!
!
!
      ! call h5ltmake_dataset_double_f(G2,"stdv_p",2,hdims2,stdv_p,h5err)
      ! call h5ltmake_dataset_double_f(G2,"mean_p",2,hdims2,mean_p,h5err)
!
      ! 
   !
!
      !call h5ltmake_dataset_double_f(derivative_id,"duzdrsq",2,hdims2,duzdrsq,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"dutdrsq",2,hdims2,dutdrsq,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"durdrsq",2,hdims2,durdrsq,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"duzdzsq",2,hdims2,duzdzsq,h5err)       
      !call h5ltmake_dataset_double_f(derivative_id,"dutdzsq",2,hdims2,dutdzsq,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"durdzsq",2,hdims2,durdzsq,h5err)
      !  
!
!
      !call h5ltmake_dataset_double_f(derivative_id,"uuDT3",2,hdims2,uuDT3,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"uuPST1",2,hdims2,uuPST1,h5err)
!
!
      !call h5ltmake_dataset_double_f(derivative_id,"ttDT31",2,hdims2,ttDT31,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"ttPST1",2,hdims2,ttPST1,h5err)
!
!
!
      !call h5ltmake_dataset_double_f(derivative_id,"rrDT31",2,hdims2,rrDT31,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"rrPST1",2,hdims2,rrPST1,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"rrPDT1",2,hdims2,rrPDT1,h5err)
!
!
!
      !call h5ltmake_dataset_double_f(derivative_id,"kTDT2",2,hdims2,kTDT2,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"kVDT1",2,hdims2,kVDT1,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"kDT4",2,hdims2,kDT4,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"kDT5",2,hdims2,kDT5,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"kDT61",2,hdims2,kDT61,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"kDT62",2,hdims2,kDT62,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"kDT63",2,hdims2,kDT63,h5err)
      !call h5ltmake_dataset_double_f(derivative_id,"kDT73",2,hdims2,kDT73,h5err)
!


         


       
       !call h5gclose_f(G1,h5err)
       
       !call h5gclose_f(budget_id,h5err)
       !call h5gclose_f(derivative_id ,h5err)
       

  

   call initialiseSTD()

end subroutine saveStats


subroutine h5dump_parallel2D(fid, name, ndims, dims, strow, rank, size, comm, info, data, ierr)
   use h5lt
   use hdf5
   use mpif
   implicit none

   integer(hid_t), intent(in) :: fid
   character(len=*), intent(in) :: name
   integer, intent(in) :: ndims
   integer(hsize_t), dimension(ndims), intent(in) :: dims
   integer, intent(in) :: rank, size
   integer, intent(in) :: comm, info, strow
   real(kind=8), dimension(:,:), intent(in) :: data
   integer :: ierr

   integer(hid_t) :: dset
   integer(hid_t) :: dspace, mspace
   integer(hid_t) :: plist_id
   integer(hsize_t), dimension(ndims) :: start, nooffset, totaldims
   integer :: mpierr
   integer :: local_rows, global_rows, rows_per_rank

   ! Obtener el número total de filas
   global_rows = dims(1)
   ! Calcular el número de filas por proceso
   rows_per_rank = global_rows / size
   ! Calcular el número total de filas después de la distribución
   totaldims(1) = global_rows

   ! Calcular el rango de filas para el proceso actual
   start(1) = rank * rows_per_rank

   ! Calcular el número de filas para este proceso
   if (rank == size - 1) then
       ! El último proceso maneja las filas restantes
       local_rows = global_rows - start(1)
   else
       local_rows = rows_per_rank
   end if

   ! Create the global dataspace
   call h5screate_simple_f(ndims, totaldims, dspace, ierr)

   ! Create the global dataset
   call h5dcreate_f(fid, name, H5T_IEEE_F64BE, dspace, dset, ierr)

   ! Create the local dataspace
   call h5screate_simple_f(ndims, (/local_rows, dims(2)/), mspace, ierr)
   ! Select the hyperslab in the global dataset
   call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, start, (/local_rows, dims(2)/), ierr)

   ! Create data transfer mode property list
   call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
   call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)

   ! Commit the memspace to the disk
   call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, data, (/local_rows, dims(2)/), ierr, mspace, dspace, plist_id)

   ! Close property list
   call h5pclose_f(plist_id, ierr)

   ! Close datasets and dataspaces
   call h5sclose_f(mspace, ierr)
   call h5dclose_f(dset, ierr)
   call h5sclose_f(dspace, ierr)
end subroutine h5dump_parallel2D

end module sta
