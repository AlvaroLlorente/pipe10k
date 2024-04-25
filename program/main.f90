!*************************************************************************
! main.f90 (executable)
!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
 PROGRAM MAIN
!*************************************************************************
   ! use type
   use wksp
   use velocity
   use io
   use sta, only:compute_sta, initialisestd, savestats
   use parameters
   implicit none
   real :: d_start, d_stop
   double precision :: steptimer = 0d0, retau

   call h5open_f(h5err)
   call initialise()

  
   do while(.not.terminate())
      

       if(mpi_rnk==0) steptimer=MPI_Wtime()

      call vel_imposesym()
      call vel_transform()
      call vel_nonlinear()
      call var_null(2)
      
      if(d_timestep<0d0) then
         call vel_maxtstep()
         call tim_new_tstep()
         if(tim_new_dt)  &
            call vel_matrices()
      end if        
      call var_null(1)

      ! Estadisticas on the fly aquí!!!!
      if (tim_step>1) then
         if (mod(tim_step,s_step)==1) then
            call compute_sta()
         endif
         call io_write2files()
         csta = csta + 1
      endif
      
      call vel_predictor()
      tim_it = 1
      
      ! Cae en este bucle, que es el predictor corrector
      do while(tim_it/=0)
         call vel_transform()
         call vel_nonlinear()
         call var_null(2)
         call vel_corrector(r,t,z)
         call tim_check_cgce()
      end do

      tim_t    = tim_t    + tim_dt
      tim_step = tim_step + 1

      if (mod(tim_step,s_step)==1) then
         steptimer = MPI_Wtime()-steptimer
         call io_write_friction(tim_step,tim_t,steptimer,retau)
         if (mpi_rnk==0) then 
325         format(a4,i5,a14,d14.6,a8,d14.6,a7,d14.6)
            write(*,325) 'Step',tim_step,' ,t simulation',tim_t,' ,s/step',steptimer,' ,Retau',retau
         endif
      endif
      if (tim_step>f_step) exit
         

   end do
   ! end main loop   

   call cleanup()
   stop

 contains


!-------------------------------------------------------------------------
!  Termaination conditions
!-------------------------------------------------------------------------
   logical function terminate()
      logical :: file_exist
            
      if(mpi_rnk==0) then
         terminate = .false.
      
         if(tim_step==i_maxtstep) then
            terminate = .true.
            print*, 'maxtstep reached!'
         end if

         if(d_maxt>0d0 .and. tim_t>=d_maxt) then
            terminate = .true.
            print*, 'maxt reached!'
         end if

         call clk_time(d_stop)
         if(dble(d_stop-d_start)/36d9 >= d_cpuhours) then
            terminate = .true.
            print*, 'cpuhours reached!'
         end if

         if( modulo(tim_step,i_save_rate2)==0) then
            inquire(file='RUNNING', exist=file_exist)
            if(.not. file_exist) then
               terminate = .true.
               print*, 'RUNNING deleted !'
            end if
         end if
      end if

#ifdef _MPI
      call mpi_bcast(terminate,1,mpi_logical, 0,mpi_comm_world,mpi_er)
#endif

   end function terminate


!-------------------------------------------------------------------------
!  Initialisation
!-------------------------------------------------------------------------
   subroutine initialise()
      logical :: file_exist

      call mpi_precompute()
      if(mpi_rnk==0) then
         call system('touch PRECOMPUTING')
         call system('echo $HOSTNAME > HOST')
      end if

      if(mpi_rnk==0)  print*, 'precomputing function requisites...'

      call initwksp()

      call par_precompute()
      call mes_precompute()
      call var_precompute()
      call tra_precompute() !T his routine takes a lot. 
      call tim_precompute()
      call vel_precompute()
      call  io_precompute()
      call  initialiseSTD()
   
      if(mpi_rnk==0)  print*, 'loading state...'
      tim_dt = 1d99
      call io_load_state()
      call vel_matrices()
      
      if(mpi_rnk==0)  print*, 'initialising output files...'
      call io_openfiles()

      if(mpi_rnk==0) then
         open (99, file='PRECOMPUTING')
         close(99, status='delete')
         open (99, file='RUNNING')
         write(99,*) 'This file indicates a job running in this directory.'
         write(99,*) 'Delete this file to cleanly terminate the process.'
         close(99)
         print*, 'timestepping.....'
      end if
      
      call clk_time(d_start)
   
   end subroutine initialise

!-------------------------------------------------------------------------
!  Program closure
!-------------------------------------------------------------------------
   subroutine cleanup()
      logical :: file_exist
   
      if(mpi_rnk==0) then
         print*, 'cleanup...'
         call clk_time(d_stop)
         print*, ' sec/step (simulation)  = ', (d_stop-d_start)/real(tim_step)
#ifdef _CPUTIME
         print*, ' CPU time  = ', int((d_stop-d_start)/6d1), ' mins.'
#else
         print*, ' WALL time = ', int((d_stop-d_start)/6d1), ' mins.'
#endif 
      end if
      
      !!!!call io_save_state()
      !!!!call io_closefiles()
      
      call cleanwksp()

      if(mpi_rnk==0) then
         inquire(file='RUNNING', exist=file_exist)
         if(file_exist) open(99, file='RUNNING')
         if(file_exist) close(99, status='delete')
      end if      

#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
      call mpi_finalize(mpi_er)
#endif
      if(mpi_rnk==0) print*, '...done!'

   end subroutine cleanup


!-------------------------------------------------------------------------
   subroutine clk_time(t)
      real, intent(out) :: t
#ifdef _CPUTIME
      call cpu_time(t)
#else
      integer, save :: ct,ctrt,ctmx,ct_=0,wrap=0
      call system_clock(ct,ctrt,ctmx)
      if(ct<ct_) wrap = wrap + 1
      ct_ = ct
      t = (real(ct)+real(ctmx)*wrap)/real(ctrt)
#endif
   end subroutine clk_time
         
!*************************************************************************
 END PROGRAM MAIN
!*************************************************************************

