!**************************************************************************
!  IN/OUT 
!
!**************************************************************************
#include "../parallel.h"
 !module wkarray
 !implicit none
 ! use velocity
 !double precision :: wk(i_N, 0:i_pH1) 
 !end module
 
 module io
!**************************************************************************
   use type
   use velocity
   use h5lt
   use hdf5
   use sta, only: saveStats

   implicit none
   save

   character(len=256) text
   character*4 extc
   character(len = 256):: fnameima,filinp, dirinp, filstt, fnamephys, filename, index


   integer               :: io_save2,extn, utau, indice
   integer,     private  :: io_KE, io_ID, io_dt, io_pt, io_fr, io_hre, io_cf
   
   !double precision, private :: wk(i_N, 0:i_pH1) ! PAsar luego por cabecera a la rutina o a través de un modulo.
   !double precision :: v1(0:i_Th-1, i_pN) 
   
   ! type (coll), private  :: c1!,c2,c3 ! Three colls are defined here. Why! They are really big. 
                                     ! They are defined as private. They cannot be used anywhere else
                                     ! Remove in future versions. We have to pass the routine some workarray 
                                     ! 

! STATS



! ------------------------- HDF5 -------------------------------

   integer:: info,ierr
   integer(hid_t):: fid,pid, dset_id,dspace_id
   integer:: h5err
   integer(hsize_t),dimension(3):: hdims3
   integer(hsize_t),dimension(1):: hdims,maxdims
   integer(hsize_t),dimension(2):: hdims2
   
   integer(hid_t) :: header_id,sta_id ! Group identifiers
   integer(size_t)::siever
   parameter (siever = 4*1024*1024)


! --------------------- Program -------------------------------------

 contains
 
!--------------------------------------------------------------------------
!  initialiser fn
!--------------------------------------------------------------------------
   subroutine io_precompute()
     implicit none


      io_save2 = 0
      io_dt    = 0
      io_hre   = 19
      io_KE    = 0
      io_ID    = 0
      io_pt    = 0
      io_fr    = 0
      io_cf    = 89

      indice   = 1

     if(mpi_rnk.eq.0) then
        open(19,file='hre.dat',status='old')
        
        call avanza(text)
        read(text,*) extn
        call avanza(text)
        read(text,'(a256)') dirinp
        call avanza(text)
        read(text,'(a256)') filinp
        call avanza(text)
        read(text,'(a256)') filstt
        
     endif
     call MPI_BCAST(extn,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
     call MPI_BCAST(filinp,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(dirinp,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(filstt,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
     write(extc,'(i4.4)') extn-1


   end subroutine io_precompute 
 

!--------------------------------------------------------------------------
!  Open files written to every ... steps runtime
!--------------------------------------------------------------------------
   subroutine io_openfiles()
      implicit none
      character(10), save :: a = 'sequential'

      if(mpi_rnk/=0) return
      
      open(io_cf,status='unknown', file=trim(filstt)//'.'//extc//'.cf')
      write(*,*)  trim(filstt)//'.cf'
      
      if(io_dt/=0)  then 
         fnameima = trim(filstt)//'tim_step.dat'
         open(io_dt,status='unknown',access=a, file=fnameima)
      endif

      if(io_KE/=0) open(io_KE,status='unknown',access=a, file=trim(filstt)//'.vel_energy')
      if(io_ID/=0) open(io_ID,status='unknown',access=a, file=trim(filstt)//'.vel_totEID')
      if(io_pt/=0) open(io_pt,status='unknown',access=a, file=trim(filstt)//'.vel_point')
      if(io_fr/=0) open(io_fr,status='unknown',access=a, file=trim(filstt)//'.vel_friction')

      a = 'append'
   end subroutine io_openfiles

!--------------------------------------------------------------------------
!  Close files written to during runtime
!--------------------------------------------------------------------------
   subroutine io_closefiles()
      if(mpi_rnk/=0) return
      if(io_dt/=0) close(io_dt)
      if(io_KE/=0) close(io_KE)
      if(io_ID/=0) close(io_ID)
      if(io_pt/=0) close(io_pt)
      if(io_fr/=0) close(io_fr)
   end subroutine io_closefiles


!--------------------------------------------------------------------------
!  Write to files
!--------------------------------------------------------------------------
   subroutine io_write2files()
     implicit none
     

      if (mod(tim_step,i_save_rate1)==0) then
          call io_save_phys_plane()
          !call io_save_phys_field()
      endif
      if(modulo(tim_step,i_save_rate1)==0) then
         call io_save_state()
         fnameima=trim(filstt)//'.'//extc//'.'//'sth'
         call saveStats(fnameima)
         ! call io_save_spectrum()
         ! call io_save_meanprof()
         extn = extn+1
         indice = 1
      endif

      if(modulo(tim_step,i_save_rate2)==0) then
         !if(io_KE/=0) call io_write_energy()
         !if(io_ID/=0) call io_write_totEID()
         if(io_pt/=0) call io_write_pointvel()
         !if(io_fr/=0) call io_write_friction()
         if(io_dt/=0 .and. d_timestep>0d0) call io_write_timestep()
         io_save2 = io_save2+1
      end if

      if(io_dt/=0 .and. tim_new_dt) call io_write_timestep()

      if(modulo(tim_step,i_save_rate2*50)==0) then
         call io_closefiles()
         !call io_openfiles()
      end if

   end subroutine io_write2files


!--------------------------------------------------------------------------
!  Load state - start from previous solution
!--------------------------------------------------------------------------
   subroutine io_load_state()
    implicit none
      integer :: K__, M__, Mp__
      integer :: e, i, rd
      integer :: n, N_
      integer :: commu, strow
      double precision :: d
      integer, parameter :: lda = 2*i_KL+1
      double precision :: A(lda,lda,i_N)
      double precision, allocatable :: r(:)
      logical :: interp
      double precision :: writetimer
      integer :: iaux(1)
      writetimer = 0d0
      ! ------------------------------------------------------------------------
      
      commu = MPI_COMM_WORLD ! HF5 needs an integer here
      
      if(mpi_rnk==0) writetimer=MPI_Wtime()-writetimer
      filename=trim(dirinp)//'/'//trim(filinp)//'.'//extc//'.h5'

      call H5Fopen_f(trim(filename),H5F_ACC_RDONLY_F,fid,h5err) ! Salida: fid

      ! write(*,*) 'El nimbre es ',filename

      ! --------------------------------------------------------------------
      ! ------------ Checking ----------------------------------------------
      ! --------------------------------------------------------------------

      ! Does the file exists? 

      if(h5err/=0) then
         if(mpi_rnk==0) print*, 'state file not found!: '//filename
         call mpi_barrier(MPI_COMM_WORLD, mpi_er)
         call mpi_finalize(mpi_er)
         stop 'io_load_state: file not found!'
      end if

      ! In this version we are not going to interpolate the field. 
      ! This will be done during preprocessing of the field. 

      ! . 
      if (mpi_rnk==0) then
         hdims = (/1/)
        call H5LTread_dataset_int_f(fid,"N",N_,hdims,h5err)
        call H5LTread_dataset_int_f(fid,"K",K__,hdims,h5err)
        call H5LTread_dataset_int_f(fid,"M",M__,hdims,h5err)
        call H5LTread_dataset_int_f(fid,"Mp",MP__,hdims,h5err)      
        if(N_/=i_N.or.K__ /=i_K.or.M__ /=i_M.or.Mp__/=i_Mp) then
            write(*,*) ' N :', N_,  ' --> ',i_N
            write(*,*) ' K :', K__, ' --> ',i_K 
            write(*,*) ' M :', M__, ' --> ',i_M 
            write(*,*) ' Mp:', Mp__,' --> ',i_Mp 
            stop 9
        endif
    endif

      ! --------------------------------------------------------------------
      ! ------------ Checking END ------------------------------------------
      ! --------------------------------------------------------------------

      hdims=(/1/)
      if(mpi_rnk==0) then 

         call H5LTread_dataset_double_f(fid,"/time",tim_t,hdims,h5err)

         call H5LTread_dataset_double_f(fid,"/dt",tim_dt,hdims,h5err)
 
         if(d_timestep>0d0) tim_dt = d_timestep

         ! dtcor y dtcfl. 
         call H5LTread_dataset_double_f(fid,"/dtcor",tim_corr_dt,hdims,h5err)
         call H5LTread_dataset_double_f(fid,"/dtcfl",tim_cfl_dt,hdims,h5err)
         call H5Fclose_f(fid,h5err)
      endif

      call MPI_BCAST(tim_t,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tim_corr_dt,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tim_cfl_dt,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      
      ! In the previous version thre was a call to io_load_coll, with interpolation
      ! However, as we are not interpolating anymore, this call does not make sense
      ! The fields are read here. 
      
      ! Reading. Taking care of vel_ur definition 
      
      call MPI_INFO_CREATE(info,h5err)
      
      call h5pcreate_f(H5P_FILE_ACCESS_F,pid,h5err)
      call h5pset_fapl_mpio_f(pid,MPI_COMM_WORLD,info,h5err)
      call h5fopen_f(trim(filename),H5F_ACC_RDONLY_F,fid,h5err,pid)
      call h5pclose_f(pid,h5err)
      
      !Load the data to the allocated array and close the file
    
      hdims2 = (/i_N,var_H%pH1+1/) ! Este es el tamaño que tengo que recibir. 
      ! El problema es que cada uno lo pone en un sitio distinto que NO TIENE porque
      ! ser correlativo! Lo pongo en sitios definidos por PH0 
      strow = var_H%pH0 + 1

      ! Ahora mismo strow no se usa. Se manda simplemente var_H%ph1. No se como reaccionará esto a cambio de vectores


!      write(*,*) 'soy',mpi_rnki' leo ',var_H%pH1+1,' desde ',strow

!      stop
      
      call h5load_parallel(fid,"/Ur/Re",2,hdims2,strow,mpi_rnk,mpi_sze,commu,info,c1%Re,h5err)
      vel_ur%Re = c1%Re
      call h5load_parallel(fid,"/Ur/Im",2,hdims2,strow,mpi_rnk,mpi_sze,commu,info,c1%Re,h5err)
      vel_ur%Im = c1%Re
      call h5load_parallel(fid,"/Uz/Re",2,hdims2,strow,mpi_rnk,mpi_sze,commu,info,c1%RE,h5err)
      vel_uz%Re = c1%Re
      call h5load_parallel(fid,"/Uz/Im",2,hdims2,strow,mpi_rnk,mpi_sze,commu,info,c1%Re,h5err)
      vel_uz%Im = c1%Re
      call h5load_parallel(fid,"/Ut/Re",2,hdims2,strow,mpi_rnk,mpi_sze,commu,info,c1%Re,h5err)
      vel_ut%Re = c1%Re
      call h5load_parallel(fid,"/Ut/Im",2,hdims2,strow,mpi_rnk,mpi_sze,commu,info,c1%Re,h5err)
      vel_ut%Im = c1%Re

      
      call h5fclose_f(fid,h5err)

      if (mpi_rnk.eq.0) then
          writetimer = MPI_Wtime() - writetimer
          write(*,*)
          write(*,*) 'Done reading fields'
          write(*,'(a20,f10.3,a3)') 'Time spend reading:',writetimer,'sc'
      endif

      call mpi_barrier(MPI_COMM_WORLD, mpi_er)
   end subroutine io_load_state

!--------------------------------------------------------------------------
!  Save state
!--------------------------------------------------------------------------
   subroutine io_save_state()
      integer :: e, f
      integer :: rd, Hd, ReImd, strow
      integer :: r,dt,dtcor,dtcfl, Ur,Ut,Uz
      integer :: info
      integer(hid_t) :: G1, G2, G3
      

      write(extc,'(I4.4)') extn

      info = MPI_INFO_NULL
      fnameima=trim(dirinp)//trim(filinp)//'.'//extc//'.'//'h5'

      if (mpi_rnk==0) then
          write(*,*) 'Writing the file ...',trim(fnameima)
      end if

      ! Save header, only master do this. 

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)


      call h5pcreate_f(H5P_FILE_ACCESS_F,pid,h5err)
      call h5pset_fapl_mpio_f(pid,MPI_COMM_WORLD,info,h5err)
      call h5pset_sieve_buf_size_f(pid, siever, h5err)
      call h5fcreate_f(trim(fnameima),H5F_ACC_TRUNC_F,fid,h5err,H5P_DEFAULT_F,pid)
      call h5gcreate_f(fid, '/Ur', G1, h5err)
      call h5gcreate_f(fid, '/Ut', G2, h5err)
      call h5gcreate_f(fid, '/Uz', G3, h5err)
      call h5pclose_f(pid,h5err)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      strow = var_H%pH0+1
      hdims2 = (/i_N,var_H%pH1+1/) ! This a bit tricky. The actual size is no the sum of var_H%pH0(mpi_size-1) + var_H%pH1

      call h5dump_parallel(G1,"/Ur/Re",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,vel_ur%Re,h5err)
      call h5dump_parallel(G1,"/Ur/Im",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,vel_ur%Im,h5err)

      call h5dump_parallel(G2,"/Ut/Re",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,vel_ut%Re,h5err)
      call h5dump_parallel(G2,"/Ut/Im",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,vel_ut%Im,h5err)

      call h5dump_parallel(G3,"/Uz/Re",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,vel_uz%Re,h5err)
      call h5dump_parallel(G3,"/Uz/Im",2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,vel_uz%Im,h5err)

      call h5gclose_f(G1,h5err)
      call h5gclose_f(G2,h5err)
      call h5gclose_f(G3,h5err)

      call h5fclose_f(fid,h5err)

      if(mpi_rnk==0) then

        call h5fopen_f(trim(fnameima),H5F_ACC_RDWR_F,fid,h5err)
        hdims = (/1/)
        call h5ltmake_dataset_double_f(fid,"time",1,hdims,(/tim_t/),h5err)
        call h5ltmake_dataset_double_f(fid,"Re",1,hdims,(/d_Re/),h5err)
        call h5ltmake_dataset_double_f(fid,"alpha",1,hdims,(/d_alpha/),h5err)

        call h5ltmake_dataset_int_f(fid,"N" ,1,hdims,(/i_N/),h5err)
        call h5ltmake_dataset_int_f(fid,"M" ,1,hdims,(/i_M/),h5err)
        call h5ltmake_dataset_int_f(fid,"K" ,1,hdims,(/i_K/),h5err)
        call h5ltmake_dataset_int_f(fid,"Mp",1,hdims,(/i_Mp/),h5err)
        
        call h5ltmake_dataset_double_f(fid,"dt"   ,1,hdims,(/tim_dt/),h5err)
        call h5ltmake_dataset_double_f(fid,"dtcfl",1,hdims,(/tim_cfl_dt/),h5err)
        call h5ltmake_dataset_double_f(fid,"dtcor",1,hdims,(/tim_corr_dt/),h5err)

        hdims = (/i_N/)
        call h5ltmake_dataset_double_f(fid,"r"   ,1,hdims,mes_D%r(1:i_N,1),h5err)
        call h5fclose_f(fid,h5err)
      end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 


   end subroutine io_save_state

!--------------------------------------------------------------------------
!  Guardar plano de velocidades radial y axial
!--------------------------------------------------------------------------

   subroutine io_save_phys_plane()
      integer :: info
      integer :: n, n_, strow, core_rad
      integer(hid_t) :: G1, G2
      character(len=20) ::cadena, nombre_dataset1, nombre_dataset2, nombre_dataset3
      

      core_rad=(mpi_sze/_Ns)  !numero de cores radiales +1
     
     write(extc,'(I4.4)') extn
     write(index,'(I3.3)') indice

      info = MPI_INFO_NULL
      fnamephys=trim(dirinp)//trim(filinp)//'.'//extc//'.'//index
      fnamephys=trim(fnamephys)//'.'//'spp'

     ! if (mpi_rnk<core_rad) then
      ! Save header, only master do this.
     
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      write(*,*) 'Soy', mpi_rnk
      call h5pcreate_f(H5P_FILE_ACCESS_F,pid,h5err)
      call h5pset_fapl_mpio_f(pid,MPI_COMM_WORLD,info,h5err)
      call h5pset_sieve_buf_size_f(pid, siever, h5err)
      call h5fcreate_f(trim(fnamephys),H5F_ACC_TRUNC_F,fid,h5err,H5P_DEFAULT_F,pid)
      call h5gcreate_f(fid, '/radial', G1, h5err)
      !call h5gcreate_f(fid, '/axial', G2, h5err)
      call h5pclose_f(pid,h5err)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     
    
      
    
   
      hdims2=(/i_Th,i_pN/) !Dimensiones plano radial
      strow=1  !Creo que no sirve para nada, de momento
     
      
         do n = 1, mes_D%pN
            n_ = mes_D%pNi + n - 1	
               p1%Re(20,:,n)=vel_z%Re(20,:,n)+vel_U(n_)
               p3%Re(:,1,n)=vel_z%Re(:,1,n)+vel_U(n_)
         enddo

      write(*,*) 'Soy', mpi_rnk
     
      
         write(cadena, '(I1)') 1
         nombre_dataset1="/radial/vel_r_"//cadena
         nombre_dataset2="/radial/vel_t_"//cadena
         nombre_dataset3="/radial/vel_z_"//cadena

         
         if (mpi_rnk<core_rad) then  
         call h5dump_parallel2(G1,nombre_dataset1,2, hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,vel_r%Re(20,:,:),h5err)
         !call h5dump_parallel(G1,nombre_dataset2,2, hdims2,strow,mpi_rnk,core_rad,MPI_COMM_WORLD,info,vel_t%Re(20,:,:),h5err)
         !call h5dump_parallel(G1,nombre_dataset3,2, hdims2,strow,mpi_rnk,core_rad,MPI_COMM_WORLD,info,   p1%Re(20,:,:),h5err)
         end if
         write(*,*) 10
         call h5gclose_f(G1,h5err)
      
     ! hdims2=(/i_pZ,i_pN/) !Dimensiones plano axial

         write(*,*) 11
      
        ! write(cadena, '(I1)') 1
        ! nombre_dataset1="/axial/vel_r_"//cadena
        ! nombre_dataset2="/axial/vel_t_"//cadena
        ! nombre_dataset3="/axial/vel_z_"//cadena
        ! call h5dump_parallel(G2,nombre_dataset1,2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,vel_r%Re(:,1,:),h5err)
        ! call h5dump_parallel(G2,nombre_dataset2,2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,vel_t%Re(:,1,:),h5err)
        ! call h5dump_parallel(G2,nombre_dataset3,2,hdims2,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,   p3%Re(:,1,:),h5err)
         

         


      !call h5gclose_f(G2,h5err)

      call h5fclose_f(fid,h5err)
      
      

      if(mpi_rnk==0) then

        call h5fopen_f(fnamephys,H5F_ACC_RDWR_F,fid,h5err)
        hdims = (/1/)
        call h5ltmake_dataset_double_f(fid,"time",1,hdims,(/tim_t/),h5err)
        call h5ltmake_dataset_double_f(fid,"Re",1,hdims,(/d_Re/),h5err)
        call h5ltmake_dataset_double_f(fid,"alpha",1,hdims,(/d_alpha/),h5err)

        call h5ltmake_dataset_int_f(fid,"N" ,1,hdims,(/i_N/),h5err)
        call h5ltmake_dataset_int_f(fid,"M" ,1,hdims,(/i_M/),h5err)
        call h5ltmake_dataset_int_f(fid,"K" ,1,hdims,(/i_K/),h5err)
        call h5ltmake_dataset_int_f(fid,"Mp",1,hdims,(/i_Mp/),h5err)
    
        
        call h5ltmake_dataset_double_f(fid,"dt"   ,1,hdims,(/tim_dt/),h5err)

        hdims = (/i_N/)
        call h5ltmake_dataset_double_f(fid,"r"   ,1,hdims,mes_D%r(1:i_N,1),h5err)

        call h5fclose_f(fid,h5err)
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

      indice=indice+1
   
   end subroutine io_save_phys_plane

!--------------------------------------------------------------------------
!  Guardar campo de velocidades radial y axial
!--------------------------------------------------------------------------

   subroutine io_save_phys_field()
      integer :: info
      integer :: strow, n, n_
      integer(hid_t) :: G1, G2
      


      write(extc,'(I4.4)') extn

      info = MPI_INFO_NULL
      fnamephys=trim(dirinp)//trim(filinp)//'.'//extc//'.'//'sppf'


      !if (mpi_rnk==0) then
      !write(*,*) 'Writing the file ...',trim(fnamephys)
      !end if

      ! Save header, only master do this.

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call h5pcreate_f(H5P_FILE_ACCESS_F,pid,h5err)
      call h5pset_fapl_mpio_f(pid,MPI_COMM_WORLD,info,h5err)
      call h5pset_sieve_buf_size_f(pid, siever, h5err)
      call h5fcreate_f(trim(fnamephys),H5F_ACC_TRUNC_F,fid,h5err,H5P_DEFAULT_F,pid)
      call h5gcreate_f(fid, '/sta', G1, h5err)
      call h5pclose_f(pid,h5err)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      hdims3=(/i_pZ,i_Th,i_pN/) !Dimensiones campo completo
      strow=1  !Creo que no sirve para nada, de momento


      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1	
            p1%Re(:,:,n)=vel_z%Re(:,:,n)+vel_U(n_)
      enddo

         

         call h5dump_parallel(G1,"/sta/vel_r",3, hdims3,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,vel_r%Re(:,:,:),h5err)
         call h5dump_parallel(G1,"/sta/vel_t",3, hdims3,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,vel_t%Re(:,:,:),h5err)
         call h5dump_parallel(G1,"/sta/vel_z",3, hdims3,strow,mpi_rnk,mpi_sze,MPI_COMM_WORLD,info,p1%Re(:,:,:),h5err)
         

      call h5gclose_f(G1,h5err)
      

      call h5fclose_f(fid,h5err)
       
      if(mpi_rnk==0) then

        call h5fopen_f(fnamephys,H5F_ACC_RDWR_F,fid,h5err)
        hdims = (/1/)
        call h5ltmake_dataset_double_f(fid,"time",1,hdims,(/tim_t/),h5err)
        call h5ltmake_dataset_double_f(fid,"Re",1,hdims,(/d_Re/),h5err)
        call h5ltmake_dataset_double_f(fid,"alpha",1,hdims,(/d_alpha/),h5err)

        call h5ltmake_dataset_int_f(fid,"N" ,1,hdims,(/i_N/),h5err)
        call h5ltmake_dataset_int_f(fid,"M" ,1,hdims,(/i_M/),h5err)
        call h5ltmake_dataset_int_f(fid,"K" ,1,hdims,(/i_K/),h5err)
        call h5ltmake_dataset_int_f(fid,"Mp",1,hdims,(/i_Mp/),h5err)
        
        call h5ltmake_dataset_double_f(fid,"dt"   ,1,hdims,(/tim_dt/),h5err)

        hdims = (/i_N/)
        call h5ltmake_dataset_double_f(fid,"r"   ,1,hdims,mes_D%r(1:i_N,1),h5err)
        
        call h5fclose_f(fid,h5err)
        
      end if
      
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
   
   end subroutine io_save_phys_field


!--------------------------------------------------------------------------
!  save axial mean flow profile
!--------------------------------------------------------------------------
   subroutine io_save_meanprof()
      integer :: n

      if(mpi_rnk/=0) return
      write(extc,'(I4.4)') extn
      fnameima=trim(filstt)//'.'//extc//'_prof'
      
      open(11, status='unknown', file=fnameima)
      write(11,*) '% t = ', tim_t
      write(11,*) '% r  uz(r) uz(r) code'
      
      ! They are removing from the flow a kind of parabolic profile. I does not affect any other quantity

      do n = 1, i_N
         write(11,'(3e20.12)')  mes_D%r(n,1),  &
            vel_uz%Re(n,0) + 1d0-mes_D%r(n,2) ,vel_uz%Re(n,0) ! They a removing from the flow a kind of para
      end do
      close(11)

   end subroutine io_save_meanprof

!--------------------------------------------------------------------------
!  write velocity at points.  r,t,z components at 3 points.
!--------------------------------------------------------------------------
   subroutine io_write_pointvel()
      integer :: n,rt(3),r(3)
      double precision :: d(9), d_
      if(_Ns/=1) stop 'io_write_pointvel: put _Ns=1'

      d = 1d0      
      do n = 1, i_N
         d_ = dabs(mes_D%r(n,1)-0.6d0) !0.6667d0)
         if(d_<d(1)) r(1) = n 
         if(d_<d(1)) d(1) = d_
         d_ = dabs(mes_D%r(n,1)-0.6d0) !0.7550d0)
         if(d_<d(2)) r(2) = n 
         if(d_<d(2)) d(2) = d_
         d_ = dabs(mes_D%r(n,1)-0.6d0) !0.9217d0)
         if(d_<d(3)) r(3) = n 
         if(d_<d(3)) d(3) = d_
      end do
      do n = 0, mpi_sze-1
         if(mes_D%pNi_(n)<=r(1)  &
            .and. mes_D%pNi_(n)+mes_D%pN_(n)-1>=r(1)) rt(1) = n
         if(mes_D%pNi_(n)<=r(2)  &
            .and. mes_D%pNi_(n)+mes_D%pN_(n)-1>=r(2)) rt(2) = n
         if(mes_D%pNi_(n)<=r(3)  &
            .and. mes_D%pNi_(n)+mes_D%pN_(n)-1>=r(3)) rt(3) = n
      end do

      if(rt(1)==mpi_rnk) then
         d(1) = vel_r%Re(0,0,r(1)-mes_D%pNi+1)
         d(2) = vel_t%Re(0,0,r(1)-mes_D%pNi+1)
         d(3) = vel_z%Re(0,0,r(1)-mes_D%pNi+1)
      end if      
      if(rt(2)==mpi_rnk) then
         d(4) = vel_r%Re(i_pZ/3,i_Th/3,r(2)-mes_D%pNi+1)
         d(5) = vel_t%Re(i_pZ/3,i_Th/3,r(2)-mes_D%pNi+1)
         d(6) = vel_z%Re(i_pZ/3,i_Th/3,r(2)-mes_D%pNi+1)
      end if
      if(rt(3)==mpi_rnk) then
         d(7) = vel_r%Re(2*i_pZ/3,2*i_Th/3,r(3)-mes_D%pNi+1)
         d(8) = vel_t%Re(2*i_pZ/3,2*i_Th/3,r(3)-mes_D%pNi+1)
         d(9) = vel_z%Re(2*i_pZ/3,2*i_Th/3,r(3)-mes_D%pNi+1)
      end if

#ifdef _MPI
      call mpi_bcast(d(1), 3, mpi_double_precision,  &
         rt(1), MPI_COMM_WORLD, mpi_er)
      call mpi_bcast(d(4), 3, mpi_double_precision,  &
         rt(2), MPI_COMM_WORLD, mpi_er)
      call mpi_bcast(d(7), 3, mpi_double_precision,  &
         rt(3), MPI_COMM_WORLD, mpi_er)
#endif
      if(mpi_rnk/=0) return
      if(tim_step==0)  write(io_pt,*) '# r=,', real(mes_D%r(r(1),1)),  &
         real(mes_D%r(r(2),1)), real(mes_D%r(r(3),1))
      write(io_pt,'(10e16.8)') tim_t, (d(n),n=1,9)

   end subroutine io_write_pointvel


!--------------------------------------------------------------------------
!  write:  1, time;  2,  bulk vel / excess pressure fraction if fixed flux;
!          3, mean uz at r=0;  4, friction velocity. 
!--------------------------------------------------------------------------
   subroutine io_write_friction(st,t0,timr,retau)
      double precision :: Ub, Uc, Ufr, retau
      integer, intent(in) ::st
      double precision :: t0,timr

      Ub = 0.5d0 + 2d0*dot_product(vel_uz%Re(:,0),mes_D%intrdr)
      if(b_const_flux)  Ub = vel_Pr0
      Uc = 1d0 + dot_product(vel_uz%Re(1:1+i_KL,0),mes_D%dr0(:,0))
      Ufr = dot_product(vel_uz%Re(i_N-i_KL:i_N,0),mes_D%dr1(:,1))
      Ufr = dsqrt(dabs((Ufr-2d0)/d_Re))
      retau = D_Re*Ufr
      
      if(mpi_rnk/=0) return
324   format(3(i7),7(d14.6)) ! Antes estaba i5, pero superaba el lim de its en ocasiones
      write(io_cf,324) st, _Nr, _Ns, t0, timr, Ub, Uc, Ufr, retau, tim_dt
      call flush(io_cf)
      
   end subroutine io_write_friction


!--------------------------------------------------------------------------
!  write to timestep file
!--------------------------------------------------------------------------
   subroutine io_write_timestep()
   10 format(3e18.10,I2)
   11 format(1e18.10,I8,3e17.9,I2)
      if(tim_new_dt) then
         if(mpi_rnk/=0) return
         write(io_dt,11) tim_t, tim_step,  &
            tim_dt, tim_corr_dt, tim_cfl_dt, tim_cfl_dir
      else
         call vel_maxtstep()
         if(mpi_rnk/=0) return
         write(io_dt,10) tim_t, tim_corr_dt, tim_cfl_dt, tim_cfl_dir
      end if
   end subroutine io_write_timestep




   

!**************************************************************************
 end module io
!**************************************************************************

   
subroutine h5dump_parallel(fid,name,ndims,dims,strow,rank,size,comm,info,data,ierr)
  use h5lt
  use hdf5
  use mpif
  implicit none

  
  integer(hid_t), intent(in):: fid
  character(len=*), intent(in):: name
  integer, intent(in):: ndims
  integer(hsize_t), dimension(ndims), intent(in):: dims
  integer, intent(in):: rank,size
  integer, intent(in):: comm,info, strow
  real(kind = 8),intent(in):: data
  integer:: ierr

  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace
  integer(hid_t):: plist_id
  integer(hsize_t), dimension(ndims):: start,nooffset,totaldims
  integer, dimension(size):: lastdims
  integer:: mpierr
  integer:: lastdim



  start = 0
  nooffset = 0
  totaldims = dims

  lastdim = dims(ndims) ! Don't mess with ints and longs

  call MPI_ALLGATHER(lastdim,1,MPI_INTEGER,lastdims,1,MPI_INTEGER,comm,mpierr)

  totaldims(ndims) = sum(lastdims)

  !if (mpi_rnk==0) write(*,*) 'lastdims',lastdims
  ! write(*,*) mpi_rnk
  !stop


  !Create the global dataspace
  call h5screate_simple_f(ndims,totaldims,dspace,ierr)
  !Create the global dataset
  call h5dcreate_f(fid,name,H5T_IEEE_F64BE,dspace,dset,ierr)

  !Create the local dataset
  call h5screate_simple_f(ndims,dims,mspace,ierr)
  call h5sselect_hyperslab_f(mspace,H5S_SELECT_SET_F,nooffset,dims,ierr)

  !Select the hyperslab in the global dataset ! This would be  var_H%pH0_. We can send just this!
  start(ndims) = sum(lastdims(1:rank+1))-lastdims(rank+1)
  !start(ndims) = strow!!!! CHANGED!!!!!
  call h5sselect_hyperslab_f(dspace,H5S_SELECT_SET_F,start,dims,ierr)
  !Create data transfer mode property list                                                                                                                          
  call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)

  !Commit the memspace to the disk
  call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,data,dims,ierr,mspace,dspace,plist_id)

  !Close property list                                                                                                                                              
  call h5pclose_f(plist_id,ierr)

  !Close datasets and dataspaces
  call h5sclose_f(mspace,ierr)
  call h5dclose_f(dset,ierr)
  call h5sclose_f(dspace,ierr)

end subroutine h5dump_parallel

subroutine h5dump_parallel2(fid,name,ndims,dims,strow,rank,size,comm,info,data,ierr)
   use h5lt
   use hdf5
   use mpif
   implicit none
 
   
   integer(hid_t), intent(in):: fid
   character(len=*), intent(in):: name
   integer, intent(in):: ndims
   integer(hsize_t), dimension(ndims), intent(in):: dims
   integer, intent(in):: rank,size
   integer, intent(in):: comm,info, strow
   real(kind = 8),intent(in):: data
   integer:: ierr
 
   integer(hid_t):: dset
   integer(hid_t):: dspace,mspace
   integer(hid_t):: plist_id
   integer(hsize_t), dimension(ndims):: start,nooffset,totaldims
   integer, dimension(size):: lastdims
   integer:: mpierr
   integer:: lastdim
 
 
 
   start = 0
   nooffset = 0
   totaldims = dims
 
   lastdim = dims(ndims) ! Don't mess with ints and longs
 
   write(*,*) 'previo 1'

   !if (mpi_rnk<mpi_sze) then
      write(*,*) 'dentro del if'
   call MPI_ALLGATHER(lastdim,1,MPI_INTEGER,lastdims,1,MPI_INTEGER,comm,mpierr)
   !end if
   write(*,*) 'fuera del if'
   totaldims(ndims) = sum(lastdims)
 
  write(*,*) 'lastdims',lastdims
  write(*,*) mpi_rnk
   
 
   !Create the global dataspace
   call h5screate_simple_f(ndims,totaldims,dspace,ierr)

   write(*,*) '1'
   !Create the global dataset
   call h5dcreate_f(fid,name,H5T_IEEE_F64BE,dspace,dset,ierr)
   write(*,*) '2'
   !Create the local dataset
   call h5screate_simple_f(ndims,dims,mspace,ierr)
   call h5sselect_hyperslab_f(mspace,H5S_SELECT_SET_F,nooffset,dims,ierr)
   write(*,*) '3'
   !Select the hyperslab in the global dataset ! This would be  var_H%pH0_. We can send just this!
   start(ndims) = sum(lastdims(1:rank+1))-lastdims(rank+1)
   write(*,*) '4'
   !start(ndims) = strow!!!! CHANGED!!!!!
   call h5sselect_hyperslab_f(dspace,H5S_SELECT_SET_F,start,dims,ierr)
   write(*,*) '5'
   !Create data transfer mode property list                                                                                                                          
   call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
   call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
   write(*,*) '6'
   !Commit the memspace to the disk
   call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,data,dims,ierr,mspace,dspace,plist_id)
   write(*,*) '7'
 
   !Close property list                                                                                                                                              
   call h5pclose_f(plist_id,ierr)
   write(*,*) '8'
   !Close datasets and datas-paces
   call h5sclose_f(mspace,ierr)
   call h5dclose_f(dset,ierr)
   call h5sclose_f(dspace,ierr)
   write(*,*) '9'
 
 end subroutine h5dump_parallel2



! -------------------------------------------------------------------! 
! -------------------------------------------------------------------! 
! ---------------------  READING SUBROUTINES ------------------------! 
! -------------------------------------------------------------------! 
! -------------------------------------------------------------------! 
! -------------------------------------------------------------------! 

subroutine h5load_parallel(fid,name,ndims,dims,strow,rank,size,comm,info,data,ierr)
#include "../parallel.h"
  use h5lt
  use hdf5
  use mpif
  implicit none
  
  integer(hid_t), intent(in):: fid
  character(len=*), intent(in):: name
  integer, intent(in):: ndims
  integer(hsize_t), dimension(ndims), intent(in):: dims
  integer, intent(in):: rank,size, strow
  integer, intent(in):: comm,info
  real(kind = 8),intent(out):: data
  integer:: ierr

  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace
  integer(hid_t):: plist_id
  integer(hsize_t), dimension(ndims):: start,nooffset,totaldims
  integer, dimension(size):: lastdims
  integer:: mpierr

  integer:: lastdim

  start = 0
  nooffset = 0
  totaldims = dims

  lastdim = dims(ndims) ! Don't mess with ints and longs

  call MPI_ALLGATHER(lastdim,1,MPI_INTEGER,lastdims,1,MPI_INTEGER,comm,mpierr)

  totaldims(ndims) = sum(lastdims)

  !Open the global dataset and get the global dataspace
  call h5dopen_f(fid,name,dset,ierr)
  call h5dget_space_f(dset,dspace,ierr)
  call h5screate_simple_f(ndims,dims,mspace,ierr)
  call h5sselect_hyperslab_f(mspace,H5S_SELECT_SET_F,nooffset,dims,ierr)

  !Select the hyperslab in the global dataset
  start(ndims) = sum(lastdims(1:rank+1))-lastdims(rank+1)
  !start(ndims) = strow ! Changed!!!

  call h5sselect_hyperslab_f(dspace,H5S_SELECT_SET_F,start,dims,ierr)

  !Create data transfer mode property list                                                                                                                          
  call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)

  !Commit the memspace to the disk
  call h5dread_f(dset,H5T_NATIVE_DOUBLE,data,dims,ierr,mspace,dspace,plist_id)

  !Close property list                                                                                                                                              
  call h5pclose_f(plist_id,ierr)

  !Close datasets and dataspaces
  call h5sclose_f(mspace,ierr)
  call h5dclose_f(dset,ierr)
  call h5sclose_f(dspace,ierr)

end subroutine h5load_parallel


subroutine avanza(text)
   implicit none
   character(len=256) text
   do
      read(19,'(a)') text
      if(text(1:2)/='CC') exit
   enddo
end subroutine



!  
!!--------------------------------------------------------------------------
!!  save spectrum: Commented. we have to dedice how to do it. 
!!--------------------------------------------------------------------------
!   subroutine io_save_spectrum()
!      double precision :: n_(1:i_N), k_(0:i_K1), m_(0:i_M1)
!      double precision :: n__(1:i_N), k__(0:i_K1), m__(0:i_M1)
!      double precision,save :: TM(i_N,i_N), x(i_N)
!      double precision :: d(i_N), dRe(i_N), dIm(i_N)
!      logical, save :: set=.false.
!      character(4) :: cnum
!      integer :: i, n,kp
!      _loop_km_vars
!   10 format(i4,1e20.12)
!      
!      if(.not.set) then
!         set =.true.
!         do n = 0, i_N-1
!            x(n+1) = 0.5d0 * ( 1d0 + dcos(d_PI*(i_N-n)/dble(i_N)) )
!         end do
!         do n = 1, i_N
!            call cheby(n-1, 0, x, i_N, TM(1,n))
!         end do
!         call mes_mat_invert(i_N,TM,i_N)
!         TM = transpose(TM)
!      end if
!
!      n_ = 0d0
!      k_ = 0d0
!      m_ = 0d0
!      _loop_km_begin
!         dRe = matmul(vel_ur%Re(:,nh), TM)
!         dIm = matmul(vel_ur%Im(:,nh), TM)
!         d = dRe*dRe+dIm*dIm
!         dRe = matmul(vel_ut%Re(:,nh), TM)
!         dIm = matmul(vel_ut%Im(:,nh), TM)
!         d = max(d, dRe*dRe+dIm*dIm)
!         dRe = matmul(vel_uz%Re(:,nh), TM)
!         dIm = matmul(vel_uz%Im(:,nh), TM)
!         d = max(d, dRe*dRe+dIm*dIm)
!         d = dsqrt(d)
!         kp = abs(k)
!         do n = 1, i_N
!             n_(n)  = max(d(n), n_(n))
!             k_(kp) = max(d(n), k_(kp))
!             m_(m)  = max(d(n), m_(m))            
!         end do
!      _loop_km_end
!#ifdef _MPI
!      call mpi_allreduce(n_, n__, i_N, mpi_double_precision,  &
!         mpi_max, MPI_COMM_WORLD, mpi_er)
!      n_ = n__
!      call mpi_allreduce(k_, k__, i_K, mpi_double_precision,  &
!         mpi_max, MPI_COMM_WORLD, mpi_er)
!      k_ = k__
!      call mpi_allreduce(m_, m__, i_M, mpi_double_precision,  &
!         mpi_max, MPI_COMM_WORLD, mpi_er)
!      m_ = m__
!#endif
!
!      if(mpi_rnk/=0) return
!      write(cnum,'(I4.4)') io_save1
!      open(11, status='unknown', file='vel_spec'//cnum//'.dat')
!      write(11,*) '# t = ', tim_t
!      write(11,*) '# n'
!      do i = 1, i_N
!         write(11,10) i, n_(i)      
!      end do
!      write(11,*)
!      write(11,*) '# k'
!      do i = 0, i_K1
!         write(11,10) i, k_(i)      
!      end do
!      write(11,*)
!      write(11,*) '# m'
!      do i = 0, i_M1
!         write(11,10) i*i_Mp, m_(i)      
!      end do
!      close(11)
!
!   end subroutine io_save_spectrum

! I've MPI_COMM_WORLDented the wo routines below. Tehy use a HUGE quantity of memory and are not needed a this point. 
! The routine !!!      call var_coll_curl(vel_ur,vel_ut,c3, c1,c2,c3) 
! uses thre!!! coll variables just for computing the dissipation! 

!!!!--------------------------------------------------------------------------
!!!!  write to energy file
!!!!--------------------------------------------------------------------------
!!!   subroutine io_write_energy()
!!!      double precision :: E,Ek0,Em0,E_, Ek(0:i_K1), Em(0:i_M1)
!!!
!!!      call var_coll_norm(vel_ur, E,Ek,Em)
!!!      Ek0 = Ek(0)
!!!      Em0 = Em(0)
!!!      call var_coll_norm(vel_ut, E_,Ek,Em)
!!!      E   = E + E_
!!!      Ek0 = Ek0 + Ek(0)
!!!      Em0 = Em0 + Em(0)
!!!      call var_coll_norm(vel_uz, E_,Ek,Em)
!!!      E   = E + E_
!!!      Ek0 = Ek0 + Ek(0)
!!!      Em0 = Em0 + Em(0)
!!!      
!!!      if(mpi_rnk/=0) return
!!!      write(io_KE,'(4e20.12)')  tim_t, E, Ek0, Em0
!!!      
!!!      if(E-Ek0>d_minE3d .or. tim_t<20d0) return
!!!      print*, 'io_write_energy: Relaminarised!'
!!!      open(99,file='RUNNING')
!!!      close(99, status='delete')
!!!
!!!   end subroutine io_write_energy
!!!
!!!
!!!!--------------------------------------------------------------------------
!!!!  (Total) Energy/E_lam, Input/I_lam, Dissip/D_lam
!!!!--------------------------------------------------------------------------
!!!   subroutine io_write_totEID()
!!!      double precision :: E,E_,E__, Ek(0:i_K1), Em(0:i_M1)
!!!      double precision :: Inp, Diss, Enrg, Lz
!!!
!!!      call var_coll_copy(vel_uz, c3)
!!!      if(mpi_rnk==0) c3%Re(:,0) = c3%Re(:,0) + vel_U
!!!      call var_coll_norm(vel_ur, E,Ek,Em)
!!!      call var_coll_norm(vel_ut, E_,Ek,Em)
!!!      call var_coll_norm(c3,     E__,Ek,Em)
!!!      Enrg = E + E_ + E__
!!!      
!!!      call var_coll_curl(vel_ur,vel_ut,c3, c1,c2,c3)
!!!      call var_coll_norm(c1, E,Ek,Em)
!!!      call var_coll_norm(c2, E_,Ek,Em)
!!!      call var_coll_norm(c3, E__,Ek,Em)
!!!      Diss = E + E_ + E__
!!!      
!!!      Lz = 2d0*d_PI/d_alpha
!!!      Enrg = Enrg / (d_PI**2/(3d0*d_alpha))
!!!      Diss = Diss * 2d0/d_Re
!!!      Diss = Diss / abs(vel_Up(i_N)*d_PI*Lz/d_Re)
!!!                                                   !   I/I_lam = 1+beta
!!!      Inp = 1d0 + vel_Pr0
!!!      
!!!      if(mpi_rnk/=0) return
!!!      write(io_ID,'(4e20.12)')  tim_t, Enrg, Inp, Diss
!!!
!!!   end subroutine io_write_totEID
