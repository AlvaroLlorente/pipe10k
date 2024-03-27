!*************************************************************************
#include "../parallel.h"
 module type
!*************************************************************************
   use parameters
   implicit none
   save

   type harm
      integer              :: pH0,pH1, pH0_(0:_Np-1),pH1_(0:_Np-1)
   end type harm

   type spec
      double precision, allocatable, dimension(:,:) :: Re
      double precision, allocatable, dimension(:,:) :: Im
   end type spec

   type coll
      double precision, allocatable, dimension(:,:)     :: Re
      double precision, allocatable, dimension(:,:)     :: Im
   end type coll

   type phys
      double precision, allocatable, dimension(:,:,:)     :: Re
   end type phys

   type mesh
      double precision :: M(2*i_KL+1, 1-i_KL:i_N)
   end type mesh

   type lumesh
      integer,          allocatable, dimension(:)   :: ipiv
      double precision, allocatable, dimension(:,:) :: M
   end type lumesh


   type rdom
      integer          :: N, pNi,pN, pNi_(0:_Np-1),pN_(0:_Np-1)
      double precision :: r(i_N,-3:3)
      double precision :: intrdr(i_N)
      double precision :: dr0(1:1+i_KL,0:i_KL), dr1(1:1+i_KL,0:i_KL)
      type (mesh)      :: dr(i_KL)
      type (mesh)      :: radLap
   end type rdom


   type (rdom) :: mes_D

   ! type lumesh
   !    integer          :: ipiv(i_N)
   !    double precision :: M(3*i_KL+1, i_N)
   ! end type lumesh

   ! type phys
   !    double precision     :: Re(0:i_pZ-1, 0:i_Th-1, i_pN)
   ! end type phys

   ! type coll
   !    double precision     :: Re(i_N, 0:i_pH1)
   !    double precision     :: Im(i_N, 0:i_pH1)
   ! end type coll

   ! type spec
   !    double precision     :: Re(0:_Hs1, i_pN)
   !    double precision     :: Im(0:_Hs1, i_pN)
   ! end type spec

   end module