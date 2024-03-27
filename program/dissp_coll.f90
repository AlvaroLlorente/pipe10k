subroutine var_coll_dissp(comp)
   implicit none
      !double precision :: n_
      integer :: n, comp, n_
      _loop_km_vars

      ! Component and r derivative
      if (comp == 1) then
         call var_coll_copy(vel_ur,c1)
         call var_coll_meshmult(1,mes_D%dr(1),c1, c2)
      else if (comp == 2) then
         call var_coll_copy(vel_ut,c1)
         call var_coll_meshmult(1,mes_D%dr(1),c1, c2)
      else if (comp == 3) then
         call var_coll_copy(vel_uz,c1)
         call var_coll_meshmult(0,mes_D%dr(1),c1, c2)
      else
         print*, 'Dissp comp error'
      endif

      
      call tra_coll2phys(c2,vel_r) ! udr 2phys

      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
         dissr(n_) = 0
         dissr(n_) = dissr(n_) + sum(vel_r%Re(:,:,n)**2) ! diss sum
      end do


      ! Theta derivative
      _loop_km_begin
         c2%Re(:,nh) = -c1%Im(:,nh)*ad_m1r1(:,m) !ad_m1r1 lleva incluido el termino *mes_D%r(:,-1)
         c2%Im(:,nh) =  c1%Re(:,nh)*ad_m1r1(:,m) !*mes_D%r(:,-1)
      _loop_km_end

      call tra_coll2phys(c2,vel_r) ! udt 2phys

      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
         disst(n_) = 0
         disst(n_) = disst(n_) + sum(vel_r%Re(:,:,n)**2) ! diss sum
      end do


      ! Z derivative
      _loop_km_begin
         c2%Re(:,nh) = -c1%Im(:,nh)*ad_k1a1(k)
         c2%Im(:,nh) =  c1%Re(:,nh)*ad_k1a1(k)
      _loop_km_end

      call tra_coll2phys(c2,vel_r) ! udz 2phys

      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
         dissz(n_) = 0
         dissz(n_) = dissz(n_) + sum(vel_r%Re(:,:,n)**2) ! diss sum
      end do

      
      if (comp == 1) then
      dissrr = dissr + disst + dissz
      else if (comp == 2) then
      disstt = dissr + disst + dissz
      else if (comp == 3) then
      disszz = dissr + disst + dissz
      else
         print*, 'Dissp comp error'
      endif




end subroutine var_coll_dissp