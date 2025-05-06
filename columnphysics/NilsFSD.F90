!GW
module NilsFSD

      use icepack_kinds
      use icepack_parameters, only: p01, p5, c0, c1, c2, c3, c4, c10
      use icepack_parameters, only: bignum, puny, gravit, pi
      use icepack_warnings, only: warnstr, icepack_warnings_add,  icepack_warnings_aborted
      use deformation_table

      implicit none
      public :: correct_FSD


    contains

      subroutine correct_FSD(ncat, nfsd, trcrn, aice, tarea, divu, dt,d_afsd_nils, aicen, ni, nj) ! Pull in tarea from ice_grid
      integer(kind=int_kind), intent(in) :: ncat, nfsd
      real(kind=dbl_kind), dimension(ncat), intent(in) :: aicen
      integer(kind=int_kind) :: k, n
      integer(kind=int_kind), intent(in) :: ni, nj
      real(kind=dbl_kind), dimension(nfsd,ncat), intent(inout) :: trcrn
      real(kind=dbl_kind), intent(in) :: aice, divu, dt, tarea
      real(kind=dbl_kind), dimension(:), intent(inout) :: d_afsd_nils
      real(kind=dbl_kind) :: before, divu_per_day, Dx, dt_s, crit_div_b, crit_div_t, crit_model_div_b, crit_model_div_t, total_negative_change, k_negative_change, weight, aice_limit,dF_all_remaining
      real(kind=dbl_kind), dimension(nfsd) :: dP, dF,area_midpoint, dF_raw 
      real(kind=dbl_kind), dimension(nfsd,ncat) :: dF_all, int_afsd, d_afsdn_nils
      real(kind=dbl_kind), dimension(ncat) :: loss_count
      logical, dimension(ncat) :: flag_redistribute
      logical :: stop_flag, enable
      d_afsd_nils=c0
      int_afsd=trcrn
      aice_limit=0.8
      enable=.false.
      if (enable .eq. .true.) then
      if (aice>aice_limit) then
         dt_s=dt/86400 ! Convert from s to day for timestep
         Dx=sqrt(tarea)/(1000) ! Convert from m^2 to km
         divu_per_day=divu*(86400) ! convert from /s to /day 
         crit_div_b=0.011
         crit_model_div_b=crit_div_b/(((50/Dx)**-0.49) *((1/dt_s)**-0.001))
         crit_div_t=0.028
         crit_model_div_t=crit_div_t/(((50/Dx)**-0.49) *((1/dt_s)**-0.001))
         !if ((tarea<1.0001*1015856600.0) .and. (tarea>0.9999*1015856600.0)) then
         !   write(warnstr,*) 'GW sanity1', dt, divu*(86400),divu_per_day, tarea, crit_model_div_b, crit_model_div_t
         !   call icepack_warnings_add(warnstr)
         !end if
         if (nj==402 .and. ni==169) then
            write(warnstr,*) 'GW',ni,nj,divu_per_day, tarea, crit_model_div_b, crit_model_div_t
            call icepack_warnings_add(warnstr)
         end if
         if (nj.eq.419 .and. ni.eq.166) then
            write(warnstr,*) 'GW',divu_per_day, tarea, crit_model_div_b, crit_model_div_t
            call icepack_warnings_add(warnstr)
         end if 

         if (divu_per_day > crit_model_div_b) then
            if (divu_per_day > crit_model_div_t) then
               divu_per_day=crit_model_div_t
            end if
            divu_per_day = divu_per_day*((50/Dx)**-0.49) *((1/dt_s)**-0.001)
            ! Rounding to 3 decimal places
            divu_per_day = nint(divu_per_day * 1000.0) / 1000.0
            !if ((tarea<1.0001*1015856600.0) .and. (tarea>0.9999*1015856600.0)) then
            !   write(warnstr,*) 'GW sanity1', dt, divu*(86400),divu_per_day, tarea, crit_model_div_b, crit_model_div_t
             !  call icepack_warnings_add(warnstr)
            !end if
            dF_raw=c0
            call lookup_deformation(divu_per_day,nfsd,dF_raw) !dP_raw = dP*area_midpoint
           ! if ((tarea<1.0001*1015856600.0) .and. (tarea>0.9999*1015856600.0)) then
           !    write(warnstr,*) 'GW sanity2', dF_raw(1), dF_raw(2), dF_raw(10), dF_raw(11), dF_raw(12)
            !   call icepack_warnings_add(warnstr)
            !end if
            !write(warnstr,*) 'dF_raw', dF_raw(1), dF_raw(2), dF_raw(10), dF_raw(11), dF_raw(12)
            !call icepack_warnings_add(warnstr)
            dP=((Dx/50)**2)*(dt_s/1)*dF_raw
            dP=dP*((aice-aice_limit)/(1-aice_limit)) ! This includes linear aice scaling...
            ! So dP has units like per area per time
            dF=dP/(0.66*(Dx*1000)**2)
         
            
            dF_all = c0
            do n=1, ncat
               dF_all(:,n) = dF(:)
            end do
            !Now to assign (first attempt will be a basic assignment ignoring if area is left untouched)
            ! Total negative change
            total_negative_change=c0
            total_negative_change = -dF(nfsd)! Assuming nfsd=12 it is the last 1 elements that are negative, this could just be a numeber fed from deformation_table
            loss_count(:)=c0
            
            flag_redistribute=.true.

            ! Debug check
            do n=ncat,1
            do k=nfsd-1,1
               if (dF_all(k,n)<0) then
                  write(warnstr,*) 'dF_all(k) is negative, this should not happen', k, n, dF_all(k,n)
                  call icepack_warnings_add(warnstr)
               end if
            end do
            end do
            ! Loop over each category (here only one category is used; change as needed)

! First before trcrn is even effected
            do n = 1, ncat
               do k = 1, nfsd
                  if (trcrn(k, n) < 0) then
                     write(warnstr,*) 'first trcrn has negative value at (', k, ',', n, '): ', trcrn(k, n)
                     call icepack_warnings_add(warnstr)
                  end if
               end do
            end do



            !if ((tarea<1.0001*1015856600.0) .and. (tarea>0.9999*1015856600.0)) then
             !  write(warnstr,*) 'GW sanity3', trcrn(1,3), trcrn(2,3), trcrn(10,3), trcrn(11,3), trcrn(12,3), dF_all(1,3), dF_all(2,3), dF_all(10,3), dF_all(11,3), dF_all(12,3)
             !  call icepack_warnings_add(warnstr)
            !end if
! Complex redistribution scheme
            do n = 1, ncat ! Loop over each categories
               if (sum(trcrn(:, n)) > 0.1) then 
               k_negative_change = c0
               dF_all_remaining=c0
               ! Loop over k from nfsd down to 1 (Python’s range(nfsd-1, 0, -1) is equivalent)
               stop_flag = .false.
               do k = nfsd, 1, -1
                  if ( dF_all(k, n) < -puny .and. k > 7 .and. (.not. stop_flag) .and. (abs(k_negative_change)+abs(loss_count(n)))<0.99*abs(total_negative_change)) then
                     if ( abs(dF_all(k, n)) < abs(trcrn(k, n))) then
                        loss_count(n) = loss_count(n) + abs(dF_all(k, n))
                        trcrn(k, n)   = trcrn(k, n)   + dF_all(k, n)
                        dF_all(k, n)  = c0
                        stop_flag     = .true.
                     elseif ( abs(dF_all(k, n)) > abs(trcrn(k, n))) then
                        dF_all_remaining = dF_all(k, n) + trcrn(k, n)
                        loss_count(n)    = loss_count(n) + abs(trcrn(k, n))
                        dF_all(k, n)     = dF_all(k, n) + trcrn(k, n)
                        trcrn(k, n)      = c0
                        if (k>8 .and. dF_all_remaining>puny .and. (dF_all(k-1, n) + dF_all_remaining) < 0) then
                           k_negative_change=k_negative_change+dF_all(k-1, n)
                           dF_all(k-1, n)=dF_all(k-1, n)+dF_all_remaining      
                        elseif (k>8 .and. dF_all_remaining>puny .and. (dF_all(k-1, n) + dF_all_remaining) > 0) then
                           dF_all(k-1, n)=dF_all(k-1, n)+dF_all_remaining
                        end if

                        !if ( dF_all(k-1, n) < 0 .and. (dF_all(k-1, n) - dF_all_remaining) < 0 ) then
                        !   dF_all(k-1, n) = dF_all(k-1, n) - dF_all_remaining
                        !elseif ( dF_all(k-1, n) > 0 .and. (dF_all(k-1, n) - dF_all_remaining) < 0 ) then
                        !   if ( k > 8 ) then
                        !      k_negative_change = k_negative_change + dF_all(k-1, n)
                        !      dF_all(k-1, n)    = dF_all(k-1, n) - dF_all_remaining
                        !   end if
                        !   ! If k is not > 8 then do nothing in this branch
                        !else
                        !   dF_all(k-1, n) = dF_all(k-1, n) - dF_all_remaining
                        !end if
                        if (k_negative_change+abs(loss_count(n))>abs(total_negative_change)) then
                           write(warnstr,*) 'GW warning k_negative_change+abs(loss_count(n))>abs(total_negative_change)', k, n, k_negative_change, loss_count(n), total_negative_change
                           call icepack_warnings_add(warnstr)
                        end if 
                     end if
                  else
                     if ( flag_redistribute(n) ) then
                        weight = (loss_count(n) + abs(k_negative_change)) / total_negative_change
                        flag_redistribute(n) = .false.
                     end if
                     if (dF_all(k, n) * weight<0) then
                        write(warnstr,*) 'dF_all(k, n) * weight<0  (', k, ',', n, '): ', trcrn(k, n), dF_all(k,n),weight, dF_all(k, n) * weight<0
                        call icepack_warnings_add(warnstr)
                     end if 
                     trcrn(k, n) = trcrn(k, n) + dF_all(k, n) * weight
                  end if
               end do
               
               end if 
            end do
!Simple redistribution scheme
 !do n=1, ncat ! Loop over each categories
 !   k_negative_change=c0
 !   if (aicen(n)>0.1) then
 !              do k=nfsd, 1, -1
 !                 if (dF_all(k,n) .lt. 0) then
 !                    if (abs(trcrn(k,n))>abs(dF_all(k,n))) then
 !                       trcrn(k,n)=trcrn(k,n)+dF_all(k,n)
 !                       loss_count(n)=loss_count(n)+abs(dF_all(k,n))
 !                       dF_all(k,n)=c0
 !                    else
 !                       dF_all(k,n)=dF_all(k,n)+abs(trcrn(k,n))
 !                       loss_count(n)=loss_count(n)+abs(trcrn(k,n))
 !                       trcrn(k,n)=c0
 !                    end if
 !                 else if (dF_all(k,n) .gt. 0) then
 !                    if (flag_redistribute(n)==.true.) then
 !                       weight=loss_count(n)/total_negative_change
 !                       flag_redistribute(n)=.false.
 !                    end if
 !                    before=trcrn(k,n)
 !                    trcrn(k,n)=trcrn(k,n)+dF_all(k,n)*weight
 !                    if (trcrn(k,n)<c0) then
 !                       write(warnstr,*) 'trcrn has negative value at  (', k, ',', n, '): ', trcrn(k, n), before, dF_all(k,n)*weight,weight
 !                       call icepack_warnings_add(warnstr)
 !                    end if
 !                 end if
 !                 
 !              end do
 !              endif
 !           end do

            !if ((tarea<1.0001*1015856600.0) .and. (tarea>0.9999*1015856600.0)) then
             !  write(warnstr,*) 'GW sanity4', trcrn(1,3), trcrn(2,3), trcrn(10,3), trcrn(11,3), trcrn(12,3), dF_all(1,3), dF_all(2,3), dF_all(10,3), dF_all(11,3), dF_all(12,3)
            !   call icepack_warnings_add(warnstr)
            !end if



! More checks
            do n = 1, ncat
               do k = 1, nfsd
                  if (trcrn(k, n) < 0) then
                     write(warnstr,*) 'before trcrn has negative value at (', k, ',', n, '): ', trcrn(k, n)
                     call icepack_warnings_add(warnstr)
                  end if
               end do
            end do
! Numerical error fix, for slight errors 
      do n = 1, ncat
               if (sum(trcrn(:, n))> c1 .and. sum(trcrn(:, n))<1.01) then
                  trcrn(:, n) = trcrn(:, n) / sum(trcrn(:, n))
               else if (sum(trcrn(:, n))>1.01) then
                  write(warnstr,*) 'more than 1% error', sum(trcrn(:, n))
                  call icepack_warnings_add(warnstr)
               end if
            end do
   ! Numerical error fix, for slight errors around 0
   do n = 1, ncat
      do k = 1, nfsd
         if (trcrn(k, n)<c0) then
            write(warnstr,*) 'trcrn has negative value at (', k, ',', n, '): ', trcrn(k, n)
            call icepack_warnings_add(warnstr)
         end if
      end do
   end do
!Sanity check
do n = 1, ncat
   if (abs(sum(trcrn(:, n)) - sum(int_afsd(:, n))) > 0.01 * sum(int_afsd(:, n))) then
      write(warnstr,*) 'trcrn sum for category', n, 'is not within a percent of int_afsd sum', sum(trcrn(:, n)), sum(int_afsd(:, n)), trcrn(10, n), int_afsd(10, n)
      call icepack_warnings_add(warnstr)
   end if
end do
! More checks
do n = 1, ncat
   do k = 1, nfsd
      if (trcrn(k, n) < 0) then
         write(warnstr,*) 'after trcrn has negative value at (', k, ',', n, '): ', trcrn(k, n)
         call icepack_warnings_add(warnstr)
      end if
   end do
end do
!if ((tarea<1.0001*1015856600.0) .and. (tarea>0.9999*1015856600.0)) then
 !  write(warnstr,*) 'GW sanity5', trcrn(1,3), trcrn(2,3), trcrn(10,3), trcrn(11,3), trcrn(12,3)
 !  call icepack_warnings_add(warnstr)
!end if
 ! Archive the final FSD
d_afsdn_nils=c0
d_afsdn_nils=trcrn-int_afsd
do k = 1, nfsd
   d_afsd_nils(k) = c0
   do n = 1, ncat
      d_afsd_nils(k) = d_afsd_nils(k) + aicen(n)*d_afsdn_nils(k,n)
   end do ! n                                                                                                                                                                                                                              
end do 

end if ! Maybe this is the endif
else
   d_afsd_nils=c0
end if
end if

end subroutine

end module


    
