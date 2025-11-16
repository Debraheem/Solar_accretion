! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff
      
      implicit none

      real(dp),save :: Lacc, Ladd ! global for saving

      contains

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_start_step => extras_start_step
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
         s% other_adjust_mdot => other_adjust_mdot
         s% other_energy => energy_routine


        ! this part only matters for the generation of the seed, since we relax.
        s% job% initial_h2 = s% x_ctrl(1)*1d-6
        s% job% initial_he3 = s% x_ctrl(2)*1d-6
        s% job% initial_he4 = s% x_ctrl(8) !0.2703d0
        s% job% initial_h1 = 1d0 - s% job% initial_h2 - s% job% initial_he3 -s% job% initial_he4 - s% initial_z

        write (*,*) 'adopted initial abundances'
        write (*,*) 'initial_h1', s% job%initial_h1
        write (*,*) 'initial_h2', s% job%initial_h2
        write (*,*) 'initial_he3', s% job%initial_he3
        write (*,*) 'initial_he4', s% job%initial_he4
        write (*,*) 'initial_z', s% initial_z

! these routines are only called once, so we should call again in start/finish step
! for modified accretion schemes.

       ! accretion composition controls.
        s% accretion_h2 = s% x_ctrl(1)*1d-6
        s% accretion_he3 = s% x_ctrl(2)*1d-6
        s% accretion_he4 = s% x_ctrl(8) !0.2703d0
        s% accretion_h1 = 1d0 - s% job% initial_h2 - s% job% initial_he3 -s% job% initial_he4 - s% initial_z

        s% job% new_Y = s% x_ctrl(8) ! we hold he3 and h2 fixed
        s% job% new_Z = s% initial_z

      end subroutine extras_controls


        integer function extras_start_step(id)
        integer, intent(in) :: id
        integer :: ierr
        integer :: k0
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        extras_start_step = 0
        end function extras_start_step






      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

        Lacc = 0d0
        Ladd = 0d0
         
      end subroutine extras_startup

      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         use chem_def, only:  ih1
         integer, intent(in) :: id
         integer :: ierr
         real(dp) :: tol_stop
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
        
        ! tolerance on stopping condition tol_stop
        tol_stop = 1d0 + 1d-8 ! can optionally make this a control
        if (s% x_logical_ctrl(1) .and. s% star_mass > s% x_ctrl(3)*tol_stop) then
        s% dt = s% dt/2
        extras_check_model = retry
        write(*, *) 'have overshot target mass, trying again with half the timestep'
        end if

      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 1
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         real(dp) :: omega2, nu_max, dr, delta_Pg, N2, integral, r
         integer :: k
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)

        names(1) = 'Delta_Pg'
        k = 0
        r = 0
        integral = 0
        N2 = 0
        dr = 0
        delta_Pg = 0

        integral = 0
        ! we integrate using cell edges instead of center.
        do k = 2, s% nz
        N2 = s% brunt_N2(k)
        r  = s% r(k)
        dr = s% rmid(k-1) - s% rmid(k)
        if (N2 > 0d0) integral = integral + sqrt(N2)*dr/r
        end do

        delta_Pg = sqrt(2d0)*pi*pi/integral
        !write (*,*) "delta_Pg = ", delta_Pg
        vals(1) = delta_Pg

         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         names(1) = 'zbar_div_abar'

         do k=1,s% nz
            vals(k,1) = s% zbar(k)/s% abar(k)
         end do
      end subroutine data_for_extra_profile_columns
  

        ! returns either keep_going or terminate.
        ! note: cannot request retry; extras_check_model can do that.
        integer function extras_finish_step(id)
        integer, intent(in) :: id
        integer :: ierr
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        extras_finish_step = keep_going

        if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
        end function extras_finish_step


        subroutine other_adjust_mdot(id, ierr)
        use star_def
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        integer :: step, i
        real (dp) :: sigmoid
        type(star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr) ! retrieve star id

        sigmoid  = 1 /(1+exp(-(s% star_age - 500d0)/100d0))
        !write(*,*) 'sigmoid', sigmoid
        if (s% star_mass < s% x_ctrl(3)) then ! target mass
        s% mstar_dot = sigmoid * s% x_ctrl(18)*Msun/secyer ! in solar masses per year
        s% do_element_diffusion = .false. ! keep diffusion off during accretion phase
        else
        s% mstar_dot = 0d0
        s% do_element_diffusion = .true. ! turn on diffusion after accretion phase
        end if

        !write(*,*) "current Accretion rate: " ,  s% mstar_dot/Msun*secyer

        end subroutine other_adjust_mdot

        subroutine energy_routine(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type(auto_diff_real_star_order1) :: xheat
        real(dp) :: mdot_msun, alpha_heat, xpos, Mstar
        type (star_info), pointer :: s
        integer :: k
        ierr = 0
        call star_ptr(id, s, ierr)

        ! make sure variables are set.
        mdot_msun = s% mstar_dot / Msun* secyer
        Mstar=(s% star_mass *Msun) ![g]
        alpha_heat = s% x_ctrl(6)
        k = 0d0
        do k = 1, s% nz
        s% extra_heat(k) = 0d0 ! make sure extra_heat is set.
        end do

        ! Deposit heat into q set by mass accreted over timestep.
        ! 1d-8 Msun/yr * 1d2 yr = 1d-6 Msun -> q = 1d-6 [Msun] / Mstar [Msun] .
        !s% x_ctrl(5) = s% mstar_dot * s% dt / (s% star_mass * Msun ) ! fraction of total mass.
        xpos = 0.01d0 ! units of q, where 0.01 = 1% by mass.

        if (s% mstar_dot <= 0.) return
        Lacc = s% x_ctrl(7) * (standard_cgrav * (s% mstar_dot) * (s% mstar)) / s% r(1)
        Ladd = (alpha_heat/2) * standard_cgrav * &
        (s% star_mass *Msun) * (s% mstar_dot) / (s% r(1)) ![erg/s]

        !write(*,*) 'Ladd/Lacc', Ladd/Lacc
        !write(*,*) 'Ladd', Ladd
        !write(*,*) 'Lacc', Lacc
        !write(*,*) 'Mdot', s% mstar_dot

        ! relics from Thomas Steindal' setup, commented out
        !Lacc =  safe_log10((1-alpha_heat)/2*standard_cgrav * (s% star_mass *Msun) *&
        !(s% mstar_dot) / (s% r(1))  /Lsun) ![erg/s]

!$OMP PARALLEL DO PRIVATE(k)
        do k = 1, s% nz
        if (s% star_mdot > 0.0d0 .and. s% m(k)>=(1.0d0-xpos)*(s% star_mass)*Msun) then
        s% extra_heat(k)%val = 2.0d0*Ladd/(Mstar*pow2(xpos)) * &
            (s% m(k)/Mstar-(1.0d0-xpos))
        !write(*,*) 's% extra_heat(k)%val', s% extra_heat(k)%val
        !write(*,*) 'xpos', xpos
        end if
        end do
!$OMP END PARALLEL DO

        end subroutine energy_routine


      end module run_star_extras
      
