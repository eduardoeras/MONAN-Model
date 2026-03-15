module bl_taylor_kzz

    use ccpp_kinds,only: kind_phys
    use mpas_log, only: mpas_log_write

    implicit none

    public :: vertical_diffusivity

contains

    subroutine vertical_diffusivity(kzz, ustar, hpbl, sshfx, govrth, rho, cp ,karman , &
                                    its, ite, kts, kte, kpbl, za, itimestep)

        !-------------------------------------------------------------------------------
        implicit none
        !-------------------------------------------------------------------------------
        !
        !   Taylor PBL Parameterization expected inputs:
        !
        !   For Stable and Neutral PBL:
        !       u* -> Friction Velocity [ms^{-1}]
        !       h -> PBL height [m]
        !       z -> height above the groud [m]
        !       k -> Von Kármán constant g -> gravity [ms^{-2}]
        !       theta_v -> virtual potential temperature near the surface [k]
        !       w'Theta' -> Surface kinematic heat flux [k m s^{-1}]

        !   For convective PBL:
        !       w* -> Convective velocity scale [m/s]
        !       h -> PBL height [m]
        !       z -> height above the groud [m]
        !       L -> Monin-Obukhov length [m]
        !       g -> gravity [ms^{-2}]
        !       theta_0 -> Reference virtual potential temperature near the surface [k]
        !       w'Theta' -> Surface kinematic heat flux [k m s^{-1}]
        !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        !
        !   Model coordinates reference:
        !
        !   its, ite : horizontal tile bounds. Physics loops over columns using:
        !              do i = its, ite
        !              Each i corresponds to one atmospheric column.
        !  
        !   kts, kte : vertical tile bounds. Vertical loops typically use:
        !              do k = kts, kte
        !              Each k corresponds to one vertical model level.
        !
        !-------------------------------------------------------------------------------

        ! Inputs
        real(kind=kind_phys), dimension( its:ite, kts:kte ), intent(out  )   ::                   kzz ! Vertical diffusivity [m^2/s]
        real(kind=kind_phys), dimension( its:ite ),          intent(in   )   ::                 ustar ! Friction velocity
        real(kind=kind_phys), dimension( its:ite ),          intent(in   )   ::                  hpbl ! PBL height
        real(kind=kind_phys), dimension( its:ite ),          intent(in   )   ::                 sshfx ! Surface Sensible Heat Flux
        real(kind=kind_phys), DIMENSION( its:ite ),          intent(in   )   ::                govrth ! Buoyancy parameter (g/theta_0)
        real(kind=kind_phys), dimension( its:ite ),          intent(in   )   ::                   rho ! Air density
        real(kind=kind_phys),                                intent(in   )   ::                    cp ! Specific heat at constant pressure
        real(kind=kind_phys),                                intent(in   )   ::                karman ! Von Kármán constant
        integer,                                             intent(in   )   ::               its,ite ! Vertical index range for the current column
        integer,                                             intent(in   )   ::               kts,kte ! Vertical index range for the current column
        integer,              dimension( its:ite ),          intent(in   )   ::                  kpbl ! Vertical index of the PBL height
        real(kind=kind_phys), dimension( its:ite, kts:kte ), intent(in   )   ::                    za ! Height above the ground [m]

        integer,                                             intent(in   )   ::             itimestep ! Current time step for debugging purposes

        ! Local variables
        character(len=10   ), dimension( its:ite, kts:kte )                  ::             stability ! PBL stability condition 
        real(kind=kind_phys), dimension( its:ite )                           ::                     L ! Monin-Obukhov length [m]
        real(kind=kind_phys), dimension( its:ite, kts:kte )                  ::                lambda ! Local Monin-Obukhov length [m]
        real(kind=kind_phys), dimension( its:ite, kts:kte )                  ::                 zovrh ! Height above the ground normalized by PBL height
        real(kind=kind_phys), dimension( its:ite, kts:kte )                  ::                 zovrl ! Stability parameter
        real(kind=kind_phys), dimension( its:ite )                           ::                wtheta ! Surface kinematic heat flux [k m s^{-1}]
        real(kind=kind_phys), dimension( its:ite )                           ::                 wstar ! Convective velocity scale
        real(kind=kind_phys), dimension( its:ite, kts:kte )                  ::                   psi ! Nondimensional molecular dissipation rate
        real(kind=kind_phys), dimension( its:ite )                           ::           ustar_local ! Range limited friction velocity
        
        real(kind=kind_phys)                                                 ::             ustar_min ! Friction velocity threshold
        real(kind=kind_phys)                                                 ::   stability_threshold ! Stability Threshold for neutral condition
        real(kind=kind_phys)                                                 ::      alpha_1, alpha_2 ! Experimental exponents for SBL
        real(kind=kind_phys)                                                 ::              exponent ! Exponent for local Monin-Obukhov length

        integer                                                              ::                  i, k
        character(len=256)                                                   ::             debug_msg
        integer, parameter                                                   :: debug_itimestep = 256 ! Time step at which to print debug information

        !DEBUG!=========================================================================!
        !if (itimestep .eq. debug_itimestep) then

            !do i = its, ite
            !    write(debug_msg, "(F10.3)" ) ustar(i)
            !    call mpas_log_write('[TAYLOR DEBUG] ustar = '//trim(debug_msg))
            !end do

            !do i = its, ite
            !    write(debug_msg, "(F10.3)" ) hpbl(i)
            !    call mpas_log_write('[TAYLOR DEBUG] hpbl = '//trim(debug_msg))
            !end do

            !do i = its, ite
            !    if (sshfx(i) .gt. 0.0) then
            !        write(debug_msg, "(F10.3)" ) sshfx(i)
            !        call mpas_log_write('[TAYLOR DEBUG] sshfx = '//trim(debug_msg))
            !    end if
            !end do

            !do i = its, ite
            !    write(debug_msg, "(F10.3)" ) govrth(i)
            !    call mpas_log_write('[TAYLOR DEBUG] govrth = '//trim(debug_msg))
            !end do
    
            !do i = its, ite
            !    write(debug_msg, "(F10.3)" ) rho(i)
            !    call mpas_log_write('[TAYLOR DEBUG] rho = '//trim(debug_msg))
            !end do

            !write(debug_msg, "(F10.3)" ) cp
            !call mpas_log_write('[TAYLOR DEBUG] cp = '//trim(debug_msg))

            !write(debug_msg, "(F10.3)" ) karman
            !call mpas_log_write('[TAYLOR DEBUG] karman = '//trim(debug_msg))

            !i = 1
            !do k = kts, kte
            !    write(debug_msg, "(F10.3)" ) za(i,k)
            !    call mpas_log_write('[TAYLOR DEBUG] za = '//trim(debug_msg))
            !end do

        !end if
        !=========================================================================!DEBUG!

        ! Surface kinematic heat flux [k m s^{-1}]
        do i = its, ite
            wtheta(i) = sshfx(i) / (rho(i) * cp)
        end do

        ! Set a floor value for u* (Prairie Grass experiment)
        ustar_min = 0.05
        do i = its, ite
            ustar_local(i) = max(ustar(i), ustar_min)
        end do

        ! Monin-Obukhov length L [M]
        do i = its, ite
            L(i) = -(ustar_local(i) ** 3) / (karman * govrth(i) * wtheta(i))
        end do  

        ! Nondimensional vertical coordinate z / h
        do i = its, ite
            do k = kts, kpbl(i) - 1
                zovrh(i,k) = za(i,k) / hpbl(i)
            end do
        end do

        ! Stability parameter z / L
        do i = its, ite
            do k = kts, kpbl(i) - 1
                zovrl(i,k) = za(i,k) / L(i)
            end do
        end do

        ! PBL stability condition
        stability_threshold = 0.001
        do i = its, ite
            do k = kts, kpbl(i) - 1
                if (zovrl(i,k) .lt. -stability_threshold) then
                    stability(i,k) = 'UNSTABLE'
                else if (zovrl(i,k) .gt. stability_threshold) then
                    stability(i,k) = 'STABLE'
                else
                    stability(i,k) = 'NEUTRAL'
                end if
            end do
        end do

        ! Vertical diffusivity main loop
        do i = its, ite
            do k = kts, kpbl(i) - 1

                ! Stable PBL
                if (trim(stability(i,k)) == 'STABLE') then

                    ! Minnesota experiment
                    alpha_1 = 2.0
                    alpha_2 = 3.0

                    ! Local Monin-Obukhov length [m]
                    exponent = (3.0 / 2.0) * alpha_1 - alpha_2
                    lambda(i,k) = abs(L(i)) * (1.0 - zovrh(i,k)) ** exponent
                
                    ! Vertical diffusivity kzz
                    kzz(i,k) = (ustar_local(i) * hpbl(i)) * &
                            0.33 * ((1.0 - zovrh(i,k)) ** (alpha_1 / 2.0)) * zovrh(i,k) / &
                            1.0 + 3.7 * zovrh(i,k) * (hpbl(i) / lambda(i,k))

                ! Unstable PBL
                elseif (trim(stability(i,k)) == 'UNSTABLE') then

                    ! Convective velocity scale
                    wstar(i) = (govrth(i) * wtheta(i) * hpbl(i)) ** (1.0 / 3.0)

                    ! Nondimensional molecular dissipation rate
                    psi(i,k) = ((1.0 - zovrh(i,k)) ** 2.0 * &
                               (- zovrl(i,k)) ** (-2.0 / 3.0) + &
                               0.75) ** (3.0 / 2.0)

                    ! Vertical diffusivity kzz
                    kzz(i,k) = 0.15 * psi(i,k) ** (1.0/3.0) * &
                               (1 - exp(- 4 * zovrh(i,k)) - &
                               0.0003 * exp(8 * zovrh(i,k))) ** (4.0/3.0) * &
                               wstar(i) * hpbl(i)
                    
                ! Neutral PBL
                else

                    ! Cabauw experiment
                    alpha_1 = 1.5
                    alpha_2 = 1.0

                    ! Local Monin-Obukhov length [m]
                    exponent = (3.0 / 2.0) * alpha_1 - alpha_2
                    lambda(i,k) = abs(L(i)) * (1.0 - zovrh(i,k)) ** exponent
                
                    ! Vertical diffusivity kzz
                    kzz(i,k) = (ustar_local(i) * hpbl(i)) * &
                            0.33 * ((1.0 - zovrh(i,k)) ** (alpha_1 / 2.0)) * zovrh(i,k) / &
                            1.0 + 3.7 * zovrh(i,k) * (hpbl(i) / lambda(i,k))
                    
                end if

                ! DEBUG
                write(debug_msg, "(A, F14.6, A, F14.6, A, F14.6, A, F14.6, A, F14.6, A, F14.6, A, A)" ) &
                " zh = ", zovrh(i,k), &
                ", sshfx = ", sshfx(i), &
                ", ustar = ", ustar(i), &
                ", ustar_local = ", ustar_local(i), &
                ", lambda = ", lambda(i,k), &
                ", kzz = ", kzz(i,k), &
                ", stability = ", stability(i,k)
                call mpas_log_write('[TAYLOR DEBUG]'//trim(debug_msg))

            end do
            
            ! Values above PBL height must be zero
            do k = kpbl(i), kte
                kzz(i,k) = 0.0
            end do

        end do

        !DEBUG!=========================================================================!

        !do i = its, ite
        !    if (wtheta(i) .gt. 0.0) then
        !        write(debug_msg, "(F14.6)" ) wtheta(i)
        !        call mpas_log_write('[TAYLOR DEBUG] wtheta = '//trim(debug_msg))
        !    end if
        !end do
        
        !if (itimestep .eq. debug_itimestep) then

        !    do i = its, ite
        !        write(debug_msg, "(A,F14.6,A,A)" ) ' L = ', L(i), ', ', stability(i)
        !        call mpas_log_write('[TAYLOR DEBUG]' //trim(debug_msg))
        !    end do

        !    do i = its, ite
        !        do k = kts, kpbl(i) - 1
        !            write(debug_msg, "(F14.6)" ) lambda(i,k)
        !            call mpas_log_write('[TAYLOR DEBUG] final lambda = '//trim(debug_msg))
        !        end do
        !    end do

        !    do i = its, ite
        !        do k = kts, kpbl(i) - 1
        !            write(debug_msg, "(F14.6)" ) kzz(i,k)
        !            call mpas_log_write('[TAYLOR DEBUG] final kzz = '//trim(debug_msg))
        !        end do
        !    end do

        !end if

        !=========================================================================!DEBUG!

    end subroutine vertical_diffusivity

end module bl_taylor_kzz