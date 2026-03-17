module bl_taylor_kzz

    use ccpp_kinds,only: kind_phys
    use mpas_log, only: mpas_log_write

    implicit none

    public :: vertical_diffusivity

contains

    subroutine vertical_diffusivity(kzz, ustar, hpbl, sshfx, govrth, rho, cp ,karman , &
                                    its, ite, kts, kte, kpbl, za, itimestep)

        !-------------------------------------------------------------------------------!
        implicit none
        !-------------------------------------------------------------------------------!
        !                                                                               !
        !   Taylor PBL Parameterization expected inputs:                                !
        !                                                                               !
        !   For Stable and Neutral PBL:                                                 !
        !       u* -> Friction Velocity [ms^{-1}]                                       !
        !       h -> PBL height [m]                                                     !
        !       z -> height above the groud [m]                                         !
        !       k -> Von Karman constant g -> gravity [ms^{-2}]                         !
        !       theta_v -> virtual potential temperature near the surface [k]           !
        !       w'Theta' -> Surface kinematic heat flux [k m s^{-1}]                    !
        !                                                                               !
        !   For convective PBL:                                                         !
        !       w* -> Convective velocity scale [m/s]                                   !
        !       h -> PBL height [m]                                                     !
        !       z -> height above the groud [m]                                         !
        !       L -> Monin-Obukhov length [m]                                           !
        !       g -> gravity [ms^{-2}]                                                  !
        !       theta_0 -> Reference virtual potential temperature near the surface [k] !
        !       w'Theta' -> Surface kinematic heat flux [k m s^{-1}]                    !
        !                                                                               !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
        !                                                                               !
        !   Model coordinates reference:                                                !
        !                                                                               !
        !   its, ite : horizontal tile bounds. Physics loops over columns using:        !
        !              do i = its, ite                                                  !
        !              Each i corresponds to one atmospheric column.                    !
        !                                                                               !
        !   kts, kte : vertical tile bounds. Vertical loops typically use:              !
        !              do k = kts, kte                                                  !
        !              Each k corresponds to one vertical model level.                  !
        !                                                                               !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
        !                                                                               !
        !   References:                                                                 !
        !                                                                               !
        !   (1) Foken, T. (2006). 50 Years of the Monin-Obukhov Similarity Theory.      !
        !       Boundary-Layer Meteorology, 119(3), 431-447. Equation (2).              !
        !                                                                               !
        !   (2) Monin, A. S., and Obukhov, A. M. (1954). Basic laws of turbulent        !
        !       mixing in the surface layer of the atmosphere. Tr. Akad. Nauk SSSR      !
        !       Geofiz. Inst., 24(151), 163-187. Equations (25, 28, 29).                !
        !                                                                               !
        !   (3) Goulart, A., et al (2024). "Evaluation of the eddy diffusivity in a     !
        !       pollutant dispersion model in the planetary boundary layer."            !
        !       Environmental Fluid Mechanics 24.3, 489-504. Table (1).                 !
        !                                                                               !
        !   (4) Holtslag, Albert AM, and Frans TM Nieuwstadt (1986). "Scaling the       !
        !       atmospheric boundary layer." Boundary-Layer Meteorology 36.1, 201-209.  !
        !       Equation (1).                                                           !
        !                                                                               !
        !   (5) Holtslag, A. A. M., and Nieuwstadt, F. T. M. (1986). Scaling the        !
        !       atmospheric boundary layer. Boundary-Layer Meteorology, 36(1-2),        !
        !       201-209. Figure (1).                                                    !
        !                                                                               !
        !   (6) Olesen, H. R., et al. (1984). Modelling velocity spectra in the         !
        !       lower part of the planetary boundary layer. Boundary-Layer Meteorology, !
        !       29, 285-312. Figure (2).                                                !
        !                                                                               !
        !   (7) Sorbjan, Z. (1985). On similarity in the atmospheric boundary layer.    !
        !       Boundary-Layer Meteorology, 34(4), 377-397. Page (380).                 !
        !                                                                               !
        !   (8) Degrazia, G. A., and Moraes, O. L. L. (1992). A model for eddy          !
        !       diffusivity in a stable boundary layer. Boundary-Layer Meteorology,     !
        !       58(3), 205-214. Equations (3, 23).                                      !
        !                                                                               !
        !   (9) Nieuwstadt, F. T. M. (1984). The Turbulent Structure of the Stable,     !
        !       Nocturnal Boundary Layer. Journal of the Atmospheric Sciences, 41(14),  !
        !       2202-2216. Equations (30, 36, 37).                                      !
        !                                                                               !
        !  (10) Campos Velho, H. F., Degrazia, G. A., Holtslag, A. A. M., and Pielke,   !
        !       R. A. (1998). New parameterizations in RAMS for vertical turbulent      !
        !       fluxes. Technical Report, Colorado State University.                    !
        !       Equations (22, 26, 29, 49).                                             !
        !                                                                               !
        !  (11) Deardorff, J. W. (1970). Convective Velocity and Temperature Scales     !
        !       for the Unstable Planetary Boundary Layer and for Rayleigh Convection.  !
        !       Journal of the Atmospheric Sciences, 27(8), 1211-1213. Equation (1).    !
        !                                                                               !
        !  (12) Hojstrup, J. (1982). Velocity spectra in the unstable planetary         !
        !       boundary layer. Journal of the Atmospheric Sciences, 39(10), 2239-2248. !
        !       Equation (8.2).                                                         !
        !                                                                               !
        !  (13) Degrazia, G. A., Rizza, U., Mangia, C., & Tirabassi, T. (1997).         !
        !       Validation of a new turbulent parameterization for dispersion models    !
        !       in convective conditions. Boundary-Layer Meteorology, 85(2), 243-254.   !
        !       Equation (15).                                                          !
        !                                                                               !
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
        real(kind=kind_phys)                                                 ::             sshfx_min ! Surface Sensible Heat Flux lower threshold
        real(kind=kind_phys)                                                 ::            lambda_max ! Local Monin-Obukhov length top threshold
        real(kind=kind_phys)                                                 ::             ustar_min ! Friction velocity threshold
        real(kind=kind_phys)                                                 ::   stability_threshold ! Stability Threshold for neutral condition
        real(kind=kind_phys)                                                 ::      alpha_1, alpha_2 ! Experimental exponents for SBL
        real(kind=kind_phys)                                                 ::              exponent ! Exponent for local Monin-Obukhov length

        integer                                                              ::                  i, k
        character(len=256)                                                   ::             debug_msg
        integer, parameter                                                   :: debug_itimestep = 256 ! Time step at which to print debug information


        ! Surface kinematic heat flux [k m s^{-1}]. References: 1, 2
        ! Set a floor value for sshfx to avoid numerical issues.
        sshfx_min = 1.0e-3
        do i = its, ite
            if (abs(sshfx(i)) .lt. sshfx_min) then
                wtheta(i) = sign(sshfx_min, sshfx(i)) / (rho(i) * cp)
            else
                wtheta(i) = sshfx(i) / (rho(i) * cp)
            end if
        end do

        ! Set a floor value for u* (Prairie Grass experiment). Reference: 3
        ustar_min = 0.05
        do i = its, ite
            ustar_local(i) = max(ustar(i), ustar_min)
        end do

        ! Monin-Obukhov length L [M]. References: 2, 4
        do i = its, ite
            L(i) = -(ustar_local(i) ** 3) / (karman * govrth(i) * wtheta(i))
        end do  

        ! Nondimensional vertical coordinate z / h. References: 5, 6
        do i = its, ite
            do k = kts, kpbl(i) - 1
                zovrh(i,k) = za(i,k) / hpbl(i)
            end do
        end do

        ! Stability parameter z / L. Reference: 2
        do i = its, ite
            do k = kts, kpbl(i) - 1
                zovrl(i,k) = za(i,k) / L(i)
            end do
        end do

        ! PBL stability condition. Reference: 2
        ! Threshold value set the same as in the BAM model Taylor implementation.
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

                    ! Minnesota experiment. References: 7, 8
                    alpha_1 = 2.0
                    alpha_2 = 3.0

                    ! Local Monin-Obukhov length [m]. References: 7, 8, 9
                    exponent = (3.0 / 2.0) * alpha_1 - alpha_2
                    lambda(i,k) = abs(L(i)) * (1.0 - zovrh(i,k)) ** exponent
                    lambda_max = 10.0 * hpbl(i)
                    if (lambda(i,k) .gt. lambda_max) then
                        lambda(i,k) = lambda_max
                    end if
                
                    ! Vertical diffusivity kzz. References: 8, 10
                    kzz(i,k) = (ustar_local(i) * hpbl(i)) * &
                            0.33 * ((1.0 - zovrh(i,k)) ** (alpha_1 / 2.0)) * zovrh(i,k) / &
                            (1.0 + 3.7 * zovrh(i,k) * (hpbl(i) / lambda(i,k)))

                ! Unstable PBL
                elseif (trim(stability(i,k)) == 'UNSTABLE') then

                    ! Convective velocity scale. References: 10, 11
                    wstar(i) = (govrth(i) * wtheta(i) * hpbl(i)) ** (1.0 / 3.0)

                    ! Nondimensional molecular dissipation rate. Reference: 12
                    psi(i,k) = ((1.0 - zovrh(i,k)) ** 2.0 * &
                               (- zovrl(i,k)) ** (-2.0 / 3.0) + &
                               0.75) ** (3.0 / 2.0)

                    ! Vertical diffusivity kzz. References: 10, 13
                    kzz(i,k) = 0.15 * psi(i,k) ** (1.0/3.0) * &
                               (1 - exp(- 4 * zovrh(i,k)) - &
                               0.0003 * exp(8 * zovrh(i,k))) ** (4.0/3.0) * &
                               wstar(i) * hpbl(i)
                    
                ! Neutral PBL
                else

                    ! Cabauw experiment. References: 9, 10
                    alpha_1 = 1.5
                    alpha_2 = 1.0

                    ! Local Monin-Obukhov length [m]. References: 7, 8, 9
                    exponent = (3.0 / 2.0) * alpha_1 - alpha_2
                    lambda(i,k) = abs(L(i)) * (1.0 - zovrh(i,k)) ** exponent
                    lambda_max = 10.0 * hpbl(i)
                    if (lambda(i,k) .gt. lambda_max) then
                        lambda(i,k) = lambda_max
                    end if
                
                    ! Vertical diffusivity kzz. References: 8, 10
                    kzz(i,k) = (ustar_local(i) * hpbl(i)) * &
                            0.33 * ((1.0 - zovrh(i,k)) ** (alpha_1 / 2.0)) * zovrh(i,k) / &
                            (1.0 + 3.7 * zovrh(i,k) * (hpbl(i) / lambda(i,k)))
                    
                end if

                ! DEBUG
                write(debug_msg, "(A, F23.9, A, F23.9, A, F23.9, A, F23.9, A, F23.9, A, F23.9, A, A)" ) &
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

    end subroutine vertical_diffusivity

end module bl_taylor_kzz