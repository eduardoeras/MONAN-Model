module bl_taylor_kzz

    use ccpp_kinds,only: kind_phys
    use mpas_log, only: mpas_log_write

    implicit none

    private

    public :: vertical_diffusivity

contains

    subroutine vertical_diffusivity(kzz_h, kzz_q, kzz_m, gamma, stability_flag, qfx, tx, ustar, hpbl, sshfx, govrth, rho, cp ,karman , &
                                    its, ite, kts, kte, kpbl, za, cori, itimestep)

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
        !  (12) Degrazia, G. A., Rizza, U., Mangia, C., & Tirabassi, T. (1997).         !
        !       Validation of a new turbulent parameterization for dispersion models    !
        !       in convective conditions. Boundary-Layer Meteorology, 85(2), 243-254.   !
        !       Equation (18).                                                          !
        !                                                                               !
        !  (13) Basu, S., & Holtslag, A. A. (2021). Turbulent Prandtl number and        !
        !       characteristic length scales in stably stratified flows: steady-state   !
        !       analytical solutions. Environmental Fluid Mechanics, 21(6), 1273-1302.  !
        !       Page 1282. Equation (2).                                                !
        !                                                                               !
        !  (14) Foken, T. (2006). 50 years of the Monin-Obukhov similarity theory.      !
        !       Boundary-Layer Meteorology, 119(3), 431-447.                            !
        !       Equations (21a, 21b, 22a, 22b).                                         !
        !                                                                               !
        !  (15) Troen, I. B., & Mahrt, L. (1986). A simple model of the atmospheric     !
        !       boundary layer; sensitivity to surface evaporation. Boundary-Layer      !
        !       Meteorology, 37(1), 129-148. Equation (17).                             !
        !                                                                               !
        !  (16) Grachev, A. A., Andreas, E. L., Fairall, C. W., Guest, P. S., &         !
        !       Persson, P. O. G. (2007). On the turbulent Prandtl number in the        !
        !       stable atmospheric boundary layer. Boundary-layer meteorology,          !
        !       125(2), 329-341. Equation (1).                                          !
        !                                                                               !
        !  (17) Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971).   !
        !       Flux-profile relationships in the atmospheric surface layer. Journal    !
        !       of Atmospheric Sciences, 28(2), 181-189. Equation (26).                 !
        !                                                                               !
        !  (18) Degrazia, G. A., Anfossi, D., Carvalho, J. C., Mangia, C., Tirabassi,   !
        !       T., & Velho, H. C. (2000). Turbulence parameterisation for PBL          !
        !       dispersion models in all stability conditions. Atmospheric              !
        !       environment, 34(21), 3575-3583. Equation (28).                          !
        !                                                                               !   
        !-------------------------------------------------------------------------------

        ! Inputs
        real(kind=kind_phys), dimension( its:ite, kts:kte ), intent(out  )   ::                 kzz_h ! Vertical diffusivity for heat [m^2/s]
        real(kind=kind_phys), dimension( its:ite, kts:kte ), intent(out  )   ::                 kzz_q ! Vertical diffusivity for moisture [m^2/s]
        real(kind=kind_phys), dimension( its:ite, kts:kte ), intent(out  )   ::                 kzz_m ! Vertical diffusivity for momentum [m^2/s]
        real(kind=kind_phys), dimension( its:ite, kts:kte ), intent(out  )   ::                 gamma ! Counter-gradient [quantity * m/s]
        integer,              dimension( its:ite, kts:kte ), intent(out  )   ::        stability_flag ! PBL stability flag: -1=UNSTABLE, 0=NEUTRAL, 1=STABLE
        real(kind=kind_phys), dimension( its:ite ),          intent(in   )   ::                   qfx ! Surface Moisture Flux [kg m^-2 s^-1]
        real(kind=kind_phys), dimension( its:ite ),          intent(in   )   ::                    tx ! Lowest model-level air temperature [K]
        real(kind=kind_phys), dimension( its:ite ),          intent(in   )   ::                 ustar ! Friction velocity
        real(kind=kind_phys), dimension( its:ite ),          intent(in   )   ::                  hpbl ! PBL height
        real(kind=kind_phys), dimension( its:ite ),          intent(in   )   ::                 sshfx ! Surface Sensible Heat Flux
        real(kind=kind_phys), dimension( its:ite ),          intent(in   )   ::                govrth ! Buoyancy parameter (g/theta_0)
        real(kind=kind_phys), dimension( its:ite ),          intent(in   )   ::                   rho ! Air density
        real(kind=kind_phys),                                intent(in   )   ::                    cp ! Specific heat at constant pressure
        real(kind=kind_phys),                                intent(in   )   ::                karman ! Von Kármán constant
        integer,                                             intent(in   )   ::               its,ite ! Vertical index range for the current column
        integer,                                             intent(in   )   ::               kts,kte ! Vertical index range for the current column
        integer,              dimension( its:ite ),          intent(in   )   ::                  kpbl ! Vertical index of the PBL height
        real(kind=kind_phys), dimension( its:ite, kts:kte ), intent(in   )   ::                    za ! Height above the ground [m]
        real(kind=kind_phys),                                intent(in   )   ::                  cori ! Coriolis parameter

        integer,                                             intent(in   )   ::             itimestep ! Current time step for debugging purposes

        ! Local variables
        character(len=10   ), dimension( its:ite, kts:kte )                  ::             stability ! PBL stability condition 
        real(kind=kind_phys), dimension( its:ite, kts:kte )                  ::                   kzz ! Vertical diffusivity
        real(kind=kind_phys), dimension( its:ite )                           ::                     L ! Monin-Obukhov length [m]
        real(kind=kind_phys), dimension( its:ite, kts:kte )                  ::                lambda ! Local Monin-Obukhov length [m]
        real(kind=kind_phys), dimension( its:ite, kts:kte )                  ::                 zovrh ! Height above the ground normalized by PBL height
        real(kind=kind_phys), dimension( its:ite, kts:kte )                  ::                 zovrl ! Stability parameter
        real(kind=kind_phys), dimension( its:ite )                           ::                wtheta ! Surface kinematic heat flux [k m s^{-1}]
        real(kind=kind_phys), dimension( its:ite )                           ::                 wstar ! Local Convective velocity scale
        real(kind=kind_phys), dimension( its:ite )                           ::           ustar_local ! Range limited friction velocity
        real(kind=kind_phys)                                                 ::             sshfx_min ! Surface Sensible Heat Flux lower threshold
        real(kind=kind_phys)                                                 ::            lambda_max ! Local Monin-Obukhov length top threshold
        real(kind=kind_phys)                                                 ::             ustar_min ! Friction velocity threshold
        real(kind=kind_phys)                                                 ::   stability_threshold ! Stability Threshold for neutral condition
        real(kind=kind_phys)                                                 ::      alpha_1, alpha_2 ! Experimental exponents for SBL
        real(kind=kind_phys)                                                 ::              exponent ! Exponent for local Monin-Obukhov length

        real(kind=kind_phys), dimension( its:ite, kts:kte )                  ::               prandtl ! Prandtl number
        real(kind=kind_phys), dimension( its:ite, kts:kte )                  ::                 phi_h ! Dimensionless gradients of heat
        real(kind=kind_phys), dimension( its:ite, kts:kte )                  ::                 phi_m ! Dimensionless gradients of momentum

        integer                                                              ::                  i, k
        character(len=512)                                                   ::             debug_msg
        integer, parameter                                                   :: debug_itimestep = 1000! Time step at which to print debug information

        ! Flag to indicate whether moisture flux contributes to buoyancy
        integer, parameter :: moist_contribution_flag = 0


        ! DEBUG ------------------------------------------------------------------------!
        !do i = its, ite
        !    do k = kts, kpbl(i) - 1
        !        write(debug_msg, "(A, F23.9, &
        !                           A, F23.9, &
        !                           A, F23.9, &
        !                           A, F23.9, &
        !                           A, F23.9, &
        !                           A, I, &
        !                           A, F23.9)" ) &
        !        " ustar = ", ustar(i), &
        !        ", hpbl = ", hpbl(i), &
        !        ", sshfx = ", sshfx(i), &
        !        ", govrth = ", govrth(i), &
        !        ", rho = ", rho(i), &
        !        ", kpbl = ", kpbl(i), &
        !        ", za = ", za(i,k)
        !        call mpas_log_write('[TAYLOR DEBUG]'//trim(debug_msg))
        !    end do
        !end do
        !------------------------------------------------------------------------ DEBUG ! 

        ! Surface kinematic heat flux [k m s^{-1}]. References: 1, 2
        ! Set a floor value for sshfx to avoid numerical issues.
        sshfx_min = 1.0e-3
        do i = its, ite
            if (abs(sshfx(i)) .lt. sshfx_min) then
                ! To use the moisture flux contribution to buoyancy, set moist_contribution_flag to 1.
                if (moist_contribution_flag == 1) then
                    wtheta(i) = sign(sshfx_min, sshfx(i)) / (rho(i) * cp) + (0.61 * tx(i) * qfx(i) / rho(i))
                else
                    wtheta(i) = sign(sshfx_min, sshfx(i)) / (rho(i) * cp)
                end if
            else
                wtheta(i) = (sshfx(i) / (rho(i) * cp))
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

        ! Initialize stability_flag to 0 (NEUTRAL) for all levels
        do i = its, ite
            do k = kts, kte
                stability_flag(i,k) = 0
            end do
        end do

        ! PBL stability condition. Reference: 2
        ! Threshold value set by empirical studies.
        stability_threshold = 0.005
        do i = its, ite
            do k = kts, kpbl(i) - 1
                if (zovrl(i,k) .lt. -stability_threshold) then
                    stability(i,k) = 'UNSTABLE'
                    stability_flag(i,k) = -1
                else if (zovrl(i,k) .gt. stability_threshold) then
                    stability(i,k) = 'STABLE'
                    stability_flag(i,k) = 1
                else
                    stability(i,k) = 'NEUTRAL'
                    stability_flag(i,k) = 0
                end if
            end do
        end do

        ! Vertical diffusivity main loop
        do i = its, ite
            do k = kts, kpbl(i) - 1

                ! Unstable PBL
                if (trim(stability(i,k)) == 'UNSTABLE') then

                    ! Convective velocity scale. References: 10, 11
                    wstar(i) = (govrth(i) * wtheta(i) * hpbl(i)) ** (1.0 / 3.0)

                    ! Vertical diffusivity kzz. Reference: 12
                    kzz(i,k) = (0.22 * wstar(i) * hpbl(i)) * &
                               (zovrh(i,k) ** (1.0/3.0)) * &
                               ((1.0 - zovrh(i,k)) ** (1.0/3.0)) * &
                               (1.0 - exp(- 4.0 * zovrh(i,k)) - &
                               0.0003 * exp(8.0 * zovrh(i,k)))

                    ! Prandtl number for CBL. Reference: 17
                    prandtl(i,k) = 0.74 

                ! Stable PBL
                elseif (trim(stability(i,k)) == 'STABLE') then

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

                    ! Prandtl number for SBL. References: 14, 16
                    phi_h(i,k) = 0.95 + 7.8 * zovrl(i,k)
                    phi_m(i,k) = 1.0 + 6.0 * zovrl(i,k)
                    prandtl(i,k) = phi_h(i,k) / phi_m(i,k)

                    ! Prandtl number is bounded by the asymptotic limits of unity and four. Reference: 15
                    if (prandtl(i,k) .lt. 1.0) then
                        prandtl(i,k) = 1.0
                    else if (prandtl(i,k) .gt. 4.0) then
                        prandtl(i,k) = 4.0
                    end if
                    
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
                
                    ! Vertical diffusivity kzz. Reference: 18
                    kzz(i,k) = (0.4 * (1.0 - zovrh(i,k)) ** 0.85 * ustar_local(i) * za(i,k)) / &
                               (1.0 + 15.0 * cori * za(i,k) / ustar_local(i)) ** (4.0 / 3.0)

                    ! Prandtl number for small values of Ri. Reference: 14
                    prandtl(i,k) = 0.85
                    
                end if

            end do
            
            ! Values above PBL height must be zero
            do k = kpbl(i), kte
                kzz(i,k) = 0.0
                prandtl(i,k) = 0.0
            end do

        end do

        ! Heat and moisture diffusivity are assumed to be the same. Reference: 15
        ! Momentum diffusivity is given by heat diffusivity times the Prandtl number. References: 13
        do i = its,ite
            do k = kts,kte
                kzz_h(i,k) = kzz(i,k)
                kzz_q(i,k) = kzz(i,k)
                kzz_m(i,k) = kzz(i,k) * prandtl(i,k)
            end do
        end do

        ! Counter-gradient term.
        do i = its, ite
            do k = kts, kpbl(i) - 1
                if (trim(stability(i,k)) == 'UNSTABLE') then
                    call counter_gradient(gamma(i,k), zovrh(i,k), wstar(i), hpbl(i))
                end if
            end do
        end do

        ! DEBUG ------------------------------------------------------------------------!
        !do i = its, ite
        !    do k = kts, kpbl(i) - 1
        !        write(debug_msg, "(A, F23.9, &
        !                           A, F23.9, &
        !                           A, F23.9, &
        !                           A, F23.9, &
        !                           A, F23.9, &
        !                           A, F23.9, &
        !                           A, F23.9, &
        !                           A, F23.9, &
        !                           A, A)" ) &
        !        " zh = ", zovrh(i,k), &
        !        ", sshfx = ", sshfx(i), &
        !        ", ustar = ", ustar(i), &
        !        ", ustar_local = ", ustar_local(i), &
        !        ", kzz_h = ", kzz_h(i,k), &
        !        ", kzz_q = ", kzz_q(i,k), &
        !        ", kzz_m = ", kzz_m(i,k), &
        !        ", prandtl = ", prandtl(i,k), &
        !        ", stability = ", stability(i,k)
        !        call mpas_log_write('[TAYLOR DEBUG]'//trim(debug_msg))
        !    end do
        !end do
        !------------------------------------------------------------------------ DEBUG ! 

        ! DEBUG ------------------------------------------------------------------------!
        do i = its, ite
            do k = kts, kpbl(i) - 1
                if (trim(stability(i,k)) .eq. 'STABLE'.and. MOD(itimestep, debug_itimestep) .ne. 0) then
                    cycle
                end if
                write(debug_msg, "(A, F23.9, &
                                A, F23.9, &
                                A, F23.9, &
                                A, F23.9, &
                                A, F23.9, &
                                A, F23.9, &
                                A, F23.9, &
                                A, F23.9, &
                                A, F23.9, &
                                A, A)" ) &
                " zh = ", zovrh(i,k), &
                ", kzz_h = ",  kzz_h(i,k), &
                ", kzz_q = ",  kzz_q(i,k), &
                ", kzz_m = ",  kzz_m(i,k), &
                ", L = ",      L(i), &
                ", lambda = ", lambda(i,k), &
                ", zovrl = ", zovrl(i,k), &
                ", wtheta = ", wtheta(i), &
                ", wstar = ", wstar(i), &
                ", stability = ", stability(i,k)
                call mpas_log_write('[TAYLOR DEBUG]'//trim(debug_msg))
            end do
        end do
        !------------------------------------------------------------------------ DEBUG ! 

    end subroutine vertical_diffusivity

    subroutine counter_gradient(gamma, zovrh, wstar, hpbl)

        !-------------------------------------------------------------------------------!
        implicit none
        !-------------------------------------------------------------------------------!
        !                                                                               !
        !   References:                                                                 !
        !                                                                               !
        !   (1) Campos Velho, H.F., et al. (1998). New parameterizations in RAMS for    !
        !       vertical turbulent fluxes, p. 17                                        !
        !                                                                               !
        !   (2) Roberti, D. R., de Campos Velho, H. F., & Degrazia, G. A. (2004).       !
        !       Identifying counter-gradient term in atmospheric convective boundary    !
        !       layer. Inverse Problems in Science and Engineering, 12(3), 329-339.     !
        !                                                                               !
        !   (3) Welter, M. E. S. (2006). Modelagem do termo de contra gradiente na      !
        !       parametrizacao de turbulencia no modelo atmosferico brams. Mestrado em  !
        !       ciencia da computacao, Instituto Nacional de Pesquisas Espaciais,       !
        !       Sao Jose dos Campos.                                                    !
        !                                                                               !
        !   (4) Eras, E.R. (2024). Turbulence modeling based on Taylor's theory for     !
        !       INPE's BAM model. PhD Thesis, INPE, p. 46-49                            !
        !                                                                               !
        !   (5) Sorbjan, Z. (1986). Local similarity of spectral and cospectral         !
        !       characteristics in the stable-continuous boundary layer. Boundary-layer !
        !       meteorology, 35(3), 257-275.                                            !
        !                                                                               !
        !   (6) Degrazia, G. A., Rizza, U., Mangia, C., & Tirabassi, T. (1997).         !
        !       Validation of a new turbulent parameterization for dispersion models    !
        !       in convective conditions. Boundary-Layer Meteorology, 85(2), 243-254.   !
        !                                                                               !
        !   (7) Andrei Nikolaevich Kolmogorov, V. Levin, Julian Charles Roland Hunt,    !
        !       Owen Martin Phillips, David Williams; The local structure of turbulence !
        !       in incompressible viscous fluid for very large Reynolds numbers. Proc.  !
        !       A 1 July 1991; 434 (1890): 9-13.                                        !
        !                                                                               !
        !-------------------------------------------------------------------------------!

        ! Inputs
        real(kind=kind_phys),                                intent(out  )   ::   gamma ! Counter-gradient [quantity * m/s]
        real(kind=kind_phys),                                intent(in   )   ::   zovrh ! Height above the ground normalized by PBL height
        real(kind=kind_phys),                                intent(in   )   ::   wstar ! Local Convective velocity scale
        real(kind=kind_phys),                                intent(in   )   ::    hpbl ! PBL height

        ! Local variables
        real(kind=kind_phys)                                              :: kolmogorov ! kolmogorov exponent
        real(kind=kind_phys)                                              ::        q_w ! Dimensionless stability function
        real(kind=kind_phys)                                              ::    sigma_w ! Vertical velocity variance
        real(kind=kind_phys)                                              ::          l ! Mixing length

        ! Local constants
        real(kind=kind_phys), parameter ::    beta = 0.02     ! Calibration constant. Reference (3).
        real(kind=kind_phys), parameter ::     psi = 0.913    ! Dimensionless dissipation rate. Reference (2).
        real(kind=kind_phys), parameter ::     c_w = 0.4      ! Empirical constant. Reference (6, 4).
        real(kind=kind_phys), parameter ::    fm_w = 0.33     ! Frequency of the spectral peak for neutral conditions. Reference (6, 5).
   
        ! DEBUG ------------------------------------------------------------------------!
        !call mpas_log_write('[TAYLOR DEBUG] Inside GAMMA function')
        !------------------------------------------------------------------------ DEBUG ! 
        
        ! kolmogorov exponent. Reference: 7
        kolmogorov = 2.0 / 3.0

        ! Dimensionless stability function. Reference: 6
        q_w = 1.68 * zovrh / (1.0 - exp(-4.0 * zovrh) - 0.0003 * exp(8.0 * zovrh))

        ! Vertical velocity variance. Reference: 6
        sigma_w = sqrt(((0.98 * c_w) / (fm_w ** kolmogorov)) * &
                ((psi / q_w) ** kolmogorov) * &
                (zovrh ** kolmogorov) * (wstar ** 2.0))

        ! Mixing length. Reference: 1, 3
        l = 0.2 * hpbl * (1.0 - exp(-4.0 * zovrh) - 0.003 * exp(8.0 * zovrh))

        ! Mean quantity of a property. Reference: 2

    end subroutine counter_gradient

end module bl_taylor_kzz