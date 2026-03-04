module bl_taylor_kzz

    use mpas_log, only: mpas_log_write

    implicit none

    public :: vertical_diffusivity

contains

    subroutine vertical_diffusivity()

        call mpas_log_write('[TAYLOR DEBUG] inside vertical_diffusivity subroutine ')

    end subroutine vertical_diffusivity

end module bl_taylor_kzz