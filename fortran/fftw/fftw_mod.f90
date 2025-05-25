module fftw_mod
   use,intrinsic :: iso_c_binding
   include 'fftw3.f03'

   integer::Np
   integer, parameter :: Nmax = 1024
   type(C_PTR) :: plan
   real(kind=8)::factor
end module fftw_mod
