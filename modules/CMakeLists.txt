file( GLOB swmp_fortran_modules *.f90 )
add_library(swmp SHARED  ${swmp_fortran_modules})
