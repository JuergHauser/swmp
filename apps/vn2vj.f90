program vn2vj

  use my_types
  use m_inout


type(scalar_field_2d) :: vmod
type(vn2vj_conf) :: conf
  character(len=strlen)     :: arg

call get_command_argument(1,arg)
call read_vn2vj_conf(trim(arg),conf)

vmod%cn=conf%cn
call read_scalar_field_2d_vn(vmod,conf%ifn_vmod)

vmod%y0=vmod%y0-(vmod%ny-2)*vmod%dy
vmod%x0=vmod%x0-vmod%dx

call write_scalar_field_2d(vmod,conf%ofn_vmod)



end program vn2vj
