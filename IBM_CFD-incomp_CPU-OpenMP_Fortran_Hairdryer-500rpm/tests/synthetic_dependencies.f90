! Minimal dependencies for testing the UNMODIFIED reference velocity module.
! No CAD, physical flow solver, or case inputs are used.
module mod_common
 implicit none
 integer,parameter::N1M=7,N2M=7,N3M=7
 integer::NFC_INTP(3)=0,NFC_INNER(3)=0,M=1,MSUB=1
 real*8::X(7),Y(7),Z(7),XMP(7),YMP(7),ZMP(7)
end module

module mod_flowarray
 implicit none
 integer::IFC(1,3),JFC(1,3),KFC(1,3),MPI(1,3,3),INOUT(7,7,7,3)
 real*8::GFI(1,3,0:2,0:2,0:2)
end module

module mod_rotor_motion
 implicit none
 real*8::length_mm=1D0,reference_speed=1D0,omega_now=2D0,wall_x=.25D0
 integer::geometry_mode=1
 logical::moving_wall=.true.,tie_wall=.false.
contains
 real*8 function rotor_sdf(x,y,z)result(phi)
 real*8,intent(in)::x,y,z
 phi=x-wall_x
 if(geometry_mode==2)phi=min(phi,max(.6D0-x,x-.7D0))
 if(geometry_mode==3)phi=-1D0
 end function
 subroutine component_functions(x,y,z,fs,fr,fh)
 real*8,intent(in)::x,y,z
 real*8,intent(out)::fs,fr,fh
 fs=10D0;fr=10D0;fh=10D0
 if(moving_wall)then
 fr=rotor_sdf(x,y,z)
 else
 fs=rotor_sdf(x,y,z)
 endif
 if(tie_wall)then
 fs=rotor_sdf(x,y,z);fr=fs
 endif
 end subroutine
 logical function rotating_owner(x,y,z)result(moving)
 real*8,intent(in)::x,y,z
 moving=moving_wall.and..not.tie_wall
 end function
end module
