! The composite zero set is the SAME function used for IBM and exported checks.
module mod_rotor_motion
 use iso_fortran_env,only:real64
 use mod_stl_sdf,only:evaluate_sdf,evaluate_fixed
 implicit none
 real(real64),save::rotor_re,physical_diameter,rpm,kinematic_nu,inflow_ratio,ramp_time
 real(real64),save::reference_speed,reference_time,omega_max,omega_now=0D0,angle_now=0D0
 real(real64),save::ctheta=1D0,stheta=0D0,x_half,length_mm,motion_scale=1D0
 real(real64),parameter::axis_x_mm=-62.86D0,casing_r=10.745D0,seal_band=.0625D0
 logical,save::tip_sealed=.true.
contains
 subroutine load_rotor(filename)
  character(*),intent(in)::filename
  integer::u,seal_flag
  logical::has_scale
  open(newunit=u,file=filename,status='old');read(u,*)physical_diameter,rpm,kinematic_nu,inflow_ratio,ramp_time;read(u,*)x_half;close(u)
  if(physical_diameter<=0.or.abs(rpm)<1D-12.or.kinematic_nu<=0.or.ramp_time<=0)error stop 'Invalid rotor input'
  reference_speed=abs(rpm)*acos(-1D0)*physical_diameter/60D0
  reference_time=physical_diameter/reference_speed;rotor_re=reference_speed*physical_diameter/kinematic_nu
  inquire(file='rotation_scale.in',exist=has_scale)
  motion_scale=1D0
  if(has_scale)then
   open(newunit=u,file='rotation_scale.in',status='old');read(u,*)motion_scale;close(u)
  endif
  if(motion_scale<0D0.or.motion_scale>1D0)error stop 'Invalid motion control scale'
  inquire(file='tip_seal.in',exist=has_scale)
  tip_sealed=.true.
  if(has_scale)then
   open(newunit=u,file='tip_seal.in',status='old');read(u,*)seal_flag;close(u)
   if(seal_flag/=0.and.seal_flag/=1)error stop 'tip_seal must be0 or1'
   tip_sealed=seal_flag==1
  endif
  omega_max=sign(2D0,rpm)*motion_scale;length_mm=physical_diameter*1000D0
  write(*,*)'COMPACT ROTOR D,RPM,NU,UREF,RE=',physical_diameter,rpm,kinematic_nu,reference_speed,rotor_re
 end subroutine
 subroutine update_rotor(t)
  real(real64),intent(in)::t
  real(real64)::s,pi
  pi=acos(-1D0)
  if(t<ramp_time)then
   s=max(0D0,t);omega_now=omega_max*.5D0*(1D0-cos(pi*s/ramp_time))
   angle_now=omega_max*.5D0*(s-ramp_time/pi*sin(pi*s/ramp_time))
  else
   omega_now=omega_max;angle_now=omega_max*(t-.5D0*ramp_time)
  endif
  ctheta=cos(angle_now);stheta=sin(angle_now)
 end subroutine
 real(real64) function smooth(s)result(w)
  real(real64),intent(in)::s
  real(real64)::v
  v=max(0D0,min(1D0,s));w=v*v*(3D0-2D0*v)
 end function
 subroutine component_functions(x,y,z,fs,fr,fh)
  real(real64),intent(in)::x,y,z
  real(real64),intent(out)::fs,fr,fh
  real(real64)::xx,yy,zz,r,qr,qy,qz,blend,dx,dr,slab
  xx=x*length_mm+axis_x_mm;yy=y*length_mm;zz=z*length_mm;r=sqrt(yy*yy+zz*zz)
  slab=max(-163.25D0-xx,xx-10.5D0)
  if(slab>1D0)then
   fs=slab
  else
   if(r>26.5D0)then
    fs=max(25D0-r,xx-7D0)
   else
    fs=evaluate_fixed(xx,yy,zz)
    blend=smooth((r-25D0)/1D0)*smooth((xx-5D0)/1D0)
    fs=(1D0-blend)*fs+blend*max(25D0-r,xx-7D0)
   endif
   blend=smooth((xx+64.89D0)/.24D0)*smooth((-58.982D0-xx)/.268D0)*smooth((r-(casing_r-.4D0))/.15D0)
   fs=(1D0-blend)*fs+blend*(casing_r-r)
   fs=max(fs,slab)
  endif
  dx=max(-66D0-xx,xx+59.7D0)
  if(dx>.1D0.or.r>11.1D0)then
   fr=max(dx,r-11.1D0)+.1D0
  else
   qr=r
   if(tip_sealed)qr=min(r,casing_r-seal_band)
   if(r>1D-12)then
    qy=(ctheta*yy+stheta*zz)*qr/r;qz=(-stheta*yy+ctheta*zz)*qr/r
   else
    qy=0;qz=0
   endif
   fr=evaluate_sdf(xx,qy,qz)
   if(tip_sealed.and.r>casing_r-seal_band)fr=max(fr,r-casing_r)
  endif
  dx=abs(xx-(-66.255D0-38.815D0)/2D0)-(66.255D0-38.815D0)/2D0;dr=r-1.05D0
  fh=sqrt(max(dx,0D0)**2+max(dr,0D0)**2)+min(max(dx,dr),0D0)
  fs=fs/length_mm;fr=fr/length_mm;fh=fh/length_mm
 end subroutine
 real(real64) function rotor_sdf(x,y,z)result(f)
  real(real64),intent(in)::x,y,z
  real(real64)::fs,fr,fh
  call component_functions(x,y,z,fs,fr,fh);f=min(fs,fr,fh)
 end function
 logical function rotating_owner(x,y,z)result(moving)
  real(real64),intent(in)::x,y,z
  real(real64)::fs,fr,fh
  call component_functions(x,y,z,fs,fr,fh)
  moving=min(fr,fh)<fs
  if(fs<=1D-10.and.min(fr,fh)<=1D-10)moving=.false.
 end function
end module
