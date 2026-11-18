! Dedicated IBM functions, mm file coordinates / nondimensional solver interface.
module mod_stl_sdf
 use iso_fortran_env, only: int32,real32,real64
 implicit none
 type scalar_field
  integer(int32):: n(3)
  real(real64):: o(3),h(3),up(3)
  real(real32),allocatable:: a(:,:,:)
 end type
 type(scalar_field),save:: fixed_field,rotor_field
 logical,save::sdf_loaded=.false.
contains
 subroutine read_one(f,name)
  type(scalar_field),intent(inout)::f
  character(*),intent(in)::name
  character(8)::magic
  integer::u,ios
  open(newunit=u,file=trim(name),access='stream',form='unformatted',status='old',convert='little_endian',iostat=ios)
  if(ios/=0)error stop 'Cannot open IBM function field'
  read(u)magic,f%n,f%o,f%h
  if(magic/='LICAFN02'.or.any(f%n<2).or.any(f%h<=0))error stop 'Invalid IBM function header'
  allocate(f%a(f%n(1),f%n(2),f%n(3)))
  read(u,iostat=ios)f%a
  close(u)
  if(ios/=0)error stop 'Truncated IBM function field'
  f%up=f%o+f%h*real(f%n-1,real64)
  write(*,*) 'IBM FIELD ',trim(name),f%n,f%h
 end subroutine
 subroutine load_sdf(filename)
  character(*),intent(in)::filename
  character(1024)::a,b
  integer::u
  open(newunit=u,file='geometry.in',status='old');read(u,'(A)')a;read(u,'(A)')b;close(u)
  call read_one(fixed_field,trim(a));call read_one(rotor_field,trim(b))
  sdf_loaded=.true.
 end subroutine
 real(real64) function interpolate(f,x,y,z,extend_sign) result(v)
  type(scalar_field),intent(in)::f
  real(real64),intent(in)::x,y,z,extend_sign
  real(real64)::p(3),pc(3),q(3),w(3),wx,wy,wz
  integer::i(3),a,b,c
  p=[x,y,z];pc=min(max(p,f%o),f%up);q=(pc-f%o)/f%h
  i=min(max(int(floor(q)),0),int(f%n)-2);w=q-real(i,real64)
  v=0
  do c=0,1
   wz=merge(w(3),1D0-w(3),c==1)
   do b=0,1
    wy=merge(w(2),1D0-w(2),b==1)
    do a=0,1
     wx=merge(w(1),1D0-w(1),a==1)
     v=v+wx*wy*wz*real(f%a(i(1)+a+1,i(2)+b+1,i(3)+c+1),real64)
    enddo
   enddo
  enddo
  v=v+extend_sign*sqrt(sum((p-pc)**2))
 end function
 real(real64) function evaluate_sdf(x,y,z) result(v)
  real(real64),intent(in)::x,y,z
  v=interpolate(rotor_field,x,y,z,1D0)
 end function
 real(real64) function evaluate_fixed(x,y,z) result(v)
  real(real64),intent(in)::x,y,z
  v=interpolate(fixed_field,x,y,z,-1D0)
 end function
end module
