program test_actual_velocity_module
 use mod_ibm_surface_velocity
 use mod_rotor_motion,only:geometry_mode,wall_x,moving_wall,tie_wall
 implicit none
 real*8::g(3),f(3),b(3),alpha,q(3),amp,err,val,truth,saved(0:2,0:2,0:2)
 integer::i,j,k,l,changed,failed,reason,ndsolid,nchecks
 logical::ok
 character(32)::arg
 nchecks=0
 do i=1,7
 XMP(i)=.5D0*(i-4);YMP(i)=.5D0*(i-4);ZMP(i)=.5D0*(i-4)
 enddo
 X=XMP;Y=YMP;Z=ZMP
 IFC=4;JFC=4;KFC=4;MPI=0;GFI=0D0;INOUT=0
 NFC_INTP=[0,1,0];NFC_INNER=0
 call init_surface_velocity
 call get_command_argument(1,arg)
 if(trim(arg)=='failstop')then
 geometry_mode=3
 call rebuild_fluid_wall_stencils
 error stop 'Unexpected return from missing-support fail-stop'
 endif

 ! Native staggering must use different component coordinates.
 X(4)=.11D0;Y(4)=.22D0;Z(4)=.33D0
 call native_position(4,4,4,1,q)
 call require(maxval(abs(q-[.11D0,0D0,0D0]))<1D-12,'MAC u position')
 call native_position(4,4,4,2,q)
 call require(maxval(abs(q-[0D0,.22D0,0D0]))<1D-12,'MAC v position')
 call native_position(4,4,4,3,q)
 call require(maxval(abs(q-[0D0,0D0,.33D0]))<1D-12,'MAC w position')
 X=XMP;Y=YMP;Z=ZMP

 g=[0D0,0D0,0D0];f=[1D0,0D0,0D0]
 call fluid_ray(g,f,b,alpha,ok)
 call require(ok.and.abs(b(1)-.25D0)<1D-8,'first plane intersection')
 call require(abs(alpha-1D0/3D0)<1D-8,'bounded distance ratio')
 val=(1D0+alpha)*affine(b)-alpha*affine(f);truth=affine(g)
 call require(abs(val-truth)<1D-12,'affine reproduction')
 call require(abs((1D0+alpha)*2.5D0-alpha*2.5D0-2.5D0)<1D-12,'constant reproduction')
 wall_x=.8D0
 call fluid_ray(g,f,b,alpha,ok)
 call require(.not.ok,'reject excessive extrapolation')
 wall_x=.25D0;geometry_mode=2
 call fluid_ray(g,f,b,alpha,ok)
 call require(.not.ok,'reject second solid on ray')
 ok=wall_path_clear([.25D0,0D0,0D0],f)
 call require(.not.ok,'reject wall-to-fluid path through second solid')
 geometry_mode=3
 call fluid_ray(g,f,b,alpha,ok)
 call require(.not.ok,'reject no fluid intersection')
 geometry_mode=1

 ! Actual set_surface_velocity: body identity comes from B, not G.
 moving_wall=.true.;tie_wall=.false.
 call set_surface_velocity(1,2,0D0,0D0,.5D0,.25D0,0D0,1D0,1D0,0D0,2D0)
 call require(surface_owner(1,2)==1.and.abs(surface_velocity(1,2)+2D0)<1D-12,'moving wall speed at B')
 call require(abs(extension_velocity(1,2)+1D0)<1D-12,'distinct rigid extension at G')
 moving_wall=.false.
 call set_surface_velocity(1,2,0D0,0D0,.5D0,.25D0,0D0,1D0,1D0,0D0,2D0)
 call require(surface_owner(1,2)==0.and.surface_velocity(1,2)==0D0,'stationary wall owner')
 moving_wall=.true.;tie_wall=.true.
 call set_surface_velocity(1,2,0D0,0D0,.5D0,.25D0,0D0,1D0,1D0,0D0,2D0)
 call require(surface_owner(1,2)==0,'coincident owner resolves to fixed in reference')
 tie_wall=.false.

 ! Repair a real unsafe stencil. Only i=6 is an eligible native fluid donor.
 ! This forces the two-index search rather than relying on a nearer ghost.
 INOUT=0;INOUT(6,4,4,2)=1
 MPI=0;GFI=0D0;GFI(1,2,0,0,0)=4D0
 MPI(1,2,1)=1;GFI(1,2,1,0,0)=-3D0
 surface_xyz(:,1,2)=[.25D0,0D0,0D0]
 call repair_one(1,2,changed,failed,reason,ndsolid,amp,err)
 call require(changed==1.and.failed==0.and.ndsolid==1,'solid donor replaced')
 call require(abs(amp-1D0/3D0)<1D-8.and.err<1D-12,'repair alpha and affine check')
 call require(abs(GFI(1,2,2,0,0)+amp)<1D-12.and.abs(GFI(1,2,0,0,0)-1D0-amp)<1D-12,'native GFI offset two')
 call require(count(abs(GFI(1,2,:,:,:))>1D-14)==2,'wall plus one fluid donor only')
 saved=GFI(1,2,:,:,:)
 call repair_one(1,2,changed,failed,reason,ndsolid,amp,err)
 call require(changed==0.and.failed==0.and.all(GFI(1,2,:,:,:)==saved),'valid stencil retained')

 ! Unsupported repair reports failure; the separate mode checks error stop.
 INOUT=0
 call repair_one(1,2,changed,failed,reason,ndsolid,amp,err)
 call require(failed==1.and.changed==0,'unsupported donor search fails')
 XMP(4)=100D0
 call repair_one(1,2,changed,failed,reason,ndsolid,amp,err)
 call require(failed==0.and.changed==0,'outside ROI untouched')
 write(*,*)'ALL_KERNEL_CHECKS_PASS',nchecks
contains
 real*8 function affine(p)result(a)
 real*8,intent(in)::p(3)
 a=.37D0+.11D0*p(1)-.23D0*p(2)+.41D0*p(3)
 end function
 subroutine require(condition,label)
 logical,intent(in)::condition
 character(*),intent(in)::label
 if(.not.condition)then
 write(*,*)'CHECK_FAILED: ',label
 error stop 1
 endif
 nchecks=nchecks+1
 write(*,*)'PASS: ',label
 end subroutine
end program
