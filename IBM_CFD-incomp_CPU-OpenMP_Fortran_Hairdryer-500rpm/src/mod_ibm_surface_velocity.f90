! Candidate mixed-body kinematics. Geometry and original solver are preserved.
module mod_ibm_surface_velocity
 use mod_common
 use mod_flowarray,only:IFC,JFC,KFC,MPI,GFI,INOUT
 use mod_rotor_motion,only:rotor_sdf,component_functions,rotating_owner,omega_now,length_mm,reference_speed
 implicit none
 real*8,allocatable::surface_velocity(:,:),extension_velocity(:,:),surface_xyz(:,:,:)
 integer*8,allocatable::surface_owner(:,:)
contains
 subroutine init_surface_velocity
 integer::n,l,i,j,k
 real*8::x1,y1,z1
 logical::moving
 if(.not.allocated(surface_velocity))then
 allocate(surface_velocity(size(IFC,1),3),extension_velocity(size(IFC,1),3), &
 surface_xyz(3,size(IFC,1),3),surface_owner(size(IFC,1),3))
 endif
 surface_velocity=0D0;extension_velocity=0D0;surface_owner=0;surface_xyz=0D0
 do l=1,3
 do n=1,NFC_INTP(l)+NFC_INNER(l)
 i=IFC(n,l);j=JFC(n,l);k=KFC(n,l)
 x1=XMP(i);y1=YMP(j);z1=ZMP(k)
 if(l==1)x1=X(i)
 if(l==2)y1=Y(j)
 if(l==3)z1=Z(k)
 moving=rotating_owner(x1,y1,z1)
 surface_owner(n,l)=merge(1,0,moving);surface_xyz(:,n,l)=[x1,y1,z1]
 extension_velocity(n,l)=rigid_component(l,y1,z1,moving)
 surface_velocity(n,l)=extension_velocity(n,l)
 enddo
 enddo
 end subroutine

 real*8 function rigid_component(l,y,z,moving)result(v)
 integer,intent(in)::l
 real*8,intent(in)::y,z
 logical,intent(in)::moving
 v=0D0
 if(.not.moving)return
 if(l==2)v=-omega_now*z
 if(l==3)v=omega_now*y
 end function

 subroutine set_surface_velocity(n,l,x1,y1,z1,x0,y0,z0,x2,y2,z2)
 integer,intent(in)::n,l
 real*8,intent(in)::x1,y1,z1,x0,y0,z0,x2,y2,z2
 real*8::fs,fr,fh,tol
 logical::moving
 call component_functions(x0,y0,z0,fs,fr,fh)
 ! Three 20-way root subdivisions: use geometric search uncertainty only.
 tol=max(1D-10,sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)/8000D0)
 moving=min(fr,fh)<fs
 if(fs<=tol.and.min(fr,fh)<=tol)moving=.false.
 surface_owner(n,l)=merge(1,0,moving);surface_xyz(:,n,l)=[x0,y0,z0]
 surface_velocity(n,l)=rigid_component(l,y0,z0,moving)
 ! This is a fictitious rigid extension at a forcing node, not physical
 ! velocity in the overlap. Use the same body identity in the mass source.
 extension_velocity(n,l)=rigid_component(l,y1,z1,moving)
 end subroutine

 subroutine audit_surface_velocity
 integer::n,l,i,j,k,aa,bb,cc,ii,jj,kk,un,nchanged,ntotal
 real*8::x0,y0,z0,x1,y1,z1,fs,fr,fh,expect,val,err,emax
 nchanged=0;ntotal=0;emax=0D0
 open(newunit=un,file='candidate_surface_audit.dat',status='replace')
 write(un,'(A)')'# l i j k owner_node owner_surface x0mm y0mm z0mm fsmm frmm surface_v_m_s extension_v_m_s affine_error_m_s'
 do l=1,3
 do n=1,NFC_INTP(l)
 i=IFC(n,l);j=JFC(n,l);k=KFC(n,l)
 x1=XMP(i);y1=YMP(j);z1=ZMP(k)
 if(l==1)x1=X(i)
 if(l==2)y1=Y(j)
 if(l==3)z1=Z(k)
 x0=surface_xyz(1,n,l);y0=surface_xyz(2,n,l);z0=surface_xyz(3,n,l)
 call component_functions(x0,y0,z0,fs,fr,fh)
 ! Algebraic reproduction of one globally rotating linear velocity field.
 ! No claim that all donors actually belong to that body in the real flow.
 expect=rigid_component(l,y1,z1,.true.)
 val=GFI(n,l,0,0,0)*rigid_component(l,y0,z0,.true.)
 do aa=0,2
 do bb=0,2
 do cc=0,2
 if(aa+bb+cc==0)cycle
 if(abs(GFI(n,l,aa,bb,cc))<1D-14)cycle
 ii=i+aa*MPI(n,l,1);jj=j+bb*MPI(n,l,2);kk=k+cc*MPI(n,l,3)
 ! Surface points in these cases are interior, away from domain wrap.
 if(ii<1.or.ii>N1M.or.jj<1.or.jj>N2M.or.kk<1.or.kk>N3M)cycle
 if(l==2)then
 val=val+GFI(n,l,aa,bb,cc)*rigid_component(l,Y(jj),ZMP(kk),.true.)
 else if(l==3)then
 val=val+GFI(n,l,aa,bb,cc)*rigid_component(l,YMP(jj),Z(kk),.true.)
 endif
 enddo
 enddo
 enddo
 err=abs(val-expect)*reference_speed;emax=max(emax,err);ntotal=ntotal+1
 if(merge(1,0,rotating_owner(x1,y1,z1))/=surface_owner(n,l))then
 nchanged=nchanged+1
 write(un,'(6I7,8ES24.15)')l,i,j,k,merge(1,0,rotating_owner(x1,y1,z1)),surface_owner(n,l), &
 x0*length_mm-62.86D0,y0*length_mm,z0*length_mm,fs*length_mm,fr*length_mm, &
 surface_velocity(n,l)*reference_speed,extension_velocity(n,l)*reference_speed,err
 endif
 enddo
 enddo
 close(un)
 open(newunit=un,file='candidate_audit_summary.dat',status='replace')
 write(un,*)ntotal,nchanged,emax
 close(un)
 end subroutine
! Repair unsafe surface stencils using only actual wall data and native fluid.
! Constant and affine fields are exactly reproduced along each accepted ray.
! No solid donor, no velocity clipping, no silent missing-support fallback.
 subroutine native_position(i,j,k,l,q)
 integer,intent(in)::i,j,k,l
 real*8,intent(out)::q(3)
 q=[XMP(i),YMP(j),ZMP(k)]
 if(l==1)q(1)=X(i)
 if(l==2)q(2)=Y(j)
 if(l==3)q(3)=Z(k)
 end subroutine

 logical function wall_path_clear(a,b)result(ok)
 real*8,intent(in)::a(3),b(3)
 real*8::q(3)
 integer::s,ns
 ns=max(2,ceiling(sqrt(sum((b-a)**2))*length_mm/.03125D0));ok=.false.
 do s=1,ns
 q=a+(b-a)*real(s,8)/ns
 if(rotor_sdf(q(1),q(2),q(3))<=1D-10)return
 enddo
 ok=.true.
 end function

 subroutine fluid_ray(g,f,wall,alpha,ok)
 real*8,intent(in)::g(3),f(3)
 real*8,intent(out)::wall(3),alpha
 logical,intent(out)::ok
 real*8::q(3),lo,hi,t,phi
 integer::s,ns,it
 logical::crossed
 ok=.false.;alpha=huge(1D0);wall=g
 ns=max(2,ceiling(sqrt(sum((f-g)**2))*length_mm/.03125D0))
 lo=0D0;hi=1D0;crossed=.false.
 do s=1,ns
 t=real(s,8)/ns;q=g+t*(f-g);phi=rotor_sdf(q(1),q(2),q(3))
 if(phi>1D-10)then
 if(.not.crossed)then
 hi=t;lo=real(s-1,8)/ns;crossed=.true.
 endif
 else
 if(crossed)return ! do not cross an intervening second solid
 endif
 enddo
 if(.not.crossed)return
 do it=1,30
 t=.5D0*(lo+hi);q=g+t*(f-g)
 if(rotor_sdf(q(1),q(2),q(3))>1D-10)then
 hi=t
 else
 lo=t
 endif
 enddo
 t=.5D0*(lo+hi);wall=g+t*(f-g);alpha=t/max(1D-30,1D0-t)
 ! The fluid sample must be at least as far from the wall as the ghost.
 ! This bounds each new fluid coefficient by1 and wall coefficient by2.
 ok=alpha<=1D0+1D-8
 end subroutine

 subroutine repair_one(n,l,changed,failed,reason,ndsolid,amp,err)
 integer,intent(in)::n,l
 integer,intent(out)::changed,failed,reason,ndsolid
 real*8,intent(out)::amp,err
 integer::i,j,k,a,b,c,ii,jj,kk,du(3),step,kind,d1,d2,d3,ad(3),trial(3),nref
 real*8::g(3),f(3),wall(3),alpha,coef,val,truth
 logical::bad,ok,found
 changed=0;failed=0;reason=0;ndsolid=0;err=0D0;amp=0D0
 i=IFC(n,l);j=JFC(n,l);k=KFC(n,l)
 call native_position(i,j,k,l,g)
 ! Bound this correction to the rotor/casing region already under diagnosis.
 if(g(1)*length_mm-62.86D0 < -66.75D0.or.g(1)*length_mm-62.86D0 > -59D0)return
 if(sqrt(g(2)**2+g(3)**2)*length_mm > 11.1D0)return
 wall=surface_xyz(:,n,l);bad=.false.;nref=0
 if(sqrt(sum((wall-g)**2))*length_mm<1D-8.and.rotor_sdf(g(1),g(2),g(3))<-1D-9)then
 bad=.true.;reason=1 ! old no-crossing branch used a solid node as wall
 endif
 do a=0,2
 do b=0,2
 do c=0,2
 if(a+b+c==0)cycle
 coef=GFI(n,l,a,b,c)
 if(abs(coef)<1D-14)cycle
 nref=nref+1
 ii=i+a*MPI(n,l,1);jj=j+b*MPI(n,l,2);kk=k+c*MPI(n,l,3)
 if(ii<1.or.ii>N1M.or.jj<1.or.jj>N2M.or.kk<1.or.kk>N3M)then
 bad=.true.;reason=2;cycle
 endif
 if(INOUT(ii,jj,kk,l)/=1)then
 bad=.true.;ndsolid=ndsolid+1;reason=3;cycle
 endif
 call native_position(ii,jj,kk,l,f)
 if(.not.wall_path_clear(wall,f))then
 bad=.true.;reason=4
 endif
 enddo
 enddo
 enddo
 if(.not.bad)return
 found=.false.
 search: do kind=1,3 ! axis rays before diagonal rays
 do step=1,2
 do d3=-1,1
 do d2=-1,1
 do d1=-1,1
 if(abs(d1)+abs(d2)+abs(d3)/=kind)cycle
 du=[d1,d2,d3];trial=[i,j,k]+step*du
 if(any(trial<1).or.trial(1)>N1M.or.trial(2)>N2M.or.trial(3)>N3M)cycle
 if(INOUT(trial(1),trial(2),trial(3),l)/=1)cycle
 call native_position(trial(1),trial(2),trial(3),l,f)
 call fluid_ray(g,f,wall,alpha,ok)
 if(.not.ok)cycle
 found=.true.;exit search
 enddo
 enddo
 enddo
 enddo
 enddo search
 if(.not.found)then
 failed=1;return
 endif
 ! MPI remains an original +/-1 direction and GFI indexes stay within0:2.
 MPI(n,l,:)=du;ad=abs(du)*step
 GFI(n,l,:,:,:)=0D0;GFI(n,l,0,0,0)=1D0+alpha
 GFI(n,l,ad(1),ad(2),ad(3))=-alpha
 call set_surface_velocity(n,l,g(1),g(2),g(3),wall(1),wall(2),wall(3),f(1),f(2),f(3))
 changed=1;amp=alpha
 val=(1D0+alpha)*(.37D0+.11D0*wall(1)-.23D0*wall(2)+.41D0*wall(3)) &
 -alpha*(.37D0+.11D0*f(1)-.23D0*f(2)+.41D0*f(3))
 truth=.37D0+.11D0*g(1)-.23D0*g(2)+.41D0*g(3)
 err=abs(val-truth)
 end subroutine

 subroutine rebuild_fluid_wall_stencils
 integer::l,n,changed,failed,reason,ndsolid,nchanged,nfailed,nsolid,un
 real*8::amp,err,maxamp,maxerr,q(3)
 nchanged=0;nfailed=0;nsolid=0;maxamp=0D0;maxerr=0D0
 do l=1,3
 !$omp parallel do private(n,changed,failed,reason,ndsolid,amp,err,q,un) &
 !$omp reduction(+:nchanged,nfailed,nsolid) reduction(max:maxamp,maxerr)
 do n=1,NFC_INTP(l)
 call repair_one(n,l,changed,failed,reason,ndsolid,amp,err)
 nchanged=nchanged+changed;nfailed=nfailed+failed;nsolid=nsolid+ndsolid
 maxamp=max(maxamp,amp);maxerr=max(maxerr,err)
 if(failed/=0)then
 call native_position(IFC(n,l),JFC(n,l),KFC(n,l),l,q)
 !$omp critical(fluid_wall_failure_output)
 open(newunit=un,file='fluid_wall_unsupported.dat',status='unknown',position='append')
 write(un,'(7I10,3ES24.15)')M,MSUB,l,n,IFC(n,l),JFC(n,l),KFC(n,l),q*length_mm
 close(un)
 !$omp end critical(fluid_wall_failure_output)
 endif
 enddo
 enddo
 open(newunit=un,file='fluid_wall_stencils.dat',status='unknown',position='append')
 write(un,'(5I10,2ES24.15)')M,MSUB,nchanged,nfailed,nsolid,maxamp,maxerr
 close(un)
 if(nfailed>0)error stop 'Wall/fluid-only stencil lacks bounded visible fluid support; inspect unsupported list'
 if(maxerr>1D-10)error stop 'Wall/fluid-only affine reproduction failed'
 end subroutine

end module
