! Exploratory homogeneous-Neumann reconstruction of OLD physical pressure
! in fluid momentum predictors ONLY. PHI, its Poisson operator, projection,
! and the original IBM mass source are deliberately a consistent unchanged set.
! Never average pressures across an intervening solid; no global solid-P reset.
module mod_ibm_pressure_neumann
 use mod_common
 use mod_flowarray,only:P,U,V,W,INOUT,QMASS,IFC,JFC,KFC,MPI,GFI
 use mod_rotor_motion,only:rotor_sdf,length_mm,reference_time,reference_speed,omega_now
 use mod_ibm_surface_velocity,only:surface_velocity
 use, intrinsic::ieee_arithmetic
 implicit none
 real*8,allocatable::pressure_gradient(:,:,:,:),center_phi(:,:,:)
 logical,allocatable::previous_fluid(:,:,:),valid_pressure(:,:,:)
 integer,save::last_modified=0,last_failed=0
 real*8,save::last_delta=0D0,last_manufactured_error=0D0
contains
 subroutine classify_pressure
 integer::i,j,k
 !$omp parallel do private(i,j)
 do k=1,N3M
 do j=1,N2M
 do i=1,N1M
 center_phi(i,j,k)=rotor_sdf(XMP(i),YMP(j),ZMP(k))
 enddo
 enddo
 enddo
 end subroutine

 subroutine initialize_pressure_neumann
 if(allocated(previous_fluid))return
 allocate(pressure_gradient(N1M,N2M,N3M,3),center_phi(N1M,N2M,N3M), &
 previous_fluid(N1M,N2M,N3M),valid_pressure(N1M,N2M,N3M))
 call classify_pressure
 previous_fluid=center_phi>1D-10
 end subroutine

 function cross(a,b)result(c)
 real*8,intent(in)::a(3),b(3)
 real*8::c(3)
 c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
 end function

 logical function visible(a,b)result(ok)
 real*8,intent(in)::a(3),b(3)
 real*8::q(3),d
 integer::n,s
 ! Sampling at <=31.25 micrometres, the rotor-function spacing. This is
 ! a sampled visibility check, not an exact CAD intersection certificate.
 d=sqrt(sum((a-b)**2))*length_mm
 n=max(2,ceiling(d/.03125D0));ok=.false.
 do s=1,n-1
 q=a+(b-a)*(real(s,8)/n)
 if(rotor_sdf(q(1),q(2),q(3))<=1D-10)return
 enddo
 ok=.true.
 end function

 subroutine interface_normal(q,i,j,k,normal,ok)
 real*8,intent(in)::q(3)
 integer,intent(in)::i,j,k
 real*8,intent(out)::normal(3)
 logical,intent(out)::ok
 real*8::b(3),a(3),r(3),v(3),vp(3),vm(3),best,d,lo,hi,t,eps,norm
 integer::ii,jj,kk,s,l,ns
 best=huge(1D0);ok=.false.;normal=0D0
 do kk=max(1,k-2),min(N3M,k+2)
 do jj=max(1,j-2),min(N2M,j+2)
 do ii=max(1,i-2),min(N1M,i+2)
 if(center_phi(ii,jj,kk)>1D-10)cycle
 v=[XMP(ii),YMP(jj),ZMP(kk)];d=sum((v-q)**2)
 if(d<best)then
 best=d;b=v;ok=.true.
 endif
 enddo
 enddo
 enddo
 if(.not.ok)return
 ! First zero along the ray, not a root across the whole thin blade.
 ns=max(2,ceiling(sqrt(best)*length_mm/.03125D0));lo=0D0;hi=1D0
 do s=1,ns
 t=real(s,8)/ns;r=q+t*(b-q)
 if(rotor_sdf(r(1),r(2),r(3))<=1D-10)then
 hi=t;exit
 endif
 lo=t
 enddo
 do s=1,28
 t=.5D0*(lo+hi);r=q+t*(b-q)
 if(rotor_sdf(r(1),r(2),r(3))>1D-10)then
 lo=t
 else
 hi=t
 endif
 enddo
 r=q+.5D0*(lo+hi)*(b-q);eps=.001D0/length_mm
 do l=1,3
 vp=r;vm=r;vp(l)=vp(l)+eps;vm(l)=vm(l)-eps
 normal(l)=(rotor_sdf(vp(1),vp(2),vp(3))-rotor_sdf(vm(1),vm(2),vm(3)))/(2D0*eps)
 enddo
 norm=sqrt(sum(normal**2));ok=norm>1D-8
 if(ok)normal=normal/norm
 end subroutine

 subroutine solve3(a,b,x,ok)
 real*8,intent(in)::a(3,3),b(3)
 real*8,intent(out)::x(3)
 logical,intent(out)::ok
 real*8::mat(3,4),tmp(4),pivot,scale
 integer::h,r,s
 mat(:,1:3)=a;mat(:,4)=b;scale=maxval(abs(a));ok=.false.;x=0D0
 if(scale<1D-30)return
 do h=1,3
 r=h-1+maxloc(abs(mat(h:3,h)),dim=1)
 if(abs(mat(r,h))<1D-9*scale)return
 tmp=mat(h,:);mat(h,:)=mat(r,:);mat(r,:)=tmp
 pivot=mat(h,h);mat(h,:)=mat(h,:)/pivot
 do s=1,3
 if(s/=h)mat(s,:)=mat(s,:)-mat(s,h)*mat(h,:)
 enddo
 enddo
 x=mat(:,4);ok=all(ieee_is_finite(x))
 end subroutine

 subroutine reconstruct(q,i,j,k,grad,ok,ndon,err,rankout,pvalue)
 real*8,intent(in)::q(3)
 integer,intent(in)::i,j,k
 real*8,intent(out)::grad(3),err
 logical,intent(out)::ok
 integer,intent(out)::ndon,rankout
 real*8,intent(out),optional::pvalue
 real*8::normal(3),t1(3),t2(3),axis(3),d(3),point(3),v(3),mat(3,3),rhs(3), &
 test_rhs(3),coef(3),test_coef(3),w,h,p0,dp,truth(3),cov(2,2),rr(2), &
 tc(2),evec(2),lam(2),theta,aa,bb,cc,slopes(2),tslopes(2),expected(2),ev(2,2)
 real*8::fit_xy(343,2),fit_w(343),fit_p(343),fit_t(343),fit_intercept
 integer::ii,jj,kk,rad,a,b
 logical::testok
 grad=0D0;ndon=0;err=0D0;rankout=0
 if(present(pvalue))pvalue=0D0
 call interface_normal(q,i,j,k,normal,ok)
 if(.not.ok)return
 axis=0D0;axis(minloc(abs(normal),dim=1))=1D0
 t1=cross(normal,axis);t1=t1/sqrt(sum(t1*t1));t2=cross(normal,t1)
 h=min(SDX(i),SDY(j),SDZ(k))
 do rad=1,3
 mat=0D0;rhs=0D0;test_rhs=0D0;ndon=0;p0=0D0
 do kk=max(1,k-rad),min(N3M,k+rad)
 do jj=max(1,j-rad),min(N2M,j+rad)
 do ii=max(1,i-rad),min(N1M,i+rad)
 if(.not.valid_pressure(ii,jj,kk))cycle
 point=[XMP(ii),YMP(jj),ZMP(kk)]
 if(.not.visible(q,point))cycle
 d=(point-q)/h;v=[1D0,dot_product(d,t1),dot_product(d,t2)]
 w=1D0/max(.0625D0,sum(d*d));ndon=ndon+1
 if(ndon==1)p0=P(ii,jj,kk)
 dp=P(ii,jj,kk)-p0
 fit_xy(ndon,:)=v(2:3);fit_w(ndon)=w;fit_p(ndon)=dp
 fit_t(ndon)=2.3D0+.7D0*v(2)-.4D0*v(3)
 do a=1,3
 rhs(a)=rhs(a)+w*v(a)*dp
 test_rhs(a)=test_rhs(a)+w*v(a)*(2.3D0+.7D0*v(2)-.4D0*v(3))
 do b=1,3
 mat(a,b)=mat(a,b)+w*v(a)*v(b)
 enddo
 enddo
 enddo
 enddo
 enddo
 if(ndon==0)cycle
 ! Same weighted fit/cutoff; centered one-sided Jacobi avoids cancellation.
 call stable_tangent_fit(ndon,fit_xy,fit_w,fit_p,fit_t,slopes,tslopes,expected,fit_intercept,rankout,lam)
 ! Seek full tangent rank within the bounded stencil before accepting loss.
 if(rankout<2.and.rad<3)cycle
 grad=(slopes(1)*t1+slopes(2)*t2)/h
 if(present(pvalue))pvalue=p0+fit_intercept
 err=maxval(abs(tslopes-expected))
 err=max(err,abs(dot_product(normal,grad))*h/max(1D0,maxval(abs(slopes))))
 ok=err<1D-7.and.all(ieee_is_finite(grad))
 return
 enddo
 end subroutine

 subroutine prepare_pressure_neumann
 integer::i,j,k,l,ii,jj,kk,un,changed,failed,ndon,fresh,rankout,rank0,rank1,nfresh
 real*8::q(3),g(3),old,delta,err,maxerr,pvalue,pchange
 logical::ok
 call classify_pressure
 valid_pressure=previous_fluid.and.center_phi>1D-10
 fresh=count((center_phi>1D-10).and..not.previous_fluid)
 changed=0;failed=0;delta=0D0;maxerr=0D0;rank0=0;rank1=0
 nfresh=0;pchange=0D0
 ! A newly fluid pressure center must not retain its old fictitious solid
 ! P after the one-stage exclusion expires. Initialize its physical P
 ! from the same connected, old-and-current-fluid fit. Donor eligibility
 ! remains false for ALL fresh points during this reconstruction stage.
 !$omp parallel do private(i,j,q,g,pvalue,ok,ndon,err,rankout) &
 !$omp reduction(+:nfresh,failed) reduction(max:pchange,maxerr)
 do k=1,N3M
 do j=1,N2M
 do i=1,N1M
 if(center_phi(i,j,k)<=1D-10.or.previous_fluid(i,j,k))cycle
 q=[XMP(i),YMP(j),ZMP(k)]
 if(.not.in_rotor_region(q))cycle
 call reconstruct(q,i,j,k,g,ok,ndon,err,rankout,pvalue)
 if(.not.ok)then
 failed=failed+1
 !$omp critical(neumann_failure)
 open(newunit=un,file='neumann_failure.dat',status='unknown',position='append')
 write(un,*)NTIME,MSUB,0,i,j,k,ndon,q*length_mm,rankout,err
 close(un)
 !$omp end critical(neumann_failure)
 else
 pchange=max(pchange,abs(P(i,j,k)-pvalue)*reference_speed**2)
 P(i,j,k)=pvalue;nfresh=nfresh+1;maxerr=max(maxerr,err)
 endif
 enddo
 enddo
 enddo
 do l=1,3
 !$omp parallel do private(i,j,ii,jj,kk,q,g,old,ok,ndon,err,rankout) &
 !$omp reduction(+:changed,failed,rank0,rank1) reduction(max:delta,maxerr)
 do k=1,N3M
 do j=1,N2M
 do i=1,N1M
 select case(l)
 case(1);old=(P(i,j,k)-P(i-1,j,k))*VVDX(i)
 case(2);old=(P(i,j,k)-P(i,j-1,k))*VVDY(j)
 case(3);old=(P(i,j,k)-P(i,j,k-1))*VVDZ(k)
 end select
 pressure_gradient(i,j,k,l)=old
 if(INOUT(i,j,k,l)/=1)cycle
 ii=i;jj=j;kk=k
 select case(l)
 case(1);if(i==1)cycle
 ii=i-1
 case(2);if(j==1)cycle
 jj=j-1
 case(3);if(k==1)cycle
 kk=k-1
 end select
 if(valid_pressure(i,j,k).and.valid_pressure(ii,jj,kk))cycle
 q=[XMP(i),YMP(j),ZMP(k)]
 if(l==1)q(1)=X(i)
 if(l==2)q(2)=Y(j)
 if(l==3)q(3)=Z(k)
 ! Bounded mechanism trial around the rotor only. Downstream thin fixed
 ! details retain the baseline treatment; this is NOT a global wall fix.
 if(.not.in_rotor_region(q))cycle
 call reconstruct(q,i,j,k,g,ok,ndon,err,rankout)
 if(.not.ok)then
 failed=failed+1
 !$omp critical(neumann_failure)
 open(newunit=un,file='neumann_failure.dat',status='unknown',position='append')
 write(un,*)NTIME,MSUB,l,i,j,k,ndon,q*length_mm,rankout,err
 close(un)
 !$omp end critical(neumann_failure)
 else
 pressure_gradient(i,j,k,l)=g(l);changed=changed+1
 if(rankout==0)rank0=rank0+1
 if(rankout==1)rank1=rank1+1
 delta=max(delta,abs(g(l)-old)*DTCONST*reference_speed);maxerr=max(maxerr,err)
 endif
 enddo
 enddo
 enddo
 enddo
 previous_fluid=center_phi>1D-10
 last_modified=changed;last_failed=failed;last_delta=delta;last_manufactured_error=maxerr
 open(newunit=un,file='neumann_reconstruction.dat',status='unknown',position='append')
 write(un,'(8I10,4ES24.15)')NTIME,MSUB,changed,failed,fresh,8,rank0,rank1,(TIME+SUBDT)*reference_time,delta,maxerr,DTCONST*reference_time
 close(un)
 open(newunit=un,file='fresh_pressure_reconstruction.dat',status='unknown',position='append')
 write(un,'(3I10,2ES24.15)')NTIME,MSUB,nfresh,(TIME+SUBDT)*reference_time,pchange
 close(un)
 if(failed>0)error stop 'Neumann reconstruction failed: inspect donor count, rank and algebra error in neumann_failure.dat'
 end subroutine

 logical function in_rotor_region(q)result(inside)
 real*8,intent(in)::q(3)
 inside=q(1)*length_mm-62.86D0>=-66.75D0.and.q(1)*length_mm-62.86D0<=-59D0.and. &
 sqrt(q(2)**2+q(3)**2)*length_mm<=11.1D0
 end function

 subroutine pressure_stage_audit
 integer::i,j,k,l,n,a,b,c,ii,jj,kk,un,idx(4)
 real*8::d,divmax,divfluid,peak,allpeak,val,rwall,donor,wallmax,area,flux,tip
 divmax=0D0;divfluid=0D0;peak=0D0;allpeak=0D0;idx=0;wallmax=0D0
 do k=1,N3M
 do j=1,N2M
 do i=1,N1M
 d=(U(i+1,j,k)-U(i,j,k))*SSDX(i)+(V(i,j+1,k)-V(i,j,k))*SSDY(j)+ &
 (W(i,j,k+1)-W(i,j,k))*SSDZ(k)-QMASS(i,j,k)
 divmax=max(divmax,abs(d))
 if(center_phi(i,j,k)>1D-10)divfluid=max(divfluid,abs(d))
 do l=1,3
 val=node_velocity(i,j,k,l);allpeak=max(allpeak,abs(val))
 if(INOUT(i,j,k,l)==1.and.abs(val)>peak)then
 peak=abs(val);idx=[i,j,k,l]
 endif
 enddo
 enddo
 enddo
 enddo
 do l=1,3
 do n=1,NFC_INTP(l)
 i=IFC(n,l);j=JFC(n,l);k=KFC(n,l);donor=0D0
 do c=0,2
 do b=0,2
 do a=0,2
 if(a+b+c==0.or.abs(GFI(n,l,a,b,c))<1D-14)cycle
 ii=i+a*MPI(n,l,1);jj=j+b*MPI(n,l,2);kk=k+c*MPI(n,l,3)
 if(ii<1.or.ii>N1M.or.jj<1.or.jj>N2M.or.kk<1.or.kk>N3M)cycle
 donor=donor+GFI(n,l,a,b,c)*node_velocity(ii,jj,kk,l)
 enddo
 enddo
 enddo
 if(abs(GFI(n,l,0,0,0))>1D-12)then
 rwall=(node_velocity(i,j,k,l)-donor)/GFI(n,l,0,0,0)-surface_velocity(n,l)
 wallmax=max(wallmax,abs(rwall))
 endif
 enddo
 enddo
 flux=0D0
 do k=1,N3M
 do j=1,N2M
 flux=flux+(U(N1,j,k)-U(1,j,k))*SDY(j)*SDZ(k)
 enddo
 enddo
 do k=1,N3M
 do i=1,N1M
 flux=flux+(V(i,N2,k)-V(i,1,k))*SDX(i)*SDZ(k)
 enddo
 enddo
 do j=1,N2M
 do i=1,N1M
 flux=flux+(W(i,j,N3)-W(i,j,1))*SDX(i)*SDY(j)
 enddo
 enddo
 tip=abs(omega_now)*10.6752657302D0/length_mm*reference_speed
 open(newunit=un,file='pressure_stage_audit.dat',status='unknown',position='append')
 write(un,'(6I10,9ES24.15)')NTIME,MSUB,idx,(TIME+SUBDT)*reference_time,peak*reference_speed, &
 allpeak*reference_speed,tip,divmax,divfluid,wallmax*reference_speed, &
 flux*reference_speed*(length_mm*.001D0)**2,DTCONST*reference_time
 close(un)
 ! Preserve the existing 5m/s review guard, now checked each RK stage.
 ! This is a stop/report criterion, never a velocity clip or a physical limit.
 if(.not.ieee_is_finite(allpeak).or.allpeak*reference_speed>5D0)then
 open(newunit=un,file='STOP_AFTER_STEP',status='replace')
 write(un,*)'Stage velocity exceeds existing5m/s review guard',NTIME,MSUB,allpeak*reference_speed
 close(un)
 endif
 end subroutine

 real*8 function node_velocity(i,j,k,l)result(vv)
 integer,intent(in)::i,j,k,l
 select case(l)
 case(1);vv=U(i,j,k)
 case(2);vv=V(i,j,k)
 case(3);vv=W(i,j,k)
 end select
 end function
! Same weighted minimum-norm tangential fit and eigenvalue cutoff as before.
! Work on centered weighted donor columns directly to avoid subtracting
! nearly equal second moments and subtracting eigenvalues for the weak mode.
 subroutine stable_tangent_fit(n,xy,weights,values,testvalues,slopes,tslopes,expected,intercept,rankout,eigenvalues)
 integer,intent(in)::n
 real*8,intent(in)::xy(343,2),weights(343),values(343),testvalues(343)
 real*8,intent(out)::slopes(2),tslopes(2),expected(2),intercept,eigenvalues(2)
 integer,intent(out)::rankout
 real*8::meanxy(2),meanp,meant,sw,design(343,2),bp(343),bt(343),right(2,2), &
 aa,bb,cc,theta,ct,st,tmp(343),tmpr(2),cutoff
 integer::r,a,iter
 sw=sum(weights(1:n))
 meanxy(1)=sum(weights(1:n)*xy(1:n,1))/sw
 meanxy(2)=sum(weights(1:n)*xy(1:n,2))/sw
 meanp=sum(weights(1:n)*values(1:n))/sw
 meant=sum(weights(1:n)*testvalues(1:n))/sw
 do r=1,n
 design(r,:)=sqrt(weights(r))*(xy(r,:)-meanxy)
 bp(r)=sqrt(weights(r))*(values(r)-meanp)
 bt(r)=sqrt(weights(r))*(testvalues(r)-meant)
 enddo
 right=0D0;right(1,1)=1D0;right(2,2)=1D0
 ! One-sided two-column Jacobi SVD; squared column norms are eigenvalues.
 do iter=1,6
 aa=sum(design(1:n,1)**2);bb=sum(design(1:n,1)*design(1:n,2));cc=sum(design(1:n,2)**2)
 if(abs(bb)<=8D0*epsilon(1D0)*sqrt(max(0D0,aa*cc)))exit
 theta=.5D0*atan2(2D0*bb,aa-cc);ct=cos(theta);st=sin(theta)
 tmp(1:n)=design(1:n,1)
 design(1:n,1)=ct*tmp(1:n)+st*design(1:n,2)
 design(1:n,2)=-st*tmp(1:n)+ct*design(1:n,2)
 tmpr=right(:,1)
 right(:,1)=ct*tmpr+st*right(:,2)
 right(:,2)=-st*tmpr+ct*right(:,2)
 enddo
 eigenvalues=[sum(design(1:n,1)**2),sum(design(1:n,2)**2)]
 cutoff=max(1D-10,1D-8*maxval(eigenvalues))
 slopes=0D0;tslopes=0D0;expected=0D0;rankout=0
 do a=1,2
 if(eigenvalues(a)<=cutoff)cycle
 rankout=rankout+1
 slopes=slopes+right(:,a)*sum(design(1:n,a)*bp(1:n))/eigenvalues(a)
 tslopes=tslopes+right(:,a)*sum(design(1:n,a)*bt(1:n))/eigenvalues(a)
 expected=expected+right(:,a)*dot_product(right(:,a),[.7D0,-.4D0])
 enddo
 intercept=meanp-dot_product(meanxy,slopes)
 end subroutine

end module
