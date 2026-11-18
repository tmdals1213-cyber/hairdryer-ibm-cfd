! Candidate correction for solid donors in LICA's ghost interpolation.
! Geometry, GFI coefficients, fluid predictors and pressure solver are unchanged.
! Each solid donor is its simultaneous IBM target, NOT an unconstrained predictor.
module mod_ibm_coupled_targets
 use mod_common
 use mod_flowarray, only: U,V,W,IFC,JFC,KFC,MPI,GFI,INOUT
 use mod_ibm_surface_velocity, only: surface_velocity,extension_velocity,surface_xyz
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 public :: compute_coupled_targets,coupled_target,audit_coupled_targets
 real*8,allocatable,save :: coupled_target(:,:)
 integer,allocatable,save :: node_id(:,:,:)
contains
 subroutine donor_indices(n,l,ix,jy,kz)
 integer,intent(in)::n,l
 integer,intent(out)::ix(0:2),jy(0:2),kz(0:2)
 integer::a
 do a=0,2
 ix(a)=IFC(n,l)+a*MPI(n,l,1)
 jy(a)=JFC(n,l)+a*MPI(n,l,2)
 kz(a)=KFC(n,l)+a*MPI(n,l,3)
 enddo
 if(ix(2)>=N1.or.ix(2)<=0)call reindex_i(ix(1),ix(2))
 if(kz(2)>=N3.or.kz(2)<=0)call reindex_k(kz(1),kz(2))
 end subroutine

 real*8 function predictor(i,j,k,l,tmp)result(vv)
 integer,intent(in)::i,j,k,l
 real*8,intent(in)::tmp(M1M,M2M,M3M,3)
 select case(l)
 case(1);vv=U(i,j,k)+tmp(i,j,k,l)
 case(2);vv=V(i,j,k)+tmp(i,j,k,l)
 case(3);vv=W(i,j,k)+tmp(i,j,k,l)
 end select
 end function

 subroutine assemble(l,tmp,ne,col,coef,b,legacy,inner_count)
 integer,intent(in)::l
 real*8,intent(in)::tmp(M1M,M2M,M3M,3)
 integer,intent(out)::ne(:),col(:,:),inner_count
 real*8,intent(out)::coef(:,:),b(:),legacy(:)
 integer::n,a,c,d,i,j,k,ix(0:2),jy(0:2),kz(0:2),don,e
 real*8::g,wall,ext
 if(.not.allocated(node_id))allocate(node_id(N1M,N2M,N3M))
 node_id=0
 do n=1,NFC_INTP(l)+NFC_INNER(l)
 node_id(IFC(n,l),JFC(n,l),KFC(n,l))=n
 enddo
 ne=0;col=1;coef=0D0;b=0D0;legacy=0D0;inner_count=0
 do n=1,NFC_INTP(l)
 wall=0D0
 if(IMOVINGON==1)wall=surface_velocity(n,l)
 b(n)=GFI(n,l,0,0,0)*wall;legacy(n)=b(n)
 call donor_indices(n,l,ix,jy,kz)
 do a=0,2
 do c=0,2
 do d=0,2
 if(a+c+d==0)cycle
 g=GFI(n,l,a,c,d)
 if(g==0D0)cycle
 i=ix(a);j=jy(c);k=kz(d)
 if(i<1.or.i>N1M.or.j<1.or.j>N2M.or.k<1.or.k>N3M)then
 write(*,*)'IBM_COUPLING_INVALID_DONOR',l,n,i,j,k
 error stop 'IBM coupling donor outside predictor domain'
 endif
 legacy(n)=legacy(n)+g*predictor(i,j,k,l,tmp)
 if(INOUT(i,j,k,l)==1)then
 b(n)=b(n)+g*predictor(i,j,k,l,tmp)
 else
 don=node_id(i,j,k)
 if(don==0)then
 write(*,*)'IBM_COUPLING_UNASSIGNED_SOLID',l,n,i,j,k
 error stop 'Solid donor has no forcing constraint'
 endif
 if(don<=NFC_INTP(l))then
 e=ne(n)+1;ne(n)=e;col(e,n)=don;coef(e,n)=-g
 else
 ext=0D0
 if(IMOVINGON==1)ext=extension_velocity(don,l)
 b(n)=b(n)+g*ext;inner_count=inner_count+1
 endif
 endif
 enddo
 enddo
 enddo
 enddo
 end subroutine

 subroutine apply_matrix(ne,col,coef,x,y)
 integer,intent(in)::ne(:),col(:,:)
 real*8,intent(in)::coef(:,:),x(:)
 real*8,intent(out)::y(:)
 integer::n,e
 y=x
 do n=1,size(ne)
 do e=1,ne(n)
 y(n)=y(n)+coef(e,n)*x(col(e,n))
 enddo
 enddo
 end subroutine

 ! Restarted GMRES with two-pass modified Gram-Schmidt and true residual check.
 ! No assumption of diagonal dominance or bounded solid target velocity.
 subroutine solve_constraints(ne,col,coef,b,x,it,res,tol,ok)
 integer,intent(in)::ne(:),col(:,:)
 real*8,intent(in)::coef(:,:),b(:)
 real*8,intent(out)::x(:),res,tol
 integer,intent(out)::it
 logical,intent(out)::ok
 integer,parameter::nr=40,maxcycles=10
 real*8,allocatable::basis(:,:),r(:),work(:)
 real*8::h(nr+1,nr),cs(nr),sn(nr),s(nr+1),yy(nr),beta,t,q,hn
 integer::n,cycle,j,k,pass,last
 n=size(b);x=b;it=0;ok=.false.
 tol=1D-12+1D-10*maxval(abs(b))
 allocate(r(n),work(n),basis(n,nr+1))
 do cycle=1,maxcycles
 call apply_matrix(ne,col,coef,x,work);r=b-work
 res=maxval(abs(r))
 if(.not.ieee_is_finite(res).or..not.all(ieee_is_finite(x)))exit
 if(res<=tol)then
 ok=.true.;exit
 endif
 beta=sqrt(dot_product(r,r));basis(:,1)=r/beta
 h=0D0;s=0D0;s(1)=beta;last=nr
 do j=1,nr
 call apply_matrix(ne,col,coef,basis(:,j),work)
 do pass=1,2
 do k=1,j
 q=dot_product(basis(:,k),work);h(k,j)=h(k,j)+q
 work=work-q*basis(:,k)
 enddo
 enddo
 hn=sqrt(dot_product(work,work));h(j+1,j)=hn
 if(hn>1D-30)basis(:,j+1)=work/hn
 do k=1,j-1
 t=cs(k)*h(k,j)+sn(k)*h(k+1,j)
 h(k+1,j)=-sn(k)*h(k,j)+cs(k)*h(k+1,j);h(k,j)=t
 enddo
 t=sqrt(h(j,j)**2+h(j+1,j)**2)
 if(t<=1D-30.or..not.ieee_is_finite(t))then
 last=j-1;exit
 endif
 cs(j)=h(j,j)/t;sn(j)=h(j+1,j)/t
 h(j,j)=t;h(j+1,j)=0D0
 s(j+1)=-sn(j)*s(j);s(j)=cs(j)*s(j)
 it=it+1
 if(abs(s(j+1))<=tol.or.hn<=1D-30)then
 last=j;exit
 endif
 enddo
 if(last==0)exit
 yy=0D0
 do j=last,1,-1
 t=s(j)
 do k=j+1,last
 t=t-h(j,k)*yy(k)
 enddo
 yy(j)=t/h(j,j)
 enddo
 do j=1,last
 x=x+yy(j)*basis(:,j)
 enddo
 enddo
 call apply_matrix(ne,col,coef,x,work);res=maxval(abs(b-work))
 ok=all(ieee_is_finite(x)).and.ieee_is_finite(res).and.res<=tol
 end subroutine

 subroutine compute_coupled_targets(tmp)
 real*8,intent(in)::tmp(M1M,M2M,M3M,3)
 integer,allocatable::ne(:),col(:,:)
 real*8,allocatable::coef(:,:),b(:),old(:),val(:)
 integer::l,n,it,inner_count,un
 real*8::res,tol
 logical::ok
 if(.not.allocated(coupled_target))allocate(coupled_target(size(IFC,1),3))
 open(newunit=un,file='ibm_coupling.dat',status='unknown',position='append')
 do l=1,3
 n=NFC_INTP(l)
 if(n==0)cycle
 allocate(ne(n),col(26,n),coef(26,n),b(n),old(n),val(n))
 call assemble(l,tmp,ne,col,coef,b,old,inner_count)
 call solve_constraints(ne,col,coef,b,val,it,res,tol,ok)
 write(un,'(8I10,5ES24.15,L3)')M,MSUB,l,n,count(ne>0),sum(ne),inner_count,it, &
 res,tol,maxval(abs(old)),maxval(abs(val)),maxval(abs(val-old)),ok
 flush(un)
 if(.not.ok)then
 write(*,*)'IBM_COUPLING_FAILED',M,MSUB,l,it,res,tol
 error stop 'Coupled IBM constraints did not converge; no fallback or clipping'
 endif
 coupled_target(1:n,l)=val
 deallocate(ne,col,coef,b,old,val)
 enddo
 close(un)
 end subroutine

 ! Preflight only: independent affine fields and export of the actual matrix.
 ! Call after FINDFORCING; this does not advance time or modify production data.
 subroutine audit_coupled_targets(tmp)
 real*8,intent(in)::tmp(M1M,M2M,M3M,3)
 integer,allocatable::ne(:),col(:,:)
 real*8,allocatable::coef(:,:),b(:),old(:),val(:),exact(:),btest(:)
 integer::l,n,nn,it,inner_count,un,ue,test,a,c,d,i,j,k,ix(0:2),jy(0:2),kz(0:2)
 real*8::res,tol,xx,yy,zz,vv,err,algebraic
 logical::ok
 character(40)::fname
 open(newunit=ue,file='coupling_manufactured.dat',status='replace')
 do l=1,3
 nn=NFC_INTP(l)
 allocate(ne(nn),col(26,nn),coef(26,nn),b(nn),old(nn),val(nn),exact(nn),btest(nn))
 call assemble(l,tmp,ne,col,coef,b,old,inner_count)
 write(fname,'(A,I1,A)')'coupling_matrix_',l,'.bin'
 open(newunit=un,file=trim(fname),status='replace',access='stream',form='unformatted')
 write(un)nn,ne,col,coef,b,old
 write(un)IFC(1:nn,l),JFC(1:nn,l),KFC(1:nn,l),surface_xyz(:,1:nn,l),surface_velocity(1:nn,l)
 close(un)
 do test=1,2
 do n=1,nn
 xx=XMP(IFC(n,l));yy=YMP(JFC(n,l));zz=ZMP(KFC(n,l))
 if(l==1)xx=X(IFC(n,l))
 if(l==2)yy=Y(JFC(n,l))
 if(l==3)zz=Z(KFC(n,l))
 exact(n)=manufactured(xx,yy,zz,test)
 btest(n)=GFI(n,l,0,0,0)*manufactured(surface_xyz(1,n,l),surface_xyz(2,n,l),surface_xyz(3,n,l),test)
 call donor_indices(n,l,ix,jy,kz)
 do a=0,2
 do c=0,2
 do d=0,2
 if(a+c+d==0.or.GFI(n,l,a,c,d)==0D0)cycle
 i=ix(a);j=jy(c);k=kz(d)
 if(INOUT(i,j,k,l)==0)then
 if(node_id(i,j,k)<=nn)cycle
 endif
 xx=XMP(i);yy=YMP(j);zz=ZMP(k)
 if(l==1)xx=X(i)
 if(l==2)yy=Y(j)
 if(l==3)zz=Z(k)
 btest(n)=btest(n)+GFI(n,l,a,c,d)*manufactured(xx,yy,zz,test)
 enddo
 enddo
 enddo
 enddo
 call apply_matrix(ne,col,coef,exact,val);algebraic=maxval(abs(val-btest))
 call solve_constraints(ne,col,coef,btest,val,it,res,tol,ok)
 err=maxval(abs(exact-val))
 write(ue,'(6I10,5ES24.15,L3)')l,test,nn,count(ne>0),sum(ne),it,algebraic,err,res,tol,maxval(abs(val)),ok
 flush(ue)
 if(.not.ok.or.algebraic>1D-10.or.err>1D-7)error stop 'Manufactured IBM coupling test failed'
 enddo
 deallocate(ne,col,coef,b,old,val,exact,btest)
 enddo
 close(ue)
 end subroutine

 real*8 function manufactured(x,y,z,test)result(v)
 real*8,intent(in)::x,y,z
 integer,intent(in)::test
 v=.37D0
 if(test==2)v=v+.11D0*x-.23D0*y+.41D0*z
 end function
end module
