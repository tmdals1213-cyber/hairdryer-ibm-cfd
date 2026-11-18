! Far-boundary pressure opening for compact reservoir commissioning.
! Normal speed is pressure-projected on all six faces; reverse flow is allowed.
! Inflow tangential face velocity0; outflow tangential zero normal derivative.
subroutine ambient_ghosts
 use mod_common
 use mod_poiss_mgbicg,only:OPENX
 use mod_flowarray,only:U,V,W,P
 implicit none
 integer*8::i,j,k
 real*8::pb,un
 call ambient_blocked_normals
 do k=1,N3M
 do j=1,N2M
  un=-U(1,j,k);pb=-.5D0*min(un,0D0)**2;P(0,j,k)=2*pb-P(1,j,k)
  un=U(N1,j,k);pb=-.5D0*min(un,0D0)**2;P(N1,j,k)=2*pb-P(N1M,j,k)
 enddo
 enddo
 do k=1,N3M
 do i=1,N1M
  un=-V(i,1,k);pb=-.5D0*min(un,0D0)**2;P(i,0,k)=2*pb-P(i,1,k)
  un=V(i,N2,k);pb=-.5D0*min(un,0D0)**2;P(i,N2,k)=2*pb-P(i,N2M,k)
 enddo
 enddo
 do j=1,N2M
 do i=1,N1M
  un=-W(i,j,1);pb=-.5D0*min(un,0D0)**2;P(i,j,0)=2*pb-P(i,j,1)
  un=W(i,j,N3);pb=-.5D0*min(un,0D0)**2;P(i,j,N3)=2*pb-P(i,j,N3M)
 enddo
 enddo
 ! Ghost values are antisymmetric for zero face velocity on inflow.
 do k=1,N3M
 do j=1,N2M
  V(0,j,k)=merge(-V(1,j,k),V(1,j,k),U(1,j,k)>0)
  W(0,j,k)=merge(-W(1,j,k),W(1,j,k),U(1,j,k)>0)
  V(N1,j,k)=merge(-V(N1M,j,k),V(N1M,j,k),U(N1,j,k)<0)
  W(N1,j,k)=merge(-W(N1M,j,k),W(N1M,j,k),U(N1,j,k)<0)
 enddo
 enddo
 do k=1,N3M
 do i=1,N1M
  U(i,0,k)=merge(-U(i,1,k),U(i,1,k),V(i,1,k)>0)
  W(i,0,k)=merge(-W(i,1,k),W(i,1,k),V(i,1,k)>0)
  U(i,N2,k)=merge(-U(i,N2M,k),U(i,N2M,k),V(i,N2,k)<0)
  W(i,N2,k)=merge(-W(i,N2M,k),W(i,N2M,k),V(i,N2,k)<0)
 enddo
 enddo
 do j=1,N2M
 do i=1,N1M
  U(i,j,0)=merge(-U(i,j,1),U(i,j,1),W(i,j,1)>0)
  V(i,j,0)=merge(-V(i,j,1),V(i,j,1),W(i,j,1)>0)
  U(i,j,N3)=merge(-U(i,j,N3M),U(i,j,N3M),W(i,j,N3)<0)
  V(i,j,N3)=merge(-V(i,j,N3M),V(i,j,N3M),W(i,j,N3)<0)
 enddo
 enddo
end subroutine
subroutine ambient_predict
 use mod_common
 use mod_poiss_mgbicg,only:OPENX
 use mod_flowarray,only:U,V,W
 implicit none
 U(1,1:N2M,1:N3M)=U(2,1:N2M,1:N3M)
 U(N1,1:N2M,1:N3M)=U(N1M,1:N2M,1:N3M)
 V(1:N1M,1,1:N3M)=V(1:N1M,2,1:N3M)
 V(1:N1M,N2,1:N3M)=V(1:N1M,N2M,1:N3M)
 W(1:N1M,1:N2M,1)=W(1:N1M,1:N2M,2)
 W(1:N1M,1:N2M,N3)=W(1:N1M,1:N2M,N3M)
 call ambient_blocked_normals
end subroutine
subroutine ambient_project(phi)
 use mod_common
 use mod_poiss_mgbicg,only:OPENX
 use mod_flowarray,only:U,V,W
 implicit none
 real*8::phi(0:M1,0:M2,0:M3)
 U(1,1:N2M,1:N3M)=U(1,1:N2M,1:N3M)-DTCONST*2D0*phi(1,1:N2M,1:N3M)/SDX(1)
 U(N1,1:N2M,1:N3M)=U(N1,1:N2M,1:N3M)+DTCONST*2D0*phi(N1M,1:N2M,1:N3M)/SDX(N1M)
 block
 integer*8::i
 do i=1,N1M
  V(i,1,1:N3M)=V(i,1,1:N3M)-OPENX(i)*DTCONST*2D0*phi(i,1,1:N3M)/SDY(1)
  V(i,N2,1:N3M)=V(i,N2,1:N3M)+OPENX(i)*DTCONST*2D0*phi(i,N2M,1:N3M)/SDY(N2M)
  W(i,1:N2M,1)=W(i,1:N2M,1)-OPENX(i)*DTCONST*2D0*phi(i,1:N2M,1)/SDZ(1)
  W(i,1:N2M,N3)=W(i,1:N2M,N3)+OPENX(i)*DTCONST*2D0*phi(i,1:N2M,N3M)/SDZ(N3M)
 enddo
 end block
 call ambient_blocked_normals
end subroutine

subroutine ambient_blocked_normals
 use mod_common
 use mod_poiss_mgbicg,only:OPENX
 use mod_flowarray,only:V,W
 implicit none
 integer*8::i
 do i=1,N1M
  if(OPENX(i)<.5D0)then
   V(i,1,1:N3M)=0D0;V(i,N2,1:N3M)=0D0
   W(i,1:N2M,1)=0D0;W(i,1:N2M,N3)=0D0
  endif
 enddo
end subroutine
