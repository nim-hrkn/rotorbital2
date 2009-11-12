#if 0

 engler angle theta and phi  determie the new c-axis
      the new c-axis=(/ cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta) /)

      rot is the matrix to find this, i.e. rot (0 0 1)^T =  the new c-axis

      rot1(theta,phi) (1,0,0)^T -> the rotated original a-axis 

      rotate this rotated original a-axis again to the a_proj axis. This ange is psi

      now you have rot2(theta,phi,psi) 
 
      rot2 ( 0 0 1)^T = (/ cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta) /)
      rot2 ( 1 0 0)^T = the a_proj axis 

      note: 
         a_proj is given by a_tmp and a_proj is calculated by
         a_proj= a_tmp - (a_tmp, 'the rotated original a-axis') 'the rotated original a-axis'

#endif
#ifdef DEBUG_WRITE
#undef DEBUG_WRITE
#endif
      module m_eulerangle

      private:: abs3 

      contains

c return |v|
      real(8) function abs3(v)
      real(8),intent(in):: v(3)
      real(8):: r
      integer:: i
      r=0.0d0
      do i=1,3
       r=r+v(i)**2
      enddo 
      abs3= sqrt(r)
      end function abs3

c return axis defined by euler angle theta and phi 
      subroutine eulerangle2axis(theta,phi,
     o  axis)
c (rcos(θ)cos(φ),rcos(θ)sin(φ),rsin(θ))
      real(8):: theta,phi,axis(3)
      axis=(/ cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta) /)
      end subroutine eulerangle2axis

c
c find euler angle theta and phi from axis0 
c
      subroutine axis2eulerangle (axis0,
     o   theta,phi)
      implicit none
      real(8),intent(in):: axis0(3)
      real(8):: axis2(3),axis(3)
      real(8):: ta ,r2
      real(8),intent(out)::phi,theta
      integer:: i
      real(8):: eps =1.0d-5
      real(8):: v2d(2), sin_phi,cos_phi,pi

      axis=axis0
      pi=2.0d0*asin(1.0d0)
c      axis=(/1.0,1.0,1.0/)
c (rsin(θ)cos(φ),rsin(θ)sin(φ),rcos(θ))
      r2=0.0d0
      do i=1,3
      r2=r2+axis(i)**2
      enddo 
      r2=1.0d0/sqrt(r2)
      axis=axis*r2

c      write(*,*) axis

      if ( abs(axis(1)) < eps .and. abs(axis(2)) < eps ) then
           ! cos(t)=0
           if (axis(3) > 0.0d0 ) then
           theta=0.0d0
           else
           theta=pi
           endif 
           phi=0.0d0 
      else 
           v2d=(/axis(2),axis(1)/)
           r2= v2d(1)**2+v2d(2)**2
           r2=1.0d0/sqrt(r2)
           v2d=v2d*r2
           call sincos2ang(v2d(1),v2d(2), phi ) 
           sin_phi=sin(phi)
           cos_phi=cos(phi)
           if (sin_phi < eps )  then
              ! use cos(phi) 
              v2d(2)= axis(3)
              v2d(1)= axis(1)/cos_phi 
              call sincos2ang(v2d(1),v2d(2), theta)  
           else 
              ! use sin(phi) 
              v2d(2)=axis(3)
              v2d(1)=axis(2)/sin_phi
              call sincos2ang(v2d(1),v2d(2), theta)  
           endif

      endif 

c      write(*,*) theta,phi
c       call eulerangle2axis(theta,phi,axis2)

c      write(*,*) axis2

      end  subroutine axis2eulerangle



      subroutine sincos2ang(vs, vc, 
     o  angle)
      implicit none
      real(8),intent(in):: vs, vc
      real(8),intent(out):: angle 
      real(8):: angle_s, angle_c

      angle_s = asin(vs)
      angle_c = acos(vc) 
      if (vs.ge.0 .and.  vc.ge.0 ) then 
         !  [0:pi/2]
         angle = angle_s
      else if (vs.ge.0 .and. vc.le.0) then
         !  [pi/2:pi]
         angle = angle_c 
      else if (vs.le.0 .and. vc.le.0 ) then
         ! [pi:pi*3/2] 
         angle = -angle_c 
      else 
         ! [pi*3/2: 2*pi]
         angle = -angle_s
      endif
      
      end subroutine sincos2ang


c
c only axisc can determine rotated c-axis
c axisa determine new a-axis in the rotated frame
c
c find rotation matrix defind by axisa and axisc
c optional : psi0,  omit rotation defined by axisa
c
      subroutine axisac2rot(axisa,axisc,
     o  rot,
     i   psi0)
      implicit none
      real(8):: axisa(3),axisc(3)
      real(8):: rot(3,3)
      real(8),optional:: psi0

      real(8):: anewtmp(3),pi,orig_caxis(3),crot(3),theta,phi,psi
      integer:: i,j 
      real(8):: ddot

      pi=2.0d0*asin(1.0d0)

#ifdef DEBUG_WRITE
      write(*,130) 'new c-axis=', axisc / sqrt(ddot(3,axisc,1,axisc,1))
#endif

c output theta,phi
      call axis2eulerangle(axisc,theta,phi)
#ifdef DEBUG_WRITE
      write(*,130) 'theta,phi(deg)=',theta*180/pi,phi*180/pi
#endif

      anewtmp=axisa
#ifdef DEBUG_WRITE
      write(*,130) 'input temporary a-axis=',axisa
#endif

c output anewtmp,psi 
      call new_a_axis (axisc,theta,phi, anewtmp,psi )
#ifdef DEBUG_WRITE
      write(*,130) 'projected a-axis=',anewtmp
      write(*,130) 'psi(deg)=',psi/pi*180.0
#endif
 130  format(a,3F20.10)

      if (present(psi0)) then
        psi=psi0
        write(*,'(a,f10.2)') 'eulerangle::axisac2rot, comment: optinal psi found: force psi(deg)=',psi/pi*180.0 
      endif 

c output rot 
      call eulerangle2rot(theta,phi,psi,rot)

c check crot is the same
      orig_caxis=(/0.0d0,0.0d0,1.0d0/)
      crot=0.0d0
      do i=1,3
        do j=1,3
        crot(i)=rot(i,j)*orig_caxis(j)
        enddo
      enddo
#ifdef DEBUG_WRITE
      write(*,130) 'rotated c-axis=',crot
#endif
c end check

      end subroutine axisac2rot


c find rot by theta,phi,psi 
      subroutine eulerangle2rot(theta,phi,psi,
     o  rot)
      implicit none
      real(8),intent(in):: theta,phi,psi
      real(8),intent(out):: rot(3,3)

      real(8):: rotz(3,3),roty(3,3),rotz2(3,3),rottmp(3,3)
      integer:: i,j,k
c  x=r sin t cos p
c  y=r sin t sin t
c  z=r cos t

c
c    -t on y
c     c  0 -s
c     0  1 0
c     s  0  c
      roty=0.0d0
      roty(1,1)= cos(-theta) 
      roty(3,1)= sin(theta) 
      roty(2,2)= 1.0d0
      roty(1,3)= -sin(theta)
      roty(3,3)= cos(theta)

c    p on z
c    c  -s 0
c    s  c  0
c    0  0  1

      rotz=0.0d0
      rotz(1,1)= cos(phi)
      rotz(2,1)= -sin(phi)
      rotz(1,2)= sin(phi)
      rotz(2,2)= cos(phi)
      rotz(3,3)=1.0d0

c    p on z
c    c  -s 0
c    s  c  0
c    0  0  1

      rotz2=0.0d0
      rotz2(1,1)= cos(psi)
      rotz2(2,1)= -sin(psi)
      rotz2(1,2)= sin(psi)
      rotz2(2,2)= cos(psi)
      rotz2(3,3)=1.0d0

c
c  rot= rotz2*roty*rotz
c
      rottmp=0.0d0
      do i=1,3
        do j=1,3
          do k=1,3
          rottmp(k,i)= rottmp(k,i)+ roty(k,j)*rotz(j,i)
          enddo 
        enddo
      enddo

      rot=0.0d0
      do i=1,3
        do j=1,3
          do k=1,3
          rot(i,k)= rot(i,k)+ rotz2(k,j)*rottmp(j,i)
          enddo
        enddo
      enddo

c how to use? 
c       do i=1,3
c        do j=1,3
c        crot(i)=rot(i,j)*v1(j)
c        enddo
c      enddo
c
c or 
c 
c      call rotdlmm(rot, 1, nl, dlmm) 
c       do im=-il,il
c         do im2=-il,il
c           rotvalL(:,:,:,im)=  rotvalL(:,:,:,im) + 
c     .  valL(:,:,:,im2)*dlmm(im2,im,il,1)
c         enddo
c       enddo 

      end  subroutine eulerangle2rot


c  find projected a-axis and angle between the rotated a-axis defined by euler angle and the new a-axis. 
c     anewtmp0= anewtmp0 - (axisc0,anewtmp0) axisc0 
c     phi = angle betwwen anewtmp0 and a-axis defined by theta, phi 
c 
       subroutine new_a_axis(axisc0,theta,phi,
     =  anewtmp0,
     o  psi)
       implicit none
       real(8),intent(in):: axisc0(3),theta,phi
       real(8):: anewtmp0(3),psi

       real(8) :: v1(3), arot(3),r2,cos_t,angle_t,rot(3,3)
       real(8):: ctmp(3),direction,anewproj(3),f,anewtmp(3),axisc(3),
     .   crossc(3)
      integer:: i,j
      real(8):: ddot,pi
      real(8):: eps=1.0d-5

      pi=2.0d0*asin(1.0d0)
#ifdef DEBUG_WRITE
      write(*,*) 'new_a_axis'
#endif
c...normalize anewtmp
      r2=0.0d0
      do i=1,3
        r2=r2+anewtmp0(i)**2
      enddo
      anewtmp=anewtmp0/sqrt(r2 )

c...normalize axisc 
      r2=0.0d0
      do i=1,3
        r2=r2+axisc0(i)**2
      enddo
      axisc=axisc0/sqrt(r2 )

 100  format(a,3F20.10)       

#ifdef DEBUG_WRITE
      write(*,100) 'c-axis=',axisc
#endif
c...find rotation matrix by theta and phi
c...rotate (1,0,0) to find rotated a-axis
c...rota=rotated a-axis
      call eulerangle2rot(theta,phi,0.0d0,rot)
      v1=(/1.0d0,0.0d0,0.0d0/)
      arot=0.0d0
      do i=1,3
        do j=1,3
          arot(i)=arot(i)+rot(i,j)*v1(j)
        enddo
      enddo
#ifdef DEBUG_WRITE
      write(*,100) 'arot=',arot

      write(*,100) 'anewtmp=',anewtmp
      write(*,100) 'axisc=',axisc 
#endif
c... inner product (anewtmp, axisc)
      f= ddot(3,anewtmp,1,axisc,1)
#ifdef DEBUG_WRITE
      write(*,100) 'f=',f 
#endif

      if ( abs(abs(f)-1.0d0) < eps ) then 

         write(*,*) 'eulerangle::new_a_axis, warning: axis is almost parallel to anewtmp '

         write(*,*) axisc
         write(*,*) anewtmp
         write(*,100) 'f=',f

         stop 'eulerangle::new_a_axis,'

      endif

c... anewproj = anewtmp- (anewtmp, axisc) axisc 
      anewproj=anewtmp - f*axisc

      r2=0.0d0
      do i=1,3
        r2=r2+anewproj(i)**2
      enddo
      anewproj=anewproj/sqrt(r2)

#ifdef DEBUG_WRITE
      write(*,100) 'proj a-axis=',anewproj
#endif

c find angle between arot and anewproj

#ifdef DEBUG_WRITE
      write(*,100)'arot=',arot
#endif

c... angle between arot(defined by euler angle theta phi)  and anewproj
      cos_t= ddot(3,arot,1,anewproj,1)
      angle_t= acos(cos_t)

c... may be angle is -psi? 
      call cross_x(anewtmp, anewproj,crossc)

      direction=ddot(3,crossc,1,axisc,1)

      if (direction>0.0d0) then 
             psi=angle_t ! maybe minus
      else
             psi=-angle_t ! maybe minus
      endif 

c...store output 
      anewtmp0=anewproj

      end subroutine new_a_axis


      end module m_eulerangle

