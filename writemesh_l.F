#define DEBUG_WRITE
      subroutine writemesh_l(il,dlmm,nl,lmax,lmx, fileid)
      implicit none
      integer,intent(in):: il
      integer,intent(in):: nl,lmax,lmx
      integer,intent(in):: fileid
      real(8):: r(3), yl( (lmx+1)*(lmx+1) ), r2s,c((lmx+1)*(lmx+1)) 

      real(8),intent(in):: dlmm(-il:il) 

      integer,parameter:: nx=10,ny=10,nz=10
      real(8):: x(3),g,rr,x0(3)
      real(8):: val(-nx:nx,-ny:ny, -nz:nz,(lmx+1)*(lmx+1)),
     .  rotvalL(-nx:nx,-ny:ny, -nz:nz ),
     .  valL(-nx:nx,-ny:ny, -nz:nz, -lmax:lmax )

      integer:: ilmlist(2,100)

      integer::ifout=10
      character(50):: fname
      integer:: nllist(0:3)

      integer:: ix,iy,iz,id,im,im2,ilm,nilm,i,il2
      real*8 :: pi,range,g0,eps=1.0d-5


      if (il<0 .or. il>3 ) then
         write(*,*) 'writemesh_l: outofrange, il=',il 
         stop ' in writemesh_l'
      endif

      pi=2.0d0*asin(1.0d0)


      range=3.0d0
      write(*,*) lmax,lmx, (lmx+1)*(lmx+1) 
      do ix=-nx,nx
      do iy=-ny,ny
      do iz=-nz,nz
         x0(1)=range*dble(ix)/nx
         x0(2)=range*dble(iy)/ny
         if (nz.eq.0) then
           x0(3)=0.0d0
         else
         x0(3)=range*dble(iz)/nz
         endif
         rr=x0(1)**2+x0(2)**2+x0(3)**2 
         if (rr> eps) then
         x=x0/sqrt(rr) 
         endif
         call sylm(x,yl,lmx,r2s)
         call sylmnc(c,lmx)
         g0=nx/4
         g=exp(-rr/g0**2) 

         do id=1, (lmx+1)**2
         val(ix,iy,iz,id)=yl(id)*g*c(id)
         enddo

c         do id=1,(lmx+1)**2
c         write(10+id,'(3i5,10F15.5)') ix,iy,iz,x0,yl(id)*g*c(id)
c         enddo

      enddo
      enddo
c         do id=1,(lmx+1)**2
c         write(10+id,*) ' '
c         enddo
      enddo

#ifdef DEBUG_WRITE
      write(*,*) id, 'max','min'
      do id=1,(lmx+1)**2
        write(*,*) id, maxval(val(:,:,:,id)), minval(val(:,:,:,id))
      enddo
#endif

      nllist(0)=1
      nllist(1)=2
      nllist(2)=5
      nllist(3)=10
c  1 2 3 4 5
c  s x y z d 

      id=nllist(il) 
#ifdef DEBUG_WRITE
      write(*,*) 'copy',-il,il,id,id+1*2*il 
      write(*,*) lbound(vall,dim=4),ubound(vall,dim=4),lbound(val,dim=4),ubound(val,dim=4) 
#endif
      valL=0.0d0
      do i=-il,il
      valL(-nx:nx,-ny:ny,-nz:nz,i)= val(-nx:nx,-ny:ny,-nz:nz,id+il+i) 
      enddo

#ifdef DEBUG_WRITE
      write(*,'(a,10f10.5)')'dlmm=', dlmm(-il:il)
#endif

      ! rotate
      !   Ylm= \sum m' Ylm' Dm'm
      rotvalL=0.0d0
         do im2=-il,il
           rotvalL(:,:,:)=  rotvalL(:,:,:) + 
     .  valL(:,:,:,im2)*dlmm(im2)
         enddo

#ifdef DEBUG_WRITE
         write(*,*) 'max,min=',maxval(rotvalL(:,:,:)),
     .          minval(rotvalL(:,:,:)) 
#endif

      id=fileid
      ! rotate 

      ifout=10+id
      write(fname,'(a,i3.3,a)') 'out',fileid,'.cube'
      open(file=fname,unit=ifout) 
      write(*,*) 'open file=',fname
      write(ifout,*)'TEST', id
      write(ifout,*)'TEST'
      write(ifout,110) 4, -1.0d0*range, -1.0d0*range, -1.0d0*range
 110  format(I5,5F12.6)
 120  format(6E13.6)
      write(ifout,110) 2*nx+1, range*1.0d0/nx, 0.0d0, 0.0d0
      write(ifout,110) 2*ny+1, 0.0d0,range*1.0d0/ny, 0.0d0
      write(ifout,110) 2*nz+1, 0.0d0,0.0d0,range*1.0d0/nz
      write(ifout,110) 1, 0.0d0, 0.0d0, 0.0d0, 0.0d0
      write(ifout,110) 2, 0.0d0, 0.0d0, 0.0d0, 1.0d0
      write(ifout,110) 3, 0.0d0, 2.0d0, 0.0d0, 0.0d0
      write(ifout,110) 4, 0.0d0, 0.0d0, 3.0d0, 0.0d0
c      write(ifout,'(3x2i4)') 1 , 1

      do ix=-nx,nx
      do iy=-ny,ny
c         write(ifout,120) val(ix,iy,:,id)
         write(ifout,120) rotvalL(ix,iy,:)
      enddo
      enddo

      close(ifout)
      
      end subroutine writemesh_l

