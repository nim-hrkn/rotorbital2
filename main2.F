#ifdef DEBUG_WRITE
#undef DEBUG_WRITE
#endif
      program tst1
      implicit none

      call read_lmto 
      call write_lmto 
      call sub1

      stop
   
      contains
c----------------------------------------
      subroutine read_lmto
      use m_readqg
      use m_readeigen
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3

      integer:: nwin,incwfin
      real(8):: efin

      integer:: i

      call read_BZDATA()
      call pshprt(60)
      nwin=-999
      efin=-999d0
      incwfin=-1
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'

      end subroutine read_lmto
c--------------------------------
      subroutine write_lmto
      use m_genallcf_v3, only: alat,plat,natom,pos,nclass,z 
      implicit none

      integer:: i,ic
      write(*,*) '---write_lmto---'
      write(*,*) 'alat=',alat
      do i=1,3
       write(*,*) plat(:,i)
      enddo 

      write(*,*) 'natom=',natom
      do i=1,natom
         write(*,*) z(i),pos(:,i)
      enddo

      end subroutine write_lmto

c----------------------------------------
      subroutine sub1
      use m_wancenterdef 
      use m_eulerangle
      use m_search_lma 
      implicit none
      integer,parameter:: lmax=2
      integer,parameter:: nl=lmax+1
      integer,parameter:: lmx=2*(lmax-1)
c..........
      integer,parameter:: mng = 5
      real(8):: dlmm(-(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,mng)

c.........
      integer:: iitem
      integer,parameter:: mitem=100
      character(20):: item(mitem) 

      character(200):: str
      integer:: ierr

      type(t_wancenterdef):: v
c............
      integer:: nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      integer,allocatable:: m_indx(:),n_indx(:),l_indx(:),ibas_indx(:)
      integer:: l,m

      character(5):: lmstr
      integer:: ibas, n_nlmalist=3, nlmalist(3),nfound

c...........
      real(8):: axisc(3),axisa(3),rot(3,3),eps=1.0d-4
      integer:: im,im2,il,ia1,ia2

c..........
      logical:: outputwf=.false.

c...... read mnla data
      call read_hbe_d(
     o  nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg )

      write(*,*)'nlmto=',nlmtot
      allocate(m_indx(nlmtot),n_indx(nlmtot),l_indx(nlmtot),ibas_indx(nlmtot))

      call read_mnla_cphi(nlmtot,
     o m_indx,n_indx,l_indx,ibas_indx )


      do while (.true.) 

c....... read center information
       read(*,'(a)',iostat=ierr)  str 
       if (ierr.ne.0) exit
      call wancenterdef_separate(str,item,mitem,iitem) 
      call wancenterdef_key(item,iitem, v,ierr,stoponerror=.false.) 
#if DEBUG_WRITE
      call t_wancenterdef_print(v) 
#endif

c..... lmstr->lm
      lmstr= v%lmstr
      call lmstr2lm(lmstr,l,m)
      v%l=l
      v%m=m 
      ibas= v%icenter
      call search_lma(l,m,ibas,
     i    nlmtot,m_indx,n_indx,l_indx, ibas_indx,
     i    n_nlmalist,
     o    nlmalist,nfound)
#ifdef DEBUG_WRITE
      write(*,*)'nfound=',nfound
#endif

c..... axis -> rot -> dlmm
      axisc=v%axis(:,3)
      axisa=v%axis(:,1)
#if DEBUG_WRITE
      write(*,*) 'axisa=',axisa
      write(*,*) 'axisc=',axisc
#endif
      call axisac2rot(axisa,axisc,rot)

#if DEBUG_WRITE
      write(*,*)'rot'
      do im=1,3
        write(*,*) rot(:,im)
      enddo
#endif

      call rotdlmm(rot,1,nl,dlmm) 
      il= v%l
      im= v%m 
#if DEBUG_WRITE
      write(*,150) 'dlmm',il,lmstr,v%axis(:,3),v%axis(:,1)
 150  format('#',a,i5,1x,a,' c=',3f10.3,' a=',3f12.2)
      write(*,160) 'dlmm=',(dlmm(im2,im,il,1) ,im2=-il,il) 
 160  format('#',a,1x,9F20.10) 
#endif

      write(*,'(a,a)') "#",str(:len_trim(str))
      do im2=-il,il
         call search_lma(il,im2,ibas,
     i    nlmtot,m_indx,n_indx,l_indx, ibas_indx,
     i    n_nlmalist,
     o    nlmalist,nfound)
         ia1=nlmalist(1)
         ia2=nlmalist(2)
         call lm2lmstr(il,im2,lmstr)
         if (  abs( dlmm(im2,im,il,1)*v%coef ) > eps ) then
         write(*,170) ia1,ia2,v%wgauss, dlmm(im2,im,il,1)*v%coef,
     .      m_indx(ia1),n_indx(ia1),l_indx(ia1), ibas_indx(ia1),
     .      m_indx(ia2),n_indx(ia2),l_indx(ia2), ibas_indx(ia2), lmstr
         endif
      enddo
 170  format(2i5,2f10.5,5x,'# m,n,l,a=',8i3,' lm= ',a)

      if (outputwf) then
      call writemesh_l(il, dlmm(-il:il,im,il,1), nl,lmax,lmx, v%id)
      endif 

      enddo

      end subroutine sub1
c----------------------------------------
      end program tst1


