      module m_search_lma
      implicit none

      contains
c-----------------------------
      subroutine search_lmstra(lmstr,ibas,
     i    nlmto,m_indx,n_indx,l_indx, ibas_indx,
     i    nlist,
     o    list,nfound)
      integer,intent(in):: nlmto,ibas,nlist
      integer,intent(in):: m_indx(nlmto),n_indx(nlmto),l_indx(nlmto),ibas_indx(nlmto)
      character(*),intent(in):: lmstr
      integer:: list(nlist)
      integer:: nfound

      integer:: l,m

      call lmstr2lm(lmstr,l,m)
      call search_lma(l,m,ibas,
     i    nlmto,m_indx,n_indx,l_indx, ibas_indx,
     i    nlist,
     o    list,nfound)

      end subroutine search_lmstra

c-----------------------------

      subroutine search_lma(l,m,ibas,
     i    nlmto,m_indx,n_indx,l_indx, ibas_indx, 
     i    nlist,
     o    list,nfound)
      implicit none
      integer,intent(in):: nlmto,ibas,l,m,nlist
      integer,intent(in):: m_indx(nlmto),n_indx(nlmto),l_indx(nlmto),ibas_indx(nlmto)
      integer:: list(nlist)
      integer:: nfound

      
      integer:: ix,ilist
      ilist=0
      do ix=1,nlmto
         if ( ibas.eq.ibas_indx(ix) .and. m_indx(ix).eq.m .and. l_indx(ix).eq.l ) then
            ilist=ilist+1  
            if (ilist>nlist) then
                write(*,*) 'error search_ibaslm, ilist>nlist',ilist,nlist 
                stop 'search_ibaslm'
            endif
            list(ilist)= ix 
         endif
      enddo

      if (ilist.eq.0) then
         write(*,*) 'error failed to find l,m,ibas=',l,m,ibas
         stop 'search_ibaslm'
      endif
      nfound=ilist 
      end subroutine search_lma
c----------------------------

      subroutine lmstr2lm(lmstr,l,m)
      implicit none
      character(5):: lmstr
      integer:: l,m
      integer,parameter:: lmax=2
      character(5):: lmstrlist(-lmax:lmax,0:lmax)

      lmstrlist(0,0)='s'
      lmstrlist(-1,1)='py'
      lmstrlist(0,1)='pz'
      lmstrlist(1,1)='px'
      lmstrlist(-2,2)='dxy'
      lmstrlist(-1,2)='dyz'
      lmstrlist(0,2)='dz2'
      lmstrlist(1,2)='dxz'
      lmstrlist(2,2)='dx2y2'

       do l=0,lmax
         do m=-l,l
            if (lmstr.eq.lmstrlist(m,l)) then
               return
            endif
         enddo
       enddo
       write(*,*) 'failed to find lmstr'
       stop

      end subroutine lmstr2lm

      subroutine lm2lmstr(l,m,lmstr)
      implicit none
      character(5):: lmstr
      integer:: l,m
      integer,parameter:: lmax=2
      character(5):: lmstrlist(-lmax:lmax,0:lmax)

      lmstrlist(0,0)='s'
      lmstrlist(-1,1)='py'
      lmstrlist(0,1)='pz'
      lmstrlist(1,1)='px'
      lmstrlist(-2,2)='dxy'
      lmstrlist(-1,2)='dyz'
      lmstrlist(0,2)='dz2'
      lmstrlist(1,2)='dxz'
      lmstrlist(2,2)='dx2y2'

      lmstr=lmstrlist(m,l)

      end subroutine lm2lmstr

c-----------------------------
      subroutine read_mnla_cphi(ldim2,
     o m_indx,n_indx,l_indx,ibas_indx )
      implicit none
      integer,intent(in):: ldim2 ! = nlmto 
      integer:: m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2)

      integer:: ifoc ,ix,ixx,iopen 

c      ifoc = iopen('@MNLA_CPHI',1,0,0)
      ifoc= 200
      open(ifoc,file='@MNLA_CPHI',status='old')
c      ldim2 = nlmto
      read(ifoc,*)
c      allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
      do ix =1,ldim2
        read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
        if(ixx/=ix) stop  'failed to readin @MNLA_CPHI'
      enddo
      close(ifoc)
 
      end subroutine read_mnla_cphi

c-----------------------------

      subroutine read_hbe_d(
     o  nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg )
      implicit none
      integer:: ifhbe
      integer:: nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg

      integer:: iopen
c      ifhbe= iopen('hbe.d',1,0,0) 
      ifhbe      = 200
      open(ifhbe,file='hbe.d',status='old') 
      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      close(ifhbe)

      end subroutine read_hbe_d

c-----------------------------

      end  module m_search_lma
