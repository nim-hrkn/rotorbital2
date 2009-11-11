#ifdef DEBUG_WRITE
#undef DEBUG_WRITE 
#endif
      module m_wancenterdef
       implicit none

      type t_wancenterdef
      integer::  id
      character(7):: lmstr
      real(8):: wgauss
      integer:: n,l,m
      integer:: icenter 
      integer:: axisinput(3)
      real(8):: axis(3,3)
      real(8):: coef
      integer:: head,tail
      end type t_wancenterdef 

      contains

      subroutine t_wancenterdef_init(v)
      type(t_wancenterdef):: v
      v%id=0
      v%lmstr=''
      v%wgauss=2.0/.529177  ! au
      v%icenter=0
      v%n=0
      v%l=1000
      v%m=1000
      v%axisinput=0
      v%axis=0.0d0
      v%axis(1,1)=1.0d0
      v%axis(2,2)=1.0d0
      v%axis(3,3)=1.0d0
      v%coef=1.0d0
      v%head=1
      v%tail=1

      end subroutine t_wancenterdef_init 

      subroutine t_wancenterdef_print(v)
       type(t_wancenterdef):: v
      write(*,'(a,i5,2H, ,$)') 'id=',v%id
      write(*,'(a,a,2H, ,$)') 'lm=',v%lmstr
      write(*,'(a,E12.3,2H, ,$)') 'wgauss=',v%wgauss
      write(*,'(a,3i5,2H, ,$)') 'nlm=',v%n,v%l,v%m 
      write(*,'(a,i5,2H, ,$)') 'centersite=',v%icenter
      write(*,'(a,3i5,2H, ,$)') 'axisinput=',v%axisinput
      write(*,'(a,9E12.3,2H, ,$)') 'axis=',v%axis
      write(*,'(a,E12.3,2H, ,$)') 'coef=',v%coef
      write(*,'(a,2i5)') 'head,tail=',v%head,v%tail
      end subroutine t_wancenterdef_print

      subroutine do_vector(str,nid,id0,cord,vec,ierr)
      implicit none
      integer:: id0
      integer:: nid
      character(*):: str(nid)
      character(4):: cord
      real(8):: vec(3)
      integer:: ierr
      integer:: i,id

      id=id0
      ierr=0

      cord=''
#ifdef DEBUG_WRITE
      write(*,*)'str(1)=',str(id) 
#endif
      if (str(id).ne.'(') then
        ierr=1
        return
      endif
      id=id+1
#ifdef DEBUG_WRITE
      write(*,*)'step1:',str(id)
#endif
      if (str(id).eq.'cart') then
         cord='cart'
      else if (str(id).eq.'abc') then
         cord='abc'
      else 
         ierr=2
         return 
      endif
#ifdef DEBUG_WRITE
      write(*,*)'cord=',cord
#endif
      id=id+1
      do i=id,id+2
         if (i>nid) then
            ierr=3
            return 
         endif
         read(str(i),*) vec(i-id+1)
      enddo 
      id=id+3
      if (str(id).eq.')') then
         id0=id 
         ierr=0
         return
      else 
         ierr=4
         return  
      endif

      end subroutine do_vector

      subroutine do_separate(str0,item,mitem,iitem) 
      implicit  none
      character(*):: str0
      integer:: mitem
      integer:: iitem
      character(20):: item(mitem)
      character(200)::str1
      integer:: mstr,nstr

      integer:: i,n,istr,istart,iend 
      character:: c

      n=len_trim(str0)
      istr=0
      str1=''
      do i=1,n
         c=str0(i:i)
         if ( c .eq.'(' .or. c.eq.')' ) then
            istr=istr+1; str1(istr:istr)=' ' 
            istr=istr+1; str1(istr:istr)=c
            istr=istr+1; str1(istr:istr)=' ' 
         else if ( c.eq.'=' ) then
            istr=istr+1; str1(istr:istr)=c
            istr=istr+1; str1(istr:istr)=' '
         else   
            istr=istr+1; str1(istr:istr)=c
         endif
      enddo 

#ifdef DEBUG_WRITE
      write(*,*)str1(:len_trim(str1))
#endif

      n=len_trim(str1)
      istart=1
      iitem=0 
      do while (.true.) 

         do while ( str1(istart:istart).eq.' ') 
           istart=istart+1
            if (istart> n ) then 
               return 
            endif
         enddo 
         iend=istart+1
         do while ( str1(iend:iend).ne.' ' )
            iend=iend+1
            if (iend> n ) exit 
         enddo 
         iitem=iitem+1
         item(iitem)=str1(istart:iend-1)
#ifdef DEBUG_WRITE
         write(*,*) iitem,item(iitem)
#endif
         istart=iend
         if (iend.eq.n) return 
      enddo 

      end subroutine do_separate

      subroutine do_key(str,mstr,v,ierr0,stoponerror)
      integer:: mstr
       integer:: ierr0
      logical,optional:: stoponerror
      type(t_wancenterdef):: v
      character(*)::str(mstr)
      real(8):: vec(3)
      character(4):: cord
      integer:: i,ierr

      character(24),parameter:: thisprogram='m_wancenterdef::do_key, '

      ierr0=0
#ifdef DEBUG_WRITE
      write(*,*) 'mstr=',mstr
#endif
      call t_wancenterdef_init(v) 
      i=1
      do while (i<=mstr)
#ifdef DEBUG_WRITE
        write(*,*) i,mstr
#endif
        if (i>=mstr) then
          return 
        endif
#ifdef DEBUG_WRITE
        write(*,*) str(i)
#endif
        select case ( str(i) ) 
          case ('id=')
              i=i+1; read(str(i),*) v%id
              i=i+1
          case ('lm=')
             i=i+1; read(str(i),*) v%lmstr
             i=i+1
          case ('rgauss=')
             i=i+1; read(str(i),*) v%wgauss
             i=i+1
          case ('n=' )
             i=i+1; read(str(i),*) v%n
             i=i+1
          case ('centersite=')
             i=i+1; read(str(i),*) v%icenter 
             i=i+1;
          case ('axis1=' )
             i=i+1;
             call do_vector(str,mstr,i,cord, vec,ierr)
             ierr0=ierr0+ierr
             v%axis(:,1)=vec
             v%axisinput(1)=1
             i=i+1;
c          case ('axis2=' )
c             i=i+1;
c             call do_vector(str,mstr,i,cord, vec,ierr)
c             ierr0=ierr0+ierr
c             v%axis(:,2)=vec
c             v%axisinput(2)=1
c             i=i+1;
          case ('axis3=' )
             i=i+1;
             call do_vector(str,mstr,i,cord, vec,ierr)
             ierr0=ierr0+ierr
             v%axis(:,3)=vec
             v%axisinput(3)=1
             i=i+1;
          case ('coef=' , 'c=' )
             i=i+1; read(str(i),*) v%coef
             i=i+1;
          case ('truncate=')
             i=i+1
             if (str(i).eq.'head') then
                v%head=1
             else if (str(i).eq.'tail') then
                v%tail=1
             endif              
             i=i+1
          case default
             write(*,*) thisprogram,'error unknown keyword, i,str=',i,str(i)
             ierr0=ierr0+1
             i=i+1
        end select
      enddo

c error check
             call lm_accept(v%lmstr,ierr)
             if (ierr.ne.0) then
                write(*,*)thisprogram,'error: uknown lm, lm=',v%lmstr
                ierr0=ierr0+1
             endif

      if ( v%head.eq.0 .and. v%tail.eq.0 ) then
           ierr0=ierr0+1
           write(*,*) thisprogram,'error: both head and tail are off'
      endif 

      if (present(stoponerror)) then
          if ( stoponerror .and. ierr0>0) then
             write(*,*) thisprogram,'stop due to request'
             stop 'many errors in input'
          endif
       endif

      end subroutine do_key

      subroutine lm_accept(str,ierr)
      character(*):: str
      integer:: ierr
      integer,parameter:: m_a=9
      character(5):: a(m_a)=(/'s','px','pz','pz','dxy','dyz','dz2','dxz','dx2y2'/)
      integer:: i
      do i=1,m_a
         if (str.eq.a(i) ) then
             ierr=0 
             return
         endif
      enddo
      ierr=1
      end subroutine lm_accept 

      end module m_wancenterdef

