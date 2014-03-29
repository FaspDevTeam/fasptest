!====================================================================================
      subroutine specsf_DoubleCSRMatrixReadA(n,nz,filename)
!====================================================================================
      integer*4:: n,i,nzpo,nz
      character*256:: filename
!
      open(8,file=filename,form='formatted',status='unknown')
      read(8,*) n
      do i=1,n
          read(8,*)
      enddo
      read(8,*) nzpo
      nz=nzpo-1
      close(8)
!
      return
      end
!====================================================================================
      subroutine specsf_DoubleCSRMatrixReadB(n,nz,ia,ja,a,filename)
!====================================================================================
      integer*4:: ia(1),ja(1)
      integer*4:: n,i,nz
      real*8:: a(1)
      character*256:: filename
!
      open(8,file=filename,form='formatted',status='unknown')
      read(8,*)
      do i=1,n+1
          read(8,*) ia(i)
      enddo
      do i=1,nz
          read(8,*) ja(i)
      enddo
      do i=1,nz
          read(8,*) a(i)
      enddo
      close(8)
!
      return
      end
!====================================================================================

!====================================================================================
      subroutine specsf_DoubleArrayRead(n,u,filename)
!====================================================================================
      integer*4:: n,i
      real*8:: u(1)
      character*256:: filename
!
      open(8,file=filename,form='formatted',status='unknown')
      read(8,*)
      do i=1,n
          read(8,*) u(i)
      enddo
      close(8)
!
      return
      end
!====================================================================================
