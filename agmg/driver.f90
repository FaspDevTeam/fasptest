!===================================================================================!
!       SPECS (Serial and Parallel, Efficient, Cunning Solvers) (c) 2010-2011       !
!                  School of Mathematics and Computational Science                  !
!             Xiaoqiang Yue, Zhiyang Zhou, Zheng Lee and Chunsheng Feng             !
!                                Xiangtan University                                !
!===================================================================================!
!
!!
!! driver.f90
!!
!! Created by Xiaoqiang Yue on 05/05/2012
!!
!
!=======================================================================================
      program SolverAGMG
!=======================================================================================
      implicit none
      integer*4:: pid,num_rows,num_nonzeros,iter
      real :: start,finish
      real*8:: tol
      character*256:: matfile
      character*256:: rhsfile
      integer*4,dimension(:),pointer::ia,ja
      real*8,dimension(:),pointer::a,gf,u
!
!      write(*,*) 'Please enter the Problem ID:'
!      write(*,'(a$)') ' pid='
!      read 110,pid
!110   format(I5)
      pid = 3
!
      select case(pid)
          case(1)
              matfile='/home/stone/Tests/data/FDM/fdm_mat_csr_127X127X127.dat'
              rhsfile='/home/stone/Tests/data/FDM/fdm_rhs_127X127X127.dat'
          case(2)
              matfile='/home/stone/Tests/data/SPE1020.amg.dat'
              rhsfile='/home/stone/Tests/data/SPE1020.rhs.dat'
          case(3)
              matfile='/home/stone/Tests/data/SPE1040.amg.dat'
              rhsfile='/home/stone/Tests/data/SPE1040.rhs.dat'
          case(4)
              matfile='/home/stone/Tests/faspsolver/data/fdm_mat_csr_1023X1023.dat'
              rhsfile='/home/stone/Tests/faspsolver/data/fdm_rhs_1023X1023.dat'
      end select
!
      call specsf_DoubleCSRMatrixReadA(num_rows,num_nonzeros,matfile)

!
      allocate(ia(num_rows+1))
      allocate(ja(num_nonzeros))
      allocate(a(num_nonzeros))
      allocate(gf(num_rows))
      allocate(u(num_rows))
!
      call specsf_DoubleCSRMatrixReadB(num_rows,num_nonzeros,&
                                       ia,ja,a,matfile)
!
      call specsf_DoubleArrayRead(num_rows,gf,rhsfile)
!
      u=0.0d0
!
      iter=1000
      tol=1.0d-6
      call cpu_time(start)
      call dagmg(num_rows,a,ja,ia,gf,u,0,6,10,iter,tol)
      call cpu_time(finish)
      write(*,*) 'Elapsed Time',finish-start
!
      deallocate(ia)
      deallocate(ja)
      deallocate(a)
      deallocate(gf)
      deallocate(u)
!
      end
!===================================================================================!
!===================================End-Of-File=====================================!
!===================================================================================!
