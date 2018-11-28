!M3C 2018 Homework 3
!This module contains four module variables and two subroutines;
!one of these routines must be developed for this assignment.
!Module variables--
! tr_b, tr_e, tr_g: the parameters b, e, and g in the tribe competition model
! numthreads: The number of threads that should be used in parallel regions within simulate2_omp
!
!Module routines---
! simulate2_f90: Simulate tribal competition over m trials. Return: all s matrices at final time
! and fc at nt+1 times averaged across the m trials.
! simulate2_omp: Same input/output functionality as simulate2.f90 but parallelized with OpenMP

module tribes
  use omp_lib
  implicit none
  integer :: numthreads
  real(kind=8) :: tr_b,tr_e,tr_g
contains

!Simulate m trials of Cooperator vs. Mercenary competition using the parameters, tr_b and tr_e.
!Input:
! n: defines n x n grid of villages
! nt: number of time steps
! m: number of trials
!Output:
! s: status matrix at final time step for all m trials
! fc_ave: fraction of cooperators at each time step (including initial condition)
! averaged across the m trials
subroutine simulate2_f90(n,nt,m,s,fc_ave)
  implicit none
  integer, intent(in) :: n,nt,m
  integer, intent(out), dimension(n,n,m) :: s
  real(kind=8), intent(out), dimension(nt+1) :: fc_ave
  integer :: i1,j1
  real(kind=8) :: n2inv
  integer, dimension(n,n,m) :: nb,nc
  integer, dimension(n+2,n+2,m) :: s2
  real(kind=8), dimension(n,n,m) :: f,p,a,pden,nbinv
  real(kind=8), dimension(n+2,n+2,m) :: f2,f2s2
  real(kind=8), allocatable, dimension(:,:,:) :: r !random numbers
  !timing variables
  integer(kind=8) :: clock_t1,clock_t2,clock_rate
  real(kind=8) :: cpu_t1,cpu_t2,clock_time

  call system_clock(clock_t1)

  call cpu_time(cpu_t1)


  !---Problem setup----
  !Initialize arrays and problem parameters

  !initial condition
  s=1
  j1 = (n+1)/2
  s(j1,j1,:) = 0

  n2inv = 1.d0/dble(n*n)
  fc_ave(1) = sum(s)*(n2inv/m)

  s2 = 0
  f2 = 0.d0

  !Calculate number of neighbors for each point
  nb = 8
  nb(1,2:n-1,:) = 5
  nb(n,2:n-1,:) = 5
  nb(2:n-1,1,:) = 5
  nb(2:n-1,n,:) = 5
  nb(1,1,:) = 3
  nb(1,n,:) = 3
  nb(n,1,:) = 3
  nb(n,n,:) = 3

  nbinv = 1.d0/nb
  allocate(r(n,n,m))
  !---finished Problem setup---


  !----Time marching----
  do i1=1,nt

    call random_number(r) !Random numbers used to update s every time step

    !Set up coefficients for fitness calculation in matrix, a
    a = 1
    where(s==0)
      a=tr_b
    end where

    !create s2 by adding boundary of zeros to s
    s2(2:n+1,2:n+1,:) = s

    !Count number of C neighbors for each point
    nc = s2(1:n,1:n,:) + s2(1:n,2:n+1,:) + s2(1:n,3:n+2,:) + &
         s2(2:n+1,1:n,:)                  + s2(2:n+1,3:n+2,:) + &
         s2(3:n+2,1:n,:)   + s2(3:n+2,2:n+1,:)   + s2(3:n+2,3:n+2,:)

    !Calculate fitness matrix, f----
    f = nc*a
    where(s==0)
      f = f + (nb-nc)*tr_e
    end where
    f = f*nbinv
    !-----------

    !Calculate probability matrix, p----
    f2(2:n+1,2:n+1,:) = f
    f2s2 = f2*s2

    !Total fitness of cooperators in community
    p = f2s2(1:n,1:n,:) + f2s2(1:n,2:n+1,:) + f2s2(1:n,3:n+2,:) + &
           f2s2(2:n+1,1:n,:) + f2s2(2:n+1,2:n+1,:)  + f2s2(2:n+1,3:n+2,:) + &
          f2s2(3:n+2,1:n,:)   + f2s2(3:n+2,2:n+1,:)   + f2s2(3:n+2,3:n+2,:)

    !Total fitness of all members of community
    pden = f2(1:n,1:n,:) + f2(1:n,2:n+1,:) + f2(1:n,3:n+2,:) + &
           f2(2:n+1,1:n,:) + f2(2:n+1,2:n+1,:)  + f2(2:n+1,3:n+2,:) + &
          f2(3:n+2,1:n,:)   + f2(3:n+2,2:n+1,:)   + f2(3:n+2,3:n+2,:)


    p = (p/pden)*tr_g + 0.5d0*(1.d0-tr_g) !probability matrix
    !----------

    !Set new affiliations based on probability matrix and random numbers stored in R
    s = 0
    where (R<=p)
        s = 1
    end where

    fc_ave(i1+1) = sum(s)*(n2inv/m)

  end do
  call cpu_time(cpu_t2)
  print *, 'elapsed cpu time (seconds) =',cpu_t2-cpu_t1
  call system_clock(clock_t2,clock_rate)
  print *, 'elapsed wall time (seconds)= ', dble(clock_t2-clock_t1)/dble(clock_rate)

end subroutine simulate2_f90

!Simulate m trials of Cooperator vs. Mercenary competition using the parameters, tr_b and tr_e.
!Same functionality as simulate2_f90, but parallelized with OpenMP
!Parallel regions should use numthreads threads.
!Input:
! n: defines n x n grid of villages
! nt: number of time steps
! m: number of trials
!Output:
! s: status matrix at final time step for all m trials
! fc_ave: fraction of cooperators at each time step (including initial condition)
! averaged across the m trials
subroutine simulate2_omp(n,nt,m,s,fc_ave)
  implicit none
  integer, intent(in) :: n,nt,m
  integer, intent(out), dimension(n,n,m) :: s
  real(kind=8), intent(out), dimension(nt+1) :: fc_ave
  integer :: k,i1,j1,i2,i3,i4,i5,i8,i9
  real(kind=8), allocatable, dimension(:,:,:) :: r !random numbers
  !Add further variables as needed
  real(kind=8) :: n2inv
  integer, dimension(n,n,m) :: nb,nc
  integer, dimension(n+2,n+2,m) :: s2
  real(kind=8), dimension(n+2,n+2,m) :: f2, f2s2
  real(kind=8), dimension(n,n,m) :: nbinv,a,f,pden,p
  !timing variables
  integer(kind=8) :: clock_t1,clock_t2,clock_rate
  real(kind=8) :: cpu_t1,cpu_t2,clock_time

  call system_clock(clock_t1)

  call cpu_time(cpu_t1)

  !initial condition and r allocation (does not need to be parallelized)
  s=1
  j1 = (n+1)/2
  s(j1,j1,:) = 0
  allocate(r(n,n,m))
  !------------------

  !Add code here
  n2inv = 1.d0/dble(n*n)
  fc_ave(1) = sum(s)*(n2inv/m)
  s2 = 0
  f2 = 0.d0

  !---finished Problem setup---


  !---Time marching---

    do i3 = 1,nt

    !$OMP parallel do reduction(+:fc_ave)
    do k = 1,m
      do i1 = 1,n
        do i2 = 1,n
          nb(i1,i2,k) = 8
        end do
      end do
      do i4 = 2,n-1
        nb(1,i4,k) = 5
        nb(n,i4,k) = 5
        nb(i4,1,k) = 5
        nb(i4,n,k) = 5
      end do
    nb(1,1,k) = 3
    nb(1,n,k) = 3
    nb(n,1,k) = 3
    nb(n,n,k) = 3
      do i1 = 1,n
        do i2 = 1,n
          nbinv(i1,i2,k) = 1.d0/nb(i1,i2,k)
        end do
      end do
  call random_number(r) !random numbers used to update s every time step
  !Set up coefficients for fitness calculation in matrix, a
    do i1 = 1,n
      do i2 = 1,n
        a(i1,i2,k) = 1
        if(s(i1,i2,k)==0) then
          a(i1,i2,k) = tr_b
        end if
      end do
    end do
  !Create s2 by adding boundary of zeros to s
    do i4 = 2,n+1
      do i5 = 2,n+1
        s2(i4,i5,k) = s(i4-1,i5-1,k) !Look at PARALLELIZATION for this
      end do
    end do
  !Count number of C neighbours for each point
    do i1 = 1,n
      do i2 = 1,n
          nc(i1,i2,k) = s2(i1,i2,k) + s2(i1,i2+1,k) + s2(i1,i2+2,k) + &
                          s2(i1+1,i2,k)              + s2(i1+1,i2+2,k) + &
                          s2(i1+2,i2,k) + s2(i1+2,i2+1,k) + s2(i1+2,i2+2,k)
                          !s2 needs to be global
      end do
    end do
  !Calculate fitness matrix, f----
  !We will have to use reduction

    do i1 = 1,n
      do i2 = 1,n
        f(i1,i2,k) = nc(i1,i2,k)*a(i1,i2,k)
        if (s(i1,i2,k)==0) then
          f(i1,i2,k) = f(i1,i2,k) + (nb(i1,i2,k)-nc(i1,i2,k))*tr_e
        end if
      end do
    end do

    do i1 = 1,n
      do i2 = 1,n
          f(i1,i2,k) = f(i1,i2,k)*nbinv(i1,i2,k)
      end do
    end do


  !------------------------------
  !Calculate probability matrix, p----
    do i4 = 2,n+1
      do i5 = 2,n+1
        f2(i4,i5,k) = f(i4-1,i5-1,k)
      end do
    end do
    do i8 = 1,n+2
      do i9 = 1,n+2
        f2s2(i8,i9,k) = f2(i8,i9,k)*s2(i8,i9,k)
      end do
    end do

    do i1 = 1,n
      do i2 = 1,n
          p(i1,i2,k) = f2s2(i1,i2,k) + f2s2(i1,i2+1,k) + f2s2(i1,i2+2,k) + &
                      f2s2(i1+1,i2,k) + f2s2(i1+1,i2+1,k) + f2s2(i1+1,i2+2,k) + &
                      f2s2(i1+2,i2,k) + f2s2(i1+2,i2+1,k) + f2s2(i1+2,i2+2,k)
      end do
    end do

  !Total fitness of all members of community
    do i1 = 1,n
      do i2 = 1,n
          pden(i1,i2,k) = f2(i1,i2,k) + f2(i1,i2+1,k) + f2(i1,i2+2,k) + &
                          f2(i1+1,i2,k) + f2(i1+1,i2+1,k) + f2(i1+1,i2+2,k) + &
                          f2(i1+2,i2,k) + f2(i1+2,i2+1,k) + f2(i1+2,i2+2,k)
      end do
    end do

    do i1 = 1,n
      do i2 = 1,n
        p(i1,i2,k) = (p(i1,i2,k)/pden(i1,i2,k))*tr_g + 0.5d0*(1.d0-tr_g) !probability matrix
      end do
    end do

  !--------------------------
  !Set new affiliations based om probability matrix and random numbers stored in R

    do i1 = 1,n
      do i2 = 1,n
        s(i1,i2,k) = 0
        if (R(i1,i2,k)<=p(i1,i2,k)) then
          s(i1,i2,k) = 1
        end if
      end do
    end do
    do i1 = 1,n
      do i2 = 1,n
        fc_ave(i3+1) = fc_ave(i3+1) + s(i1,i2,k)*(n2inv/m)
      end do
    end do
  end do
  !$OMP end parallel do
end do
call cpu_time(cpu_t2)
print *, 'elapsed cpu time (seconds) =',cpu_t2-cpu_t1
call system_clock(clock_t2,clock_rate)
print *, 'elapsed wall time (seconds)= ', dble(clock_t2-clock_t1)/dble(clock_rate)
deallocate(r)

end subroutine simulate2_omp
end module tribes
