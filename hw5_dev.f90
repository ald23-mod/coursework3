!M3C 2018 Homework 3
!Ernest Bitkivskij
!01258803



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
  integer :: i,j,k,l,t
  real(kind=8), allocatable, dimension(:,:,:) :: r !random numbers

  real(kind=8) :: n2inv
  integer, dimension(n,n,m) :: nb,nc
  integer, dimension(n+2,n+2,m) :: s2
  real(kind=8), dimension(n+2,n+2,m) :: f2, f2s2
  real(kind=8), dimension(n,n,m):: nbinv,a,f,pden,p

  allocate(r(n,n,m))

  s=1
  l = (n+1)/2
  s(l,l,:) = 0

  n2inv = 1.d0/dble(n*n)
  fc_ave(1) = sum(s)*(n2inv/m)
  s2 = 0
  f2 = 0.d0


  do t=1,nt !Start the year loop
    !$OMP parallel do
    do k = 1,m
      do i = 1,n
        do j = 1,n
          nb(i,j,k) = 8
        end do
      end do
    end do
    !$OMP end parallel do

    !$OMP parallel do
    do k=1,m
      do i = 2,n-1
        nb(1,i,k) = 5
        nb(n,i,k) = 5
        nb(i,1,k) = 5
        nb(i,n,k) = 5
      end do

      nb(1,1,k) = 3
      nb(1,n,k) = 3
      nb(n,1,k) = 3
      nb(n,n,k) = 3

      do i = 1,n
        do j = 1,n
          nbinv(i,j,k) = 1.d0/nb(i,j,k)
        end do
      end do
    end do
    !$OMP end parallel do

    call random_number(r)

    !Set up coefficients for fitness calculation in matrix, a

    !$OMP parallel do
    do k = 1,m
      do i = 1,n
        do j = 1,n
          a(i,j,k) = 1
        end do
      end do

      do i = 1,n
        do j = 1,n
          if(s(i,j,k)==0) then
            a(i,j,k) = tr_b
          end if
        end do
      end do
    end do
    !$OMP end parallel do

    !Create s2 by adding boundary of zeros to s

    !$OMP parallel do
    do k =1,m
      do i = 2,n+1
        do j = 2,n+1
          s2(i,j,k) = s(i-1,j-1,k)
        end do
      end do
    end do
    !$OMP end parallel do

    !Count number of C neighbours for each point

    !$OMP parallel do
    do k = 1,m
      do i = 1,n
        do j = 1,n
          nc(i,j,k) = s2(i,j,k)   + s2(i,j+1,k)   + s2(i,j+2,k)   + &
                      s2(i+1,j,k)                 + s2(i+1,j+2,k) + &
                      s2(i+2,j,k) + s2(i+2,j+1,k) + s2(i+2,j+2,k)
        end do
      end do
    end do
    !$OMP end parallel do

    !Calculate fitness matrix, f
    

    !$OMP parallel do
    do k = 1,m
      do i = 1,n
        do j = 1,n
          f(i,j,k) = nc(i,j,k)*a(i,j,k)
        end do
      end do
    end do
    !$OMP end parallel do

    !$OMP parallel do
    do k = 1,m
      do i = 1,n
        do j = 1,n
          if (s(i,j,k)==0) then
            f(i,j,k) = f(i,j,k) + (nb(i,j,k)-nc(i,j,k))*tr_e
          end if
        end do
      end do
    end do
    !$OMP end parallel do

    !$OMP parallel do
    do k = 1,m
      do i = 1,n
        do j = 1,n
          f(i,j,k) = f(i,j,k)*nbinv(i,j,k)
        end do
      end do
    end do
    !$OMP end parallel do

    !Calculate probability matrix, p

    !$OMP parallel do
    do k = 1,m
      do i = 2,n+1
        do j = 2,n+1
          f2(i,j,k) = f(i-1,j-1,k)
        end do
      end do
    end do
    !$OMP end parallel do

    !$OMP parallel do
    do k = 1,m
      do i = 1,n+2
        do j = 1,n+2
          f2s2(i,j,k) = f2(i,j,k)*s2(i,j,k)
        end do
      end do
    end do
    !$OMP end parallel do

    !Total fitness of cooperators in community

    !$OMP parallel do
    do k = 1,m
      do i = 1,n
        do j = 1,n
          p(i,j,k) = f2s2(i,j,k)   + f2s2(i,j+1,k)   + f2s2(i,j+2,k)   + &
                     f2s2(i+1,j,k) + f2s2(i+1,j+1,k) + f2s2(i+1,j+2,k) + &
                     f2s2(i+2,j,k) + f2s2(i+2,j+1,k) + f2s2(i+2,j+2,k)
        end do
      end do
    end do
    !$OMP end parallel do

    !Total fitness of all members of community

    !$OMP parallel do
    do k = 1,m
      do i = 1,n
        do j = 1,n
          pden(i,j,k) = f2(i,j,k)   + f2(i,j+1,k)   + f2(i,j+2,k)   + &
                        f2(i+1,j,k) + f2(i+1,j+1,k) + f2(i+1,j+2,k) + &
                        f2(i+2,j,k) + f2(i+2,j+1,k) + f2(i+2,j+2,k)
        end do
      end do
    end do
    !$OMP end parallel do

    !$OMP parallel do
    do k = 1,m
      do i = 1,n
        do j = 1,n
          p(i,j,k) = (p(i,j,k)/pden(i,j,k))*tr_g + 0.5d0*(1.d0-tr_g) !probability matrix
        end do
      end do
    end do
    !$OMP end parallel do

    !Set new affiliations based om probability matrix and random numbers stored in R

    !$OMP parallel do
    do k = 1,m
      do i = 1,n
        do j = 1,n
          s(i,j,k) = 0
        end do
      end do
    end do
    !$OMP end parallel do

    !$OMP parallel do
    do k = 1,m
      do i = 1,n
        do j = 1,n
          if (R(i,j,k)<=p(i,j,k)) then
            s(i,j,k) = 1
          end if
        end do
      end do
    end do
    !$OMP end parallel do

    fc_ave(t+1) = sum(s)*(n2inv/m)

  end do

  deallocate(r)

end subroutine simulate2_omp

end module tribes
