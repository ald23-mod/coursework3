!Anas Lasri Doukkali
!CID: 01209387
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
  !As asked by question one, I will now proceed to explain my approach to
  !the parallelization of the simulate_f90 subroutine above.
  !After carrefuly reading the question and analyzing the code we realize that the
  !variable 'm' will be the variable subject to parallelization. However at first sight
  !we also note that the loop concerning the year variable 'nt' will produce
  !inneficiencies if left as the outside loop. We hence need to think of a way
  !of having the 'm' loop outside, this way optimizing our problem as best possible
  !The way to tackle the issue is as follows:
  !We will first change 'r' to be two-dimensional instead of 3-dimensional as in
  !the previous version of the code. This is so that we can reduce massively the
  !storage space to be used in our machine and hence free more memory for our
  !calculations.
  !Secondly, now that we have the 'm' loop outside and that r is two-dimensional
  !we can reuse 'r' in our calculations and hence reduce the dimension of all the
  !variables used in the vectorization except the one we are considering as returns
  !inside both loops, and hence, with the exception of 's' all the other variables
  !will now also be two-dimensional. This is highly efficient and will allow us
  !to have a much faster code with the parallelization.
  !After having set the code as explained above I will proceed now to parallelize.
  !This is done by first setting the number of threads with the command that was
  !given in the lectures.
  !We know open the parallel region by calling OMP with firstprivate, private, and
  !a reduction as well. As these methods were called in lectures, firstprivates
  !will be used for the variables that were defined before the parallel region and
  !initiallized also. As for private, we pick the variable for which each thread
  !will require an own copy of the variable. Lastly we notice that we need to use
  !a 'plus' reduction for one of our variables given the way it was defined 'fc_ave'
  !This function was tested using the edited version of the hw3_main.f90 that runs
  !said function and provided the correct and expected results whether accuracy
  !or performance-wise.
  !You can also see I added the time clock controls to be able to track the
  !performance of the codes.
  implicit none
  integer, intent(in) :: n,nt,m
  integer, intent(out), dimension(n,n,m) :: s
  real(kind=8), intent(out), dimension(nt+1) :: fc_ave
  real(kind=8), dimension(nt+1) ::fc
  integer :: j1,k,t
  real(kind=8), allocatable, dimension(:,:) :: r !random numbers
  real(kind=8) :: n2inv
  integer, dimension(n,n) :: nb,nc
  integer, dimension(n+2,n+2) :: s2
  real(kind=8), dimension(n+2,n+2) :: f2,f2s2
  real(kind=8), dimension(n,n) :: nbinv,a,f,pden,p
  !timing variables
  integer(kind=8) :: clock_t1,clock_t2,clock_rate
  real(kind=8) :: cpu_t1,cpu_t2,clock_time
  call system_clock(clock_t1)
  call cpu_time(cpu_t1)
  !$ call omp_set_num_threads(numthreads)
  allocate(r(n,n))

  s=1
  j1 = (n+1)/2
  s(j1,j1,:) = 0

  n2inv = 1.d0/dble(n*n)
  fc_ave(1) = sum(s)*(n2inv/m)
  s2 = 0
  f2 = 0.d0

  nb = 8
  nb(1,2:n-1) = 5
  nb(n,2:n-1) = 5
  nb(2:n-1,1) = 5
  nb(2:n-1,n) = 5
  nb(1,1) = 3
  nb(1,n) = 3
  nb(n,1) = 3
  nb(n,n) = 3
  nbinv = 1.d0/nb
  !$OMP parallel do firstprivate(s2,f2),private(t,a,nc,f,f2s2,p,pden),reduction(+:fc)
  do k = 1,m

    do t = 1,nt
      call random_number(r)
      !Set up coefficients for fitness calculation in matrix, a
      a = 1
      where(s(:,:,k)==0)
        a = tr_b
      end where

      !Create s2 by adding boundary of zeros to s
      s2(2:n+1,2:n+1) = s(:,:,k)

      !Count number of C neighbours for each point
      nc = s2(1:n,1:n)   + s2(1:n,2:n+1)   + s2(1:n,3:n+2)   + &
           s2(2:n+1,1:n)                   + s2(2:n+1,3:n+2) + &
           s2(3:n+2,1:n) + s2(3:n+2,2:n+1) + s2(3:n+2,3:n+2)

      !Calculate fitness matrix, f
      f = nc*a
      where (s(:,:,k)==0)
        f = f + (nb-nc)*tr_e
      end where
      f = f*nbinv

      !Calculate probability matrix, p
      f2(2:n+1,2:n+1) = f
      f2s2 = f2*s2

      !Total fitness of cooperators in community
      p =    f2s2(1:n,1:n)   + f2s2(1:n,2:n+1)   + f2s2(1:n,3:n+2)   + &
             f2s2(2:n+1,1:n) + f2s2(2:n+1,2:n+1) + f2s2(2:n+1,3:n+2) + &
             f2s2(3:n+2,1:n) + f2s2(3:n+2,2:n+1) + f2s2(3:n+2,3:n+2)

      !Total fitness of all members of community
      pden = f2(1:n,1:n)     + f2(1:n,2:n+1)     + f2(1:n,3:n+2)     + &
             f2(2:n+1,1:n)   + f2(2:n+1,2:n+1)   + f2(2:n+1,3:n+2)   + &
             f2(3:n+2,1:n)   + f2(3:n+2,2:n+1)   + f2(3:n+2,3:n+2)

      p = (p/pden)*tr_g + 0.5d0*(1.d0-tr_g)

      !Set new affiliations based on probability matrix and random numbers stored in R
      s(:,:,k) = 0
      where (r(:,:)<=p)
        s(:,:,k) = 1
      end where

      fc(t+1) = sum(s(:,:,k))

    end do

    fc_ave = fc_ave + (fc)*(n2inv/m)

  end do
!$OMP end parallel do
call cpu_time(cpu_t2)
call system_clock(clock_t2,clock_rate)
print *, 'elapsed cpu time (seconds) =',cpu_t2-cpu_t1
print *, 'elapsed wall time (seconds)= ', dble(clock_t2-clock_t1)/dble(clock_rate)
deallocate(r)
end subroutine simulate2_omp
end module tribes
