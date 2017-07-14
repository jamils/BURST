
       program main


!------linkages.

!      called by - none
!      calls     - [subroutine] bbn, put_numass, nu_assign_equil,
!                               sdratio


!-----remarks.

!     Program to call bbn and sdratio.


!------modules.

       use mainvar
       !use bbnvar
       use bbnvar_v2
       use massvar
       use intermod
       use sdratiovar
       use nu_assign

#ifdef prllel
       use mpi
#endif

       implicit none


       save


!------local variables.

       integer m,n !indicies

       real(dl) barlow !lowest value of obh2loop

       real(dl) obh2loop !ob*h**2
       real(dl) tneff !\tilde{N}_{eff}
       real(dl) neffparm !parameter to change Neff

       type(nuvar), dimension(:,:), allocatable :: numain !\nu variable

       real(dl) ypri !4He Mass fraction
       real(dl) yd !D relative abundance
       real(dl) yhe3 !3He relative abundance
       real(dl) yli7 !7Li relative abundance

       real(dl), dimension(:), allocatable :: lepnum

       real(dl), dimension(:), allocatable :: stmass !sterile nu masses
       real(dl), dimension(:), allocatable :: stmix !sterile nu \sin(2\theta) relative to nu_e
       real(dl), dimension(:), allocatable :: numass !all nu masses

       real(dl) parfact

       !character (len = 6) numchar1
       !character (len = 29 + len(numchar1)) form1

       !need to change len when using a different directory (cannot
       !allocate character length with gfortran 4.7.2)

       real(dl) time1, time2 !clock times for beginning and ending run

#ifdef prllel
       integer ierr !mpi variable for error instance
#endif


!------procedure.

#ifdef prllel
       !initialize mpi
       call mpi_init(ierr)
       call mpi_comm_size(mpi_comm_world,nprocs,ierr)
       call mpi_comm_rank(mpi_comm_world,rank,ierr)
#endif


!------set values from main_params.ini.

       open (unit=10, file='main_params.ini', status='unknown')

       read (10,*) barlow

       read (10,*) mntau
       read (10,*) h
       read (10,*) och2

       read (10,*) nbins !must be 4n + 1
       read (10,*) nnu
       if (nnu.lt.3) nnu = 3
       allocate(lepnum(nnu)) !assign dimension of lepnum
       allocate(xi(nnu)) !assign dimension of xi
       read (10,*) lepnum(:)
       xi = 0._dl
       if (nnu.gt.3) then
         allocate(stmix(nnu-3)) !assign dimension of stmix, if applicable
         allocate(stmass(nnu-3)) !assign dimension of stmass, if applicable
         read (10,*) stmix(:)
         read (10,*) stmass(:)
       else
         read (10,'(/)')
       end if

       read (10,*) temp0k

       read (10,*) bbnflag
       read (10,*) transflag
       read (10,*) equilflag
       read (10,*) yeflag
       read (10,*) saveflag
       read (10,*) mswflag
       read (10,*) frepostr

       close (unit=10)


!------allocate arrays.

       allocate(numain(nnu,2)) !assign dimension of numain
       allocate(numass(nnu)) !assign dimension of numass


!------open output files.

#ifdef prllel
       if (rank.eq.0) then
#endif
       if (bbnflag) open (unit=21, file =frepostr // &
                         'bbn.dat', status='unknown')  !output file.

       if (equilflag) open (unit=22, file = frepostr // &
                           'equil.dat', status='unknown')  !output file.

       if (yeflag) open (unit=23, file = frepostr // &
                        'ye.dat', status='unknown')  !output file.
#ifdef prllel
       end if
#endif


#ifdef prllel
       if (rank.eq.0) then
#endif
       if (transflag) then
         open (unit=52, file = frepostr // &
              'sumrule_52.dat', status='unknown')  !output file.
       end if
#ifdef prllel
       end if
#endif

#ifdef prllel
       if (rank.eq.0) then
#endif
       if (transflag) then
         open (unit=53, file = frepostr // &
              'sumrule_53.dat', status='unknown')  !output file.
       end if
#ifdef prllel
       end if
#endif

#ifdef prllel
       if (rank.eq.0) then
#endif
       if (transflag) then
         open (unit=72, file = frepostr // &
              'nu_72.dat', status='unknown')  !output file.
       end if !transflag
#ifdef prllel
       end if
#endif

#ifdef prllel
       if (rank.eq.0) then
#endif
       if (transflag) then
         open (unit=73, file = frepostr // &
              'nu_73.dat', status='unknown')  !output file.
       end if !transflag
#ifdef prllel
       end if
#endif

#ifdef prllel
       if (rank.eq.0) then
#endif
       if (transflag) then
         open (unit=74, file = frepostr // &
              'nu_74.dat', status='unknown')  !output file.
       end if !transflag
#ifdef prllel
       end if
#endif

#ifdef prllel
       if (rank.eq.0) then
#endif
       if (transflag) then
         open (unit=75, file = frepostr // &
              'nu_75.dat', status='unknown')  !output file.
       end if !transflag
#ifdef prllel
       end if
#endif


#ifdef prllel
       if (rank.eq.0) then
#endif
       open (unit=82, file = frepostr // &
            'params_82.dat', status='unknown')  !output file.
#ifdef prllel
       end if
#endif

#ifdef prllel
       if (rank.eq.0) then
#endif
       open (unit=87, file = frepostr // &
            'abunds_87.dat', status='unknown')  !output file.
#ifdef prllel
       end if
#endif


       call main_init


       !15Apr2015 EG: Do not change nbins during parameter runs.



!------get time at start of computation.

       call cpu_time(time1) 



       obh2loop = barlow

       parfact = 1._dl
       do n=1,2
         do m=1,nnu
           allocate(numain(m,n)%occfrac(nbins))
           if (m.le.3) then
             if (n.eq.1) xi(m) = lep_to_xi(lepnum(m))
             call nu_assign_equil(numain(m,n),parfact*xi(m),0._dl,1._dl)
           else
             call nu_assign_null(numain(m,n),stmass(m-3),stmix(m-3))
             !call nu_assign_equil(numain(m,n),xi(m),stmass(m-3),1._dl)
           end if
         end do
         parfact = -1._dl
       end do
       !if (rank.eq.0) write (*,*) numain(1,1)%occfrac(1:3)
       !if (rank.eq.0) write (*,*) numain(1,2)%occfrac(1:3)


       neffparm = 0._dl
       if (bbncompflag) then
         !call bbn(obh2loop,neffparm,ypri,yd,yhe3,yli7,numain)
         call bbn_yes(obh2loop,neffparm,ypri,yd,yhe3,yli7,numain)
       else
         ypri = 0._dl
       end if
#ifdef prllel
       if (rank.eq.0) then
#endif
       !write (*,*) obh2loop, tneff, ypri, yd
#ifdef prllel
       end if !(rank.eq.0)
#endif


       call numass_assign(numain(1:3,:))


       tneff = 3._dl
       call sdratio(obh2loop,neffparm,ypri,tneff,numain)

#ifdef prllel
       if (rank.eq.0) then
#endif
       write (*,*) '\omega_b: ', obh2loop
       if (size(numain, 1).eq.4) write (*,*) 'm_{\nu_s}: ', numain(4,1)%mass
       if (size(numain, 1).eq.4) write (*,*) '\sin(2\theta): ', numain(4,1)%mix
       write (*,*) '\tneff: ', tneff
       write (*,*) 'Y_P: ', ypri
       write (*,*) 'D/H: ', yd
#ifdef prllel
       end if !(rank.eq.0)
#endif



       deallocate(lepnum)
       deallocate(xi)
       if (nnu.gt.3) then
         deallocate(stmix)
         deallocate(stmass)
       end if
       do n=1,2
         do m=1,nnu
           deallocate(numain(m,n)%occfrac)
           !nullify(numain(m,n)%occfrac)
         end do
       end do
       deallocate(numain)
       deallocate(numass)
       call trans_dealloc

!------close data files.

#ifdef prllel
       if (rank.eq.0) then
#endif

       if (bbnflag) close(unit=21, status = 'keep')
       if (equilflag) close(unit=22, status = 'keep')
       if (yeflag) close(unit=23, status = 'keep')
       if (transflag) then
         close(unit=52, status = 'keep')
         close(unit=53, status = 'keep')
         close(unit=72, status = 'keep')
         close(unit=73, status = 'keep')
         close(unit=74, status = 'keep')
         close(unit=75, status = 'keep')
       end if !transflag
       close(unit=82, status = 'keep')
       close(unit=87, status = 'keep')

#ifdef prllel
       end if !(rank.eq.0)
#endif


!------get time at end of computation and print duration.

       call cpu_time(time2) 

#ifdef prllel
       if (rank.eq.0) then
#endif

       write (*,50) time2 - time1

50     format('Duration (in s): ',1pe12.3)

#ifdef prllel
       end if !(rank.eq.0)
#endif


#ifdef prllel
       call mpi_barrier(mpi_comm_world,ierr)
       call mpi_finalize(ierr)
#endif


       stop


       end program main

