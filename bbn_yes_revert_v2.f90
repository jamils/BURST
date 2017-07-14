
       subroutine bbn_yes_revert_v2(t,dt,bbnvalues)


!------linkages.

!      called by - [subroutine] bbn_yes
!      calls     - none


!------remarks.

!      Reverts to a bbn computation.


!------modules.

       use bbnvar_v2
#ifdef prllel
       use mpi
#endif


       implicit none


!------throughput variables.

       real(dl), intent(out) :: t !time
       real(dl), intent(out) :: dt !time step
       type(bbnevolvar), intent(inout) :: bbnvalues !bbn evolution variable


!------local variables.

#ifdef prllel
       integer ierr
#endif


!------procedure.

#ifdef prllel
       if (rank.eq.0) then
#endif

       !open(unit=85,file=frepostr // 'savestate',form='unformatted', &
       !         status='unknown')
       open(unit=85,file= 'savestate',form='unformatted', &
                !status='unknown')
                status='old')

       read(85) t
       read(85) dt
       read(85) bbnvalues%tpl
       read(85) bbnvalues%asf
       read(85) bbnvalues%phie
       read(85) bbnvalues%sen
       read(85) bbnvalues%hv
       read(85) bbnvalues%y(:)
       read(85) bbnvalues%nuoccprob(1,1)%occfrac(:)
       read(85) bbnvalues%nuoccprob(1,2)%occfrac(:)
       read(85) bbnvalues%nuoccprob(2,1)%occfrac(:)
       read(85) bbnvalues%nuoccprob(2,2)%occfrac(:)
       read(85) bbnvalues%nuoccprob(3,1)%occfrac(:)
       read(85) bbnvalues%nuoccprob(3,2)%occfrac(:)
       read(85) thm(:)
       read(85) xout(:,:)
       read(85) thmout(:,:)
       read(85) t9out(:)
       read(85) tout(:)
       read(85) dtout(:)
       read(85) senout(:)
       read(85) hubout(:)
       read(85) tcmevout(:)
       read(85) tcmdec
       read(85) asfdec
       read(85) is
       read(85) ip
       read(85) it
       read(85) nurhoflag
       read(85) snu, stot
       read(85) dydt(:)


       close(85,status='keep')

       !write(40,*) t
       !write(40,*) dt
       !write(40,*) bbnvalues%tpl
       !write(40,*) bbnvalues%asf
       !write(40,*) bbnvalues%phie
       !write(40,*) bbnvalues%sen
       !write(40,*) bbnvalues%hv
       !write(40,*) bbnvalues%y(:)
       !write(40,*) bbnvalues%nuoccprob(1,1)%occfrac(:)
       !write(40,*) bbnvalues%nuoccprob(1,2)%occfrac(:)
       !write(40,*) bbnvalues%nuoccprob(2,1)%occfrac(:)
       !write(40,*) bbnvalues%nuoccprob(2,2)%occfrac(:)
       !write(40,*) bbnvalues%nuoccprob(3,1)%occfrac(:)
       !write(40,*) bbnvalues%nuoccprob(3,2)%occfrac(:)
       !write(40,*) thm(:)
       !write(40,*) xout(:,:)
       !write(40,*) thmout(:,:)
       !write(40,*) t9out(:)
       !write(40,*) tout(:)
       !write(40,*) dtout(:)
       !write(40,*) senout(:)
       !write(40,*) hubout(:)
       !write(40,*) tcmevout(:)
       !write(40,*) tcmdec
       !write(40,*) asfdec
       !write(40,*) is
       !write(40,*) ip
       !write(40,*) it
       !write(40,*) nurhoflag
       !write(40,*) snu, stot
       !write(40,*) dydt(:)
       !call flush(40)
#ifdef prllel
       end if !(rank.eq.0)
       call mpi_barrier(mpi_comm_world,ierr)

       call mpi_bcast(t,1,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(dt,1,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(bbnvalues%tpl,1,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(bbnvalues%asf,1,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(bbnvalues%phie,1,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(bbnvalues%sen,1,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(bbnvalues%hv,1,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(bbnvalues%y,nnuc,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(bbnvalues%nuoccprob(1,1)%occfrac,nbins,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(bbnvalues%nuoccprob(1,2)%occfrac,nbins,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(bbnvalues%nuoccprob(2,1)%occfrac,nbins,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(bbnvalues%nuoccprob(2,2)%occfrac,nbins,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(bbnvalues%nuoccprob(3,1)%occfrac,nbins,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(bbnvalues%nuoccprob(3,2)%occfrac,nbins,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(thm,14,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(xout,nnuc*itmax,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(thmout,8*itmax,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(t9out,itmax,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(dtout,itmax,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(senout,itmax,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(hubout,itmax,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(tcmevout,itmax,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(tcmdec,1,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(asfdec,1,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(is,1,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(ip,1,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(it,1,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(nurhoflag,1,mpi_logical,0,mpi_comm_world,ierr)
       call mpi_bcast(snu,1,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(stot,1,mpi_double_precision,0,mpi_comm_world,ierr)
       call mpi_bcast(dydt,nnuc,mpi_double_precision,0,mpi_comm_world,ierr)
#endif

       return


       end subroutine bbn_yes_revert_v2

