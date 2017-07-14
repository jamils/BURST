
       subroutine bbn_yes_save_state_v2(t,bbnvalues,dt)


!------linkages.

!      called by - [subroutine] bbn_yes_save_out_v2
!      calls     - none


!------remarks.

!      output accumulator for save state.


!------modules.

       use bbnvar_v2


       implicit none


!------throughput variables.

       real(dl), intent(in) :: t !time
       type(bbnevolvar), intent(in) :: bbnvalues !bbn evolution variable
       real(dl), intent(in) :: dt !time step


!------local variables.

       integer i,j !indicies
       real(dl) tcmev, tplmev !temperatures
       real(dl) alln, allp, allnp


!------procedure.


#ifdef prllel
       if (rank.eq.0) then
#endif

       !open(unit=85,file=frepostr // 'savestate',form='unformatted', &
       !         status='unknown')
       open(unit=85,file= 'savestate',form='unformatted', &
                status='unknown')
       write (85) t
       write (85) dt
       write (85) bbnvalues%tpl
       write (85) bbnvalues%asf
       write (85) bbnvalues%phie
       write (85) bbnvalues%sen
       write (85) bbnvalues%hv
       write (85) bbnvalues%y(:)
       write (85) bbnvalues%nuoccprob(1,1)%occfrac(:)
       write (85) bbnvalues%nuoccprob(1,2)%occfrac(:)
       write (85) bbnvalues%nuoccprob(2,1)%occfrac(:)
       write (85) bbnvalues%nuoccprob(2,2)%occfrac(:)
       write (85) bbnvalues%nuoccprob(3,1)%occfrac(:)
       write (85) bbnvalues%nuoccprob(3,2)%occfrac(:)
       write (85) thm(:)
       write (85) xout(:,:)
       write (85) thmout(:,:)
       write (85) t9out(:)
       write (85) tout(:)
       write (85) dtout(:)
       write (85) senout(:)
       write (85) hubout(:)
       write (85) tcmevout(:)
       write (85) tcmdec
       write (85) asfdec
       write (85) is
       write (85) ip
       write (85) it-1
       write (85) nurhoflag
       write (85) snu, stot
       write (85) dydt(:)

       close(85,status='keep')

#ifdef prllel
       end if
#endif


       return


       end subroutine bbn_yes_save_state_v2

