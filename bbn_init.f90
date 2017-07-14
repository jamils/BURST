
       subroutine bbn_init
 
!------linkages.

!      called by - [subroutine] main_init
!      calls     - none


!------remarks.

!      Initiates the BBN variables.


!------modules.

       !use bbnvar
       use bbnvar_v2
       use ratenp
       use gauss


       implicit none


!------local variables.

       integer nplgndrbins
       integer nplagurbins
       logical successflag
       integer glagbesselbins


!------procedure.


!20----input initialization information

       !read in reaction parameters:
       iform(:) = int(reacpr(:,2))!reaction type.
       ii(:)    = int(reacpr(:,3))!incoming nuclide type.
       jj(:)    = int(reacpr(:,4))!incoming nuclide type.
       kk(:)    = int(reacpr(:,5))!outgoing nuclide type.
       ll(:)    = int(reacpr(:,6))!outgoing nuclide type.
       rev(:)   = reacpr(:,7)     !reverse reaction coefficient.
       q9(:)    = reacpr(:,8)     !energy released.

       irun = 1 !do full run.
       isize = nnuc
       jsize = nrec

!------set values from bbn_params.ini.

       open (unit=20, file='bbn_params.ini', status='unknown')

       read (20,*) bbncompflag
       read (20,*) cy
       read (20,*) ct
       read (20,*) t9i
       read (20,*) t9f
       read (20,*) ytmin
       read (20,*) inc
       read (20,*) dt1
       read (20,*) dtlow
       read (20,*) epstol
       read (20,*) nplgndrbins
       read (20,*) nplagurbins
       read (20,*) glagbesselbins

       close (unit=20)


!------allocate arrays:

       allocate(npnuoccarray(1,1,nbins)) !assign dimension of npnuoccarray

       allocate(nplgndr%abscissas(nplgndrbins)) !assign dimension of nplgndr%abscissas
       allocate(nplgndr%weights(nplgndrbins)) !assign dimension of nplgndr%weights

       allocate(sclgndra(nplgndrbins)) !assign dimension of sclgndra
       allocate(sclgndrw(nplgndrbins)) !assign dimension of sclgndrw

       allocate(nplagur%abscissas(nplagurbins)) !assign dimension of nplagur%abscissas
       allocate(nplagur%weights(nplagurbins)) !assign dimension of nplagur%weights

       allocate(sclagura(nplagurbins)) !assign dimension of sclagura
       allocate(sclagurw(nplagurbins)) !assign dimension of sclagurw

       allocate(glagbessel%abscissas(glagbesselbins)) !assign dimension of glag%besselabscissas
       allocate(glagbessel%weights(glagbesselbins)) !assign dimension of glagbessel%weights

!------initialize arrays:

       call cgqf(1,nplgndr%abscissas,nplgndr%weights,successflag)
       call cgqf(5,nplagur%abscissas,nplagur%weights,successflag,0._dl)

       call cgqf(5,glagbessel%abscissas,glagbessel%weights,successflag,0._dl)
       glagbessel%weights(:) = exp(glagbessel%abscissas(:))*glagbessel%weights(:)

       return


       end subroutine bbn_init
