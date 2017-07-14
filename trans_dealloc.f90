
       subroutine trans_dealloc


!------linkages.

!      called by - [program] main
!      calls     - none


!------remarks.

!      Subroutine to deallocate transport arrays.


!------modules.

       use transvar


       implicit none



!------local variables.


!------procedure.


!------deallocate arrays:


       deallocate(eps_bins%abscissas)
       deallocate(eps_bins%weights)
       deallocate(glgndr%abscissas)
       deallocate(glgndr%weights)
       deallocate(glagur%abscissas)
       deallocate(glagur%weights)
       deallocate(prec_trans)
       deallocate(woccfrac)
       deallocate(ratescatt)
       deallocate(frsratescatt)
       deallocate(rateannih)
       deallocate(frsrateannih)
       deallocate(scatt)
       deallocate(frsscatt)
       deallocate(annih)
       deallocate(frsannih)
       deallocate(glegai)
       deallocate(glegwi)
       deallocate(glegai2)
       deallocate(glegwi2)
       deallocate(glagai2)
       deallocate(glagwi2)
       deallocate(glegai3)
       deallocate(glegwi3)
       deallocate(e2a)
       deallocate(e2w)
       deallocate(e3a)
       deallocate(e3w)
       deallocate(glagai3)
       deallocate(glagwi3)
       deallocate(ggenai2)
       deallocate(ggenwi2)
       deallocate(e3r2a)
       deallocate(e3r2w)
       deallocate(algndrout)
       deallocate(wlgndrout)
       deallocate(alagurout)
       deallocate(wlagurout)
       deallocate(algndrin)
       deallocate(wlgndrin)
       deallocate(alagurin)
       deallocate(wlagurin)
       deallocate(agenout)
       deallocate(wgenout)
       deallocate(agenin)
       deallocate(wgenin)
       deallocate(eouta)
       deallocate(eoutw)
#ifdef prllel
       deallocate(rsprllel)
       deallocate(frsrsprllel)
       deallocate(raprllel)
       deallocate(frsraprllel)
#endif


       return


       end subroutine trans_dealloc
