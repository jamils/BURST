

       module intermod


!------linkages.

!      used by main

!------remarks.

!      Module for interfaces to allow non-explicit shaped arrays.


       implicit none




       interface bbn_yes_interface
         subroutine bbn_yes(obh2,neffparm,ypri,yd,yhe3,yli7,nubbn)
           use bbnvar_v2
           implicit none
           real(dl), intent(in) :: obh2 !Baryon number
           real(dl), intent(in) :: neffparm !Parameter for changing Neff
           real(dl), intent(out) :: ypri !4He Mass Fraction
           real(dl), intent(out) :: yd !D relative abundance
           real(dl), intent(out) :: yhe3 !3He relative abundance
           real(dl), intent(out) :: yli7 !7Li relative abundance
           type(nuvar), dimension(:,:), intent(inout), target :: nubbn !output neutrino occupation probabilities
         end subroutine bbn_yes
       end interface bbn_yes_interface

       interface sdratio_interface
         subroutine sdratio(obh2in,neffparm,ypri,tneff,nuvarin)
           use sdratiovar
           implicit none
           real(dl), intent(in) :: obh2in !input baryon number
           real(dl), intent(in) :: neffparm !parameter for Neff
           real(dl), intent(in) :: ypri !primordial 4He mass fraction
           real(dl), intent(inout) :: tneff !Neff target in; Neff calculated out
           type (nuvar), intent(inout), dimension(:,:), optional, target :: nuvarin !nu occupation fractions and info.
         end subroutine sdratio
       end interface sdratio_interface

       interface assign_interface
         subroutine assign_recom(nurecom,moccfrac,numass,neffparm)
           use mainvar
           implicit none
           type(nuvar), intent(out) :: nurecom !nu parameters for recom
           !15Apr2015 EG: Use nuvar instead of real below.
           real(dl), dimension(:), intent(in) :: moccfrac !nu occfrac in mass eigenbasis
           real(dl), intent(in) :: numass !nu mass values
           real(dl), intent(in) :: neffparm !Parameter used for dilution factor
         end subroutine assign_recom
       end interface assign_interface


       end module intermod
