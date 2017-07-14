
       module nse


!------modules.

       use mainvar


       implicit none


!------data.

!    nuclide and corresponding number
!    --------------------------------
!    1) D
!    2) T
!    3) 3He
!    4) 4He
!    5) 6Li
!    6) 7Li
!    7) 7Be

       real(dl), dimension(4) :: atomz = (/1._dl,1._dl,2._dl,2._dl/) !atomic numbers
       real(dl), dimension(4) :: atomn = (/1._dl,2._dl,1._dl,2._dl/) !neutron numbers
       real(dl), dimension(4) :: deltaq = (/2.2246_dl,8.4818_dl,7.7180_dl,28.296_dl/) !binding energies (in MeV)
       real(dl), dimension(4) :: atomm = (/1875.6_dl,2808.9_dl,2809.4_dl,3727.4_dl/) !atomic masses (in MeV)
       real(dl), dimension(4) :: spina =(/3._dl,2._dl,2._dl,1._dl/) !spin partition functions
       !real(dl), dimension(7) :: z = (/1._dl,1._dl,2._dl,2._dl,3._dl,3._dl,4._dl/) !atomic numbers
       !real(dl), dimension(7) :: n = (/1._dl,2._dl,1._dl,2._dl,3._dl,4._dl,3._dl/) !neutron numbers
       !real(dl), dimension(7) :: deltaq = (/2.2246_dl,8.4818_dl,7.7180_dl,28.296_dl,30.461_dl,37.711_dl,35.556_dl/) !binding energies (in MeV)
       !real(dl), dimension(7) :: m = (/1875.6_dl,2808.9_dl,2809.4_dl,3727.4_dl,5603.1_dl,6535.4_dl,6536.2_dl/) !atomic masses (in MeV)
       !real(dl), dimension(7) :: spina =(/3._dl,2._dl,2._dl,1._dl,3._dl,4._dl,4._dl/) !spin partition functions


       contains


!--------------------------------------------------------------------

       subroutine nse_det(tpl,yneu,ypro,rb)


!------linkages.

!      called by - [subroutine] save_out
!      calls     - none




!------remarks.

!      Calculates NSE abundances.

       real(dl), intent(in) :: tpl !plasma temperature in MeV
       real(dl), intent(in) :: yneu !free neutron abundance
       real(dl), intent(in) :: ypro !free proton abundance
       real(dl), intent(in) :: rb !baryon energy density in cgs

!------local variables.

       integer j
       real(dl) nb !baryon number density in MeV
       real(dl), dimension(4) :: ya !nse abundances
       real(dl), dimension(4) :: equil_y !nse relative abundances

       real(dl) :: amu = 931.494061_dl !atomic mass unit in MeV
       real(dl) :: spinp = 2._dl !proton spin
       real(dl) :: spinn = 2._dl !neutron spin

!------procedure.

       !all units MeV.

       !baryon number density:
       nb = rb/amu/232011.568_dl       

       do j=1,4
         ya(j) = (2._dl*pi)**(1.5_dl*(atomz(j) + atomn(j) - 1._dl))* &
              spina(j)/spinp**atomz(j)/spinn**atomn(j)* &
              ((atomz(j) + atomn(j))/amu**(atomz(j) + atomn(j) - 1._dl))**1.5_dl* &
              tpl**(1.5_dl*(1._dl - atomz(j) - atomn(j)))*exp(deltaq(j)/tpl)* &
              nb**(atomz(j) + atomn(j) - 1._dl)* &
              ypro**atomz(j)*yneu**atomn(j)
       end do


       do j=1,4
         equil_y(j) = ya(j)/ypro
         if (equil_y(j).gt.2._dl) then
           equil_y(j) = 2._dl
         end if
       end do

       write (22,25) tpl,(equil_y(j),j=1,4)

25     format (1p10e12.3)


       return


       end subroutine nse_det


!--------------------------------------------------------------------

       subroutine nse_get(tpl,yneu,ypro,rb,y)


!------linkages.

!      called by - [subroutine] save_out
!      calls     - none




!------remarks.

!      Calculates NSE abundances.

       real(dl), intent(in) :: tpl !plasma temperature in MeV
       real(dl), intent(in) :: yneu !free neutron abundance
       real(dl), intent(in) :: ypro !free proton abundance
       real(dl), intent(in) :: rb !baryon energy density in cgs
       real(dl), dimension(:), intent(out) :: y !abundances

!------local variables.

       integer j
       real(dl) nb !baryon number density in MeV
       real(dl), dimension(4) :: ya !nse abundances
       real(dl), dimension(4) :: equil_y !nse relative abundances

       real(dl) :: amu = 931.494061_dl !atomic mass unit in MeV
       real(dl) :: spinp = 2._dl !proton spin
       real(dl) :: spinn = 2._dl !neutron spin

!------procedure.

       !all units MeV.

       !baryon number density:
       nb = rb/amu/232011.568_dl       

       do j=1,4
         y(j+2) = (2._dl*pi)**(1.5_dl*(atomz(j) + atomn(j) - 1._dl))* &
              spina(j)/spinp**atomz(j)/spinn**atomn(j)* &
              ((atomz(j) + atomn(j))/amu**(atomz(j) + atomn(j) - 1._dl))**1.5_dl* &
              tpl**(1.5_dl*(1._dl - atomz(j) - atomn(j)))*exp(deltaq(j)/tpl)* &
              nb**(atomz(j) + atomn(j) - 1._dl)* &
              ypro**atomz(j)*yneu**atomn(j)
       end do


       return


       end subroutine nse_get


       end module nse
