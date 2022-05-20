!========================================================================================
   subroutine cal_hroi ( instrument, grid_id, iroi, ngxn ) 

   implicit none

   character (len=8), intent(in)          :: instrument
   integer, intent(in)                    :: grid_id, iroi
   integer, intent(inout)                 :: ngxn

!   if ( instrument == 'Radar   ' .or. instrument == 'aircft  ' .or.  &
!        instrument == 'metar   ' .or. instrument == 'satwnd  ' ) then
!!Successive covariance localization (SCL) technique:
!!http://hfip.psu.edu/fuz4/2011/WengZhang2011MWR.pdf
   if ( instrument == 'Radar   ' ) then
        if ( grid_id == 1 ) then
             ngxn = 1
        else if ( grid_id == 2 ) then
             ngxn = 1
             if( mod(iroi,3) == 1 ) ngxn = 3
        else if ( grid_id == 3 ) then
             ngxn = 1
             if( mod(iroi,3) == 1 ) ngxn = 3
             if( mod(iroi,9) == 1 ) ngxn = 9
        else if ( grid_id == 4 ) then
             ngxn = 1
             if( mod(iroi,3) == 1 ) ngxn = 3
             if( mod(iroi,9) == 1 ) ngxn = 9
             !if( mod(iroi,27) == 1) ngxn = 27
        endif
   else
!        ngxn = 3**(grid_id-1)
! need to adjust HROI in namelist.enkf to match the (hroi_n_grid)*(grid spacing)=hroi_km
! namelist.enkf sets hroi_n_grid
       ngxn = 1
   endif

   end subroutine cal_hroi
!========================================================================================
   subroutine corr(dx,dy,dz,ngx,ngz,corr_coef)
! This is alternative routine to calculate corr_coef, if schur_matrix
! requires too much memory to be calculated before hand.

  implicit none

  real, intent(in)    :: dx,dy,dz        !dx: x distance(grids) from model grid to obs
  integer, intent(in) :: ngx, ngz        !ngx: horrizontal cutting off distances
  integer :: i,j,k
  real, intent(out)   :: corr_coef

  integer :: horradi
  real :: k1
  real :: comp_cov_factor
  real :: horrad, distance, distanceh

! calculate horizontal radius at height k
     horrad = (real(ngx)/real(ngz))*sqrt(real(ngz**2. - dz**2.))
     horradi = int(horrad)   ! grid points within the radius

! equivalence of k in terms of dx  ! added FZ 2004/09/09
     k1 = dz * real(ngx) / real(ngz)

        distanceh = sqrt(real(dx**2. + dy**2.))   ! hor. distance from z-axis

        if ((dx==0.) .and. (dy==0.) .and. (dz==0.)) then
           corr_coef = 1.

        else if (distanceh<=horrad) then

           distance  = sqrt( real(dx**2. + dy**2. + k1**2.))   ! 3-d distance from obs

           corr_coef  = comp_cov_factor(dble(distance),dble(ngx/2.))

        else
           corr_coef  = 0.

        end if

   end subroutine corr
!========================================================================================
   subroutine corr_matrix(nx,nz,cmatr)

! written by Altug Aksoy, 05/20/2003

! this subroutine computes the coefficient matrix to be used for
! compact correlation function calculations. the matrix is computed
! once and stored so that it can be refered to when calculating the gain.

  implicit none

  integer, intent(in) :: nx,nz
  integer :: i,j,k
  integer :: centerx, centerz, horradi, totradi

  real :: k1
  real , intent(out), dimension(2*nx+1,2*nx+1,2*nz+1) :: cmatr
  real :: comp_cov_factor
  real :: horrad, distance, distanceh, term1, term2, totrad, schur

  cmatr = 0.

  centerx = nx+1   ! location of origin (obs) within the matrix
  centerz = nz+1

  do k = -nz,nz

! calculate horizontal radius at height k
     horrad = (real(nx)/real(nz))*sqrt(real(nz**2. - k**2.))
     horradi = int(horrad)   ! grid points within the radius

! equivalence of k in terms of dx  ! added FZ 2004/09/09
     k1 = real(k) * real(nx) / real(nz)

     do j = -horradi,horradi
     do i = -horradi,horradi

        distanceh = sqrt(real(i**2. + j**2.))   ! hor. distance from z-axis

        if ((i==0) .and. (j==0) .and. (k==0)) then
           cmatr(centerx,centerx,centerz) = 1.

        else if (distanceh<=horrad) then

           distance  = sqrt( real(i**2. + j**2. + k1**2.))   ! 3-d distance from obs
!           distance  = sqrt( real(i**2. + j**2. ))   ! 2-d distance from obs

           schur  = comp_cov_factor(dble(distance),dble(nx/2.))
           cmatr(centerx+i, centerx+j,centerz+k) = schur

        end if

      enddo
      enddo

   enddo

   end subroutine corr_matrix
!========================================================================================
   subroutine corr_matrix_h(nx,cmatr)
! this subroutine computes the coefficient matrix to be used for
! compact correlation function calculations. the matrix is computed
! once and stored so that it can be refered to when calculating the gain.

  implicit none

  integer, intent(in) :: nx
  integer :: i,j,k
  integer :: centerx, horradi, totradi

  real , intent(out), dimension(2*nx+1,2*nx+1) :: cmatr
  real :: comp_cov_factor
  real :: horrad, distance, distanceh, term1, term2, totrad, schur

  cmatr = 0.

  centerx = nx+1   ! location of origin (obs) within the matrix
     do j = -nx, nx
     do i = -nx, nx
       distanceh = sqrt(real(i**2. + j**2.))   ! hor. distance from z-axis

        if ((i==0) .and. (j==0)) then
           cmatr(centerx,centerx) = 1.

        else if (distanceh<=real(nx)) then
           distance  = sqrt( real(i**2. + j**2. ))   ! 2-d distance from obs

           schur  = comp_cov_factor(dble(distance),dble(nx/2.))
           cmatr(centerx+i, centerx+j) = schur

        end if

     enddo
     enddo

  end subroutine corr_matrix_h
!==============================================================================
   function comp_cov_factor(z_in, c)

   implicit none

   real comp_cov_factor
   double precision z_in, c
   double precision z, r

!  Computes a covariance cutoff function from Gaspari and Cohn
!  (their eqn. 4.10) QJRMS, 125, 723-757.

!  z_in is the distance while c is the cutoff distance.
!  For distances greater than 2c, the cov_factor returned goes to 0.

   z = dabs(z_in)
   r = z / c

   if(z >= 2*c) then
      comp_cov_factor = 0.0
   else if(z >= c .and. z < 2*c) then
      comp_cov_factor =                                               &
          ( ( ( ( r/12.  -0.5 )*r  +0.625 )*r +5./3. )*r  -5. )*r     &
                                                 + 4. - 2./(3.*r)
   else
      comp_cov_factor =                                               &
          ( ( ( -0.25*r +0.5 )*r +0.625 )*r  -5./3. )*r**2 + 1.
   endif

   end function comp_cov_factor
!==============================================================================
   subroutine corr_elf(varname, satid, ch, ca, ngx, ngz, kk)
! This is routine to give correlation coefficient from ELF

   implicit none

   character (len=10), intent(in)  :: varname
   character (len=12), intent(in)  :: satid
   integer, intent(in)    :: ch
   real, intent(in)       :: ca
   integer, intent(inout) :: ngx, ngz
   real, intent(inout)    :: kk
   ! input for ELFs computed offline
   integer, parameter     :: n_ca = 17
   real, parameter        :: clevs = 2.0
   integer, dimension(n_ca) :: ngx_list, ngz_list, kk_list
 
   ngx_list = ngx
   ngz_list = ngz
   kk_list  = int(kk)
 
   if (trim(adjustl(satid)) == 'abi_gr' .or. trim(adjustl(satid)) == 'ahi_h8') then
    if (ch == 8) then
      select case (trim(adjustl(varname)))
        case('T', 'TSK')
          ngx_list = (/66, 72, 74, 90, 94, 84, 88, 96, 96, 80, 84, 92, 88, 80, 94, 98, 78 /)
          ngz_list = (/87, 137, 145, 164, 160, 169, 183, 208, 208, 214, 220, 220, 214, 220, 214, 203, 198 /)
          kk_list = (/45, 44, 45, 40, 41, 39, 36, 31, 31, 30, 29, 29, 30, 29, 30, 32, 31 /)
        case('QVAPOR')
          ngx_list = (/68, 72, 72, 72, 84, 80, 80, 80, 82, 82, 80, 84, 72, 76, 76, 84, 96 /)
          ngz_list = (/82, 104, 124, 160, 169, 174, 178, 188, 198, 198, 183, 193, 188, 178, 178, 178, 203 /)
          kk_list = (/40, 41, 42, 40, 39, 38, 37, 35, 32, 33, 35, 34, 34, 37, 37, 36, 32 /)
        case('QCLOUD')
          ngx_list = (/68, 70, 82, 80, 76, 76, 68, 72, 62, 72, 76, 68, 80, 78, 64, 78, 74 /)
          ngz_list = (/193, 178, 188, 188, 178, 193, 188, 193, 193, 198, 208, 214, 209, 214, 198, 208, 203 /)
          kk_list = (/33, 34, 35, 34, 37, 33, 33, 34, 33, 32, 31, 30, 30, 29, 31, 31, 30 /)
        case('QRAIN')
          ngx_list = (/78, 88, 64, 76, 80, 68, 64, 76, 80, 68, 58, 66, 86, 90, 82, 66, 72 /)
          ngz_list = (/247, 233, 239, 241, 239, 255, 264, 273, 273, 266, 277, 284, 273, 273, 273, 273, 266 /)
          kk_list = (/25, 27, 26, 25, 26, 24, 23, 22, 22, 22, 21, 21, 22, 22, 22, 22, 22 /)
        case('QICE')
          ngx_list = (/74, 88, 86, 92, 76, 80, 84, 76, 82, 88, 86, 74, 88, 82, 74, 72, 70 /)
          ngz_list = (/148, 183, 178, 174, 174, 169, 165, 169, 169, 178, 178, 155, 155, 155, 144, 152, 152 /)
          kk_list = (/33, 36, 37, 38, 36, 38, 39, 38, 39, 37, 37, 42, 42, 42, 40, 38, 39 /)
        case('QSNOW')
          ngx_list = (/80, 92, 80, 88, 84, 104, 94, 76, 72, 78, 84, 72, 80, 82, 72, 68, 76 /)
          ngz_list = (/226, 203, 193, 188, 220, 214, 208, 183, 174, 188, 188, 165, 169, 164, 132, 123, 169 /)
          kk_list = (/28, 32, 34, 34, 29, 30, 31, 33, 36, 34, 35, 38, 38, 40, 41, 40, 39 /)
        case('QGRAUP')
          ngx_list = (/120, 88, 68, 74, 96, 86, 72, 76, 64, 58, 76, 68, 64, 74, 80, 64, 70 /)
          ngz_list = (/255, 226, 193, 188, 208, 203, 188, 193, 208, 198, 188, 174, 160, 174, 169, 169, 160 /)
          kk_list = (/24, 28, 34, 35, 31, 32, 33, 33, 31, 32, 34, 36, 37, 37, 36, 38, 39 /)
        case('PH')
          ngx_list = (/32, 0, 74, 84, 96, 74, 52, 38, 78, 74, 72, 74, 76, 32, 32, 76, 68 /)
          ngz_list = (/99, 0, 214, 255, 226, 155, 96, 93, 220, 125, 114, 178, 188, 73, 83, 209, 214 /)
          kk_list = (/29, 1, 20, 24, 28, 30, 30, 29, 29, 31, 31, 32, 33, 31, 31, 30, 29 /)
        case('MU')
          ngx_list = (/0, 0, 32, 32, 32, 52, 72, 62, 60, 72, 62, 52, 48, 60, 48, 48, 60 /)
          ngz_list = (/0, 0, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999 /)
          kk_list = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
        case('P', 'PSFC')
          ngx_list = (/32, 0, 32, 32, 32, 64, 114, 122, 100, 100, 84, 64, 72, 56, 78, 62, 96 /)
          ngz_list = (/75, 0, 72, 72, 72, 214, 220, 220, 214, 214, 214, 214, 214, 203, 214, 220, 214 /)
          kk_list = (/50, 1, 48, 49, 49, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 29, 30 /)
        case('U', 'U10')
          ngx_list = (/84, 72, 72, 76, 74, 92, 104, 128, 114, 88, 68, 68, 80, 104, 70, 104, 32 /)
          ngz_list = (/183, 151, 128, 125, 151, 193, 208, 226, 226, 220, 214, 214, 208, 208, 208, 203, 157 /)
          kk_list = (/35, 42, 43, 44, 43, 34, 31, 28, 28, 29, 30, 30, 31, 31, 31, 32, 32 /)
        case('V', 'V10')
          ngx_list = (/98, 82, 76, 80, 88, 96, 90, 120, 88, 76, 76, 76, 72, 88, 84, 72, 92 /)
          ngz_list = (/239, 198, 174, 183, 188, 233, 233, 226, 220, 215, 226, 220, 187, 214, 208, 193, 188 /)
          kk_list = (/26, 33, 36, 35, 34, 27, 27, 28, 29, 28, 28, 28, 29, 30, 31, 34, 35 /)
      end select
    else
       write(*,*) 'ELFs ERROR!!! Please prepare ELFs for channel ',ch
    endif
   else
      write(*,*) 'ELFs ERROR!!! Please prepare ELFs for ',trim(adjustl(satid))
   endif
   ngx = ngx_list(min(int(ca/clevs)+1,n_ca))
   ngz = ngz_list(min(int(ca/clevs)+1,n_ca))
   kk = real(kk_list(min(int(ca/clevs)+1,n_ca)))
   end subroutine corr_elf
!========================================================================================

