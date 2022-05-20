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
             if( mod(iroi,27) == 1) ngxn = 27
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
   integer, parameter     :: n_ca = 19
   real, parameter        :: clevs = 2.0
   integer, dimension(n_ca) :: ngx_list, ngz_list, kk_list
 
   ngx_list = ngx
   ngz_list = ngz
   kk_list  = int(kk)
 
   if (trim(adjustl(satid)) == 'abi_gr' .or. trim(adjustl(satid)) == 'ahi_h8') then
    if (ch == 8) then
      select case (trim(adjustl(varname)))
        case('T')
          ngx_list = (/102, 36, 30, 68, 66, 50, 40, 26, 40, 36, 64, 60, 58, 56, 68, 60, 32, 102, 94 /)
          ngz_list = (/6, 10, 10, 25, 31, 39, 31, 29, 30, 120, 92, 101, 93, 120, 120, 20, 120, 101, 85 /)
          kk_list = (/45, 44, 42, 47, 49, 44, 41, 41, 42, 6, 22, 19, 21, 7, 6, 6, 9, 16, 6 /)
        case('QVAPOR')
          ngx_list = (/32, 44, 36, 28, 28, 32, 36, 40, 46, 36, 32, 30, 28, 46, 52, 46, 102, 102, 102 /)
          ngz_list = (/20, 25, 27, 29, 33, 35, 38, 45, 59, 47, 39, 40, 43, 37, 36, 50, 87, 120, 120 /)
          kk_list = (/41, 42, 42, 41, 40, 40, 38, 36, 32, 35, 37, 37, 37, 37, 36, 34, 23, 11, 13 /)
        case('QCLOUD')
          ngx_list = (/10, 2, 30, 86, 8, 10, 10, 24, 16, 14, 16, 18, 8, 18, 22, 24, 22, 22, 42 /)
          ngz_list = (/2, 2, 8, 3, 15, 11, 14, 17, 21, 26, 23, 18, 14, 18, 15, 19, 24, 10, 68 /)
          kk_list = (/40, 46, 6, 46, 39, 41, 43, 41, 40, 38, 40, 42, 44, 40, 38, 37, 36, 40, 20 /)
        case('QRAIN')
          ngx_list = (/8, 2, 12, 16, 10, 22, 24, 24, 28, 32, 32, 32, 30, 28, 30, 30, 40, 50, 44 /)
          ngz_list = (/3, 2, 35, 33, 43, 42, 48, 69, 71, 74, 81, 79, 66, 72, 72, 58, 60, 62, 77 /)
          kk_list = (/40, 44, 21, 22, 22, 18, 18, 12, 14, 16, 16, 13, 13, 15, 15, 17, 17, 15, 15 /)
        case('QICE')
          ngx_list = (/38, 52, 48, 42, 34, 38, 34, 38, 44, 34, 30, 28, 28, 32, 62, 6, 12, 48, 76 /)
          ngz_list = (/14, 21, 22, 20, 19, 21, 21, 21, 21, 23, 19, 19, 18, 18, 38, 36, 13, 7, 6 /)
          kk_list = (/39, 43, 43, 44, 45, 46, 46, 46, 47, 49, 49, 50, 50, 52, 44, 40, 35, 55, 43 /)
        case('QSNOW')
          ngx_list = (/52, 50, 44, 54, 44, 38, 38, 36, 36, 36, 34, 36, 36, 40, 40, 44, 40, 44, 42 /)
          ngz_list = (/26, 29, 38, 37, 53, 52, 58, 46, 45, 52, 49, 41, 40, 40, 33, 35, 40, 20, 47 /)
          kk_list = (/29, 32, 30, 34, 30, 31, 30, 36, 37, 39, 37, 40, 40, 42, 43, 43, 43, 50, 43 /)
        case('QGRAUP')
          ngx_list = (/52, 22, 14, 44, 16, 24, 26, 26, 26, 30, 30, 30, 32, 32, 34, 30, 40, 50, 102 /)
          ngz_list = (/19, 5, 15, 25, 44, 46, 47, 50, 51, 58, 55, 57, 51, 45, 43, 44, 51, 49, 49 /)
          kk_list = (/21, 15, 15, 36, 31, 35, 35, 37, 38, 37, 39, 40, 39, 42, 41, 42, 41, 42, 42 /)
        case('PH')
          ngx_list = (/26, 0, 16, 92, 94, 94, 60, 52, 52, 36, 32, 56, 56, 68, 70, 88, 102, 102, 98 /)
          ngz_list = (/1, 1, 5, 95, 113, 120, 120, 120, 120, 120, 120, 120, 117, 120, 120, 120, 120, 120, 118 /)
          kk_list = (/29, 6, 46, 6, 6, 6, 10, 6, 6, 8, 23, 27, 29, 32, 25, 22, 24, 24, 24 /)
        case('MU')
          ngx_list = (/0, 0, 2, 8, 6, 24, 74, 70, 40, 80, 78, 90, 64, 102, 80, 72, 102, 102, 102 /)
          ngz_list = (/999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999 /)
          kk_list = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
        case('P')
          ngx_list = (/16, 0, 2, 8, 6, 24, 74, 70, 40, 76, 72, 72, 60, 102, 78, 76, 102, 102, 102 /)
          ngz_list = (/2, 1, 2, 15, 11, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 99, 113, 120, 120 /)
          kk_list = (/51, 6, 49, 49, 50, 30, 30, 30, 30, 30, 30, 30, 30, 30, 21, 20, 20, 24, 24 /)
        case('U')
          ngx_list = (/86, 84, 46, 62, 102, 102, 34, 60, 52, 60, 50, 76, 96, 102, 68, 54, 58, 62, 102 /)
          ngz_list = (/7, 11, 27, 9, 9, 8, 12, 16, 22, 26, 15, 15, 25, 25, 42, 30, 16, 18, 16 /)
          kk_list = (/50, 40, 49, 54, 55, 55, 10, 10, 8, 7, 42, 43, 44, 44, 44, 42, 42, 54, 54 /)
        case('V')
          ngx_list = (/102, 20, 2, 100, 50, 42, 38, 38, 38, 52, 50, 54, 46, 44, 42, 34, 32, 102, 102 /)
          ngz_list = (/120, 47, 64, 11, 8, 26, 83, 81, 109, 116, 102, 102, 78, 100, 86, 51, 64, 73, 93 /)
          kk_list = (/9, 35, 55, 52, 53, 18, 6, 6, 6, 6, 6, 6, 6, 6, 6, 33, 33, 30, 18 /)
        case('PSFC')
          ngx_list = (/0, 0, 2, 2, 2, 24, 74, 70, 44, 80, 76, 74, 62, 102, 80, 76, 102, 102, 102 /)
          ngz_list = (/999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999 /)
          kk_list = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
        case('U10')
          ngx_list = (/44, 34, 18, 2, 2, 102, 62, 60, 54, 60, 72, 74, 56, 56, 88, 102, 102, 102, 102 /)
          ngz_list = (/999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999 /)
          kk_list = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
        case('V10')
          ngx_list = (/102, 90, 62, 44, 36, 38, 40, 38, 40, 50, 48, 52, 46, 44, 42, 36, 102, 102, 102 /)
          ngz_list = (/999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999 /)
          kk_list = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
        case('TSK')
          ngx_list = (/12, 2, 0, 24, 42, 50, 18, 42, 50, 64, 68, 102, 2, 90, 76, 16, 56, 62, 102 /)
          ngz_list = (/999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999 /)
          kk_list = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
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

