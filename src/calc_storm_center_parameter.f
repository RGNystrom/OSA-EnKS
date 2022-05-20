  Program calc_storm_center_parameter

!----------------------------------------------------------------
! This program is used to caluclate the average parameter values
! within distance of TC cented

! Input files:  atcf_file:  tcvitals.dat

!-----work flow 
! 1, readin observed center from TCVitials file
! 2, find the i,j by mapping the lat, lon of TCVitals
! 3, Calculate mean:
!    while: r<=100

!----------------------------------------------------------------
  use netcdf
  implicit none

  integer              :: iost,ilat,ilon,cdfid,rcode,ix,jx,kx,i,j,k,tcind1_i,tcind1_j
  integer              :: varnum, ii,jj,kk, m, fid
  character(len= 1)    :: ilatc, ilonc
  real                 :: dx, tclat, tclon, dis0, dis1, r, a, Rmin, Rmax

  real, allocatable, dimension(:,:) :: xlat, xlong

  character (len=10), allocatable, dimension(:) :: varname
  character (len=10)                            :: var

  real, allocatable, dimension(:,:,:)   :: dat

  real                 :: p_tot, p_cnt, p_mean

  Rmin = 0
  Rmax = 100
!----------------------------------------------------------------
! Get storm center from "tcvitals.dat"
  open(10, file="tcvitals.dat", status='old', form = 'formatted', iostat = iost )
  if ( iost .ne. 0 ) then
     stop 'NO a-deck file found!'
  endif
  read(10,'(32x, i4, a1, i5, a1)')ilat,ilatc,ilon,ilonc
  tclat=ilat/10.
  if (ilatc.eq.'S')tclat=-ilat/10.
  tclon=ilon/10.
  if (ilonc.eq.'W')tclon=-ilon/10.
  !write(*,*)'TC center:', tclat, tclon

!----------------------------------------------------------------
! Map the TC center to the WRF domain
  rcode = nf_open ('wrfinput', NF_NOWRITE, cdfid)
  rcode = nf_get_att_int (cdfid, nf_global, 'WEST-EAST_PATCH_END_UNSTAG', ix)
  rcode = nf_get_att_int (cdfid, nf_global, 'SOUTH-NORTH_PATCH_END_UNSTAG', jx)
  rcode = nf_get_att_int (cdfid, nf_global, 'BOTTOM-TOP_PATCH_END_UNSTAG', kx)
  rcode = nf_get_att_real(cdfid, nf_global, 'DX', dx)
  dx=dx/1000.
  rcode = nf_close(cdfid)
  !write(*,*)'ix,jx,kx=',ix,jx,kx
  allocate( xlat (ix, jx) ) 
  allocate( xlong (ix, jx) ) 
  call get_variable2d( 'wrfinput   ', 'XLAT      ', ix, jx, 1, xlat )
  call get_variable2d( 'wrfinput   ', 'XLONG     ', ix, jx, 1, xlong )
  if ( tclat.gt.xlat(1,1) .and. tclat.lt.xlat(1,jx) .and.      &
       tclon.gt.xlong(1,1) .and. tclon.lt.xlong(ix,1) ) then
     dis0=999999. 
     do j = 1, jx 
     do i = 1, ix
        dis1 = (tclat-xlat(i,j))**2 + (tclon-xlong(i,j))**2
        if ( dis1 <= dis0 ) then
           dis0 = dis1
           tcind1_i = i
           tcind1_j = j
        end if
     end do
     end do
  else
     tcind1_i = -9999
     tcind1_j = -9999
  end if
 !write(*,*)'TC center in WRFout:', tcind1_i, tcind1_j

!----------------------------------------------------------------
! variable Mean to calculate
  varnum = 3 
  allocate( varname (varnum) )
  varname = (/'SFC_ALPHA ', 'SFC_VC    ', 'SFC_BETA  '/)

!----------------------------------------------------------------
  do_wrf_var  : do m = 1, varnum

!.... get dimensions
      var = varname(m)
      call wrf_var_dimension ( var, ix, jx, kx, ii, jj, kk )
      allocate( dat   ( ii, jj, kk ) )
      !write(*,*)var, ii, jj, kk

!....... get data
      if ( kk > 1 ) then
         call get_variable3d('wrfinput',     var, ii, jj, kk, 1, dat)
      else if ( kk == 1 ) then
         call get_variable2d('wrfinput',     var, ii, jj, 1, dat)
      endif

!....... calc mean for dat
      p_cnt = 0.
      p_tot = 0.
      p_mean = 0.
      do j = 1, jj
      do i = 1, ii
         !r = sqrt((i-tcind1_i)**2+(j-tcind1_j)**2)*dx
         r = (((i-tcind1_i)**2+(j-tcind1_j)**2)**0.5)*dx
         if ( r < Rmax ) then
            !write(*,*)i,j
            do k = 1,kk
                p_cnt = p_cnt + 1.
                p_tot = p_tot + dat(i,j,k)
                !write(*,*)p_tot
            end do
         endif
      end do
      end do

      !write(*,*)'dat:',dat(i,j,k)
      !write(*,*)'cnt:',p_cnt
      !write(*,*)'tot:',p_tot
      p_mean = p_tot/p_cnt

      write(*,*)p_mean

!....... put back to wrfinput
      call open_file('wrfinput   ', nf_write, fid)
      if ( kk > 1 ) then
         call write_variable3d(fid, var, ii, jj, kk, 1, dat)
      else if ( kk == 1 ) then
         call write_variable2d(fid, var, ii, jj, 1, dat)
      endif 
      call close_file( fid )

      deallocate( dat  )
  end do do_wrf_var

  !write(*,*)'!!! Successful completion of calc_storm_center_parameter.exe!!!'


end Program calc_storm_center_parameter

!==============================================================================
subroutine wrf_var_dimension ( var, ix, jx, kx, ii, jj, kk )

   character(len=10), intent(in)   :: var
   integer, intent(in)             :: ix, jx, kx
   integer, intent(out)            :: ii, jj, kk

   ii = ix
   jj = jx
   kk = kx
   if      ( var == 'U         ' ) then
      ii = ix + 1
   else if ( var == 'V         ' ) then
      jj = jx + 1
   else if ( var == 'W         ' .or. var == 'PH        ' .or. var == 'PHB       ' ) then
      kk = kx + 1
   else if ( var == 'MU        ' .or. var == 'MUB       ' .or. var == 'Q2        '  &
        .or. var == 'T2        ' .or. var == 'TH2       ' .or. var == 'PSFC      '  &
        .or. var == 'SST       ' .or. var == 'TSK       ' .or. var == 'XICE      '  &
        .or. var == 'SFROFF    ' .or. var == 'UDROFF    ' .or. var == 'IVGTYP    '  &
        .or. var == 'ISLTYP    ' .or. var == 'VEGFRA    ' .or. var == 'GRDFLX    '  &
        .or. var == 'SNOW      ' .or. var == 'SNOWH     ' .or. var == 'CANWAT    '  &
        .or. var == 'MAPFAC_M  ' .or. var == 'F         '  &
        .or. var == 'E         ' .or. var == 'SINALPHA  ' .or. var == 'COSALPHA  '  &
        .or. var == 'HGT       ' .or. var == 'TSK       ' .or. var == 'RAINC     '  &
        .or. var == 'RAINNC    ' .or. var == 'SWDOWN    ' .or. var == 'GLW       '  &
        .or. var == 'XLAT      ' .or. var == 'XLONG     ' .or. var == 'TMN       '  &
        .or. var == 'XLAND     ' .or. var == 'PBLH      ' .or. var == 'HFX       '  &
        .or. var == 'QFX       ' .or. var == 'LH        ' .or. var == 'SNOWC     '  &
        .or. var == 'SR        ' .or. var == 'POTEVP    ' .or. var == 'U10       '  &
        .or. var == 'SFC_ALPHA ' .or. var == 'SFC_VC    ' .or. var == 'SFC_BETA  '  &
        .or. var == 'V10       '  ) then
      kk = 1
   else if ( var == 'MAPFAC_U  ' ) then
      kk = 1
      ii = ix + 1
   else if ( var == 'MAPFAC_V  ' ) then
      kk = 1
      jj = jx + 1
   else if ( var == 'FNM       ' .or. var == 'FNP       '  &
        .or. var == 'RDNW      ' .or. var == 'RDN       '  &
        .or. var == 'DNW       ' .or. var == 'DN        '  &
        .or. var == 'ZNU       '                          ) then
      ii = 1
      jj = 1
   else if ( var == 'ZNW       '                          ) then
      ii = 1
      jj = 1
      kk = kx + 1
   endif

end subroutine wrf_var_dimension

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
