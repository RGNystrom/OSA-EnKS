  Program relocate_vortex_to_HPI

!----------------------------------------------------------------
! This program is used to replace environment by GFS

! Input files:  atcf_file:  tcvitals.dat
! wrfinput (or output): wrfinput_gfs and wrfinput

!-----work flow 
! 1, readin observed center from TCVitials file
! 2, find the i,j by mapping the lat, lon of TCVitals
! 3, replace the environemt:
!    Xnew = a*wrfinput + (1-a)*wrfinput_gfs
!    while: r<=300, a=1.0, r is the distance from the grid to center
!           r> 600, a=0.0,
!           300<r<=600, a=(r-300)/(600-300) 

!----------------------------------------------------------------
  use netcdf
  use map_utils
  use constants
  implicit none

  integer              :: iost,ilat,ilon,cdfid,rcode,ix,jx,kx,i,j,k,tcind1_i,tcind1_j
  integer              :: varnum, ii,jj,kk, m, fid, ens_ind1_i, ens_ind1_j,i2,j2
  character(len= 1)    :: ilatc, ilonc
  real                 :: dx, tclat, tclon, dis0, dis1, r, a, Rmin, Rmax
  real                 :: ens_lon, ens_lat, Rmin2

  real, allocatable, dimension(:,:) :: xlat, xlong

  character (len=10), allocatable, dimension(:) :: varname
  character (len=10)                            :: var
  character (len=10) :: inputfile
  character (len=80) :: times

  real, allocatable, dimension(:,:,:)   :: gfs     !!! initial field
  real, allocatable, dimension(:,:,:)   :: dat
  real, allocatable, dimension(:,:,:)   :: dat2
  real, allocatable, dimension(:,:,:)   :: pert
  real, allocatable, dimension(:,:,:)   :: pmean

  real, allocatable, dimension(:,:,:)   :: z
  real, allocatable, dimension(:,:,:)   :: t
  real, allocatable, dimension(:,:,:)   :: ph
  real, allocatable, dimension(:,:,:)   :: phb
  real, allocatable, dimension(:,:,:)   :: qc
  real, allocatable, dimension(:,:,:)   :: qr
  real, allocatable, dimension(:,:,:)   :: qv
  real, allocatable, dimension(:,:,:)   :: p
  real, allocatable, dimension(:,:,:)   :: pb
  real, allocatable, dimension(:,:)   :: slp
  real, allocatable, dimension(:,:)   :: slpsmth
  real, allocatable, dimension(:)   :: znw
  real, allocatable, dimension(:)   :: znu
  real, dimension(3)   :: center

  type(proj_info)    :: proj

  Rmin = 300
  Rmax = 600
  Rmin2 = 600
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
  write(*,*)'HPI TC center:', tclat, tclon

  inputfile='wrfinput'
!----------------------------------------------------------------
! Map the TC center to the WRF domain
  rcode = nf_open ('wrfinput', NF_NOWRITE, cdfid)
  rcode = nf_get_att_int (cdfid, nf_global, 'WEST-EAST_PATCH_END_UNSTAG', ix)
  rcode = nf_get_att_int (cdfid, nf_global, 'SOUTH-NORTH_PATCH_END_UNSTAG', jx)
  rcode = nf_get_att_int (cdfid, nf_global, 'BOTTOM-TOP_PATCH_END_UNSTAG', kx)
  rcode = nf_get_att_real(cdfid, nf_global, 'DX', dx)
  dx=dx/1000.
  rcode = nf_close(cdfid)
  write(*,*)'ix,jx,kx',ix,jx,kx
  allocate( xlat (ix, jx) )
  allocate( xlong (ix, jx) )
  call get_variable2d( inputfile, 'XLONG     ', ix, jx, 1, xlong )
  call get_variable2d( inputfile, 'XLAT      ', ix, jx, 1, xlat )
  write(*,*)'Check!'
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
  write(*,*)'HPI TC center in WRFout:', tcind1_i, tcind1_j

!----------------------------------------------------------------
! Get storm center from ens
  ! Get domain info
  !call hurricane_center_assimilation(wrf_file,ix,jx,kx,int(tcind1_i),int(tcind1_j),&
  !                                          center,znu,znw,xlong,xlat,proj)
  allocate( znw (kx+1) )
  call get_variable1d( inputfile, 'ZNW       ', kx+1, 1, znw )
  allocate( znu (kx) )
  call get_variable1d( inputfile, 'ZNU       ', kx, 1, znu )
  allocate( t (ix,jx,kx) )
  call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, t)
  allocate( ph (ix,jx,kx+1) )
  call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
  allocate( phb (ix,jx,kx+1) )
  call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
  allocate( qv (ix,jx,kx) )
  call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
  allocate( qc (ix,jx,kx) )
  call get_variable3d(inputfile, 'QCLOUD    ', ix, jx, kx,   1, qc )
  allocate( qr (ix,jx,kx) )
  call get_variable3d(inputfile, 'QRAIN     ', ix, jx, kx,   1, qr )
  allocate( p (ix,jx,kx) )
  call get_variable3d(inputfile, 'P         ', ix, jx, kx,   1, p  )
  allocate( pb (ix,jx,kx) )
  call get_variable3d(inputfile, 'PB        ', ix, jx, kx,   1, pb )

  ph = (ph + phb)/g
  p  = p + pb

  write(*,*)'Check2!'
  allocate( z (ix,jx,kx) )
!.... interp PH to unstaggered eta level
  do j = 1, jx
  do i = 1, ix
      z(i,j,1:kx) = (ph(i,j,1:kx)*(znw(2:kx+1)-znu(1:kx))+ph(i,j,2:kx+1)*(znu(1:kx)-znw(1:kx)))/(znw(2:kx+1)-znw(1:kx))
      do k=1,kx
      if( p(i,j,k)/100. >= 1200. ) then
          write(*,'(a,3i4,2f8.0)')'P error at:',i,j,k,p(i,j,k)/100.,pb(i,j,k)/100.
      endif
      enddo
  enddo
  enddo
  write(*,*)'Check3!'
!.... calculate slp
   allocate( slp (ix,jx) )
   call compute_seaprs(ix, jx, kx, z, t, p, qv, slp, 0)
  write(*,*)'Check4!'

   slpsmth = slp
   call smooth(slpsmth, ix, jx, 305)

!.... find the minumum slp
   center(3) = 200000.    !pa
   do j = 1, jx
   do i = 1, ix
         if( slpsmth(i,j) < center(3) )then
             center(3) = slpsmth(i,j)
             center(1) = i*1.
             center(2) = j*1.
         endif
   enddo
   enddo
   center(3) = slp(int(center(1)), int(center(2)))
   
  ens_ind1_i=center(1)
  ens_ind1_j=center(2)
  write(*,*)'Ens TC center in WRFout:', ens_ind1_i, ens_ind1_j

!----------------------------------------------------------------
! variables will be replaced
  varnum = 26 
  allocate( varname (varnum) )
  varname = (/'U         ', 'V         ', 'W         ', 'PH        ', 'PHB       ', &
              'T         ', 'MU        ', 'MUB       ', 'P         ', &
              'Q2        ', 'T2        ', 'TH2       ', 'PSFC      ', 'PB        ', &
              'U10       ', 'V10       ', 'QVAPOR    ', 'QCLOUD    ', 'QRAIN     ', &
              'QICE      ', 'QSNOW     ', 'QGRAUP    ', 'TSK       ', 'HGT       ', &
              'TMN       ', 'SST       '/)

!----------------------------------------------------------------
  do_wrf_var  : do m = 1, varnum

!.... get dimensions
      var = varname(m)
      write(*,*)'Check5!'
      write(*,*)'ix,jx,kx',ix,jx,kx
      call wrf_var_dimension ( inputfile, var, ix, jx, kx, ii, jj, kk )
      write(*,*)'Check6!'
      allocate( gfs   ( ii, jj, kk ) )
      allocate( dat   ( ii, jj, kk ) )
      allocate( dat2   ( ii, jj, kk ) )
      !allocate( pert   ( ii, jj, kk ) )
      !allocate( pmean   ( ii, jj, kk ) )
      write(*,*)var, ii, jj, kk

!....... get data
      if ( kk > 1 ) then
         call get_variable3d('wrfinput_gfs', var, ii, jj, kk, 1, gfs)
         call get_variable3d('wrfinput',     var, ii, jj, kk, 1, dat)
         call get_variable3d('wrfinput',     var, ii, jj, kk, 1, dat2)
         !call get_variable3d('wrfinput_pmean',     var, ii, jj, kk, 1, pmean)
         !call get_variable3d('wrfinput_pert',     var, ii, jj, kk, 1, pert)
      else if ( kk == 1 ) then
         call get_variable2d('wrfinput_gfs', var, ii, jj, 1, gfs)
         call get_variable2d('wrfinput',     var, ii, jj, 1, dat)
         call get_variable2d('wrfinput',     var, ii, jj, 1, dat2)
         !call get_variable2d('wrfinput_pmean',     var, ii, jj, 1, pmean)
         !call get_variable2d('wrfinput_pert',     var, ii, jj, 1, pert)
      endif

!....... perform relocation
     do j = 1, jj
      do i = 1, ii
         !r = sqrt((i-tcind1_i)**2+(j-tcind1_j)**2)*dx
         r = (((i-tcind1_i)**2+(j-tcind1_j)**2)**0.5)*dx
         if ( r > Rmin2 ) then
            i2 = i
            j2 = j
         else if ( r <= Rmin2 ) then
            i2 = tcind1_i+(i-ens_ind1_i)
            j2 = tcind1_j+(j-ens_ind1_j)
         endif
         do k = 1, kk
            if (i2 <= ii .AND. i2 >= 1 .AND. j2 <= jj .AND. j2 >= 1 ) then
               dat2(i2,j2,k)=dat(i,j,k)
            endif
         end do
      end do
      end do

!....... combine gfs and work to dat
      do j = 1, jj
      do i = 1, ii
         !r = sqrt((i-tcind1_i)**2+(j-tcind1_j)**2)*dx
         r = (((i-tcind1_i)**2+(j-tcind1_j)**2)**0.5)*dx
         if ( r > Rmax ) then
            a = 0.0
         else if ( r < Rmin ) then
            a = 1.0
         else if ( r <= Rmax .and. r >= Rmin )then
            a = (Rmax-r)/(Rmax-Rmin)
         endif
         do k = 1, kk
            dat(i,j,k)=dat2(i,j,k)*a+gfs(i,j,k)*(1.-a)
            !dat(i,j,k)=dat(i,j,k)*a+gfs(i,j,k)*(1.-a)+ (1.-a)*(pert(i,j,k)-pmean(i,j,k)) 
         end do
      end do
      end do

!....... put back to wrfinput
      call open_file('wrfinput   ', nf_write, fid)
      if ( kk > 1 ) then
         call write_variable3d(fid, var, ii, jj, kk, 1, dat)
      else if ( kk == 1 ) then
         call write_variable2d(fid, var, ii, jj, 1, dat)
      endif 
      call close_file( fid )

      deallocate( gfs  )
      deallocate( dat  )
      deallocate( dat2  )
      !deallocate( pert  )
      !deallocate( pmean  )
  end do do_wrf_var

  write(*,*)'!!! Successful completion of relocate_vortex_to_HPI.exe!!!'


end Program relocate_vortex_to_HPI

!==============================================================================
!subroutine wrf_var_dimension ( var, ix, jx, kx, ii, jj, kk )
!
!   character(len=10), intent(in)   :: var
!   integer, intent(in)             :: ix, jx, kx
!   integer, intent(out)            :: ii, jj, kk
!
!   ii = ix
!   jj = jx
!   kk = kx
!   if      ( var == 'U         ' ) then
!      ii = ix + 1
!   else if ( var == 'V         ' ) then
!      jj = jx + 1
!   else if ( var == 'W         ' .or. var == 'PH        ' .or. var == 'PHB       ' ) then
!      kk = kx + 1
!   else if ( var == 'MU        ' .or. var == 'MUB       ' .or. var == 'Q2        '  &
!        .or. var == 'T2        ' .or. var == 'TH2       ' .or. var == 'PSFC      '  &
!        .or. var == 'SST       ' .or. var == 'TSK       ' .or. var == 'XICE      '  &
!        .or. var == 'SFROFF    ' .or. var == 'UDROFF    ' .or. var == 'IVGTYP    '  &
!        .or. var == 'ISLTYP    ' .or. var == 'VEGFRA    ' .or. var == 'GRDFLX    '  &
!        .or. var == 'SNOW      ' .or. var == 'SNOWH     ' .or. var == 'CANWAT    '  &
!        .or. var == 'MAPFAC_M  ' .or. var == 'F         '  &
!        .or. var == 'E         ' .or. var == 'SINALPHA  ' .or. var == 'COSALPHA  '  &
!        .or. var == 'HGT       ' .or. var == 'TSK       ' .or. var == 'RAINC     '  &
!        .or. var == 'RAINNC    ' .or. var == 'SWDOWN    ' .or. var == 'GLW       '  &
!        .or. var == 'XLAT      ' .or. var == 'XLONG     ' .or. var == 'TMN       '  &
!        .or. var == 'XLAND     ' .or. var == 'PBLH      ' .or. var == 'HFX       '  &
!        .or. var == 'QFX       ' .or. var == 'LH        ' .or. var == 'SNOWC     '  &
!        .or. var == 'SR        ' .or. var == 'POTEVP    ' .or. var == 'U10       '  &
!        .or. var == 'V10       '  ) then
!      kk = 1
!   else if ( var == 'MAPFAC_U  ' ) then
!      kk = 1
!      ii = ix + 1
!   else if ( var == 'MAPFAC_V  ' ) then
!      kk = 1
!      jj = jx + 1
!   else if ( var == 'FNM       ' .or. var == 'FNP       '  &
!        .or. var == 'RDNW      ' .or. var == 'RDN       '  &
!        .or. var == 'DNW       ' .or. var == 'DN        '  &
!        .or. var == 'ZNU       '                          ) then
!      ii = 1
!      jj = 1
!   else if ( var == 'ZNW       '                          ) then
!      ii = 1
!      jj = 1
!      kk = kx + 1
!   endif
!
!end subroutine wrf_var_dimension
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
