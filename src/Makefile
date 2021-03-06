FC=mpif90
FCFLAG=-f90=ifort -i4 -w -ftz -align all -fno-alias -fp-model precise -assume realloc_lhs -convert big_endian -traceback

all: enkf.mpi enkf_hydro.mpi replace_perturbationmean_by_initial.exe replace_geo_by_initial.exe replace_xfmean_with_gfs.exe enkf_gfs_hybrid_inflation.exe soairborne_to_3dvar.exe so_to_3dvar.exe ensemble_mean.exe replace_mean.exe replace_environment_by_gfs.exe replace_mean_outside_site.exe replace_environment_by_gfs.exe cal_inflate.mpi calc_storm_center_parameter.exe relocate_vortex_to_HPI.exe 

clean:
	rm -f *.o *.mod *.mpi *.exe

mpi_module.o:
	$(FC) -c $(FCFLAG) -FR mpi_module.f
module_structure.o:
	$(FC) -c $(FCFLAG) -FR module_structure.f
module_netcdf.o:
	$(FC) -c $(FCFLAG) -I${NETCDF}/include module_netcdf.f
module_map_utils.o:
	$(FC) -c $(FCFLAG) module_map_utils.f
module_wrf_tools.o:
	$(FC) -c $(FCFLAG) -FR module_wrf_tools.f
module_radar.o:
	$(FC) -c $(FCFLAG) -FR module_radar.f
cal_roi.o:
	$(FC) -c $(FCFLAG) -FR cal_roi.f
cal_abei.o:
	$(FC) -c $(FCFLAG) -FR cal_abei.f
hurricane_center.o:
	$(FC) -c $(FCFLAG) -FR hurricane_center.f
obs_io.o:
	$(FC) -c $(FCFLAG) -FR obs_io.f
xb.o:
	$(FC) -c $(FCFLAG) -FR -I$(NETCDF)/include -I$(CRTM_DIR)/include xb.f
sub_enkf_util.o:
	$(FC) -c $(FCFLAG) -FR sub_enkf_util.f
enkf.o:
	$(FC) -c $(FCFLAG) -FR enkf.f
enkf_hydro.o:
	$(FC) -c $(FCFLAG) -FR enkf_hydro.f
main.o:
	$(FC) -c $(FCFLAG) -FR main.f
ensemble_mean.o:
	$(FC) -c $(FCFLAG) -FR ensemble_mean.f
replace_mean.o:
	$(FC) -c $(FCFLAG) -FR replace_mean.f
replace_environment_by_gfs.o:
	$(FC) -c $(FCFLAG) -FR replace_environment_by_gfs.f
cal_inflate_bt.o:
	$(FC) -c $(FCFLAG) -FR cal_inflate_bt.f
calc_storm_center_parameter.o:
	$(FC) -c $(FCFLAG) -FR calc_storm_center_parameter.f
relocate_vortex_to_HPI.o:
	$(FC) -c $(FCFLAG) -FR relocate_vortex_to_HPI.f

enkf.mpi: module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o cal_abei.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o enkf.o
	$(FC) -o enkf.mpi $(FCFLAG) -FR -I$(NETCDF)/include  module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o cal_abei.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o enkf.o -L$(CRTM_DIR)/lib -lcrtm -L$(NETCDF)/lib -lnetcdff main.f

enkf_hydro.mpi: module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o cal_abei.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o enkf_hydro.o
	$(FC) -o enkf_hydro.mpi $(FCFLAG) -FR -I$(NETCDF)/include  module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o cal_abei.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o enkf_hydro.o -L$(CRTM_DIR)/lib -lcrtm -L$(NETCDF)/lib -lnetcdff main.f

replace_perturbationmean_by_initial.exe: module_netcdf.o
	$(FC) -o replace_perturbationmean_by_initial.exe $(FCFLAG) -FR -I$(NETCDF)/include module_netcdf.o -L$(NETCDF)/lib -lnetcdff replace_perturbationmean_by_initial.f

replace_geo_by_initial.exe: module_netcdf.o
	$(FC) -o replace_geo_by_initial.exe $(FCFLAG) -FR -I$(NETCDF)/include module_netcdf.o -L$(NETCDF)/lib -lnetcdff replace_geo_by_initial.f

replace_mean_outside_site.exe: module_structure.o module_netcdf.o module_wrf_tools.o module_map_utils.o
	$(FC) -o replace_mean_outside_site.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o module_netcdf.o module_wrf_tools.o module_map_utils.o -L$(NETCDF)/lib -lnetcdff replace_mean_outside_site.f

replace_xfmean_with_gfs.exe: module_structure.o module_netcdf.o module_wrf_tools.o module_map_utils.o
	$(FC) -o replace_xfmean_with_gfs.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o module_netcdf.o module_wrf_tools.o module_map_utils.o -L$(NETCDF)/lib -lnetcdff replace_xfmean_with_gfs.f

enkf_gfs_hybrid_inflation.exe: module_structure.o module_netcdf.o module_wrf_tools.o module_map_utils.o
	$(FC) -o enkf_gfs_hybrid_inflation.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o module_netcdf.o module_wrf_tools.o module_map_utils.o -L$(NETCDF)/lib -lnetcdff enkf_gfs_hybrid_inflation.f

soairborne_to_3dvar.exe: module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o
	$(FC) -o soairborne_to_3dvar.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o -L$(CRTM_DIR)/lib -lcrtm -L$(NETCDF)/lib -lnetcdff soairborne_to_3dvar.f

so_to_3dvar.exe: module_structure.o mpi_module.o module_netcdf.o module_wrf_tools.o module_map_utils.o module_radar.o cal_roi.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o
	$(FC) -o so_to_3dvar.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o module_radar.o cal_roi.o hurricane_center.o obs_io.o xb.o sub_enkf_util.o -L$(CRTM_DIR)/lib -lcrtm -L$(NETCDF)/lib -lnetcdff so_to_3dvar.f

ensemble_mean.exe: module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o sub_enkf_util.o
	$(FC) -o ensemble_mean.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o sub_enkf_util.o -L$(CRTM_DIR)/lib -lcrtm -L$(NETCDF)/lib -lnetcdff ensemble_mean.f

replace_mean.exe: module_netcdf.o mpi_module.o
	$(FC) -o replace_mean.exe $(FCFLAG) -FR -I$(NETCDF)/include module_netcdf.o mpi_module.o -L$(NETCDF)/lib -lnetcdff replace_mean.f

replace_environment_by_gfs.exe: module_netcdf.o mpi_module.o
	$(FC) -o replace_environment_by_gfs.exe $(FCFLAG) -FR -I$(NETCDF)/include module_netcdf.o mpi_module.o -L$(NETCDF)/lib -lnetcdff replace_environment_by_gfs.f

calc_storm_center_parameter.exe: module_netcdf.o mpi_module.o
	$(FC) -o calc_storm_center_parameter.exe $(FCFLAG) -FR -I$(NETCDF)/include module_netcdf.o mpi_module.o -L$(NETCDF)/lib -lnetcdff calc_storm_center_parameter.f

wrf2obs3dvar: module_netcdf.o module_map_utils.o module_structure.o module_wrf_tools.o
	rm -f wrf2obs3dvar.exe
	$(FC) -o wrf2obs3dvar.exe $(FCFLAG) -FR -I$(NETCDF)/include module_netcdf.o module_map_utils.o module_structure.o module_wrf_tools.o -L$(NETCDF)/lib -lnetcdff wrf2obs3dvar.f #time_calculation.f

cal_inflate.mpi: module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o sub_enkf_util.o
	$(FC) -o cal_inflate.mpi $(FCFLAG) -FR -I$(NETCDF)/include  module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o sub_enkf_util.o -L$(CRTM_DIR)/lib -lcrtm -L$(NETCDF)/lib -lnetcdff cal_inflate.f

relocate_vortex_to_HPI.exe: module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o sub_enkf_util.o
	$(FC) -o relocate_vortex_to_HPI.exe $(FCFLAG) -FR -I$(NETCDF)/include module_structure.o mpi_module.o module_netcdf.o module_map_utils.o module_wrf_tools.o sub_enkf_util.o -L$(NETCDF)/lib -lnetcdff relocate_vortex_to_HPI.f
