# this script is customized for chrysalis

# ----- if you need to clone the model -----
# the current working branches for mpas production runs
branch=evasinha/iac-gcam-james-paper   # e3sm branch to check out
giac_branch=master  # iac submodule branch
gcam_branch=aldivi/gcam_gcam6_merge_co2     # gcam submodule branch

# set these two variables to your desired paths
topdir=$HOME/e3sm  # home directory or other build area top
clone=E3SM   # Subdirectory name to clone e3sm into (the default is E3SM if you don't give the $clone argument below to git)

# If clone dir doesn't exist, get the clone
if [[ ! -d $topdir/$clone ]]; then
    # Clone and check out branch
    git clone git@github.com:E3SM-Project/E3SM.git $clone
    cd $clone
    git checkout $branch
    git submodule update --init --recursive

    # Check out submodule development branches
    # these are set correctly in the branch above
    #cd components/gcam/src
    #git checkout $giac_branch
    #cd iac/gcam
    #git checkout $gcam_branch
fi

# ----- Directory paths -----
export USERID=$USER
export E3SM_DIR=$HOME/$clone
export E3SM_SCRIPTS_DIR=${E3SM_DIR}/cime/scripts
export E3SM_OUTPUT_DIR=/lcrc/group/e3sm/${USERID}/e3sm_scratch
export SURFACE_DATA_DIR=/lcrc/group/e3sm/data/inputdata/lnd/clm2/surfdata_map

# ------ Create new case -----
# you can customize your CASE_NAME by adding to string below
export RES=ne30pg2_f09_oEC60to30v3 # non-default grids are: atm:ne30np4.pg2  lnd:0.9x1.25  ocnice:oEC60to30v3  rof:null  glc:null  wav:null   mask is: oEC60to30v3
#export COMPSET=SSP245_EAM%CMIP6_ELM%CNPRDCTCBC_MPASSI%PRES_DOCN%DOM_SROF_SGLC_SWAV_GCAM_BGC%LNDATM
export COMPSET=SSP245_ZATM_BGC
export CASEID=$(date '+%Y%m%d')
export CASE_NAME=${CASEID}_${COMPSET}_${RES}
export CASE_ARCHIVE_DIR=${E3SM_OUTPUT_DIR}/${CASE_NAME}/archive

export E3SM_CASE_DIR=${E3SM_SCRIPTS_DIR}/${CASE_NAME}

export iesm_dyn_source=landuse.timeseries_0.9x1.25_HIST_simyr2015_c201021.nc

# this is gcam scenario and machine dependent
gcam_rdir=/lcrc/group/e3sm/data/inputdata/iac/giac/gcam/gcam_6_0/restart/ssp2rcp45/chrysalis/

ldir=/lcrc/group/e3sm/data/inputdata/lnd/clm2/rawdata/LUT_input_files_current
gcam_idir=/lcrc/group/e3sm/data/inputdata/iac/giac/gcam/gcam_6_0
glm_idir=/lcrc/group/e3sm/data/inputdata/iac/giac/glm
glm2iacdir=/lcrc/group/e3sm/data/inputdata/iac/giac/glm2iac

# this used only if read_scalars is .true. below in user_nl_gcam
# this points to the location of scalar values calculated from a control run to be used for one-way impact simulations
scalardir=/lcrc/group/e3sm/ac.dalei.hao/E3SM_GCAM_simulations/20241226_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3_CONTROL/run

# Scratch directory and subdirectories for this build.  
scratchdir=${E3SM_OUTPUT_DIR}/${CASE_NAME}
rundir=$scratchdir/run

# Delete old case and run directory
#rm -rf ${E3SM_CASE_DIR}
#rm -rf ${E3SM_OUTPUT_DIR}/${CASE_NAME}

cd ${E3SM_SCRIPTS_DIR}

./create_newcase \
 --case ${CASE_NAME} \
 --compset ${COMPSET} \
 --res ${RES} \
 --output-root ${E3SM_OUTPUT_DIR} \
 --project e3sm \
 --machine chrysalis

# ----- Modify user namelist fles -----

cd ${E3SM_CASE_DIR}

export finidat_COMPSET_alias=I20TREAMELMCNPRDCTCBCBGC
export finidat_CASEID=20241204
export finidat_case=${finidat_CASEID}_${finidat_COMPSET_alias}_${RES}
export RUN_REFDATE=2015-01-01

export domainpath=/lcrc/group/e3sm/data/inputdata/share/domains
export lnd_domainfile=domain.lnd.0.9x1.25_oEC60to30v3.231108.nc
export atm_domainfile=domain.lnd.ne30pg2_oEC60to30v3.200220.nc

# suplphos = 'ALL' sets supplemental phosphorus as active for all vegetation types

# Setting do_harvest == .false. because the iac takes care of this
# transient pft flag automatically set to false when flanduse.timeseries is not set 

cat >> user_nl_elm << EOF
&elm_inparm
 hist_mfilt = 1, 365, 1
 hist_nhtfrq = 0, -24, 0
 hist_dov2xy = .true., .true., .false.
 hist_fincl2 = 'TBOT', 'TREFMXAV', 'TREFMNAV', 'RAIN', 'SNOW', 'SNOWDP'
 hist_fincl3 = 'GPP', 'ER', 'HR', 'NPP'
 model_year_align_pdep = 2000
 stream_year_first_pdep = 2000
 stream_year_last_pdep = 2000
 stream_fldfilename_ndep = '/lcrc/group/e3sm/data/inputdata/lnd/clm2/ndepdata/fndep_elm_cbgc_exp_simyr1849-2101_1.9x2.5_ssp245_c240903.nc'
 model_year_align_ndep = 2015
 stream_year_first_ndep = 2015
 stream_year_last_ndep = 2100
 stream_fldfilename_popdens = '/lcrc/group/e3sm/data/inputdata/lnd/clm2/firedata/elmforc.ssp2_hdm_0.5x0.5_simyr1850-2101_c20200623.nc'
 model_year_align_popdens = 2015
 stream_year_first_popdens = 2015
 stream_year_last_popdens = 2100
 check_finidat_fsurdat_consistency = .false.
 check_finidat_year_consistency = .false.
 check_dynpft_consistency = .false.
 do_budgets = .false.
 do_harvest = .false.
EOF

# with co2_flag = .true. we also have to add -co_cycle to CAM_CONFIG_OPTS; this is done below via /.xmlchange to add the four co2 flux tracers
cat >> user_nl_eam << EOF
&co2_cycle_nl
 co2_flag = .true.
 co2_readFlux_fuel = .false.
 co2_readFlux_aircraft = .false.
 co2_readFlux_ocn = .true.
 co2flux_ocn_file = '/lcrc/group/acme/public_html/e3sm_support/compset_generation/ssp245/ocn/co2_flux/merged_regridded_data/fgco2_CESM2_SSP245_ne30pg2_2015-2100.nc'
&chem_surfvals_nl
 co2vmr = 0.000001e-6
 scenario_ghg = 'RAMPED'
 bndtvghg = '/lcrc/group/e3sm/ac.dalei.hao/GCAM_E3SM/ghg_concentration_data_SSP245_GCAM/GHG_GCAM_SSP245_Annual_Global_2015-2102_c20240712.nc'
 co2_conserv_error_tol_per_year = 1.e-5
EOF

# some options for user_nl_gcam
# to run historical period, add run_gcam = .false.
# to spinup gcam for restart and baseline files:
#  add gcam_spinup = .true.
# to run one-way impacts use read_scalars = .true. and provide scalar files from a control simulation

# to run a control simulation without two-way human feedbacks set:
#    elm_ehc_carbon_scaling = .false.
#    elm_ehc_agyield_scaling = .false.
cat >> user_nl_gcam << EOF
&gcam_inparm
 ehc_eam_co2_emissions = .true.
 elm_ehc_carbon_scaling = .true.
 elm_ehc_agyield_scaling = .true.
 gcam_config = '/lcrc/group/e3sm/data/inputdata/iac/giac/gcam/gcam_6_0/configuration/configuration_ssp245_in_E3SM.xml'
 base_hr_file = '/lcrc/group/e3sm/data/inputdata/iac/giac/gcam/gcam_6_0/data/base_f09_20241204_ZATM_annAvgMonthly_2010-2014_hr.csv'
 base_npp_file = '/lcrc/group/e3sm/data/inputdata/iac/giac/gcam/gcam_6_0/data/base_f09_20241204_ZATM_annAvgMonthly_2010-2014_npp.csv'
 base_pft_file = '/lcrc/group/e3sm/data/inputdata/iac/giac/gcam/gcam_6_0/data/base_f09_20241204_ZATM_annAvgMonthly_2010-2014_pft_wt.csv'
EOF

# ----- Case setup -----
# IAC runs once a year, so we have to set the NCPL
# options accordingly (we can't give fractional IAC_NCPL).
# want to run 85 years (through 2099), get a little more than 1 year per hour
#    so run in two batches; 45 years, then continue for 40 years
./xmlchange -append CAM_CONFIG_OPTS="-co2_cycle"
./xmlchange SAVE_TIMING=TRUE
./xmlchange RUN_TYPE=hybrid
./xmlchange GET_REFCASE=TRUE
./xmlchange RUN_REFDIR=/lcrc/group/e3sm/ac.eva.sinha/E3SM_GCAM_lnd_init/${finidat_case}
./xmlchange RUN_REFCASE=${finidat_case}
./xmlchange RUN_REFDATE=${RUN_REFDATE}
./xmlchange RUN_STARTDATE=2015-01-01
./xmlchange ATM_DOMAIN_PATH=${domainpath}
./xmlchange LND_DOMAIN_PATH=${domainpath}
./xmlchange ATM_DOMAIN_FILE=${atm_domainfile}
./xmlchange LND_DOMAIN_FILE=${lnd_domainfile}
./xmlchange NCPL_BASE_PERIOD=year
./xmlchange ATM_NCPL=17520
./xmlchange IAC_NCPL=1
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=6
./xmlchange REST_N=5
./xmlchange JOB_QUEUE=slurm
./xmlchange JOB_WALLCLOCK_TIME=08:00:00
./xmlchange MAX_MPITASKS_PER_NODE=64
./xmlchange MAX_TASKS_PER_NODE=64
./xmlchange NTASKS_ATM=3328
./xmlchange NTASKS_CPL=3328
./xmlchange NTASKS_OCN=3328
./xmlchange NTASKS_ICE=3328
./xmlchange NTASKS_LND=3328
./xmlchange SSTICE_DATA_FILENAME=/lcrc/group/acme/public_html/e3sm_support/compset_generation/ssp245/ocn/sst/merged_regridded_data/sst_ice_GFDL-ESM4_ssp245_r2i1p1f1_gr_201501-210012_land_interpolated.nc
./xmlchange SSTICE_YEAR_START=2015
./xmlchange SSTICE_YEAR_END=2100
./xmlchange SSTICE_YEAR_ALIGN=2015
./xmlchange ROOTPE=0
./xmlchange NTHRDS=1

#./xmlchange DOUT_S_ROOT=${CASE_ARCHIVE_DIR}

./case.setup

# GCAM and GLM currently read some configuration and input files from
# the current case/run directories, as we haven't yet put them in a standard
# place or modified code to look for them there.  Thus, we manually
# copy them over for now.
cp -r $gcam_idir/input $scratchdir
cp $gcam_idir/configuration/log_conf.xml $rundir
cp $glm_idir/glm.fut.conf.chrysalis  $rundir/glm.fut.conf

# separate restarts for each machine and gcam scenario
cp $gcam_rdir/restart.* $rundir

cp $ldir/iESM_Ref_CropPast2015_c10142019.nc $rundir/iESM_Init_CropPast.nc
cp $ldir/surfdata_360x720_mcrop2015_c07082020.nc $rundir/surfdata_360x720_mcrop_init.nc
cp $glm2iacdir/$iesm_dyn_source $rundir
cp $glm2iacdir/surfdata_360x720_potveg.nc $rundir
cp $glm2iacdir/mksurf_landuse_iESM_720x360.nc $rundir
cp $glm2iacdir/iac_in_chrysalis $rundir/iac_in

#### remember to re-do these three copies before each run
cp $rundir/iESM_Init_CropPast.nc $rundir/iESM_Dyn_CropPast.nc
cp $rundir/surfdata_360x720_mcrop_init.nc  $rundir/surfdata_360x720_mcrop_dyn.nc
cp $rundir/$iesm_dyn_source  $rundir/surfdata_iESM_dyn.nc

# this copies the impact scalars to the run directory
# use this only if read_scalars = .true. in user_nl_gcam
#cp $scalardir/scalars* $rundir

# ----- Case build -----
./case.build

# ----- Run model -----
./case.submit
