###########
#
# read and plot global cam carbon fluxes and stocks
# standard cam h0 history files for two cases
# also plot the differences between cases
#  these are monthly files
#  read resolution from files - both cases need to have same resolution and extent
#   time = 1
#
# input values are monthly averages
#  fluxes: kgC/m^2/s
#  stocks: kgC/m^2
#  concentration, volume mixing ratio: mol/molair (store and output in ppmv)
#
# output values are:
#  cumulative per month or year for carbon fluxes (Pg/mm or Pg/yy) (just the carbon)
#  monthly or annual average for carbon stocks (Pg) (just the carbon, also for the input vmr variables)
#  ppmv for concentration values
#
# can store only ~50 monthly grids, so store only one at a time
#
# also output a csv file of the plotted annual data
#
# cam landfrac data is the same as clm landfrac data
# do not need the land area of each grid cell
#  area(lat, lon) area of whole grid cell, km^2
#  landfrac(lat, lon) fraction of grid cell that is land
#  landmask(lat, lon) 1=land, 0=ocean
#  missing/fill value = 1.0E36f
#

# this takes ~20 seconds per year on edison login node

cat("started plot_cam_carbon.r", date(), "\n\n")

library(ncdf4)

##### in file

#inname = "chronref_ft_bcase"
inname = "ftrans_maxfor"

#base_name = "b.e11.B20TRBPRP.f09_g16.iESM_chronref_fulltransient.001.cam.h0."
#base_name = "b.e11.B20TRBPRP.f09_g16.iESM_chronref_ftrans_maxfor.001.cam.h0."
base_name = "b.e11.B20TRBPRP.f09_g16.iESM_chronref_ftrans_maxfor.001.cam.h0."

#casearch1 = "/glade/p/cesm/sdwg_dev/iESM_histlulcc/archive/b.e11.B20TRBPRP.f09_g16.iESM_chronref_fulltransient.001/"

# the maxfor case is split into two directories, but the files names within are consistent with each other
# the non-cori directory includes 1850-1914, and the cori directory includes 1950-2004
# use MAXFOR = TRUE for the maxfor case, and MAXFOR = FALSE for any other cases
MAXFOR = TRUE
casearch1 = "/global/project/projectdirs/m1204/shix/b.e11.B20TRBPRP.f09_g16.iESM_chronref_ftrans_maxfor.001/"
casearch2 = "/global/project/projectdirs/m1204/shix/b.e11.B20TRBPRP.f09_g16.iESM_chronref_ftrans_maxfor.001_Cori/"
last_year1 = 1914

datadir1 = paste0(casearch1, "atm/hist/")
datadir2 = paste0(casearch2, "atm/hist/")

#datadir1 = "/Users/adivi/projects/iesm/diagnostics/hist_lulcc/chronref_fulltransient4_bcase/"

#outdir = "/glade/u/home/avdivitt/iesm/cesm1_1_iesm27/scripts/b.e11.B20TRBPRP.f09_g16.iESM_chronref_fulltransient.001/"
outdir = "/global/u1/a/avdivitt/iesm/iesm_hist_ctrl2/"
#outdir = "/Users/adivi/projects/iesm/diagnostics/hist_lulcc/chronref_fulltransient4_bcase/"

ntag = ".nc"

start_year = 1850
end_year = 2004

##### reference file

#refname = "chronref_ft_bcase"
#refname = "iesm_hist_ctrl2"
refname = "ftrans_minfor"

#ref_name = "b.e11.B20TRBPRP.f09_g16.iESM_chronref_fulltransient.001.cam.h0."
#ref_name = "b.e11.B20TRBPRP.f09_g16.iESM_exp12_ctrl.002.cam.h0."
ref_name = "b.e11.B20TRBPRP.f09_g16.iESM_chronref_ftrans_minfor.001.cam.h0."

#refarch = "/glade/p/cesm/sdwg_dev/iESM_histlulcc/archive/b.e11.B20TRBPRP.f09_g16.iESM_chronref_fulltransient.001/"
#refarch = "/global/project/projectdirs/m1204/shix/b.e11.B20TRBPRP.f09_g16.iESM_exp12_ctrl.002/"
refarch = "/scratch2/scratchdirs/avdivitt/b.e11.B20TRBPRP.f09_g16.iESM_chronref_ftrans_minfor.001/"
refdir = paste0(refarch, "atm/hist/")
#refdir = "/Users/adivi/projects/iesm/diagnostics/hist_lulcc/chronref_fulltransient4_bcase/"

pdfout = paste(outdir, base_name, start_year, "-", end_year, "_vs_", refname, "_c.pdf", sep = "")
csvout = paste(outdir, base_name, start_year, "-", end_year, "_vs_", refname, "_c.csv", sep = "")

# this script accomodates only whole years

# variables

# surface co2 fluxes; lon,lat,time; kg/m2/s:
# SFCO2, SFCO2_FFF, SFCO2_LND, SFCO2_OCN

# co2 stocks; lon,lat,time; kg/m2:
# TMCO2, TMCO2_FFF, TMCO2_LND, TMCO2_OCN
# not sure what the fff, lnd, and ocn represent because they do not add up to the total co2

# burden stocks; lon,lat,time; kg/m2:
# BURDENBC (black carbon)

# observed vmr variables; time; mol/molair:
# co2vmr, ch4vmr

# the variables have to be grouped
start_month = 1
end_month = 12

# co2 fluxes, lon,lat,time; kg/m2/s
var_names_flux = c("SFCO2", "SFCO2_FFF", "SFCO2_LND", "SFCO2_OCN")
num_flux_vars = 4
ind_flux_vars = c(1:num_flux_vars)

# column stocks, lon,lat,time; kg/m2
var_names_col = c("TMCO2", "TMCO2_FFF", "TMCO2_LND", "TMCO2_OCN", "BURDENBC")
num_col_vars = 5
ind_col_vars = c((num_flux_vars+1):(num_flux_vars+num_col_vars))
tmco2_ind = num_flux_vars + 1
bc_ind = num_flux_vars + 5

# input co2 and methane volume mixing ratios; global avg; mol/molair (store and output in ppmv)
var_names_vmr = c("co2vmr", "ch4vmr")
num_vmr_vars = 2
ind_vmr_vars = c((num_flux_vars+num_col_vars+1):(num_flux_vars+num_col_vars+num_vmr_vars))
co2vmr_ind = num_flux_vars + num_col_vars + 1
ch4vmr_ind = num_flux_vars + num_col_vars + 2

# the number of variables read in
num_in_vars = num_flux_vars + num_col_vars + num_vmr_vars

# extra C output vars:
var_names_extra_c = c("co2vmr_PgC", "ch4vmr_PgC")
num_extra_c_vars = 2
ind_extra_c_vars = c((num_flux_vars+num_col_vars+num_vmr_vars+1):(num_flux_vars+num_col_vars+num_vmr_vars+num_extra_c_vars))
co2vmr_pgc_ind = num_flux_vars + num_col_vars + num_vmr_vars + 1
ch4vmr_pgc_ind = num_flux_vars + num_col_vars + num_vmr_vars + 2

# extra ppmv output vars:
var_names_extra_ppmv = c("TMCO2_ppmv")
num_extra_ppmv_vars = 1
ind_extra_ppmv_vars = c((num_flux_vars+num_col_vars+num_vmr_vars+num_extra_c_vars+1):(num_flux_vars+num_col_vars+num_vmr_vars+num_extra_c_vars+num_extra_ppmv_vars))
tmco2_ppmv_ind = num_flux_vars + num_col_vars + num_vmr_vars + num_extra_c_vars + 1

num_vars = num_flux_vars + num_col_vars + num_vmr_vars + num_extra_c_vars + num_extra_ppmv_vars
var_names = c(var_names_flux, var_names_col, var_names_vmr, var_names_extra_c, var_names_extra_ppmv)

num_plot_vars = num_vars
ind_plot_vars = c(ind_flux_vars, ind_col_vars, ind_vmr_vars, ind_extra_c_vars, ind_extra_ppmv_vars)
plot_var_names = c(var_names_flux, var_names_col, "co2vmr_ppmv", "ch4vmr_ppmv", var_names_extra_c, var_names_extra_ppmv)

# most are stock outputs
ylabs = c(rep("PgC flux per",num_flux_vars), rep("PgC avg for",num_col_vars), rep("ppmv avg for",num_vmr_vars), rep("PgC avg for",num_extra_c_vars), rep("ppmv avg for",num_extra_ppmv_vars))

num_years = end_year - start_year + 1
tot_months = num_years * 12

# days per month
dayspermonth = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# days per year
daysperyear = sum(dayspermonth)

# seconds per month
secpermonth = 86400 * dayspermonth

# convert grams to Pg
g2Pg = 1E-15

# convert kilograms to Pg
kg2Pg = 1E-12

# convert from km^2 to m^2
kmsq2msq = 1E6

# convert from Kelvin to Celsius
K2C = -273.15

# convert from vmr to ppmv
vmr2ppmv = 1E6

# surface area of earth in m^2
# based on authalic radius, which is spherical radius that gives same surface area as the reference ellipsoid, which is wgs84
surfarea = 4 * (4.0 * atan(1.0)) * 6371007 * 6371007

# total avg mass of atmosphere (kg)
mass_atm_avg = 5.1480E18

# molecular mass of carbon (g/mol)
mm_c = 12.01

# molecular mass of CO2 (g/mol)
mm_co2 = 44.01

# molecualar mass of CH4 (g/mol)
mm_ch4 = 16.04

# molar mass of atmosphere avg (g/mol) (dry)
mm_atm = 28.97

# get some necessary info from the first file 
fname = paste(refdir, ref_name, start_year, "-01", ntag, sep = "")	
nc_id = nc_open(fname)
num_lon = nc_id$dim$lon$len
num_lat = nc_id$dim$lat$len
num_time = nc_id$dim$time$len
landfrac_cam = ncvar_get(nc_id, varid = "LANDFRAC", start = c(1,1,1), count = c(num_lon, num_lat, num_time))
gw = ncvar_get(nc_id, varid = "gw", start = c(1), count = c(num_lat))
gw_sum = sum(gw)
gw_glb_sum = gw_sum * num_lon
nc_close(nc_id)

# arrays for storing 1 month of time series of (time, lat, lon) variables
grid_vardata = array(dim=c(num_in_vars, num_lon, num_lat, 1))
grid_refdata = array(dim=c(num_in_vars, num_lon, num_lat, 1))
# arrays for storing monthly time series of global sums of (time, lat, lon) variables
glb_vardata = array(dim = c(num_vars, tot_months))
glb_refdata = array(dim = c(num_vars, tot_months))
glb_diffdata = array(dim = c(num_vars, tot_months))
glb_vardata[,] = 0.0
glb_refdata[,] = 0.0

# arrays for storing annual time series of global sums
glbann_vardata = array(dim = c(num_vars, num_years))
glbann_refdata = array(dim = c(num_vars, num_years))
glbann_diffdata = array(dim = c(num_vars, num_years))
glbann_vardata[,] = 0.0
glbann_refdata[,] = 0.0


# loop over year
yind = 0
for (y in start_year:end_year) {
	yind = yind + 1
	cat("starting year", y, date(), "\n")
	# loop over month
	for (mind in start_month:end_month) {
		if(mind > 9) {
			mtag = mind
		}else {
			mtag = paste(0, mind, sep = "")
		}
		# need to read the maxfor files from different directories, based on year
		if(MAXFOR){
			if(y <= last_year1){
				fname = paste(datadir1, base_name, y, "-", mtag, ntag, sep = "")
			} else {
				fname = paste(datadir2, base_name, y, "-", mtag, ntag, sep = "")
			}
		} else {
			fname = paste(datadir1, base_name, y, "-", mtag, ntag, sep = "")
		}
		rname = paste(refdir, ref_name, y, "-", mtag, ntag, sep = "")		
		nc_id = nc_open(fname)
		nc_id_ref = nc_open(rname)
		totmind = (y - start_year) * 12 + mind
		
		# loop over flux variables, km/m2/s
		for (vind in ind_flux_vars) {
			grid_vardata[vind,,,1] = ncvar_get(nc_id, varid = var_names[vind], start = c(1,1,1), count = c(num_lon, num_lat, num_time))
			grid_refdata[vind,,,1] = ncvar_get(nc_id_ref, varid = var_names[vind], start = c(1,1,1), count = c(num_lon, num_lat, num_time))
			
			# multiply by gauss weights and sum the globe to prepare for calculating the global values
			for (lonind in 1:num_lon) {
				glb_vardata[vind, totmind] = glb_vardata[vind, totmind] + sum(grid_vardata[vind,lonind,,1] * gw, na.rm = TRUE)
				glb_refdata[vind, totmind] = glb_refdata[vind, totmind] + sum(grid_refdata[vind,lonind,,1] * gw, na.rm = TRUE)
			}
			
			# convert to Pg C per month
			glb_vardata[vind, totmind] = glb_vardata[vind, totmind] / gw_glb_sum * surfarea * secpermonth[mind] * kg2Pg * mm_c / mm_co2
			glb_refdata[vind, totmind] = glb_refdata[vind, totmind] / gw_glb_sum * surfarea * secpermonth[mind] * kg2Pg * mm_c / mm_co2
			
			# calculate monthly differences
			glb_diffdata[vind, totmind] = glb_vardata[vind, totmind] - glb_refdata[vind, totmind]

			# sum to get total annual fluxes
			glbann_vardata[vind, yind] = glbann_vardata[vind, yind] + glb_vardata[vind, totmind] 
			glbann_refdata[vind, yind] = glbann_refdata[vind, yind] + glb_refdata[vind, totmind]
			
			# get the annual difference at the end of each year
			if (mind == end_month) { 
				glbann_diffdata[vind, yind] = glbann_vardata[vind, yind] - glbann_refdata[vind, yind]
			}
			
		} # end loop over flux vars
		
		# loop over column stock variables, kg/m2
		for (vind in ind_col_vars) {
			grid_vardata[vind,,,1] = ncvar_get(nc_id, varid = var_names[vind], start = c(1,1,1), count = c(num_lon, num_lat, num_time))
			grid_refdata[vind,,,1] = ncvar_get(nc_id_ref, varid = var_names[vind], start = c(1,1,1), count = c(num_lon, num_lat, num_time))
			
			# multiply by gauss weights and sum the globe to prepare for calculating the global values
			for (lonind in 1:num_lon) {
				glb_vardata[vind, totmind] = glb_vardata[vind, totmind] + sum(grid_vardata[vind,lonind,,1] * gw, na.rm = TRUE)
				glb_refdata[vind, totmind] = glb_refdata[vind, totmind] + sum(grid_refdata[vind,lonind,,1] * gw, na.rm = TRUE)
			}
			
			# convert to Pg C
			if (vind == bc_ind) { # input is carbon mass/area
				glb_vardata[vind, totmind] = glb_vardata[vind, totmind] / gw_glb_sum * surfarea * kg2Pg
				glb_refdata[vind, totmind] = glb_refdata[vind, totmind] / gw_glb_sum * surfarea * kg2Pg
			} else { # input is CO2 mass/area
				glb_vardata[vind, totmind] = glb_vardata[vind, totmind] / gw_glb_sum * surfarea * kg2Pg * mm_c / mm_co2
				glb_refdata[vind, totmind] = glb_refdata[vind, totmind] / gw_glb_sum * surfarea * kg2Pg * mm_c / mm_co2
				
				# also convert the total co2 column stock to ppmv
				if (vind == tmco2_ind) {
					glb_vardata[tmco2_ppmv_ind, totmind] = glb_vardata[vind, totmind] * mm_atm / mm_c / mass_atm_avg / kg2Pg * vmr2ppmv
					glb_refdata[tmco2_ppmv_ind, totmind] = glb_refdata[vind, totmind] * mm_atm / mm_c / mass_atm_avg / kg2Pg * vmr2ppmv
				}
			}
			
			# calculate monthly differences
			glb_diffdata[vind, totmind] = glb_vardata[vind, totmind] - glb_refdata[vind, totmind]

			# sum monthly stock values for averaging
			glbann_vardata[vind, yind] = glbann_vardata[vind, yind] + glb_vardata[vind, totmind] * dayspermonth[mind]
			glbann_refdata[vind, yind] = glbann_refdata[vind, yind] + glb_refdata[vind, totmind] * dayspermonth[mind]
			
			# calc the annual stock averages at the end of each year
			# and get the annual difference at the end of each year
			if (mind == end_month) {
				glbann_vardata[vind, yind] = glbann_vardata[vind, yind] / daysperyear
				glbann_refdata[vind, yind] = glbann_refdata[vind, yind] / daysperyear
				glbann_diffdata[vind, yind] = glbann_vardata[vind, yind] - glbann_refdata[vind, yind]
			}
			
		} # end loop over column vars
		
		# loop over input vmr variables, mol/molair; store them in ppmv
		for (vind in ind_vmr_vars) {
			grid_vardata[vind,1,1,1] = ncvar_get(nc_id, varid = var_names[vind], start = c(1), count = c(num_time))
			grid_refdata[vind,1,1,1] = ncvar_get(nc_id_ref, varid = var_names[vind], start = c(1), count = c(num_time))
			
			# convert to ppmv
			glb_vardata[vind, totmind] = vmr2ppmv * grid_vardata[vind,1,1,1]
			glb_refdata[vind, totmind] = vmr2ppmv * grid_refdata[vind,1,1,1]
				
			# convert to carbon mass
			if (vind == co2vmr_ind) {
				glb_vardata[co2vmr_pgc_ind, totmind] = grid_vardata[vind,1,1,1] * mm_c / mm_atm * mass_atm_avg * kg2Pg
				glb_refdata[co2vmr_pgc_ind, totmind] = grid_refdata[vind,1,1,1] * mm_c / mm_atm * mass_atm_avg * kg2Pg
			}
			if (vind == ch4vmr_ind) {
				glb_vardata[ch4vmr_pgc_ind, totmind] = grid_vardata[vind,1,1,1] * mm_c / mm_atm * mass_atm_avg * kg2Pg
				glb_refdata[ch4vmr_pgc_ind, totmind] = grid_refdata[vind,1,1,1] * mm_c / mm_atm * mass_atm_avg * kg2Pg
			}

			# monthly difference
			glb_diffdata[vind, totmind] = glb_vardata[vind, totmind] - glb_refdata[vind, totmind]
				
			# sum the values for averaging
			glbann_vardata[vind, yind] = glbann_vardata[vind, yind] + glb_vardata[vind, totmind] * dayspermonth[mind]
			glbann_refdata[vind, yind] = glbann_refdata[vind, yind] + glb_refdata[vind, totmind] * dayspermonth[mind]

			# calc the annual average at the end of each year,
			# and get the annual difference at the end of each year
			if (mind == end_month) { 
				glbann_vardata[vind, yind] = glbann_vardata[vind, yind] / daysperyear
				glbann_refdata[vind, yind] = glbann_refdata[vind, yind] / daysperyear
				glbann_diffdata[vind, yind] = glbann_vardata[vind, yind] - glbann_refdata[vind, yind]
			}
			
		} # end loop over vmr vars
		
		# loop over extra vars; one loop for both types because they aggregate the same way
		# the values have been filled above, but diffs and aggregation have to happen
		for (vind in c(ind_extra_c_vars,ind_extra_ppmv_vars)) {
			# monthly difference
			glb_diffdata[vind, totmind] = glb_vardata[vind, totmind] - glb_refdata[vind, totmind]
				
			# sum the values for averaging
			glbann_vardata[vind, yind] = glbann_vardata[vind, yind] + glb_vardata[vind, totmind] * dayspermonth[mind]
			glbann_refdata[vind, yind] = glbann_refdata[vind, yind] + glb_refdata[vind, totmind] * dayspermonth[mind]

			# calc the annual average at the end of each year,
			# and get the annual difference at the end of each year
			if (mind == end_month) { 
				glbann_vardata[vind, yind] = glbann_vardata[vind, yind] / daysperyear
				glbann_refdata[vind, yind] = glbann_refdata[vind, yind] / daysperyear
				glbann_diffdata[vind, yind] = glbann_vardata[vind, yind] - glbann_refdata[vind, yind]
			}
			
		} # end loop over extra vars
		
		nc_close(nc_id)
		nc_close(nc_id_ref)
	} # end for mind loop over months
	
} # end for y loop over years

# plot the variables in one file
pdf(file = pdfout, paper = "letter")

stitle = paste(inname, "-", refname)

plotind = 0
for(vind in ind_plot_vars) {
	
	plotind = plotind + 1

	# months
	ymin = min(glb_vardata[vind,], glb_refdata[vind,], na.rm = TRUE)
	ymax = max(glb_vardata[vind,], glb_refdata[vind,], na.rm = TRUE)
	ymin = ymin - 0.15 * (ymax - ymin)
	plot.default(x = c(1:tot_months), y = glb_vardata[vind,], main = plot_var_names[plotind],
		xlab = "Month", ylab = paste(ylabs[plotind], "month"),
		ylim = c(ymin,ymax), mgp = c(2, 1, 0),
		type = "n")
	lines(x = c(1:tot_months), y = glb_vardata[vind,], lty = 1, col = "black")
	lines(x = c(1:tot_months), y = glb_refdata[vind,], lty = 2, col = "gray50")
	legend(x = "bottomleft", lty = c(1,2), col = c("black", "gray50"), legend = c(inname, refname))
	# month diff
	ptitle = paste(plot_var_names[plotind], "diff")
	plot(x = c(1:tot_months), y = glb_diffdata[vind,], xlab = "Month", ylab = paste(ylabs[plotind], "month"), main = ptitle, sub = stitle, type = "l")

	# years
	ymin = min(glbann_vardata[vind,], glbann_refdata[vind,], na.rm = TRUE)
    ymax = max(glbann_vardata[vind,], glbann_refdata[vind,], na.rm = TRUE)
	ymin = ymin - 0.15 * (ymax - ymin)
	plot.default(x = c(start_year:end_year), y = glbann_vardata[vind,], main = plot_var_names[plotind],
		xlab = "Year", ylab = paste(ylabs[plotind], "year"),
		ylim = c(ymin,ymax), mgp = c(2, 1, 0),
		type = "n")
	lines(x = c(start_year:end_year), y = glbann_vardata[vind,], lty = 1, col = "black")
	lines(x = c(start_year:end_year), y = glbann_refdata[vind,], lty = 2, col = "gray50")
	legend(x = "bottomleft", lty = c(1,2), col = c("black", "gray50"), legend = c(inname, refname))	
	# year diff
	ptitle = paste(plot_var_names[plotind], "diff")
	plot(x = c(start_year:end_year), y = glbann_diffdata[vind,], xlab = "Year", ylab = paste(ylabs[plotind], "year"), main = ptitle, sub = stitle, type = "l")
		
} # end for vind loop over plots

dev.off()

# write the global annual values to a csv file, via a data frame table
case = c(rep(inname, num_years), rep(refname, num_years))
year = rep(c(start_year:end_year), 2)
SFCO2_PgCperyr = c(glbann_vardata[1,], glbann_refdata[1,])
SFCO2_FFF_PgCperyr = c(glbann_vardata[2,], glbann_refdata[2,])
SFCO2_LND_PgCperyr = c(glbann_vardata[3,], glbann_refdata[3,])
SFCO2_OCN_PgCperyr = c(glbann_vardata[4,], glbann_refdata[4,])
TMCO2_PgC_avg = c(glbann_vardata[5,], glbann_refdata[5,])
TMCO2_FFF_PgC_avg = c(glbann_vardata[6,], glbann_refdata[6,])
TMCO2_LND_PgC_avg = c(glbann_vardata[7,], glbann_refdata[7,])
TMCO2_OCN_PgC_avg = c(glbann_vardata[8,], glbann_refdata[8,])
BURDENBC_PgC_avg = c(glbann_vardata[9,], glbann_refdata[9,])
co2vmr_ppmv_avg = c(glbann_vardata[10,], glbann_refdata[10,])
ch4vmr_ppmv_avg = c(glbann_vardata[11,], glbann_refdata[11,])
co2vmr_PgC_avg = c(glbann_vardata[12,], glbann_refdata[12,])
ch4vmr_PgC_avg = c(glbann_vardata[13,], glbann_refdata[13,])
TMCO2_ppmv_avg = c(glbann_vardata[14,], glbann_refdata[14,])

df = data.frame(case, year, SFCO2_PgCperyr, SFCO2_FFF_PgCperyr, SFCO2_LND_PgCperyr, SFCO2_OCN_PgCperyr, TMCO2_PgC_avg, TMCO2_FFF_PgC_avg, TMCO2_LND_PgC_avg, TMCO2_OCN_PgC_avg, BURDENBC_PgC_avg, co2vmr_ppmv_avg, ch4vmr_ppmv_avg, co2vmr_PgC_avg, ch4vmr_PgC_avg, TMCO2_ppmv_avg, stringsAsFactors=FALSE)
write.csv(df, csvout, row.names=FALSE)

cat("finished plot_cam_carbon.r", date(), "\n\n")
