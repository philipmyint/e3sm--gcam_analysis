###########
#
# read and plot standard elm 2 h0 history files
#  these are monthly files
#  read resolution from files
#   time = 1
#
# note that the dimension order listed below is that in the files
#
# need the land area of each grid cell (might eventually need pftmask)
#  area(lon, lat) area of whole grid cell, km^2
#  landfrac(lon, lat) fraction of grid cell that is land
#  landmask(lon, lat) 1=land, 0=ocean
#  missing/fill value = 1.0E36f
#

# C flux variables:
#  gC/m^2/s
#  mean over the time period (month)
#  missing/fill value = 1.0E36f
#  need to convert to Pg C per month
#  land area coverage appears to be correct
#  only one time record
#  float DWT_CONV_CFLUX_GRC(lon, lat, time); only on first timestep of year - compare with dribbled
#  float DWT_CONV_CFLUX_DRIBBLED(lon, lat, time)
#  float ER(lon, lat, time); total ecosystem respiration, autotrophic plus heterotrophic
#  float GPP(lon, lat, time)
#  float HR(lon, lat, time)
#  float LAND_UPTAKE(lon, lat, time); NEE minus LAND_USE_FLUX, negative for uptake
#  float LAND_USE_FLUX(lon, lat, time); total C emitted from conversion and wood products
#  float NBP(lon, lat, time); includes fir, landuse, and harvest flux, positive for sink
#  float NEE(lon, lat, time); includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source
#  float NEP(lon, lat, time); excludes fire, landuse, harvest, and hrv_xsmrpool flux, positive for sink
#  float NPP(lon, lat, time)
#  float PFT_FIRE_CLOSS(lon, lat, time); patch-level fire c loss for non-peat fires outsiede land-type converted region
#  float WOOD_HARVESTC(lon, lat, time)
#  float WOODC_ALLOC(lon, lat, time)
#  float WOODC_LOSS(lon, lat, time)

# area variables:
#  float FAREA_BURNED(lon, lat, time); fraction of column, but not sure how to normalize because different columns are added, and in the code this is per second?
#  float PCT_LANDUNIT(lon, lat, ltype, time); percent of land unit on topounit (not implemented)
#  float PCT_NAT_PFT(lon, lat, natpft, time); percent of pft on veg landunit (not implemented)
#  float TLAI(lon, lat, time); total projected leaf area index, unitless

# carbon state variables:
#  float TOTECOSYSC(lon, lat, time); gC/m^2, includes veg and excludes cpool and product pools
#  float TOTSOMC(lon, lat, time); gC/m^2
#  float TOTVEGC(lon, lat, time); gC/m^2, excludes cpool
#  float WOODC(lon, lat, time); gC/m^2

cat("started plot_elm_hist0", date(), "\n\n")

library(ncdf4)

#datadir = "/Users/adivi/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd/"
#outdir = "/Users/adivi/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd/"
#datadir = "~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1deg_14feb2023/"
#outdir = "~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1deg_14feb2023/"
#datadir = "~/projects/e3sm/iesmv2/test_diagnostics/gcam6/"
#outdir = "~/projects/e3sm/iesmv2/test_diagnostics/gcam6/"
datadir = "~/projects/e3sm/iesmv2/test_diagnostics/gcam6_1deg/"
outdir = "~/projects/e3sm/iesmv2/test_diagnostics/gcam6_1deg/"

#base_name_hist = "I20TRQIAELMCNPRDCTCBC_f19_g16_trans.elm.h0"
#base_name_fut = "20230110_zlnd_gcam_test.elm.h0"
#base_name_hist = "20230214_I20TRQIAELMCNPRDCTCBC_f09_g16_trans.elm.h0"
#base_name_fut = "20230214_zlnd_gcam_test.elm.h0"
#base_name_hist = "20230412_I20TRQIAELMCNPRDCTCBC_f19_g16_trans.elm.h0"
#base_name_fut = "18January2023_zlnd_gcam6.elm.h0"
base_name_hist = "20230214_I20TRQIAELMCNPRDCTCBC_f09_g16_trans.elm.h0"
base_name_fut = "27March2023_zlnd_gcam6_f09_g16.elm.h0"

ntag = ".nc"

start_year_hist = 2013
end_year_hist = 2014
months_hist = c(1,6)
num_months_hist = length(months_hist)
num_years_hist = end_year_hist - start_year_hist + 1
tot_months_hist = num_years_hist * num_months_hist

start_year_fut = 2015
end_year_fut = 2015
months_fut = c(1:12)
num_months_fut = length(months_fut)
num_years_fut = end_year_fut - start_year_fut + 1
tot_months_fut = num_years_fut * num_months_fut

tot_months = (num_months_hist*num_years_hist) + (num_months_fut*num_years_fut)

plot_months = array(dim=tot_months)

pdfout = paste(outdir, "check_h0_", start_year_hist, "-", end_year_fut, ".pdf", sep = "")
csvout = paste(outdir, "check_h0_", start_year_hist, "-", end_year_fut, ".csv", sep = "")

var_flux_names = c("DWT_CONV_CFLUX_GRC", "DWT_CONV_CFLUX_DRIBBLED", "ER", "GPP", "HR", "LAND_UPTAKE", "LAND_USE_FLUX",
	"NBP", "NEE", "NEP", "NPP", "PFT_FIRE_CLOSS", "WOOD_HARVESTC", "WOODC_ALLOC", "WOODC_LOSS")
num_flux_vars = length(var_flux_names)

var_dimarea_names = c("FAREA_BURNED", "TLAI")
num_dimarea_vars = length(var_dimarea_names)

var_cstate_names = c("TOTECOSYSC", "TOTSOMC", "TOTVEGC", "WOODC")
num_cstate_vars = length(var_cstate_names)

# seconds per month
secpermonth = 86400 * c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# convert grams to Pg
g2Pg = 1E-15

# convert from km^2 to m^2
kmsq2msq = 1E6

# get some necessary info from the first file 
fname = paste(datadir, base_name_hist, ".", start_year_hist, "-01", ntag, sep = "")	
nc_id = nc_open(fname)
num_lon = nc_id$dim$lon$len
num_lat = nc_id$dim$lat$len
num_gridcell = nc_id$dim$gridcell$len
num_pft = nc_id$dim$natpft$len
num_lunit = nc_id$dim$ltype$len
num_time = nc_id$dim$time$len

area = ncvar_get(nc_id, varid = "area", start = c(1,1), count = c(num_lon, num_lat))  # km^2
landfrac = ncvar_get(nc_id, varid = "landfrac", start = c(1,1), count = c(num_lon, num_lat))
landmask = ncvar_get(nc_id, varid = "landmask", start = c(1,1), count = c(num_lon, num_lat))
pftmask = ncvar_get(nc_id, varid = "pftmask", start = c(1,1), count = c(num_lon, num_lat))
nc_close(nc_id)

# convert this to m^2
landarea = area * landfrac * landmask * pftmask * kmsq2msq

# arrays for storing time series of variables
grid_flux_data = array(dim=c(num_flux_vars, num_lon, num_lat, tot_months))
grid_dimarea_data = array(dim=c(num_dimarea_vars, num_lon, num_lat, tot_months))
grid_cstate_data = array(dim=c(num_cstate_vars, num_lon, num_lat, tot_months))

# arrays for storing time series of global sums of variables
glb_flux_data = array(dim = c(num_flux_vars, tot_months))
glb_dimarea_data = array(dim = c(num_dimarea_vars, tot_months))
glb_cstate_data = array(dim = c(num_cstate_vars, tot_months))

totmind=0
# loop over year
for (yr in c(c(start_year_hist:end_year_hist),c(start_year_fut:end_year_fut))) {
	cat("starting year", yr, date(), "\n")
	
	if(yr <= end_year_hist & yr >= start_year_hist ){
		months = months_hist
		base_name = base_name_hist
	} else if(yr <= end_year_fut & yr >= start_year_fut ){
		months = months_fut
		base_name = base_name_fut
	}
	
	# loop over month
	for (mind in months) {
		if(mind > 9) {
			dtag = paste0(yr, "-", mind)
		}else {
			dtag = paste0(yr, "-", 0, mind)
		}
		fname = paste(datadir, base_name, ".", dtag, ntag, sep = "")
		nc_id = nc_open(fname)
		totmind = totmind + 1
		plot_months[totmind] = dtag
		
		# loop over flux variables
		for (vind in 1:num_flux_vars) {
			grid_flux_data[vind,,,totmind] = ncvar_get(nc_id, varid = var_flux_names[vind], start = c(1,1,1), count = c(num_lon, num_lat, num_time))
			# convert to Pg C per month
			grid_flux_data[vind,,,totmind] = grid_flux_data[vind,,,totmind] * landarea * secpermonth[mind] * g2Pg
			glb_flux_data[vind, totmind] = sum(grid_flux_data[vind,,,totmind], na.rm = TRUE)
		}
		
		# loop over dimensionless area variables
		# note that global sums are of the dimensionless data, not normalized
		for (vind in 1:num_dimarea_vars) {
			grid_dimarea_data[vind,,,totmind] = ncvar_get(nc_id, varid = var_dimarea_names[vind], start = c(1,1,1), count = c(num_lon, num_lat, num_time))
			# no conversion
			glb_dimarea_data[vind, totmind] = sum(grid_dimarea_data[vind,,,totmind], na.rm = TRUE)
		}
		
		# loop over cstate variables
		for (vind in 1:num_cstate_vars) {
			grid_cstate_data[vind,,,totmind] = ncvar_get(nc_id, varid = var_cstate_names[vind], start = c(1,1,1), count = c(num_lon, num_lat, num_time))
			# convert to Pg C
			grid_cstate_data[vind,,,totmind] = grid_cstate_data[vind,,,totmind] * landarea * g2Pg
			glb_cstate_data[vind, totmind] = sum(grid_cstate_data[vind,,,totmind], na.rm = TRUE)
		}
		
		nc_close(nc_id)
	} # end mind loop
} # end yr loop

data_df = data.frame(month=c(1:tot_months))

# plot the global sum variables in one file
pdf(file = pdfout, paper = "letter")

for(vind in 1:num_flux_vars) {

	plot(x = c(1:tot_months), y = glb_flux_data[vind,], xlab = "Months", ylab = "Pg C per month", main = var_flux_names[vind], type = "l", xaxt="n")
	axis(1, at = c(1:tot_months), labels = plot_months)
	data_df$ncol = glb_flux_data[vind,]
	names(data_df)[ncol(data_df)] = var_flux_names[vind]
}

for(vind in 1:num_dimarea_vars) {

	plot(x = c(1:tot_months), y = glb_dimarea_data[vind,], xlab = "Months", ylab = "dimensionless sum", main = var_dimarea_names[vind], type = "l", xaxt="n")
	axis(1, at = c(1:tot_months), labels = plot_months)
	data_df$ncol = glb_dimarea_data[vind,]
	names(data_df)[ncol(data_df)] = var_dimarea_names[vind]
}

for(vind in 1:num_cstate_vars) {

	plot(x = c(1:tot_months), y = glb_cstate_data[vind,], xlab = "Months", ylab = "Pg C", main = var_cstate_names[vind], type = "l", xaxt="n")
	axis(1, at = c(1:tot_months), labels = plot_months)
	data_df$ncol = glb_cstate_data[vind,]
	names(data_df)[ncol(data_df)] = var_cstate_names[vind]
}

dev.off()

write.csv(data_df, file=csvout,row.names = FALSE)

cat("finished plot_elm_hist0", date(), "\n\n")