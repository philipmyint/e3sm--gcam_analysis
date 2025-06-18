# plot_elm_h0_pfts.r

# plots time series of h0 pfts across historical and future boundary

# uses the landfrac in the h0 file, which is the land frac that the atmosphere sees (based on the domain file)
# also try using the static land info from the historical land surface file

library(raster)
library(rasterVis)
library(ncdf4)
library(parallel)
library(plyr)

#datadir = "~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_20dec2022/"
#outdir = "~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_20dec2022/"
#datadir = "~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1deg_14feb2023/"
#outdir = "~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1deg_14feb2023/"
datadir = "~/projects/e3sm/iesmv2/test_diagnostics/ensemble_restarts/"
outdir = "~/projects/e3sm/iesmv2/test_diagnostics/ensemble_restarts/"

#hist_base = "20221220_I20TRQIAELMCNPRDCTCBC_f19_g16_trans.elm.h0"
hist_base = "20230214_I20TRQIAELMCNPRDCTCBC_f09_g16_trans.elm.h0"

#fut_base = "20230110_zlnd_gcam_test.elm.h0"
fut_base = "20230214_zlnd_gcam_test.elm.h0"
fut_base = "20250110_I20TREAMELMCNPRDCTCBCBGC_ne30pg2_f09_oEC60to30v3.elm.h0"

nctag = ".nc"

#hist_lt_file = "~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_/landuse.timeseries_1.9x2.5_hist_simyr1850-2015_c221220.nc"
hist_lt_file = "~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1deg_14feb2023/landuse.timeseries_0.9x1.25_hist_simyr1850-2015_c201021.nc"

#h0_dom_file = "~/projects/e3sm/iesmv2/test_diagnostics/domain.lnd.fv1.9x2.5_gx1v6.090206.nc"
h0_dom_file = "~/projects/e3sm/iesmv2/inputs/domain.lnd.fv0.9x1.25_gx1v6.090309.nc"

start_year_hist = 2013
end_year_hist = 2014
months_hist = c(1,6)
num_months_hist = length(months_hist)
num_years_hist = end_year_hist - start_year_hist + 1
tot_months_hist = num_years_hist * num_months_hist

start_year_fut = 2015
end_year_fut = 2022
months_fut = c(1:12)
num_months_fut = length(months_fut)
num_years_fut = end_year_fut - start_year_fut + 1
tot_months_fut = num_years_fut * num_months_fut

tot_months = (num_months_hist*num_years_hist) + (num_months_fut*num_years_fut)

plot_months = array(dim=tot_months)

PROJ4_STRING = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# only the 16 valid ones are included in processing; these are in order of storage
pft_names = c("Bare", "NdlEgTmpTr", "NdlEgBorTr", "NdlDcBorTr", "BrdEgTrpTr", "BrdEgTmpTr", "BrdDcTrpTr", "BrdDcTmpTr", "BrdDcBorTr",
	"BrdEgTmpSh", "BrdDcTmpSh", "BrdDcBorSh", "C3ArcGr", "C3NonArcGr", "C4Gr","Crop")
num_pfts = length(pft_names)

h0_pft_areas_glb = array(dim=c(tot_months, num_pfts))

harvest_names = c("HARVEST_VH1", "HARVEST_VH2", "HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3")
num_harv = length(harvest_names)

pdfout = paste(outdir, "check_h0_pfts_", start_year_hist, "-", end_year_fut, ".pdf", sep = "")

#### get some h0 info that doesn't change
if(months_hist[1] > 9) {
	dtag = months_hist[1]
}else {
	dtag = paste0(0, months_hist[1])
}
h0file = paste0(datadir,hist_base, ".",start_year_hist,"-",dtag,nctag)
cat("processing", h0file, "\n")
h0 = nc_open(h0file)

# get the dimensions and the cell centers
# use the given lat lon order
nlon_out = h0$dim$lon$len
nlat_out = h0$dim$lat$len
ncells_out = nlon_out * nlat_out
lon_out = ncvar_get(h0,varid="lon")
lat_out = ncvar_get(h0,varid="lat")

# grid cell area
h0_cellarea = raster(h0file, values = TRUE, varname = "area")

# read in the fraction of grid cell that is not ocean (fraction)
# this is the land portion of a grid cell
h0_landfrac = raster(h0file, values = TRUE, varname = "landfrac")

# landmask matches landfrac, but read it in anyway
h0_landmask = raster(h0file, values = TRUE, varname = "landmask")

# this determines if pfts are actually there or not
h0_pftmask = raster(h0file, values = TRUE, varname = "pftmask")

# read in the land unit percentages of land portion of grid cell
# 9 land units in order: veg, crop, ice, multiple ice, lake, wetland, urban tbd, urban hd, urban md
h0_veg = raster(h0file, values = TRUE, varname = "PCT_LANDUNIT", band = 1)
h0_crop = raster(h0file, values = TRUE, varname = "PCT_LANDUNIT", band = 2)
h0_ice = raster(h0file, values = TRUE, varname = "PCT_LANDUNIT", band = 3)
h0_multice = raster(h0file, values = TRUE, varname = "PCT_LANDUNIT", band = 4)
h0_lake = raster(h0file, values = TRUE, varname = "PCT_LANDUNIT", band = 5)
h0_wetland = raster(h0file, values = TRUE, varname = "PCT_LANDUNIT", band = 6)
h0_urban_tbd = raster(h0file, values = TRUE, varname = "PCT_LANDUNIT", band = 7)
h0_urban_hd = raster(h0file, values = TRUE, varname = "PCT_LANDUNIT", band = 8)
h0_urban_md = raster(h0file, values = TRUE, varname = "PCT_LANDUNIT", band = 9)

# calculate the area of the land (km^2)
# may need to drop the pftmask because it cuts a few land cells
h0_land_area = h0_cellarea * h0_landfrac * h0_pftmask * h0_landmask
h0_land_area_glb = cellStats(h0_land_area, stat = "sum")
cat("land area (km^2)", h0_land_area_glb, "\n")
levelplot(h0_land_area, margin=FALSE, main = "h0 land area (km^2)")

# check the vegetated land unit as percent of land portion of a grid cell
# the max and min of the sum should be 100
h0_lu_sum = h0_veg + h0_crop + h0_ice + h0_multice + h0_lake + h0_wetland + h0_urban_tbd + h0_urban_hd + h0_urban_md
cat("h0 land unit max is",  cellStats(h0_lu_sum, stat="max"), "min is", cellStats(h0_lu_sum, stat="min"), "\n")
h0_veglu = 100 - h0_crop - h0_ice - h0_multice - h0_lake - h0_wetland - h0_urban_tbd - h0_urban_hd - h0_urban_md
dhv = h0_veglu - h0_veg
nzh_inds = Which(dhv != 0, cells = TRUE)
cat("number of h0 non-zero veg lu diffs is",  length(nzh_inds), "\n")
cat("max diffs is", max(dhv[nzh_inds]), "min diffs is", min(dhv[nzh_inds]), "\n")
levelplot(h0_veg, margin = FALSE, main = "h0 veg land unit, percent of land portion")

# calculate the area of the vegetated land unit (km^2)
h0_veglu_area = h0_cellarea * h0_landfrac * h0_pftmask * h0_landmask * h0_veg / 100
h0_veglu_area_glb = cellStats(h0_veglu_area, stat = "sum")
cat("veg land unit area (km^2)", h0_veglu_area_glb, "\n")
levelplot(h0_veglu_area, margin=FALSE, main = "h0 veg land unit area (km^2)")

nc_close(h0)

### get the lt land info instead

hist_lt = nc_open(hist_lt_file)

num_lon_hist = hist_lt$dim$lsmlon$len
num_lat_hist = hist_lt$dim$lsmlat$len
num_time_hist = hist_lt$dim$time$len

# note that landfrac and landmask coincide in the th lt file
# need to set the geo parameters

area_hist = raster(hist_lt_file, values = TRUE, varname = "AREA")
extent(area_hist) <- h0_cellarea
projection(area_hist) <- h0_cellarea
lf_hist = raster(hist_lt_file, values = TRUE, varname = "LANDFRAC_PFT")
extent(lf_hist) <- h0_cellarea
projection(lf_hist) <- h0_cellarea
mask_hist = raster(hist_lt_file, values = TRUE, varname = "PFTDATA_MASK")
extent(mask_hist) <- h0_cellarea
projection(mask_hist) <- h0_cellarea
veg_hist = raster(hist_lt_file, values = TRUE, varname = "PCT_NATVEG")
extent(veg_hist) <- h0_cellarea
projection(veg_hist) <- h0_cellarea

veg_area_hist = veg_hist / 100 * lf_hist * area_hist * mask_hist

nc_close(hist_lt)

totmind=0
start_fut_mind = -1
# loop over year
for (yr in c(c(start_year_hist:end_year_hist),c(start_year_fut:end_year_fut))) {
	cat("starting year", yr, date(), "\n")
	
	if(yr <= end_year_hist & yr >= start_year_hist ){
		months = months_hist
		base_name = hist_base
		HIST_FLAG = TRUE
	} else if(yr <= end_year_fut & yr >= start_year_fut ){
		months = months_fut
		base_name = fut_base
		HIST_FLAG = FALSE
	}

	# loop over month
	for (mind in months) {
		if(mind > 9) {
			dtag = paste0(yr, "-", mind)
		}else {
			dtag = paste0(yr, "-", 0, mind)
		}
		fname = paste(datadir, base_name, ".", dtag, nctag, sep = "")
		cat("processing", fname, "\n")
		nc_id = nc_open(fname)
		totmind = totmind + 1
		plot_months[totmind] = dtag
		
		# keep track of the boundary values
		if (HIST_FLAG) { end_hist_mind = totmind
		}else {
			if (start_fut_mind == -1 ) { start_fut_mind = totmind }
		} 
		
		# read in the pft percentages of veg land unit
		# band is the pft in this case
		# 17 pfts in order: 1 bare, 8 tree, 3 shrub, 3 grass, 1 crop, 1 empty
		# use a loop here to read into a stack; don't include the empty pft
		h0_pft_stack = stack()
		for (i in 1:num_pfts) {
			rin = raster(fname, values = TRUE, varname = "PCT_NAT_PFT", band = i)
			h0_pft_stack <- addLayer(h0_pft_stack, rin)
		}

		# calculate the pft areas (km^2)
		#h0_pft_areas = h0_pft_stack / 100 * h0_veglu_area
		h0_pft_areas = h0_pft_stack / 100 * veg_area_hist
		h0_pft_areas_glb[totmind,] = cellStats(h0_pft_areas, stat = "sum")
		
		nc_close(nc_id)
	} # end mind loop

	
} # end yr loop

forest = array(dim=c(tot_months))
shrub = array(dim=c(tot_months))
grass = array(dim=c(tot_months))
crop = array(dim=c(tot_months))
bare = array(dim=c(tot_months))
forest[] = 0
shrub[] = 0
grass[] = 0
crop[] = 0
bare[] = 0

# aggregatge some pfts
for(i in 1: num_pfts) {

	if (i >=2 & i <=9) {
		forest = forest + h0_pft_areas_glb[,i]
	}
	
	if (i >=10 & i <=12) {
		shrub = shrub + h0_pft_areas_glb[,i]
	}
	
	if (i >=13 & i <=15) {
		grass = grass + h0_pft_areas_glb[,i]
	}
	
	if (i == 16) {
		crop = crop + h0_pft_areas_glb[,i]
	}
	
	if (i == 1) {
		bare = bare + h0_pft_areas_glb[,i]
	}
	
}

# plot the global timeseries of each pft

pdf(pdfout)

for ( p in 1:num_pfts) {

	plot(x = c(1:tot_months), y = h0_pft_areas_glb[,p], xlab = "Months", ylab = "km^2", main = paste(pft_names[p], "h0 global area"),
		sub = paste("hist last =", h0_pft_areas_glb[end_hist_mind,p], "; fut first =", h0_pft_areas_glb[start_fut_mind,p]), type = "l", xaxt="n")
	axis(1, at = c(1:tot_months), labels = plot_months)
}

plot(x = c(1:tot_months), y = forest, xlab = "Months", ylab = "km^2", main = paste("forest h0 global area"),
	sub = paste("hist last =", forest[end_hist_mind], "; fut first =", forest[start_fut_mind]), type = "l", xaxt="n")
axis(1, at = c(1:tot_months), labels = plot_months)

plot(x = c(1:tot_months), y = shrub, xlab = "Months", ylab = "km^2", main = paste("shrub h0 global area"),
	sub = paste("hist last =", shrub[end_hist_mind], "; fut first =", shrub[start_fut_mind]), type = "l", xaxt="n")
axis(1, at = c(1:tot_months), labels = plot_months)

plot(x = c(1:tot_months), y = grass, xlab = "Months", ylab = "km^2", main = paste("grass h0 global area"),
	sub = paste("hist last =", grass[end_hist_mind], "; fut first =", grass[start_fut_mind]), type = "l", xaxt="n")
axis(1, at = c(1:tot_months), labels = plot_months)

#plot(x = c(1:tot_months), y = crop, xlab = "Months", ylab = "km^2", main = paste("crop h0 global area"),
#	sub = paste("hist last =", crop[end_hist_mind], "; fut first =", crop[start_fut_mind]), type = "l", xaxt="n")
#axis(1, at = c(1:tot_months), labels = plot_months)

plot(x = c(1:tot_months), y = bare, xlab = "Months", ylab = "km^2", main = paste("bare h0 global area"),
	sub = paste("hist last =", bare[end_hist_mind], "; fut first =", bare[start_fut_mind]), type = "l", xaxt="n")
axis(1, at = c(1:tot_months), labels = plot_months)

dev.off()






