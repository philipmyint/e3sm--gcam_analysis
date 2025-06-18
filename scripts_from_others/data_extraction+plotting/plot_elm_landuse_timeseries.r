### code to make simple plots of e3sm landuse timeseries

library(ncdf4)

pftnames = c("bare", "ndlevrtmptree", "ndlevrbortree", "ndldecbortree", "brdevrtrptree", "brdevrtmptree", "brddectrptree", "brddectmptree", "brddecbortree", "brdevrtmpshrub", "brddectmpshrub", "brddecborshrub", "c3arcgrass", "c3nonarcgrass", "c4grass", "crop", "na")

# only do the valid ones
# crop is pft 16 and so don't need to do it again to group the generics at the end
npft = 16

# calculate only the harvested forest area, as per ELM, rather than the total harvested area in the veg unit, as per the input
harvest_names = c("HARVEST_VH1", "HARVEST_VH2", "HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3")
num_harv = length(harvest_names)


# areas are in km^2

#### future file

#future = nc_open("~/projects/e3sm/diagnostics/landuse.timeseries_0.5x0.5_SSP1_RCP26_simyr2015-2100_c220317.nc")
#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd/surfdata_iESM_dyn.nc")
#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_hybrid/surfdata_1.9x2.5_SSP5_RCP85_simyr2015_c210916.nc")
#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_hybrid/surfdata_iESM.nc")
#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_hybrid/landuse.timeseries_1.9x2.5_SSP5_RCP85_simyr2015-2100_c210916.nc")
#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1feb2023/surfdata_iESM_dyn.nc")
#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1feb2023/landuse.timeseries_0.9x1.25_HIST_simyr2015_c201021.nc")
#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/surfdata_iESM_dyn_2deg.nc")
#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/surfdata_iESM_2deg.nc")
#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1deg_14feb2023/surfdata_iESM.nc")
#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1deg_14feb2023/surfdata_iESM_dyn.nc")
#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/gcam6/surfdata_iESM_dyn.nc")
#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_20231010_zlnd_gcam6_SSP2/surfdata_iESM_dyn.nc")
#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/20240327_SSP245_ZATM_ne30pg2_f09_oEC60to30v3/surfdata_iESM_dyn.nc")

#future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_20240718_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3/surfdata_iESM_dyn.nc")

future = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_20240730_SSP245_ZATM/feedbacks/surfdata_iESM_dyn.nc")

num_lon_future = future$dim$lsmlon$len
num_lat_future = future$dim$lsmlat$len

num_time_future = future$dim$time$len
#num_time_future = 1

pfts_future = ncvar_get(future, varid="PCT_NAT_PFT")
area_future = ncvar_get(future, varid="AREA")
lf_future = ncvar_get(future, varid="LANDFRAC_PFT")
mask_future = ncvar_get(future, varid="PFTDATA_MASK")
veg_future = ncvar_get(future, varid="PCT_NATVEG")

year_future = ncvar_get(future, varid="YEAR")
nyear_future=length(year_future)
#year_future = 2015
#nyear_future = 1

# read in the harvest fractions of veg land unit and add them together
harvest_frac_future = array(dim=c(num_lon_future, num_lat_future, num_harv, num_time_future))
harvest_frac_future[,,,] = 0
for (i in 1:num_harv) {
	harvest_frac_future[,,i,] = ncvar_get(future, varid=harvest_names[i])
}

nc_close(future)

# get the indices of the unique years, by taking the last index of each
list_inds = split(seq_along(year_future), year_future)
tinds_future = unname(unlist(lapply(list_inds, FUN=max)))

#tinds_future= c(1:nyear_future)
year_future = year_future[tinds_future]
nyear_future = length(tinds_future)

pa_future = array(dim=c(nyear_future))
forest_future = array(dim=c(nyear_future))
shrub_future = array(dim=c(nyear_future))
grass_future = array(dim=c(nyear_future))
crop_future = array(dim=c(nyear_future))
bare_future = array(dim=c(nyear_future))
harvest_area_future = array(dim=c(nyear_future))
forest_frac_future = array(dim=c(num_lon_future, num_lat_future, nyear_future))

veg_area_future = veg_future / 100 * lf_future * area_future * mask_future

#pdf("~/projects/e3sm/diagnostics/check_ssp1_rcp26.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd/check_surfdata_iESM_dyn.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_hybrid/check_fsurdat.pdf"
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_hybrid/check_ssp5_rcp85.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1feb2023/check_surfdata_iESM_dyn.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1feb2023/check_lt15.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/check_surfdata_iESM_dyn_2deg.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/check_surfdata_iESM_2deg.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1deg_14feb2023/check_surfdata_iESM.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1deg_14feb2023/check_surfdata_iESM_dyn.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/gcam6/check_surfdata_iESM_dyn.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_20231010_zlnd_gcam6_SSP2/check_surfdata_iESM_dyn.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/20240327_SSP245_ZATM_ne30pg2_f09_oEC60to30v3/check_surfdata_iESM_dyn.pdf")

#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_20240718_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3/check_surfdata_iESM_dyn.pdf")

pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_20240730_SSP245_ZATM/feedbacks/check_surfdata_iESM_dyn.pdf")


forest_future[] = 0
shrub_future[] = 0
grass_future[] = 0
crop_future[] = 0
bare_future[] = 0
harvest_area_future[] = 0
forest_frac_future[,,] = 0

for(i in 1: npft) {
	pa_future[] = 0
	for(t in 1: nyear_future) {
		pa_future[t] = sum(pfts_future[,,i,tinds_future[t]] / 100 * veg_area_future, na.rm=TRUE)
		#pa_future[t] = sum(pfts_future[,,i] / 100 * veg_area_future, na.rm=TRUE)
		
		# get the gridded forest frac total
		if (i >=2 & i <=9) {
			forest_frac_future[,,t] = forest_frac_future[,,t] + pfts_future[,,i,tinds_future[t]] / 100
		}
		
	}
	plot(year_future,pa_future,main=pftnames[i], type="l", sub=paste(year_future[1],pa_future[1]))
	
	if (i >=2 & i <=9) {
		forest_future = forest_future + pa_future
	}
	
	if (i >=10 & i <=12) {
		shrub_future = shrub_future + pa_future
	}
	
	if (i >=13 & i <=15) {
		grass_future = grass_future + pa_future
	}
	
	if (i == 16) {
		crop_future = crop_future + pa_future
	}
	
	if (i == 1) {
		bare_future = bare_future + pa_future
	}
	
}

for (i in 1:num_harv) {
	for(t in 1: nyear_future) {
		harvest_area_future[t] = harvest_area_future[t] + sum(harvest_frac_future[,,i,tinds_future[t]] * veg_area_future * forest_frac_future[,,t], na.rm=TRUE)
	}
}

plot(year_future,forest_future,main="forest km^2", type="l", sub=paste(year_future[1],forest_future[1]))
plot(year_future,shrub_future,main="shrub km^2", type="l", sub=paste(year_future[1],shrub_future[1]))
plot(year_future,grass_future,main="grass km^2", type="l", sub=paste(year_future[1],grass_future[1]))
#plot(year_future,crop_future,main="crop km^2", type="l", sub=paste(year_future[1],crop_future[1]))
plot(year_future,bare_future,main="bare km^2", type="l", sub=paste(year_future[1],bare_future[1]))

plot(year_future,harvest_area_future,main="harvest area km^2", type="l", sub=paste(year_future[1],harvest_area_future[1]))

dev.off()



#### historical file

#hist = nc_open("~/projects/e3sm/diagnostics/landuse.timeseries_0.5x0.5_HIST_simyr1850-2015_c200924.nc")
#hist = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd/landuse.timeseries_1.9x2.5_hist_simyr1850-2015_c221220.nc")
#hist = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1deg_14feb2023/landuse.timeseries_0.9x1.25_hist_simyr1850-2015_c201021.nc")
#hist = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/gcam6_1deg/surfdata_iESM_dyn_zlnd_hist.nc")
#hist = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_20230509_2deg_zlnd_hist/surfdata_iESM_dyn.nc")
#hist = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_20230927_zlnd_gcam6_SSP2_no_carbon_scaling/surfdata_iESM_dyn.nc")
#hist = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/20240327_SSP245_ZATM_ne30pg2_f09_oEC60to30v3_MOD/surfdata_iESM_dyn.nc")

#hist = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_20240718_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3/no_feedbacks/surfdata_iESM_dyn.nc")

hist = nc_open("~/projects/e3sm/iesmv2/test_diagnostics/eva_20240730_SSP245_ZATM/no_feedbacks/surfdata_iESM_dyn.nc")

num_lon_hist = hist$dim$lsmlon$len
num_lat_hist = hist$dim$lsmlat$len
num_time_hist = hist$dim$time$len

pfts_hist = ncvar_get(hist, varid="PCT_NAT_PFT")
area_hist = ncvar_get(hist, varid="AREA")
lf_hist = ncvar_get(hist, varid="LANDFRAC_PFT")
mask_hist = ncvar_get(hist, varid="PFTDATA_MASK")
veg_hist = ncvar_get(hist, varid="PCT_NATVEG")
year_hist = ncvar_get(hist, varid="YEAR")
nyear_hist=length(year_hist)

# read in the harvest fractions of veg land unit and add them together
harvest_frac_hist = array(dim=c(num_lon_hist, num_lat_hist, num_harv, num_time_hist))
harvest_frac_hist[,,,] = 0
for (i in 1:num_harv) {
	harvest_frac_hist[,,i,] = ncvar_get(hist, varid=harvest_names[i])
}

nc_close(hist)

# get the indices of the unique years, by taking the last index of each
list_inds = split(seq_along(year_hist), year_hist)
tinds_hist = unname(unlist(lapply(list_inds, FUN=max)))

#tinds_hist= c(1:nyear_hist)
year_hist = year_hist[tinds_hist]
nyear_hist = length(tinds_hist)

pa_hist = array(dim=c(nyear_hist))
forest_hist = array(dim=c(nyear_hist))
shrub_hist = array(dim=c(nyear_hist))
grass_hist = array(dim=c(nyear_hist))
crop_hist = array(dim=c(nyear_hist))
bare_hist = array(dim=c(nyear_hist))
harvest_area_hist = array(dim=c(nyear_hist))
forest_frac_hist = array(dim=c(num_lon_future, num_lat_future, nyear_hist))

veg_area_hist = veg_hist / 100 * lf_hist * area_hist * mask_hist

# get the static land info from the h0 file
#h0_file = "~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_20dec2022/20221220_I20TRQIAELMCNPRDCTCBC_f19_g16_trans.elm.h0.2010-01.nc"
#h0 = nc_open(h0_file)
#h0_cellarea = ncvar_get(h0, varid="area")
#h0_landfrac = ncvar_get(h0, varid="landfrac")
#h0_pftmask = ncvar_get(h0, varid="pftmask")
#h0_lu = ncvar_get(h0, varid="PCT_LANDUNIT")
#h0_veg = h0_lu[,,1]

#h0_veg_area = h0_veg / 100 * h0_landfrac * h0_cellarea * h0_pftmask

#pdf("~/projects/e3sm/iesmv2/test_diagnostics/check_hist.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_20dec2022/check_hist_h0norm.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1deg_14feb2023/check_hist.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/gcam6_1deg/check_hist_zlnd.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_20230509_2deg_zlnd_hist/hist_zlnd.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_20230927_zlnd_gcam6_SSP2_no_carbon_scaling/check_surfdata_iesm_dyn.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/20240327_SSP245_ZATM_ne30pg2_f09_oEC60to30v3_MOD/check_surfdata_iesm_dyn.pdf")

#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_20240718_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3/no_feedbacks/check_surfdata_iesm_dyn.pdf")

pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_20240730_SSP245_ZATM/no_feedbacks/check_surfdata_iESM_dyn.pdf")

forest_hist[] = 0
shrub_hist[] = 0
grass_hist[] = 0
crop_hist[] = 0
bare_hist[] = 0
harvest_area_hist[] = 0
forest_frac_hist[,,] = 0

for(i in 1:npft) {
	pa_hist[] = 0
	for(t in 1: nyear_hist) {
		pa_hist[t] = sum(pfts_hist[,,i,tinds_hist[t]] / 100 * veg_area_hist, na.rm=TRUE)
		#pa_hist[t] = sum(pfts_hist[,,i,t] / 100 * h0_veg_area, na.rm=TRUE)
		
		# get the gridded forest frac total
		if (i >=2 & i <=9) {
			forest_frac_hist[,,t] = forest_frac_hist[,,t] + pfts_hist[,,i,tinds_hist[t]] / 100
		}
	}
	plot(year_hist,pa_hist,main=pftnames[i], type="l", sub=paste(year_hist[nyear_hist],pa_hist[nyear_hist]))
	
	if (i >=2 & i <=9) {
		forest_hist = forest_hist + pa_hist
	}
	
	if (i >=10 & i <=12) {
		shrub_hist = shrub_hist + pa_hist
	}
	
	if (i >=13 & i <=15) {
		grass_hist = grass_hist + pa_hist
	}
	
	if (i == 16) {
		crop_hist = crop_hist + pa_hist
	}
	
	if (i == 1) {
		bare_hist = bare_hist + pa_hist
	}
	
}
	
for (i in 1:num_harv) {
	for(t in 1: nyear_hist) {
		harvest_area_hist[t] = harvest_area_hist[t] + sum(harvest_frac_hist[,,i,tinds_hist[t]] * veg_area_hist * forest_frac_hist[,,t], na.rm=TRUE)
		#harvest_area_hist[t] = harvest_area_hist[t] + sum(harvest_frac_hist[,,i,t] * h0_veg_area, na.rm=TRUE)
	}
}

plot(year_hist,forest_hist,main="forest km^2", type="l", sub=paste(year_hist[nyear_hist],forest_hist[nyear_hist]))
plot(year_hist,shrub_hist,main="shrub km^2", type="l", sub=paste(year_hist[nyear_hist],shrub_hist[nyear_hist]))
plot(year_hist,grass_hist,main="grass km^2", type="l", sub=paste(year_hist[nyear_hist],grass_hist[nyear_hist]))
#plot(year_hist,crop_hist,main="crop km^2", type="l", sub=paste(year_hist[nyear_hist],crop_hist[nyear_hist]))
plot(year_hist,bare_hist,main="bare km^2", type="l", sub=paste(year_hist[nyear_hist],bare_hist[nyear_hist]))

plot(year_hist,harvest_area_hist,main="harvest area km^2", type="l", sub=paste(year_hist[nyear_hist],harvest_area_hist[nyear_hist]))

dev.off()



if(FALSE){

##### plot the historical and future in the same time series

#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1feb2023/check_all.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/check_all_single_2deg.pdf")
#pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_zlnd_1deg_14feb2023/check_all.pdf")
pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_20231010_zlnd_gcam6_SSP2/check_hist_future_surfdata_iESM_dyn.pdf")

plot(c(year_hist, year_future[2:nyear_future]), c(crop_hist, crop_future[2:nyear_future]), main="crop km^2", type="n",
	sub=paste(year_hist[nyear_hist], "hist =", crop_hist[nyear_hist], "future =", crop_future[1]))
lines(year_hist, crop_hist, type="l")
lines(year_future, crop_future, type="l")

plot(c(year_hist, year_future[2:nyear_future]), c(forest_hist, forest_future[2:nyear_future]), main="forest km^2", type="n",
	sub=paste(year_hist[nyear_hist], "hist =", forest_hist[nyear_hist], "future =", forest_future[1]))
lines(year_hist, forest_hist, type="l")
lines(year_future, forest_future, type="l")

plot(c(year_hist, year_future[2:nyear_future]), c(shrub_hist, shrub_future[2:nyear_future]), main="shrub km^2", type="n",
	sub=paste(year_hist[nyear_hist], "hist =", shrub_hist[nyear_hist], "future =", shrub_future[1]))
lines(year_hist, shrub_hist, type="l")
lines(year_future,shrub_future, type="l", sub=paste(year_future[1],shrub_future[1]))

plot(c(year_hist, year_future[2:nyear_future]), c(grass_hist, grass_future[2:nyear_future]), main="grass km^2", type="n",
	sub=paste(year_hist[nyear_hist], "hist =", grass_hist[nyear_hist], "future =", grass_future[1]))
lines(year_hist, grass_hist, type="l")
lines(year_future,grass_future, type="l")

plot(c(year_hist, year_future[2:nyear_future]), c(bare_hist, bare_future[2:nyear_future]), main="bare km^2", type="n",
	sub=paste(year_hist[nyear_hist], "hist =", bare_hist[nyear_hist], "future =", bare_future[1]))
lines(year_hist, bare_hist, type="l")
lines(year_future, bare_future, type="l")

plot(c(year_hist, year_future[2:nyear_future]), c(harvest_area_hist, harvest_area_future[2:nyear_future]), main="harvest area km^2", type="n",
	sub=paste(year_hist[nyear_hist], "hist =", harvest_area_hist[nyear_hist], "future =", harvest_area_future[1]))
lines(year_hist, harvest_area_hist, type="l")
lines(year_future, harvest_area_future, type="l")

dev.off()


################ this is for another specific file, but generalized the code above

s2r45 = nc_open("~/projects/e3sm/giac_v2/check_previous_version/landuse.timeseries_ne30np4_SSP2_RCP45_simyr2015-2100_200728.nc")

pfts = ncvar_get(s2r45, varid="PCT_NAT_PFT")
area = ncvar_get(s2r45, varid="AREA")
lf = ncvar_get(s2r45, varid="LANDFRAC_PFT")
mask = ncvar_get(s2r45, varid="PFTDATA_MASK")
veg = ncvar_get(s2r45, varid="PCT_NATVEG")
year = ncvar_get(s2r45, varid="YEAR")


nc_close(s2r45)

nyear=length(year)

pdf("~/projects/e3sm/diagnostics/check_ssp2_rcp45.pdf")
pa = array(dim=c(nyear))
forest = array(dim=c(nyear))
shrub = array(dim=c(nyear))
grass = array(dim=c(nyear))
forest[] = 0
shrub[] = 0
grass[] = 0
for(i in 1:npft) {
	pa[] = 0
	veg_area = veg / 100 * lf * area * mask
	for(t in 1: nyear) {
		pa[t] = sum(pfts[,i,t] / 100 * veg_area, na.rm=TRUE)
	}
	plot(year,pa,main=pftnames[i], type="l", sub=paste(year[1],pa[1]))
	
	if (i >=2 & i <=9) {
		forest = forest + pa
	}
	
	if (i >=10 & i <=12) {
		shrub = shrub + pa
	}
	
	if (i >=13 & i <=15) {
		grass = grass + pa
	}
	
}

plot(year,forest,main="forest", type="l", sub=paste(year[1],forest[1]))
plot(year,shrub,main="shrub", type="l", sub=paste(year[1],shrub[1]))
plot(year,grass,main="grass", type="l", sub=paste(year[1],grass[1]))


dev.off()


} # end if(FALSE)


##############################################################################
# this is for comparing the same time series

####

futlab = "with feedbacks"
futcol = "with_fdbks"

histlab = "no feedbacks"


# this truncates the 'future' data to the length of the 'hist' data, or vice-versa
if(nyear_future >= nyear_hist) {
	nyear_future=nyear_hist
	year_future=year_future[1:nyear_future]
	pa_future = pa_future[1:nyear_future]
	forest_future = forest_future[1:nyear_future]
	shrub_future = shrub_future[1:nyear_future]
	grass_future = grass_future[1:nyear_future]
	crop_future = crop_future[1:nyear_future]
	bare_future = bare_future[1:nyear_future]
	harvest_area_future = harvest_area_future[1:nyear_future]
} else {
	nyear_hist=nyear_future
	year_hist=year_hist[1:nyear_hist]
	pa_hist = pa_hist[1:nyear_hist]
	forest_hist = forest_hist[1:nyear_hist]
	shrub_hist = shrub_hist[1:nyear_hist]
	grass_hist = grass_hist[1:nyear_hist]
	crop_hist = crop_hist[1:nyear_hist]
	bare_hist = bare_hist[1:nyear_hist]
	harvest_area_hist = harvest_area_hist[1:nyear_hist]
}

pdf("~/projects/e3sm/iesmv2/test_diagnostics/eva_20240730_SSP245_ZATM/compare_surfdata_iESM_dyn_fdbk_nofdbk.pdf")

# this enables the plot below even if the above plots have been run
forest_future[] = 0
shrub_future[] = 0
grass_future[] = 0
crop_future[] = 0
bare_future[] = 0
harvest_area_future[] = 0
harvest_area_temp_future = harvest_area_future
forest_hist[] = 0
shrub_hist[] = 0
grass_hist[] = 0
crop_hist[] = 0
bare_hist[] = 0
harvest_area_hist[] = 0
harvest_area_temp_hist = harvest_area_hist

for(i in 1: npft) {
	pa_future[] = 0
	pa_hist[] = 0
	for(t in 1: nyear_future) {
		pa_future[t] = sum(pfts_future[,,i,tinds_future[t]] / 100 * veg_area_future, na.rm=TRUE)
		pa_hist[t] = sum(pfts_hist[,,i,tinds_hist[t]] / 100 * veg_area_hist, na.rm=TRUE)
	}
	plot(year_future,pa_future,main=pftnames[i], type="n", sub=paste(year_future[nyear_future], futlab, pa_future[nyear_future], histlab, pa_hist[nyear_future]),
	     xlab="year", ylab="global km^2",
	     ylim=c(min(c(pa_future,pa_hist),na.rm=TRUE),max(c(pa_future,pa_hist),na.rm=TRUE)))
	#lines(year_future, pa_future, type="l", col="blue", lty="solid")     
	#lines(year_hist,pa_hist, type="l", col="red", lty="solid")
	#legend(x="bottom", legend = c(futlab, histlab), col=c("blue", "red"), lty=c("solid", "solid"))
	lines(year_future, pa_future, type="l", col="black", lty="solid")     
	lines(year_hist,pa_hist, type="l", col="orange", lty="solid")
	legend(x="bottom", legend = c(futlab, histlab), col=c("black", "orange"), lty=c("solid", "solid"))
	
	diff = pa_future - pa_hist
	diff_percent = 100 * (diff / pa_hist)
	cat("pft", pftnames[i], "max ( min ) diff area (", futlab, "-", histlab, ") is", max(diff, na.rm=TRUE), "(", min(diff, na.rm=TRUE), ") km^2\n")
	cat("pft", pftnames[i], "max ( min )", futlab, "percent diff area from", histlab, "is", max(diff_percent, na.rm=TRUE), "(", min(diff_percent, na.rm=TRUE), ")\n\n")
	
	if (i >=2 & i <=9) {
		forest_future = forest_future + pa_future
		forest_hist = forest_hist + pa_hist
	}
	
	if (i >=10 & i <=12) {
		shrub_future = shrub_future + pa_future
		shrub_hist = shrub_hist + pa_hist
	}
	
	if (i >=13 & i <=15) {
		grass_future = grass_future + pa_future
		grass_hist = grass_hist + pa_hist
	}
	
	if (i == 16) {
		crop_future = crop_future + pa_future
		crop_hist = crop_hist + pa_hist
	}
	
	if (i == 1) {
		bare_future = bare_future + pa_future
		bare_hist = bare_hist + pa_hist
	}
	
}

for (i in 1:num_harv) {
	harvest_area_temp_hist[] = 0
	harvest_area_temp_future[] = 0
	for(t in 1: nyear_future) {
		harvest_area_temp_future[t] = sum(harvest_frac_future[,,i,tinds_future[t]] * veg_area_future * forest_frac_future[,,t], na.rm=TRUE)
		harvest_area_future[t] = harvest_area_future[t] + harvest_area_temp_future[t]
		harvest_area_temp_hist[t] = sum(harvest_frac_hist[,,i,tinds_hist[t]] * veg_area_hist * forest_frac_hist[,,t], na.rm=TRUE)
		harvest_area_hist[t] = harvest_area_hist[t] + harvest_area_temp_hist[t]
		cat(futlab, harvest_names[i],"year", t, "area is", sum(harvest_frac_future[,,i,tinds_future[t]] * veg_area_future, na.rm=TRUE), "km^2 \n")
		cat(histlab, harvest_names[i],"year", t, "area is", sum(harvest_frac_hist[,,i,tinds_hist[t]] * veg_area_hist, na.rm=TRUE), "km^2 \n")
		
		diff = harvest_frac_future[,,i,tinds_future[t]] - harvest_frac_hist[,,i,tinds_hist[t]]
		diff_percent = 100 * (diff / harvest_frac_hist[,,i,tinds_hist[t]])
		cat(harvest_names[i], futlab, "-", histlab,"year", t, "max frac diff is", max(diff,na.rm=TRUE), "\n")
		cat(harvest_names[i], futlab, "-", histlab,"year", t, "max frac diff percent is", max(diff_percent,na.rm=TRUE), "\n")
		cat(harvest_names[i], futlab, "-", histlab,"year", t, "min frac diff is", min(diff,na.rm=TRUE), "\n")
		cat(harvest_names[i], futlab, "-", histlab,"year", t, "min frac diff percent is", min(diff_percent,na.rm=TRUE), "\n")
	}
	plot(year_future, harvest_area_temp_future,main=harvest_names[i], type="n", sub=paste(year_future[nyear_future], futlab, harvest_area_temp_future[nyear_future],
	     histlab, harvest_area_temp_hist[nyear_future]), xlab="year", ylab="global km^2",
	     ylim=c(min(c(harvest_area_temp_future,harvest_area_temp_hist),na.rm=TRUE),max(c(harvest_area_temp_future,harvest_area_temp_hist),na.rm=TRUE)))
	#lines(year_future, harvest_area_temp_future, type="l", col="blue", lty="solid")
	#lines(year_hist, harvest_area_temp_hist, type="l", col="red", lty="solid")
	#legend(x="bottom", legend = c(futlab, histlab), col=c("blue", "red"), lty=c("solid", "solid"))
	lines(year_future, harvest_area_temp_future, type="l", col="black", lty="solid")
	lines(year_hist, harvest_area_temp_hist, type="l", col="orange", lty="solid")
	legend(x="bottom", legend = c(futlab, histlab), col=c("black", "orange"), lty=c("solid", "solid"))
	
	cat("hafut", harvest_area_temp_future, "\nhahist", harvest_area_temp_hist,"\n")
	
	diff = harvest_area_temp_future - harvest_area_temp_hist
	diff_percent = 100 * (diff / harvest_area_temp_hist)
	cat(harvest_names[i], "max ( min ) diff area (", futlab, "-", histlab, ") is", max(diff, na.rm=TRUE), "(", min(diff, na.rm=TRUE), ") km^2\n")
	cat(harvest_names[i], "max ( min )", futlab, "percent diff area from", histlab, "is", max(diff_percent, na.rm=TRUE), "(", min(diff_percent, na.rm=TRUE), ")\n\n")
}

plot(year_future,forest_future,main="forest km^2", type="n", sub=paste(year_future[nyear_future], futlab, forest_future[nyear_future],
     histlab, forest_hist[nyear_future]), xlab="year", ylab="global km^2",
     ylim=c(min(c(forest_future,forest_hist),na.rm=TRUE),max(c(forest_future,forest_hist),na.rm=TRUE)))
#lines(year_future, forest_future, type="l", col="blue", lty="solid")
#lines(year_hist, forest_hist, type="l", col="red", lty="solid")
#legend(x="bottom", legend = c(futlab, histlab), col=c("blue", "red"), lty=c("solid", "solid"))
lines(year_future, forest_future, type="l", col="black", lty="solid")
lines(year_hist, forest_hist, type="l", col="orange", lty="solid")
legend(x="bottom", legend = c(futlab, histlab), col=c("black", "orange"), lty=c("solid", "solid"))
diff = forest_future - forest_hist
diff_percent = 100 * (diff / forest_hist)
cat("forest max ( min ) diff area (", futlab, "-", histlab, ") is", max(diff, na.rm=TRUE), "(", min(diff, na.rm=TRUE), ") km^2\n")
cat("forest max ( min )", futlab, "percent diff area from", histlab, "is", max(diff_percent, na.rm=TRUE), "(", min(diff_percent, na.rm=TRUE), ")\n\n")

plot(year_future,shrub_future,main="shrub km^2", type="n", sub=paste(year_future[nyear_future], futlab, shrub_future[nyear_future],
     histlab, shrub_hist[nyear_future]), xlab="year", ylab="global km^2",
     ylim=c(min(c(shrub_future,shrub_hist),na.rm=TRUE),max(c(shrub_future,shrub_hist),na.rm=TRUE)))
#lines(year_future, shrub_future, type="l", col="blue", lty="solid")
#lines(year_hist, shrub_hist, type="l", col="red", lty="solid")
#legend(x="bottom", legend = c(futlab, histlab), col=c("blue", "red"), lty=c("solid", "solid"))
lines(year_future, shrub_future, type="l", col="black", lty="solid")
lines(year_hist, shrub_hist, type="l", col="orange", lty="solid")
legend(x="bottom", legend = c(futlab, histlab), col=c("black", "orange"), lty=c("solid", "solid"))
diff = shrub_future - shrub_hist
diff_percent = 100 * (diff / shrub_hist)
cat("shrub max ( min ) diff area (", futlab, "-", histlab, ") is", max(diff, na.rm=TRUE), "(", min(diff, na.rm=TRUE), ") km^2\n")
cat("shrub max ( min )", futlab, "percent diff area from", histlab, "is", max(diff_percent, na.rm=TRUE), "(", min(diff_percent, na.rm=TRUE), ")\n\n")

plot(year_future,grass_future,main="grass km^2", type="n", sub=paste(year_future[nyear_future], futlab, grass_future[nyear_future],
     histlab, grass_hist[nyear_future]), xlab="year", ylab="global km^2",
     ylim=c(min(c(grass_future,grass_hist),na.rm=TRUE),max(c(grass_future,grass_hist),na.rm=TRUE)))
#lines(year_future, grass_future, type="l", col="blue", lty="solid")
#lines(year_hist, grass_hist, type="l", col="red", lty="solid")
#legend(x="bottom", legend = c(futlab, histlab), col=c("blue", "red"), lty=c("solid", "solid"))
lines(year_future, grass_future, type="l", col="black", lty="solid")
lines(year_hist, grass_hist, type="l", col="orange", lty="solid")
legend(x="bottom", legend = c(futlab, histlab), col=c("black", "orange"), lty=c("solid", "solid"))
diff = grass_future - grass_hist
diff_percent = 100 * (diff / grass_hist)
cat("grass max ( min ) diff area (", futlab, "-", histlab, ") is", max(diff, na.rm=TRUE), "(", min(diff, na.rm=TRUE), ") km^2\n")
cat("grass max ( min )", futlab, "percent diff area from", histlab, "is", max(diff_percent, na.rm=TRUE), "(", min(diff_percent, na.rm=TRUE), ")\n\n")


plot(year_future,crop_future,main="crop km^2", type="n", sub=paste(year_future[nyear_future], futlab, crop_future[nyear_future],
     histlab, crop_hist[nyear_future]), xlab="year", ylab="global km^2",
     ylim=c(min(c(crop_future,crop_hist),na.rm=TRUE),max(c(crop_future,crop_hist),na.rm=TRUE)))
#lines(year_future, crop_future, type="l", col="blue", lty="solid")
#lines(year_hist, crop_hist, type="l", col="red", lty="solid")
#legend(x="bottom", legend = c(futlab, histlab), col=c("blue", "red"), lty=c("solid", "sollid"))
lines(year_future, crop_future, type="l", col="black", lty="solid")
lines(year_hist, crop_hist, type="l", col="orange", lty="solid")
legend(x="bottom", legend = c(futlab, histlab), col=c("black", "orange"), lty=c("solid", "solid"))
diff = crop_future - crop_hist
diff_percent = 100 * (diff / crop_hist)
cat("crop max ( min ) diff area (", futlab, "-", histlab, ") is", max(diff, na.rm=TRUE), "(", min(diff, na.rm=TRUE), ") km^2\n")
cat("crop max ( min )", futlab, "percent diff area from", histlab, "is", max(diff_percent, na.rm=TRUE), "(", min(diff_percent, na.rm=TRUE), ")\n\n")

plot(year_future,bare_future,main="bare km^2", type="n", sub=paste(year_future[nyear_future], futlab, bare_future[nyear_future],
     histlab, bare_hist[nyear_future]), xlab="year", ylab="global km^2",
     ylim=c(min(c(bare_future,bare_hist),na.rm=TRUE),max(c(bare_future,bare_hist),na.rm=TRUE)))
#lines(year_future, bare_future, type="l", col="blue", lty="solid")
#lines(year_hist, bare_hist, type="l", col="red", lty="solid")
#legend(x="bottom", legend = c(futlab, histlab), col=c("blue", "red"), lty=c("solid", "solid"))
lines(year_future, bare_future, type="l", col="black", lty="solid")
lines(year_hist, bare_hist, type="l", col="orange", lty="solid")
legend(x="bottom", legend = c(futlab, histlab), col=c("black", "orange"), lty=c("solid", "solid"))
diff = bare_future - bare_hist
diff_percent = 100 * (diff / bare_hist)
cat("bare max ( min ) diff area (", futlab, "-", histlab, ") is", max(diff, na.rm=TRUE), "(", min(diff, na.rm=TRUE), ") km^2\n")
cat("bare max ( min )", futlab, "percent diff area from", histlab, "is", max(diff_percent, na.rm=TRUE), "(", min(diff_percent, na.rm=TRUE), ")\n\n")

plot(year_future, harvest_area_future,main="total forest harvest km^2", type="n", sub=paste(year_future[nyear_future], futlab, harvest_area_future[nyear_future],
     histlab, harvest_area_hist[nyear_future]), xlab="year", ylab="global km^2",
     ylim=c(min(c(harvest_area_future,harvest_area_hist),na.rm=TRUE),max(c(harvest_area_future,harvest_area_hist),na.rm=TRUE)))
#lines(year_future, harvest_area_future, type="l", col="blue", lty="solid")
#lines(year_hist, harvest_area_hist, type="l", col="red", lty="solid")
#legend(x="bottom", legend = c(futlab, histlab), col=c("blue", "red"), lty=c("solid", "solid"))
lines(year_future, harvest_area_future, type="l", col="black", lty="solid")
lines(year_hist, harvest_area_hist, type="l", col="orange", lty="solid")
legend(x="bottom", legend = c(futlab, histlab), col=c("black", "orange"), lty=c("solid", "solid"))
diff = harvest_area_future - harvest_area_hist
diff_percent = 100 * (diff / harvest_area_hist)
cat("total harvest max ( min ) diff area (", futlab, "-", histlab, ") is", max(diff, na.rm=TRUE), "(", min(diff, na.rm=TRUE), ") km^2\n")
cat("total harvest max ( min )", futlab, "percent diff area from", histlab, "is", max(diff_percent, na.rm=TRUE), "(", min(diff_percent, na.rm=TRUE), ")\n\n")

dev.off()

# write the aggregate time series data

out_agg_data = data.frame(year=year_future, forest_future, forest_hist, shrub_future, shrub_hist, grass_future, grass_hist, crop_future,
							crop_hist, bare_future, bare_hist, harvest_area_future, harvest_area_hist)
names(out_agg_data) = c("year", paste0("forest_", futcol), "forest", paste0("shrub_", futcol), "shrub", paste0("grass_", futcol), "grass", paste0("crop_", futcol), "crop",
                         paste0("bare_", futcol), "bare", paste0("harvest_", futcol), "harvest")
out_agg_data$units="area_km2"

write.csv(out_agg_data, file="~/projects/e3sm/iesmv2/test_diagnostics/eva_20240730_SSP245_ZATM/compare_surfdata_iESM_dyn_fdbk_nofdbk.csv")

