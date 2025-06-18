# example of extracting GCAM output data from xml file
# see plot_gcam_prices.r in the metarepo for additional examples

# first install the rgcam library (https://github.com/JGCRI/rgcam)
# you can do this locally or on the HPC
# in R: require(devtools)
#       install_github('JGCRI/rgcam', build_vignettes=TRUE)

library(rgcam)

# a query file is required to extract data from the xml output
#    this file defines which data to extract
qfile = "land_allocation_yieldscaler_cdensity_co2price_query.xml"

# set working directory and include paths in names below as needed
gcam_xml_path = "./"
gcam_xml_out = "GCAMDBOutput_E3SM_gcam.xml"
# you can name this whatever you want
gcam_proj_file = "s2r4p5_e3sm.proj"
outdir = "./"

# create the project data file

# open the xml file - java needs lots of memory
conn <- localDBConn(dbPath=gcam_xml_path, dbFile=gcam_xml_out, maxMemory="32g")
# extract the scenario (assume there is only one so don't have to name it) and create the project data file
#    but scenario name is defined in gcam config file by, e.g.: <Value name="scenarioName">4p5</Value>

# the file argument is taken as a literal string now, so:
# clobber=TRUE overwrites the existing file
prj <- addScenario(conn, gcam_proj_file, scenario=NULL, queryFile=readChar(qfile, file.info(qfile)$size), clobber=TRUE)

# if the project data file already exists, just read it in:
# prj <- loadProject(gcam_proj_file)

# show summary of data in project file
str(prj)


# plot co2 emissions by region

gcam_co2_reg <-getQuery(prj0, query="CO2 emissions by region")
co2_reg = as.data.frame(gcam_co2_reg)
co2_reg$scenario = "gcam_example"

# aggregate to globe and then add the global row
temp_globe = aggregate(value ~ Units + scenario + year, data = co2_reg, FUN = sum)
region = rep("Global", nrow(temp_globe))
temp_globe = cbind(region, temp_globe, stringsAsFactors=FALSE)
co2_reg = rbind(co2_reg, temp_globe, stringsAsFactors=FALSE, make.row.names = FALSE)

pdf(file=paste0(outdir,"gcam_co2_reg_emissions",".pdf"))

for( r in c(unique(co2_reg$region)) ){
	
	plot_data = mc[mc$year >= 2015 & mc$region == r,]
	
	# plot the data

	p1 <- ggplot(data=plot_data) +
		geom_line(aes(x=year, y=value.x, color="gcam_example"), linewidth=1) +
		labs(y=paste("MT CO2"), title=paste(r,"CO2 emissions")) +
		theme_bw() +
		scale_color_manual(values=c( gcam_example = "orange"))
	
	print(p1)

} # end r loop over regions

dev.off()

write.csv(mc, file=paste0(outdir,"gcam_co2_reg_emissions",".csv"))


