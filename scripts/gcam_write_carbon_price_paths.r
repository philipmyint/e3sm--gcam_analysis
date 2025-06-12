# write_carbon_price_paths.r

# extract the co2 price info from the GCAM xml output and put it into two carbon tax xml gcam configuration files:
#    carbon_tax_<scen>_nearterm.xml
#    carbon_tax_<scen>_LTG_high.xml
# value units are 1990$/tC

# actually the near-term values are not adjusted, so this file does not need to be created
# only the global carbon tax trajectory is calculated by the target finder
# the near-term and long-term taxes are complementary and should not have overlapping years,
#    but the phase-in file determines which values to use in the sim

# these are the near and long term carbon prices for a scenario when running e3sm-gcam
# in the gcam config file:
#		<Value name = "long-term-co2">../input/policy/carbon_tax_<scen>_LTG_high.xml
#		<Value name = "near-term-co2">../input/policy/carbon_tax_<scen>_nearterm.xml

# note that this phase in file does not change:
#		<Value name = "co2-link">../input/policy/2025_target_finder_phasein.xml</Value>
# and has to match the policy file in terms of start year for carbon tax. e.g.:
#       <Value name="policy-target-file">../input/policy/forcing_target_4p5.xml</Value>
#       has a 2025 start year for carbon pricing

# also note that when running the target finder mode, this tax file has initial values of zero:
#       <Value name = "long-term-co2">../input/policy/carbon_tax_0.xml</Value>

# the near term market is regional and is not recalculated
# the long term market is global and is filled out for just the first region and then copied to the rest

# prior to loading R on compy:
#    install the java jars:
#       download jars.zip from https://github.com/JGCRI/modelinterface/releases
#       move the zip file to ${HOME}/libs then unzip
#    cd <e3sm-rundir> (or where the .xml file is)
#    module load java/1.8.0_31
#    cp <clone>/components/gcam/src/iac/gcam/exe/XMLDBDriver.* .
#    export CLASSPATH='${$HOME}/libs/jars/*:XMLDBDriver.jar'
# i don't think the following is necessary for two reasons:
#    rgcam can read the xml file
#    when trying to create the basexdb it expects there to be an .../output/inf.basex file, but no output folder or database exists for e3sm-gcam runs
#    java XMLDBDriver --db-path=../output/database_basexdb_<scen> --doc-name=NAME --xml-GCAMDBOutput_<scen>.xml

# set the working directory to the run directory (or where the .xml file is)

# this script assumes that there is only one scenario in the gcam output

# Note that the proj file is much, much small than the xml
# to transfer data off the machines, it is better to create the proj data file first, then transfer it instead of the xml file

# this is diagnostic as it is not necessary to create a new near term file
NEARTERM = FALSE

######## set the gcam xml output file name and the desired project data file name

gcam_xml_out="/Users/adivi/projects/e3sm/iesmv2/inputs/carbon_price_paths/ssp2rcp45/GCAMDBOutput_SSP2_4p5_for_E3SM.xml"
gcam_proj_file="s2r4p5_e3sm.proj"
data_dir = "ssp2rcp45"

# set this to true if the project data file needs to be created
CREATE_PROJ_FILE = FALSE

########

# this gets some other info as well
qfile = "land_allocation_yieldscaler_cdensity_co2price_query.xml"

nt_outfile_base = "carbon_tax_nearterm_"
lt_outfile_base = "carbon_tax_LTG_"

# rather than use template files, set up the required text here
xmlversion_text = '<?xml version="1.0" encoding="UTF-8"?>'
xml_struct = c("scenario", "world", "region", "ghgpolicy")
xml_ghgpolicy = c("market", "fixedTax")
name_tag = " name="
year_tag = " year="
xml_ext=".xml"
num_reg=32
global_tag = "global"
nt_year = 2020
# long term years are 2025 - 2100 in 5-year intervals
lt_years = rep(2025, 16)
for(i in 1:length(lt_years)){
   lt_years[i] = lt_years[i] + (i-1)*5
}
tab="\t"

# note that currently near term values are by region and long term values are global
near_term_tag = "CO2_NearTerm"
long_term_tag = "CO2_LTG"

require(devtools)
# if rgcam is not installed
# install_github('JGCRI/rgcam', build_vignettes=TRUE)
library(rgcam)

###### if there is no proj data file yet:

if(CREATE_PROJ_FILE) {

	# open the xml file - java needs lots of memory
	conn <- localDBConn(dbPath= '', dbFile=gcam_xml_out, maxMemory="32g")
	# extract the scenario (assume there is only one so don't have to name it) and create the project data file
	#    but scenario name is defined in gcam config file by, e.g.: <Value name="scenarioName">4p5</Value>
	#prj <- addScenario(conn, gcam_proj_file, scenario=NULL, queryFile=qfile, clobber=TRUE)
	# the file argument is taken as a literal string now, so:
	prj <- addScenario(conn, gcam_proj_file, scenario=NULL, queryFile=readChar(qfile, file.info(qfile)$size), clobber=TRUE)

} else {

	##### if the proj file has been created do this:
	prj <- loadProject(gcam_proj_file)
	#####

}


# get the tibble of co2 price data (the query name is defined in the qfile)
co2<-getQuery(prj, query="CO2 prices")
scen_name = co2$scenario[1]

# take the global market out of this region vector
regions = unique(sub("CO2.*", "", co2$market))[1:num_reg]

# start the xml output files
# add the double quotes as necessary

if(NEARTERM){
   nt_outfile = paste0(data_dir,"/",nt_outfile_base,scen_name,xml_ext)
   cat(xmlversion_text,"\n", file= nt_outfile)
   cat(paste0("<",xml_struct[1],name_tag,'"',scen_name,'">\n'), file= nt_outfile, append=TRUE)
   cat(paste0(tab,"<",xml_struct[2],">\n"), file= nt_outfile, append=TRUE)
}

lt_outfile = paste0(data_dir,"/",lt_outfile_base,scen_name,xml_ext)

cat(xmlversion_text,"\n", file= lt_outfile)
cat(paste0("<",xml_struct[1],name_tag,'"',scen_name,'">\n'), file= lt_outfile, append=TRUE)
cat(paste0(tab,"<",xml_struct[2],">\n"), file= lt_outfile, append=TRUE)

# loop over regions to fill out the carbon tax values
for(r in 1:num_reg){

    if(NEARTERM){	
	   # near term
	   cat(paste0(tab,tab,"<",xml_struct[3],name_tag,'"',regions[r],'">\n'), file= nt_outfile, append=TRUE)
	   cat(paste0(tab,tab,tab,"<",xml_struct[4],name_tag,'"',near_term_tag,'">\n'), file= nt_outfile, append=TRUE)
	   cat(paste0(tab,tab,tab,tab,"<", xml_ghgpolicy[1],">",regions[r],"</",xml_ghgpolicy[1],">\n"), file= nt_outfile, append=TRUE)
	   nt_value=format( round(co2$value[co2$market==paste0(regions[r],near_term_tag) & co2$year==nt_year],2), nsmall=2 )
	   cat(paste0(tab,tab,tab,tab,"<", xml_ghgpolicy[2],year_tag,'"',nt_year,'">',nt_value,"</",xml_ghgpolicy[2],">\n"), file= nt_outfile, append=TRUE)
	   cat(paste0(tab,tab,tab,"</",xml_struct[4],">\n"), file= nt_outfile, append=TRUE)
	   cat(paste0(tab,tab,"</",xml_struct[3],">\n"), file= nt_outfile, append=TRUE)
	}
	
	# long term - global only
	cat(paste0(tab,tab,"<",xml_struct[3],name_tag,'"',regions[r],'">\n'), file= lt_outfile, append=TRUE)
	cat(paste0(tab,tab,tab,"<",xml_struct[4],name_tag,'"',long_term_tag,'">\n'), file= lt_outfile, append=TRUE)
	cat(paste0(tab,tab,tab,tab,"<", xml_ghgpolicy[1],">",global_tag,"</",xml_ghgpolicy[1],">\n"), file= lt_outfile, append=TRUE)
	if(r==1){
		for(l in lt_years){
			lt_value=format( round(co2$value[co2$market==paste0(global_tag,long_term_tag) & co2$year==l],3), nsmall=3)
			cat(paste0(tab,tab,tab,tab,"<", xml_ghgpolicy[2],year_tag,'"',l,'">',lt_value,"</",xml_ghgpolicy[2],">\n"), file= lt_outfile, append=TRUE)
		}
	}
	cat(paste0(tab,tab,tab,"</",xml_struct[4],">\n"), file= lt_outfile, append=TRUE)
	cat(paste0(tab,tab,"</",xml_struct[3],">\n"), file= lt_outfile, append=TRUE)
	
	
} # end r loop over regions

# close world and scenario

if(NEARTERM){
   cat(paste0(tab,"</",xml_struct[2],">\n"), file= nt_outfile, append=TRUE)
   cat(paste0("</",xml_struct[1],">\n"), file= nt_outfile, append=TRUE)
}

cat(paste0(tab,"</",xml_struct[2],">\n"), file= lt_outfile, append=TRUE)
cat(paste0("</",xml_struct[1],">\n"), file= lt_outfile, append=TRUE)








