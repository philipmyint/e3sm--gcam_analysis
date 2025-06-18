require(devtools)
library(dplyr)
library(rgcam)
library(rjson)

# Create a vector of command line arguments, each specifying a JSON file with user inputs.
args = commandArgs(trailingOnly=TRUE)

# Stop the script execution if no JSON files were specified.
num_json_files = length(args)
if(num_json_files == 0)
{
    stop("Usage: Rscript gcam_extract_csv_from_project_files.r `path/to/json/input/file(s)'")
    geterrmessage()
}

# Read all JSON files entered on the command line and produce a .csv file for each block of each JSON file.
for (json_index in 1:num_json_files)
{
    # Read the current JSON file and put the inputs specified in each block of the file into a list.
    list_of_inputs = fromJSON(file=args[json_index])
    for (inputs_index in 1:length(list_of_inputs))
    {
        inputs = list_of_inputs[[inputs_index]]
        # Here, variable refers to the specific quantity we want to extract from the project files produced from the GCAM-generated xml output files.
        variable = inputs$variable

        # Each scenario for the current block will be read in as a dataframe from the project file associated with the scenario.
        # Load this dataframe into a list, then combine all the dataframes in the list into a single dataframe containing data on all the scenarios.
        list_of_dataframes = list()
        for (scenario_index in 1:length(inputs$scenarios))
        {
            scenario = inputs$scenarios[scenario_index]
            project_file = inputs$project_files[scenario_index]

            prj = loadProject(project_file)
            data = getQuery(prj, query=variable)
            df = as.data.frame(data)
            df$scenario = scenario
            list_of_dataframes[[scenario_index]] = df
        }
        df = bind_rows(list_of_dataframes)

        # Write the dataframe to a output file in .csv format.
        output_file = inputs$output_file
        if (endsWith(output_file, ".csv"))
        {
            write.csv(df, output_file)
        }
        else
        {
            write.csv(df, paste0(output_file, ".csv"))
        }
    }
}