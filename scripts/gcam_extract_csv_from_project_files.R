require(devtools)
library(dplyr)
library(rgcam)
library(rjson)

# Create a vector of command line arguments, each specifying a JSON file with user inputs.
args = commandArgs(trailingOnly=TRUE)

# Stop the script execution if no JSON files were specified.
numJSON = length(args)
if(numJSON == 0)
{
    stop("Usage: Rscript gcam_extract_csv_from_project_files.R `path/to/json/input/file(s)'")
    geterrmessage()
}

# Read all JSON files entered on the command line and produce a .csv file for each block of each JSON file.
for (indexJSON in 1:numJSON)
{
    # Read the current JSON file and put the inputs specified in each block of the file into a list.
    listOfInputs = fromJSON(file=args[indexJSON])
    for (inputsIndex in 1:length(listOfInputs))
    {
        inputs = listOfInputs[[inputsIndex]]
        # Here, variable refers to the specific quantity we want to extract from the project files produced from the GCAM-generated XML output files.
        variable = inputs$variable

        # Each scenario for the current block will be read in as a dataframe from the project file associated with the scenario.
        # Load this dataframe into a list, then combine all the dataframes in the list into a single dataframe containing data on all the scenarios.
        listOfDataframes = list()
        for (indexScenario in 1:length(inputs$scenarios))
        {
            scenario = inputs$scenarios[indexScenario]
            projectFile = inputs$projectFiles[indexScenario]

            project = loadProject(projectFile)
            data = getQuery(project, query=variable)
            df = as.data.frame(data)
            df$scenario = scenario
            listOfDataframes[[indexScenario]] = df
        }
        df = bind_rows(listOfDataframes)

        # Write the dataframe to a output file in .csv format.
        outputFile = inputs$outputFile
        if (endsWith(outputFile, ".csv"))
        {
            write.csv(df, outputFile)
        }
        else
        {
            write.csv(df, paste0(outputFile, ".csv"))
        }
    }
}