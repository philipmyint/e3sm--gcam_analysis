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
}

# Cache for loaded project files to avoid loading the same file multiple times across all variables.
# Each unique project file is loaded at most once regardless of how many variables reference it.
projectCache = list()

# Read all JSON files entered on the command line and produce a .csv file for each block of each JSON file.
for (indexJSON in 1:numJSON)
{
    # Check that the JSON file exists before trying to read it.
    if (!file.exists(args[indexJSON]))
    {
        stop(paste0("Error: JSON file not found: ", args[indexJSON]))
    }

    # Read the current JSON file and put the inputs specified in each block of the file into a list.
    listOfInputs = fromJSON(file=args[indexJSON])
    for (indexInputs in 1:length(listOfInputs))
    {
        inputs = listOfInputs[[indexInputs]]
        # Here, variable refers to the specific quantity we want to extract from the project files produced from the GCAM-generated XML output files.
        variable = inputs$variable

        # Check that scenarios and projectFiles arrays are the same length.
        if (length(inputs$scenarios) != length(inputs$projectFiles))
        {
            stop(paste0("Error: 'scenarios' (length ", length(inputs$scenarios), ") and 'projectFiles' (length ",
                        length(inputs$projectFiles), ") must be the same length for variable '", variable, "'."))
        }

        # Check that the output directory exists before attempting to write.
        outputFile = inputs$outputFile
        outputDir = dirname(outputFile)
        if (!dir.exists(outputDir))
        {
            stop(paste0("Error: output directory does not exist: ", outputDir))
        }

        # Each scenario for the current block will be read in as a dataframe from the project file associated with the scenario.
        # Load this dataframe into a list, then combine all the dataframes in the list into a single dataframe containing data on all the scenarios.
        listOfDataframes = list()
        for (indexScenario in 1:length(inputs$scenarios))
        {
            scenario = inputs$scenarios[indexScenario]
            projectFile = inputs$projectFiles[indexScenario]

            # Check that the project file exists before trying to load it.
            if (!file.exists(projectFile))
            {
                stop(paste0("Error: project file not found: ", projectFile, " (scenario '", scenario, "', variable '", variable, "')"))
            }

            # Load the project file only if it has not been loaded before; otherwise use the cached version.
            if (is.null(projectCache[[projectFile]]))
            {
                message(paste0("Loading project file: ", projectFile, "..."))
                projectCache[[projectFile]] = tryCatch(
                    loadProject(projectFile),
                    error = function(e) stop(paste0("Error loading project file '", projectFile, "': ", e$message))
                )
            }
            project = projectCache[[projectFile]]

            message(paste0("[", indexInputs, "/", length(listOfInputs), "] Extracting '", variable,
                           "' for scenario '", scenario, "'..."))
            data = tryCatch(
                getQuery(project, query=variable),
                error = function(e) stop(paste0("Error: variable '", variable, "' not found in '", projectFile, "': ", e$message))
            )
            df = as.data.frame(data)
            df$scenario = scenario
            listOfDataframes[[indexScenario]] = df
        }
        df = bind_rows(listOfDataframes)

        # Write the dataframe to a output file in .csv format.
        if (endsWith(outputFile, ".csv"))
        {
            write.csv(df, outputFile)
        }
        else
        {
            write.csv(df, paste0(outputFile, ".csv"))
        }
        message(paste0("  Written to: ", outputFile))
    }
}
