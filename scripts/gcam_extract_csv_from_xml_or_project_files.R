# do not run this in e3sm unified
# need to load the R module first, then use Rscript to run this script from the command line

message(paste0("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

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
    stop("Usage: Rscript gcam_extract_csv_from_xml_or_project_files.R `path/to/json/input/file(s)'")
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

        # Extract required inputs.
        variable   = inputs$variable
        outputFile = inputs$outputFile
        scenarios  = inputs$scenarios
        projectFiles = inputs$projectFiles

        # Extract optional inputs.
        # xmlFiles: raw GCAM XML output files from which project files are created (one per scenario).
        #           Only used when createProjectFiles is TRUE.
        xmlFiles = if (!is.null(inputs$xmlFiles)) inputs$xmlFiles else NULL
        # createProjectFiles: if TRUE, create project files from xmlFiles before extracting data.
        #                     If FALSE (default), projectFiles must already exist and queryFile is not needed.
        createProjectFiles = if (!is.null(inputs$createProjectFiles)) inputs$createProjectFiles else FALSE
        # queryFile: path to the XML query file containing the queries to run against the XML output.
        #            Only required when createProjectFiles is TRUE. Defaults to ../workflow/e3sm_gcam_queries.xml.
        queryFile = if (!is.null(inputs$queryFile)) inputs$queryFile else "../workflow/e3sm_gcam_queries.xml"
        # maxMemory: Java heap memory limit for createProjectFiles step. Defaults to '32g'.
        maxMemory = if (!is.null(inputs$maxMemory)) inputs$maxMemory else "32g"

        # Check that scenarios and projectFiles arrays are the same length.
        if (length(scenarios) != length(projectFiles))
        {
            stop(paste0("Error: 'scenarios' (length ", length(scenarios), ") and 'projectFiles' (length ",
                        length(projectFiles), ") must be the same length for variable '", variable, "'."))
        }

        # If creating project files, xmlFiles must also be specified and the same length as scenarios.
        if (createProjectFiles)
        {
            if (is.null(xmlFiles))
            {
                stop(paste0("Error: 'createProjectFiles' is true but 'xmlFiles' is not specified for variable '", variable, "'."))
            }
            if (length(xmlFiles) != length(scenarios))
            {
                stop(paste0("Error: 'xmlFiles' (length ", length(xmlFiles), ") and 'scenarios' (length ",
                            length(scenarios), ") must be the same length for variable '", variable, "'."))
            }
        }

        # If createProjectFiles is TRUE but all project files already exist, warn and skip creation.
        if (createProjectFiles && all(file.exists(projectFiles)))
        {
            message("Warning: project files already exist. Not creating new project files from model output. ",
                    "To create new project files either delete/rename existing files or rename project files ",
                    "in input json file.")
            createProjectFiles = FALSE
        }

        # Check that the output directory exists before attempting to write; create it if not.
        outputDir = dirname(outputFile)
        if (!dir.exists(outputDir))
        {
            message(paste0("Output directory does not exist: '", outputDir, "'. Creating it now..."))
            dir.create(outputDir, recursive=TRUE)
        }

        # If creating project files, read the query file contents now.
        # queryFile is only needed for the addScenario step; it is not required when reading
        # from existing project files since getQuery looks up pre-computed results by name.
        queryFileContents = NULL
        if (createProjectFiles)
        {
            if (!file.exists(queryFile))
            {
                stop(paste0("Error: query file not found: '", queryFile, "'. ",
                            "'queryFile' is required when 'createProjectFiles' is true."))
            }
            queryFileContents = readChar(queryFile, file.info(queryFile)$size)
        }

        # Each scenario for the current block will be read into a dataframe from its project file.
        # Load all dataframes into a list, then combine into a single dataframe containing all scenarios.
        listOfDataframes = list()
        for (indexScenario in 1:length(scenarios))
        {
            scenario    = scenarios[indexScenario]
            projectFile = projectFiles[indexScenario]

            # If requested, create the project file from the raw GCAM XML output before loading it.
            # This uses Java via rgcam's localDBConn and addScenario, same pattern as write_carbon_price_paths.r.
            if (createProjectFiles)
            {
                xmlFile = xmlFiles[indexScenario]
                if (!file.exists(xmlFile))
                {
                    stop(paste0("Error: XML file not found: '", xmlFile, "' (scenario '", scenario, "')"))
                }
                projectFileDir = dirname(projectFile)
                if (!dir.exists(projectFileDir))
                {
                    message(paste0("Directory for project file does not exist: '", projectFileDir, "'. Creating it now..."))
                    dir.create(projectFileDir, recursive=TRUE)
                }
                message(paste0("Creating project file '", projectFile, "' from '", xmlFile, "'..."))
                conn = tryCatch(
                    localDBConn(dbPath='', dbFile=xmlFile, maxMemory=maxMemory),
                    error = function(e) stop(paste0("Error opening XML file '", xmlFile, "': ", e$message))
                )
                tryCatch(
                    addScenario(conn, projectFile, scenario=NULL, queryFile=queryFileContents, clobber=TRUE),
                    error = function(e) stop(paste0("Error creating project file '", projectFile, "': ", e$message))
                )
                # Invalidate any cached version of this project file since it has just been (re)created.
                projectCache[[projectFile]] = NULL
                message(paste0("  Project file created: ", projectFile))
            }

            # Check that the project file exists before trying to load it.
            if (!file.exists(projectFile))
            {
                stop(paste0("Error: project file not found: '", projectFile,
                            "' (scenario '", scenario, "', variable '", variable, "'). ",
                            "Set 'createProjectFiles': true and provide 'xmlFiles' to create it automatically."))
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
                           "' for scenario '", scenario, "' from '", projectFile, "'..."))
            data = tryCatch(
                getQuery(project, query=variable),
                error = function(e) {
                    message(paste0("Warning: variable '", variable, "' not found in '", projectFile,
                                   "' - skipping. This may mean the query returned no data for this ",
                                   "simulation configuration. Error: ", e$message))
                    return(NULL)
                }
            )
            if (is.null(data)) next
            df = as.data.frame(data)
            df$scenario = scenario
            listOfDataframes[[indexScenario]] = df
        }
        df = bind_rows(listOfDataframes)

        # Write the dataframe to an output file in .csv format.
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
message(paste0("End time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
