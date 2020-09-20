data(camtrapsMultiSeason)
data(recordTableSampleMultiSeason)

recordTableSampleMultiSeason <- recordTableSampleMultiSeason[, c("Station", "Species", "DateTimeOriginal")]


dateFormat <- "%d/%m/%Y"

camop_season <- cameraOperation(CTtable         = camtrapsMultiSeason,
                                stationCol   = "Station",
                                setupCol     = "Setup_date",
                                sessionCol   = "session",
                                retrievalCol = "Retrieval_date",
                                hasProblems  = TRUE,
                                dateFormat   = dateFormat
)

View(camop_season)

DetHist_multi <- detectionHistory(recordTable      = recordTableSampleMultiSeason,
camOp                  = camop_season,
stationCol             = "Station",
speciesCol             = "Species",
species                = "VTA",
occasionLength         = 10,
day1                   = "station",
recordDateTimeCol      = "DateTimeOriginal",
includeEffort          = TRUE,
scaleEffort            = FALSE,
timeZone               = "UTC",
unmarkedMultFrameInput = FALSE
)

