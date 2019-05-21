library(RMySQL)
library(GenomicRanges)
library(dplyr)
library(stringr)
library(tibble)
library(rtracklayer)

options(useFancyQuotes = FALSE)

# Define database connection groups.
specimenManagement.DBgroup <- 'specimen_management'
currentIntSites.DBgroup    <- 'intsites_miseq'
prev454IntSites.DBgroup    <- '454intsites'


# Retrieve all sample data.
sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn1    <- dbConnect(MySQL(), group = specimenManagement.DBgroup)
sampleData <- dbGetQuery(dbConn1, "select * from gtsp where Trial = 'SCID1_Paris_Cavazzana'")


# Retreive all samples with current intSite data.
dbConn2        <- dbConnect(MySQL(), group = currentIntSites.DBgroup)
samplesClause  <- paste0('sampleName like ', paste0(sapply(sampleData$SpecimenAccNum, function(x){ sQuote(paste0(x, '-%'))  }), collapse=' or sampleName like '))
intSiteSamples <- unique(gsub('\\-\\d+$', '', unname(unlist(dbGetQuery(dbConn2, sprintf("select sampleName from samples where %s", samplesClause))))))
sapply(dbListConnections(MySQL()), dbDisconnect)


dbConn  <- dbConnect(MySQL(), group = prev454IntSites.DBgroup)
sets <- unname(unlist(dbGetQuery(dbConn, "select distinct setname from genetherapy_samples where setname like '%XSCID%'")))
sets <- unique(sets[grepl('\\-XSCIDp\\d\\d?\\-', sets, ignore.case = TRUE)])
sapply(dbListConnections(MySQL()), dbDisconnect)
  
remove_NA <- function(x) x[! is.na(x)]
sets.ids <- remove_NA(str_extract(sets, 'GTSP\\d+'))


# Stop if any of the final 454 sets are in the illuminaSets vector.
illuminaSets <- readLines('IlluminaSetsIn454DB')
if(any(toupper(sets) %in% toupper(illuminaSets))) stop()


# Retreive all data from the sets table.
dbConn   <- dbConnect(MySQL(), group = prev454IntSites.DBgroup)
sets.tbl <- dplyr::select(dbGetQuery(dbConn, 'select * from sets'), name, freeze, cellType, setDesc)


# For each set, retrieve the associated alignments and the number of times each alignment was observed (read count).
r <- dplyr::bind_rows(lapply(sets, function(setName){
      sql <- paste0("select psl, cloneCount from sets_psl where cloneCount >= 1 and name = '", setName, "'")
      psls <- dbGetQuery(dbConn, sql)
      
      if(nrow(psls) > 0){
        sql <- paste0("select id, strand, restriction, freeze, blockCount, qSize, tName, tStart, tEnd from psl where id in (", paste0(sQuote(psls$psl), collapse = ','), ")")
        d   <- suppressWarnings(dbGetQuery(dbConn, sql))
        
        if(nrow(d) > 0){ 
          d$reads    <- psls[match(d$id, psls$psl),]$cloneCount
          d$setName  <- setName
          d$posid    <- paste0(d$tName, d$strand, ifelse(d$strand == '+', d$tStart, d$tEnd))
          return(d)
        } else {
          return(tibble())
        }
      } else {
        return(tibble())
      }
}))
    

# Join meta data from the sets table and attempt to create time points from set names.
r <- dplyr::left_join(r, sets.tbl, by = c('setName' = 'name'))
r$patient <- str_match(toupper(r$setName), 'XSCID(P[P]?\\d+)\\-')[,2]
r$tp1 <- str_match(toupper(r$setName), '\\-([\\d\\_]+M)')[,2]
r$tp2 <- str_match(toupper(r$setName), '\\-([^\\-]+)\\-GTSP')[,2]

# Stop if time points were created with more than one method.
if(length(base::intersect(which(is.na(r$tp1)), which(is.na(r$tp2)))) != 0) stop('One or more time points could not be uniquely determined.')


# Decide upon a time point source, format it and rename data columns.
r$timePoint <- ifelse(is.na(r$tp1), r$tp2, r$tp1)
r$timePoint <- paste0('m', sub('M', '', sub('_', '.', r$timePoint)))

r <- dplyr::select(r, tName, tStart, tEnd, strand, restriction, reads, patient, posid, freeze.x, cellType, timePoint, setName, setDesc)
names(r) <- c('seqnames', 'start', 'end', 'strand', 'restriction', 'reads', 'patient', 'posid', 'freeze', 'cellType', 'timePoint', 'setName', 'setDesc')


# Extract GTSP ids from set names, default to NA if not possible.
r$GTSP <- str_extract(r$setName, 'GTSP\\d+')


# Average reads for fragments with duplicate restriction digests 
# then sum reads for each non-standardized fragment.
r2 <- dplyr::mutate(r, group = group_indices(r, patient, cellType, timePoint),
                GTSP  =  ifelse(is.na(GTSP), paste0('LS-', group), GTSP)) %>%
      dplyr::group_by(GTSP, seqnames, start, end, strand, restriction) %>%
      dplyr::mutate(reads = mean(reads)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()  %>%
      dplyr::group_by(GTSP, seqnames, start, end, strand) %>%
      dplyr::mutate(reads = sum(reads),
            restrictions = paste(restriction, collapse = ', ')) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(-restriction, -group, -setDesc, -freeze, -posid, -setName, -restrictions, seqnames, start, end, strand, reads, patient, cellType, timePoint)


# Lift data over from hg18 to hg38.
chain <- import.chain('hg18ToHg38.over.chain')
r3    <- as_tibble(data.frame(unlist(liftOver(makeGRangesFromDataFrame(r2, keep.extra.columns = TRUE), chain))))


# Add additional metadata so that data can be incorperated into patient report maker software.
r3$patient    <- paste0('p', r3$patient)
r3$trial      <- 'SCID1_Paris_Cavazzana'
r3$dataSource <- '454'
r3$timePoint  <- toupper(r3$timePoint)
r3$timePoint  <- gsub('_', '.', r3$timePoint)
r3$sampleName <- paste0(r3$GTSP, '-0')


# Convert time point to days and months elapsed.
r3$timePointType<- stringr::str_extract(r3$timePoint, '[DMY]')
r3$timePointType[which(is.na(r3$timePointType))] <- 'X'

r4 <- do.call(rbind, lapply(split(r3, r3$timePointType), function(x){
  n <- as.numeric(stringr::str_match(x$timePoint, '[\\d\\.]+')) * ifelse(grepl('\\-', x$timePoint), -1, 1)
  
  if(x$timePointType[1] == 'D'){
    x$timePointMonths <- n / 30.4167
    x$timePointDays   <- n
  } else if(x$timePointType[1] == 'M'){
    x$timePointMonths <- n
    x$timePointDays   <- n * 30.4167
  } else if(x$timePointType[1] == 'Y'){
    x$timePointMonths <- n * 12
    x$timePointDays   <- n * 365
  } else {
    message('Warning - could not determine date unit for: ', x$timePointType[1])
    x$timePointMonths <- n
    x$timePointDays   <- n 
  }
  x
}))
r4$timePointType <- NULL

saveRDS(r4, file = 'legacyData.rds')