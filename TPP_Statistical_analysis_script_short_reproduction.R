
### TPP data analysis reproduction

### Load raw data
raw_data <- read.csv("https://osf.io/24qkr/download")

### Sessions conducted with the test accounts or without lab_IDs are excluded (these IDs were preregistered to be excluded):
lab_IDs_to_exclude <- c("", "18155ef201564afbb81f6a8b74aa9a033eac51ec6595510eca9606938ffaced3", "ece83ceb8611d1926746e5bb3597ed1e8cb5d336521331b31961d5c0348883cf", "bd2dd15be34863e9efb77fbddfe744382a9c62c6a497e8bcf3097a47905b905b", "fff9cb9dcc3ac735fc25a59f424e98278a731c23ccd57276d292996c2ba7784f")
data_nontest <- raw_data[!(raw_data[,"laboratory_ID_code"] %in% lab_IDs_to_exclude), ]

### We exclude all rows which are not erotic trials. (Including rows that don't contain trial data, like demographics etc.)
data_nontest_trials = data_nontest[!is.na(data_nontest[, "trial_number"]),]
data_nontest_trials_erotic = data_nontest_trials[data_nontest_trials[, "reward_type"] == "erotic", ]

### we restrict the dataset to the only extend until the preregistered stopping point. So any data that was recorded after trial 37836 are not taken into account in the analysis. 
current_stopping_point = 37836
data_BF = data_nontest_trials_erotic[1:current_stopping_point,]

### Results
data_BF[,"sides_match"] = as.factor(tolower(as.logical(data_BF[,"sides_match"])))
successes = sum(as.logical(data_BF[,"sides_match"]))

# percentage of correct guesses
successes/nrow(data_BF)
