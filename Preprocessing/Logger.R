list.of.packages <- c("log4r", "docstring")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

library(log4r)
library(docstring)

get_logger <- function(file_name = "", threshold = "DEBUG", subDir = "Logging") {
    #' A simple logger initiator
    #'
    #' @description A simple logger initiator, using `log4r` package.
    #' the file and directory will be automatically created if not exists.
    #' @param file_name character. the actual logger file,
    #' `*.log` file is recommended.'
    #' @param threshold character. the logger level, default: "DEBUG"
    #' level options (from low to high): ["DEBUG, "INFO", "ERROR"]
    #' @param subDir character. log file directory, default: "Logging"
    #' @return It returns a logger object for logging messages.
    #' @note you may log the message by [log4r::info, log4r::debug, log4r::error]
    #' depends on the logger's threshold, the message should be recorded at `./Logging/file_name`.
    #' @examples
    #' r$> logger <- get_logger("my_logger.txt")
    #' r$> log4r::info(logger, "Info_message.")
    #' # INFO  [2021-07-02 17:22:58] Info_message.
    #' r$> readLines("Logging/my_logger.txt")
    #' [1] "INFO  [2021-07-02 17:22:58] Info_message."
    #' # changing threhold
    #' logger <- get_logger("my_logger.txt", threshold = "INFO")
    #' r$> log4r::error(logger, "Error_message.")
    #' ERROR [2021-07-02 17:50:09] Error_message.
    #' r$> log4r::debug(logger, "Debug_message.")
    #' r$> readLines("Logging/my_logger.txt")
    #' [1] "INFO  [2021-07-02 17:22:58] Info_message."
    #' [2] "ERROR [2021-07-02 17:50:09] Error_message."

    logging_dir <- file.path(getwd(), subDir)
    if (!dir.exists(logging_dir)) {
        dir.create(logging_dir)
    }

    logfile_path <- file.path(logging_dir, paste(Sys.Date(), file_name))
    if (!file.exists(logfile_path)) {
        file.create(logfile_path)
    }

    my_console_appender <- log4r::console_appender(
        layout = default_log_layout()
    )

    my_file_appender <- log4r::file_appender(logfile_path,
        append = TRUE,
        layout = default_log_layout()
    )

    return(log4r::logger(
        threshold = threshold, appenders =
            list(my_console_appender, my_file_appender)
    ))
}

sanity_check <- function(condition, case,
                         logger, caller = match.call()[[1]], warning = FALSE) {
    log4r::debug(
        logger = logger,
        message = paste("[", caller, "] checking: ", case)
    )

    if (condition) {
        log4r::debug(logger = logger, message = "...passed")
        return(TRUE)
    }

    message <- paste("::: [", caller, "] failed. :::\n", "\t- REASON: ", case)
    if (warning) {
        log4r::warn(logger = logger, message = message)
        return(FALSE)
    } else {
        log4r::error(logger = logger, message = message)
        stop(message)
    }
}

logger_init_msg <- function(curr_func, logger, msg = "") {
    log4r::debug(
        logger = logger,
        message = paste("*** [", curr_func, "] initialised.", msg, "***")
    )
}

logger_complete_msg <- function(curr_func, logger, msg = "") {
    log4r::debug(
        logger = logger,
        message = paste("*** [", curr_func, "] completed.", msg, "***")
    )
}