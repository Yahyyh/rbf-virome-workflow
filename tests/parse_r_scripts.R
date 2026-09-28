#!/usr/bin/env Rscript

# Base-R syntax gate for continuous integration or a clean Conda environment.

args <- commandArgs(trailingOnly = TRUE)
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
default_root <- normalizePath(
  file.path(dirname(sub("^--file=", "", script_argument[1L])), ".."),
  winslash = "/",
  mustWork = TRUE
)
repository_root <- if (length(args) == 1L) {
  normalizePath(args[1L], winslash = "/", mustWork = TRUE)
} else default_root

r_directory <- file.path(repository_root, "scripts", "r")
r_files <- list.files(r_directory, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
if (length(r_files) == 0L) {
  stop(sprintf("No R scripts found under: %s", r_directory), call. = FALSE)
}

failures <- character(0)
for (path in sort(r_files)) {
  result <- tryCatch(
    {
      parse(file = path, keep.source = TRUE)
      NULL
    },
    error = function(error) conditionMessage(error)
  )
  if (is.null(result)) {
    message(sprintf("PARSE OK: %s", path))
  } else {
    failures <- c(failures, sprintf("%s: %s", path, result))
  }
}

if (length(failures) > 0L) {
  stop(paste(c("R parse failures:", failures), collapse = "\n"), call. = FALSE)
}
