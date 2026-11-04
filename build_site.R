# Build the static Shinylive website from the app folder.
#
# Folder layout (names can be changed below):
#   <project folder>/
#     build_site.R
#     shinylive_app/app.R      <- the app (no data files in here!)
#
# Run with source("build_site.R") or the Source button in RStudio. The script
# finds its own folder, so the working directory does not matter.

# ---- Settings ---------------------------------------------------------------
# Leave PROJECT_DIR empty to use the folder this script is in, or give a full
# path, e.g. "C:/Users/me/Documents/my_site" (use / not \).
PROJECT_DIR <- ""
APP_DIR     <- "shinylive_app"   # app folder, inside PROJECT_DIR
SITE_DIR    <- "docs"            # output website, inside PROJECT_DIR
SERVE_AFTER <- TRUE              # open a local test server when done

# ---- Locate the project folder ----------------------------------------------
script_dir <- function() {
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)        # source()
  if (!is.null(of)) return(dirname(normalizePath(of)))
  if (requireNamespace("rstudioapi", quietly = TRUE) &&               # RStudio
      rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()                                                             # fallback
}
root <- if (nzchar(PROJECT_DIR)) normalizePath(PROJECT_DIR, mustWork = TRUE) else script_dir()
app  <- file.path(root, APP_DIR)
site <- file.path(root, SITE_DIR)

message("Project folder: ", root)
if (!dir.exists(app))
  stop("App folder not found: ", app,
       "\nCreate it (with app.R inside) or set PROJECT_DIR / APP_DIR at the top of this script.")
if (!file.exists(file.path(app, "app.R")))
  stop("No app.R in ", app)

# ---- Packages ---------------------------------------------------------------
for (pkg in c("shinylive", "httpuv"))
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)

# ---- Safety check: never publish data files with the site --------------------
other <- setdiff(list.files(app, recursive = TRUE), "app.R")
data_files <- grep("\\.(bed|tsv|txt|csv|gz|rds|rdata)$", other, ignore.case = TRUE, value = TRUE)
if (length(data_files))
  stop("Remove data files from ", app, " before exporting: ", paste(data_files, collapse = ", "))
if (length(other))
  message("Note: these files will also be published: ", paste(other, collapse = ", "))

# ---- Export -----------------------------------------------------------------
message("Exporting (the first run downloads the Shinylive runtime and packages; ",
        "this can take 5-20 minutes)...")
shinylive::export(app, site)
message("Website written to: ", site)

# ---- Local test ---------------------------------------------------------------
# Opens the site in your browser. The console stays busy while the server runs:
# press Esc (RStudio) or Ctrl+C to stop it. Then push the docs/ folder to GitHub
# and enable Pages (Settings > Pages > Deploy from branch, folder /docs).
if (SERVE_AFTER) httpuv::runStaticServer(site)