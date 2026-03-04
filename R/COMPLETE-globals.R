# This creates the environment when the package is loaded/sourced.
# It is internal to the package namespace.
COMPLETE_env <- new.env(parent = emptyenv())

.tc_env <- new.env(parent = emptyenv())
.tc_env$stack <- list()

`%notin%` <- Negate(`%in%`)

.merge_style <- function(parent, child) {
  
  modifyList(parent, child)
}

#' gen_hex_color implementation
#'
#' Generates N random hex colors using the rbg() function 
#' 
#' @author Google Gemini
#' @param n (integer) Number of hex colors to generate
#' @param bright? (bool) Should the hex colors be bright? (Default: TRUE)
#' @return N hex colors
#' @export
gen_hex_color <- function(n = 1, bright = TRUE) {
  min_val <- if(bright) 0.4 else 0.0
  rgb(runif(n, min_val, 1), runif(n, min_val, 1), runif(n, min_val, 1))
}

#' get_contrast_color implementation
#'
#' Get contrast color of a hex color code
#' 
#' @author Google Gemini
#' @param n (integer) Number of hex colors to generate
#' @param bright? (bool) Should the hex colors be bright? (Default: TRUE)
#' @return N hex colors
#' @export
get_contrast_color <- function(hex_color) {
  # Convert Hex to RGB and normalize to 0-1 range
  rgb_val <- col2rgb(hex_color) / 255
  
  # Calculate perceived brightness
  luminance <- 0.299 * rgb_val[1] + 0.587 * rgb_val[2] + 0.114 * rgb_val[3]
  
  # Return Black for bright backgrounds, White for dark backgrounds
  if (luminance > 0.5) {
    return("#000000") 
  } else {
    return("#FFFFFF")
  }
}

#' dim_color implementation
#'
#' Dim color of a hex color code
#' 
#' @author Google Gemini
#' @param hex_color (string) hex color code
#' @return Dim hhex color
#' @seealso [COMPLETE::gen_hex_color()]
#' @export
dim_color <- function(hex_color) {
  # Convert Hex to RGB
  rgb_val <- col2rgb(hex_color) / 255
  # Reduce intensity by 30%
  dimmed_rgb <- rgb_val * 0.7
  # Convert back to Hex
  rgb(dimmed_rgb[1], dimmed_rgb[2], dimmed_rgb[3])
}

#' termcolor implementation
#'
#' Prints the output of enclosed expr in color
#' 
#' @author Google Gemini
#' @param expr R Expression
#' @param color FG color
#' @param type Foreground/Background (<fg>/bg). Leave empty for fg and bg
#' @param bg_color Background color if type == "bg"
#' @param bold? Boldface?
#' @param italic? Italicface?
#' @param underline? Underline?
#' @param strikethrough? Strike-through?
#' @param url? URL?
#' @param auto_contrast? Auto adjust output colors based on contrast? (Default: TRUE)
#' @return NULL
#' @export
termcolor <- function(expr, color = NULL, bg_color = NULL, type = NULL, 
                      bold = NULL, italic = NULL, underline = NULL, 
                      strikethrough = NULL, url = NULL, auto_contrast = TRUE) {
  
  options(cli.num_colors = 256)
  expr_sub <- substitute(expr)
  parent_env <- parent.frame()
  
  is_nested <- length(.tc_env$stack) > 0
  
  parent_style <- if(is_nested) tail(.tc_env$stack, 1)[[1]] else {
    # DEFAULT RECOVERY: If not nested, start with NO colors/styles
    list(color = NULL, bg_color = NULL, type = "fg", bold = FALSE, 
         italic = FALSE, underline = FALSE, strikethrough = FALSE, url = NULL)
  }
  # # 1. Improved Inheritance Logic
  # parent_style <- if(length(.tc_env$stack) > 0) tail(.tc_env$stack, 1)[[1]] else {
  #   list(color = NULL, bg_color = NULL, type = "fg", bold = FALSE, italic = FALSE,
  #        underline = FALSE, strikethrough = FALSE, url = NULL)
  # }
  
  # Inherit from parent ONLY if the argument is missing from the current call
  # current_style <- list(
  #   color         = if(!missing(color)) color else parent_style$color,
  #   bg_color      = if(!missing(bg_color)) bg_color else parent_style$bg_color,
  #   type          = if(!missing(type)) type else parent_style$type,
  #   bold          = if(!missing(bold)) bold else parent_style$bold,
  #   italic        = if(!missing(italic)) italic else parent_style$italic,
  #   underline     = if(!missing(underline)) underline else parent_style$underline,
  #   strikethrough = if(!missing(strikethrough)) strikethrough else parent_style$strikethrough,
  #   url           = if(!missing(url)) url else parent_style$url
  # )
  current_style <- list(
    color         = if(!missing(color)) color else (if(is_nested) parent_style$color else NULL),
    italic        = if(!missing(italic)) italic else (if(is_nested) parent_style$italic else FALSE),
    bg_color         = if(!missing(bg_color)) bg_color else (if(is_nested) parent_style$bg_color else NULL),
    type         = if(!missing(type)) type else (if(is_nested) parent_style$type else NULL),
    bold         = if(!missing(bold)) bold else (if(is_nested) parent_style$bold else NULL),
    italic         = if(!missing(italic)) italic else (if(is_nested) parent_style$italic else NULL),
    underline         = if(!missing(underline)) underline else (if(is_nested) parent_style$underline else NULL),
    strikethrough         = if(!missing(strikethrough)) strikethrough else (if(is_nested) parent_style$strikethrough else NULL),
    url         = if(!missing(url)) url else (if(is_nested) parent_style$url else NULL)
  )
  
  # 2. Logic for rendering and auto-contrast
  if (!is.null(current_style$bg_color)) {
    # If BG exists and user didn't specify a type, default to "both" to ensure it shows
    if (missing(type)) current_style$type <- "both"
    
    # Only auto-contrast if a specific color wasn't provided in THIS call
    if (auto_contrast && missing(color)) {
      current_style$color <- get_contrast_color(current_style$bg_color)
    }
  } else {
    # If BG is NULL, ensure we aren't trying to render it
    if (missing(type) || current_style$type == "both") current_style$type <- "fg"
  }
  
  # ... [Push state to stack and capture output as before] ...
  .tc_env$stack <- c(.tc_env$stack, list(current_style))
  on.exit({ .tc_env$stack <- head(.tc_env$stack, -1) }, add = TRUE)
  
  tmp <- tempfile(); tcon <- file(tmp, open = "wt")
  
  # ONLY sink output. This is stackable and safe for nesting.
  sink(tcon, type = "output")
  
  res <- tryCatch({
    withCallingHandlers(
      withVisible(eval(expr_sub, envir = parent_env)),
      
      # 1. Capture and Format Warnings
      warning = function(w) {
        warn_msg <- cli::style_bold(paste("Warning:", w$message))
        warn_msg <- cli::make_ansi_style("#FFA500")(warn_msg)
        # We cat() this into the current output sink
        cat("\001", warn_msg, "\002\n", sep = "")
        invokeRestart("muffleWarning")
      },
      
      # 2. Capture standard messages (The fix for the "Grey Box")
      message = function(m) {
        # biomartr messages like "Dataset not found..." arrive here.
        # We route them to cat() so they hit our styled buffer instead of stderr.
        cat(m$message) 
        # Muffling prevents the message from reaching the actual console/stderr
        invokeRestart("muffleMessage")
      }
    )
  }, error = function(e) {
    err_msg <- cli::style_bold(paste("Error:", e$message))
    err_msg <- cli::make_ansi_style("#FF0000")(err_msg)
    cat("\001", err_msg, "\002\n", sep = "")
    return(NULL)
  }, finally = {
    # Close ONLY the output sink
    sink(type = "output")
    if(exists("tcon") && isOpen(tcon)) close(tcon)
  })
  
  full_output <- readLines(tmp, warn = FALSE)
  unlink(tmp)
  
  if (length(full_output) > 0 && (cli::is_ansi_tty() || cli::is_dynamic_tty()) ) {
    style_out <- paste(full_output, collapse = "\n")
    parts <- strsplit(style_out, "(?<=\\002)|(?=\\001)", perl = TRUE)[[1]]
    res_str <- ""
    
    for (p in parts) {
      if (p == "" || p == "\n") {
        res_str <- paste0(res_str, p); next
      }
      
      if (startsWith(p, "\001")) {
        res_str <- paste0(res_str, p)
      } else {
        clean_p <- sub("\n$", "", p)
        had_newline <- grepl("\n$", p)
        styled_p <- clean_p
        
        if (isTRUE(current_style$bold))          styled_p <- cli::style_bold(styled_p)
        if (isTRUE(current_style$italic))        styled_p <- cli::style_italic(styled_p)
        if (isTRUE(current_style$underline))     styled_p <- cli::style_underline(styled_p)
        if (isTRUE(current_style$strikethrough)) styled_p <- cli::style_strikethrough(styled_p)
        
        # Color Logic: ONLY apply BG if current_style$bg_color is not NULL
        if (current_style$type %in% c("fg", "both") && !is.null(current_style$color)) {
          styled_p <- cli::make_ansi_style(current_style$color)(styled_p)
        }
        if (current_style$type %in% c("bg", "both") && !is.null(current_style$bg_color)) {
          styled_p <- cli::make_ansi_style(current_style$bg_color, bg = TRUE)(styled_p)
        }
        
        if (!is.null(current_style$url) && nzchar(current_style$url)) {
          styled_p <- cli::style_hyperlink(styled_p, url = current_style$url)
        }
        
        if (had_newline) styled_p <- paste0(styled_p, "\n")
        res_str <- paste0(res_str, "\001", styled_p, "\002")
      }
    }
    
    if (length(.tc_env$stack) == 1) {
      res_str <- gsub("\001|\002", "", res_str)
      res_str <- paste0(res_str, "\033[0m") 
    }
    cat(res_str, "\n")
  }else{
    cat(cli::ansi_strip(full_output))
  }
  
  if (!is.null(res) && res$visible) return(res$value) else return(invisible(res$value))
}

# termcolor <- function(expr, color = NULL, bg_color = NULL, type = NULL, 
#                       bold = NULL, italic = NULL, auto_contrast = TRUE) {
#   
#   # CRITICAL FIX 1: Force CLI to generate colors even when sink() is active
#   options(cli.num_colors = 256)
#   
#   expr_sub <- substitute(expr)
#   parent_env <- parent.frame()
#   
#   # 1. Inheritance Logic
#   parent_style <- if(length(.tc_env$stack) > 0) tail(.tc_env$stack, 1)[[1]] else {
#     list(color = NULL, bg_color = NULL, type = "fg", bold = FALSE, italic = FALSE)
#   }
#   
#   current_style <- list(
#     color    = if(!is.null(color))    color    else parent_style$color,
#     bg_color = if(!is.null(bg_color)) bg_color else parent_style$bg_color,
#     type     = if(!is.null(type))     type     else parent_style$type,
#     bold     = if(!is.null(bold))     bold     else parent_style$bold,
#     italic   = if(!is.null(italic))   italic   else parent_style$italic
#   )
#   
#   if (auto_contrast && !is.null(current_style$bg_color)) {
#     current_style$color <- get_contrast_color(current_style$bg_color)
#     current_style$type <- "both"
#   } else if (!is.null(bg_color) && is.null(type)) {
#     current_style$type <- "both"
#   }
#   
#   # print(parent_style)
#   # print(current_style)
#   
#   # 2. Push state to stack
#   .tc_env$stack <- c(.tc_env$stack, list(current_style))
#   on.exit({ .tc_env$stack <- head(.tc_env$stack, -1) }, add = TRUE)
#   
#   # 3. Capture Output (Warnings/Errors are shielded immediately)
#   tmp <- tempfile(); tcon <- file(tmp, open = "wt")
#   sink(tcon, type = "output"); sink(tcon, type = "message")
#   
#   res <- tryCatch({
#     withCallingHandlers(
#       withVisible(eval(expr_sub, envir = parent_env)),
#       warning = function(w) {
#         warn_msg <- cli::style_bold(paste("Warning:", w$message))
#         warn_msg <- cli::make_ansi_style("#FFA500")(warn_msg)
#         message("\001", warn_msg, "\002") # Shield from parent
#         invokeRestart("muffleWarning")
#       }
#     )
#   }, error = function(e) {
#     err_msg <- cli::style_bold(paste("Error:", e$message))
#     err_msg <- cli::make_ansi_style("#FF0000")(err_msg)
#     message("\001", err_msg, "\002") # Shield from parent
#     return(NULL)
#   }, finally = {
#     sink(type = "message"); sink(type = "output")
#     close(tcon)
#   })
#   
#   # 4. Process Output with Shielding Logic
#   full_output <- readLines(tmp, warn = FALSE)
#   unlink(tmp)
#   
#   if (length(full_output) > 0) {
#     style_out <- paste(full_output, collapse = "\n")
#     
#     # CRITICAL FIX 2: Split text into protected and unprotected chunks
#     parts <- strsplit(style_out, "(?<=\\002)|(?=\\001)", perl = TRUE)[[1]]
#     res_str <- ""
#     
#     for (p in parts) {
#       if (p == "") next
#       
#       if (startsWith(p, "\001")) {
#         # This chunk was styled by a child; DO NOT TOUCH IT!
#         res_str <- paste0(res_str, p)
#       } else {
#         # Unprotected chunk; apply THIS level's styling
#         styled_p <- p
#         if (current_style$bold)   styled_p <- cli::style_bold(styled_p)
#         if (current_style$italic) styled_p <- cli::style_italic(styled_p)
#         
#         if (current_style$type %in% c("fg", "both") && !is.null(current_style$color)) {
#           styled_p <- cli::make_ansi_style(current_style$color)(styled_p)
#         }
#         if (current_style$type %in% c("bg", "both")) {
#           bg_val <- if (!is.null(current_style$bg_color)) current_style$bg_color else current_style$color
#           if (!is.null(bg_val)) styled_p <- cli::make_ansi_style(bg_val, bg = TRUE)(styled_p)
#         }
#         
#         # Wrap the freshly styled chunk in shields for the parent layer above
#         res_str <- paste0(res_str, "\001", styled_p, "\002")
#       }
#     }
#     
#     # If we are at the outermost block, strip shields and print to console
#     if (length(.tc_env$stack) == 1) {
#       res_str <- gsub("\001|\002", "", res_str)
#     }
#     
#     cat(res_str, "\n")
#   }
#   
#   # 5. Return Assignment compatibility
#   if (!is.null(res) && res$visible) return(res$value) else return(invisible(res$value))
# }

termcolor <- compiler::cmpfun(termcolor)
dim_color <- compiler::cmpfun(dim_color)
get_contrast_color <- compiler::cmpfun(get_contrast_color)
gen_hex_color <- compiler::cmpfun(gen_hex_color)

# termcolor <- function(expr,
#                       color = "#FFFFFF",
#                       bg_color = NULL,
#                       type = NULL,
#                       bold = FALSE,
#                       italic = FALSE,
#                       auto_contrast = FALSE) {
# 
#   expr_sub <- substitute(expr)
# 
#   # Determine type automatically
#   if (auto_contrast && !is.null(bg_color)) {
#     color <- get_contrast_color(bg_color)
#     type <- "both"
#   } else if (is.null(type)) {
#     type <- if (!is.null(bg_color)) "both" else "fg"
#   }
# 
#   tmp <- tempfile()
#   tcon <- file(tmp, open = "wt")
# 
#   sink(tcon, type = "output")
#   sink(tcon, type = "message")
# 
#   res <- NULL
#   error_occured <- FALSE
# 
#   tryCatch({
#     # res <- withVisible(eval(expr_sub, envir = parent.frame()))
# 
#     wrapper <- function() {
#       eval(expr_sub, envir = parent.frame())
#     }
# 
#     res <- withVisible(wrapper())
# 
#   },
#   error = function(e) {
# 
#     error_occured <<- TRUE
# 
#     # Close sinks before styling error
#     sink(type = "message")
#     sink(type = "output")
#     close(tcon)
# 
#     full_output <- readLines(tmp, warn = FALSE)
#     unlink(tmp)
# 
#     # Print captured output first
#     if (length(full_output) > 0) {
#       cat(paste(full_output, collapse = "\n"), "\n")
#     }
# 
#     # Styled bold red error
#     err_msg <- paste("Error:", e$message)
# 
#     err_msg <- cli::style_bold(err_msg)
#     err_msg <- cli::make_ansi_style("#FF0000")(err_msg)
# 
#     cat(err_msg, "\n")
# 
#     return(invisible(NULL))
#   })
# 
#   if (!error_occured) {
#     sink(type = "message")
#     sink(type = "output")
#     close(tcon)
# 
#     full_output <- readLines(tmp, warn = FALSE)
#     unlink(tmp)
# 
#     style_out <- paste(full_output, collapse = "\n")
# 
#     if (length(full_output) > 0) {
# 
#       if (bold) style_out <- cli::style_bold(style_out)
#       if (italic) style_out <- cli::style_italic(style_out)
# 
#       if (type %in% c("fg", "both")) {
#         style_out <- cli::make_ansi_style(color)(style_out)
#       }
# 
#       if (type %in% c("bg", "both")) {
#         bg_val <- if (!is.null(bg_color)) bg_color else color
#         style_out <- cli::make_ansi_style(bg_val, bg = TRUE)(style_out)
#       }
# 
#       cat(style_out, "\n")
#     }
# 
#     if (!is.null(res) && res$visible) {
#       # return(invisible(res$value))
#       return(res$value)
#     }
#   }
# }

# termcolor <- function(expr, color = "#FFFFFF", bg_color = NULL, 
#                       type = NULL, bold = FALSE, italic = FALSE,
#                       auto_contrast = FALSE) { 
#   # Automatically determine the best text color if requested
#   if (auto_contrast && !is.null(bg_color)) {
#     color <- get_contrast_color(bg_color)
#     type <- "both"
#   } else if (is.null(type)) {
#     type <- if (!is.null(bg_color)) "both" else "fg"
#   }
#   
#   tmp <- tempfile()
#   tcon <- file(tmp, open = "wt")
#   sink(tcon, type = "output")
#   sink(tcon, type = "message")
#   
#   res <- tryCatch({
#     withVisible(eval(expr, envir = parent.frame()))
#   }, error = function(e) {
#     cat("Error: ", e$message, "\n")
#     return(NULL)
#   }, finally = {
#     sink(type = "message")
#     sink(type = "output")
#     close(tcon)
#   })
#   
#   full_output <- readLines(tmp, warn = FALSE)
#   unlink(tmp)
#   
#   style_out <- paste(full_output, collapse = "\n")
#   
#   if (length(full_output) > 0) {
#     if (bold) style_out <- cli::style_bold(style_out)
#     if (italic) style_out <- cli::style_italic(style_out)
#     
#     if (type %in% c("fg", "both")) {
#       style_out <- cli::make_ansi_style(color)(style_out)
#     }
#     
#     if (type %in% c("bg", "both")) {
#       bg_val <- if (!is.null(bg_color)) bg_color else color
#       style_out <- cli::make_ansi_style(bg_val, bg = TRUE)(style_out)
#     }
#     
#     cat(style_out, "\n")
#   }
#   
#   if (!is.null(res) && res$visible) return(invisible(res$value))
# }

# termcolor <- function(expr, color = "#FFFFFF", bg_color = NULL, type = NULL) {
#   # 1. Smart Type Detection
#   if (is.null(type)) {
#     type <- if (!is.null(bg_color)) "both" else "fg"
#   }
#   
#   # 2. Capture ALL output (stdout + stderr) simultaneously
#   tmp <- tempfile()
#   tcon <- file(tmp, open = "wt")
#   sink(tcon, type = "output")
#   sink(tcon, type = "message")
#   
#   res <- tryCatch({
#     # Evaluate expression in the calling environment
#     withVisible(eval(expr, envir = parent.frame()))
#   }, error = function(e) {
#     message("Error in expression: ", e$message)
#     return(NULL)
#   }, finally = {
#     sink(type = "message")
#     sink(type = "output")
#     close(tcon)
#   })
#   
#   full_output <- readLines(tmp, warn = FALSE)
#   unlink(tmp)
#   
#   # 3. Build the Style Function using cli (consistent with Tidyverse)
#   style_fn <- function(x) x
#   
#   if (type %in% c("fg", "both")) {
#     style_fn <- cli::combine_ansi_styles(style_fn, cli::make_ansi_style(color))
#   }
#   
#   if (type %in% c("bg", "both")) {
#     target_bg <- if (!is.null(bg_color)) bg_color else color
#     style_fn <- cli::combine_ansi_styles(style_fn, cli::make_ansi_style(target_bg, bg = TRUE))
#   }
#   
#   # 4. Print the styled output block
#   if (length(full_output) > 0) {
#     cat(style_fn(paste(full_output, collapse = "\n")), "\n")
#   }
#   
#   # Return result invisibly
#   if (!is.null(res) && res$visible) return(invisible(res$value))
# }

# gen_hex_color <- function(n = 1) {
#   # runif(n) generates n random numbers between 0 and 1
#   # rgb() converts these decimal fractions into #RRGGBB format
#   rgb(runif(n), runif(n), runif(n))
# }

# termcolor <- function(expr, color = "#FFFFFF", type = "fg", bg_color = NULL) {
#   # 1. Capture stdout and stderr simultaneously to a temp file
#   tmp <- tempfile()
#   tcon <- file(tmp, open = "wt")
#   sink(tcon, type = "output")
#   sink(tcon, type = "message")
#   
#   res <- tryCatch({
#     withVisible(eval(expr, envir = parent.frame()))
#   }, finally = {
#     sink(type = "message")
#     sink(type = "output")
#     close(tcon)
#   })
#   
#   full_output <- readLines(tmp, warn = FALSE)
#   unlink(tmp)
#   
#   # 2. Define the styling logic using cli
#   style_fn <- switch(type,
#                      "fg"   = cli::make_ansi_style(color),
#                      "bg"   = cli::make_ansi_style(color, bg = TRUE),
#                      "both" = cli::combine_ansi_styles(
#                        cli::make_ansi_style(color), 
#                        cli::make_ansi_style(if(!is.null(bg_color)) bg_color else color, bg = TRUE)
#                      ),
#                      cli::make_ansi_style(color) # Default to fg if type is unknown
#   )
#   
#   # 3. Print the styled output
#   if (length(full_output) > 0) {
#     cat(style_fn(paste(full_output, collapse = "\n")), "\n")
#   }
#   
#   # Return original result invisibly
#   if (res$visible) return(invisible(res$value))
# }

# force_reset_future <- function() {
#   message("Executing Nuclear Future Reset...")
#   
#   # 1. Access the internal registry safely
#   # We use future::: because it's not exported to the user
#   reg <- future:::ClusterRegistry
#   
#   try({
#     # Get the current cluster
#     cl <- reg("get")
#     if (!is.null(cl)) {
#       message("Stopping active cluster nodes...")
#       parallel::stopCluster(cl)
#     }
#   }, silent = TRUE)
#   
#   # 2. Wipe the registry clean
#   reg("set", NULL)
#   
#   # 3. Revert to sequential plan
#   plan(sequential)
#   
#   # 4. Final garbage collection
#   gc()
#   
#   message("Reset complete. Backend is now sequential.")
# }

# force_reset_future <- function() {
#   message("Executing Nuclear Future Reset...")
#   
#   # 1. Force R to stop all cluster workers without trying to 'talk' to them
#   try(parallel::stopCluster(future:::ClusterRegistry("get")), silent = TRUE)
#   
#   # 2. Clear the internal registry
#   future:::ClusterRegistry("set", NULL)
#   
#   # 3. Set plan to sequential to clear the stack
#   plan(sequential)
#   
#   # 4. Garbage collection
#   gc()
#   
#   message("Reset complete.")
# }

# Loads the parameters from package level env
#
# @param param parameter
# @param ... dot_args from or for a function
# @return Value on successful fetch
# check_global(param, ...){
#   dot_args <- list(...)
#   value <- NULL
#   if(param %in% names(dot_args)){
#     value <- dot_args[[param]]
#   }else{
#     value <- COMPLETE_env[[param]]
#   }
#   if(isTRUE(is.null(value))){
#     stop(paste("[check_global()] Error:",param,"is empty. value is NULL.\n"))
#   }
# }