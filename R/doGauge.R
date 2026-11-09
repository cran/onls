doGauge <- function(output = c("console", "html", "pdf"), outdir = tempdir(), fontsize = 14,
                    silent = TRUE, file = system.file("demo", "gauge.R", package = "onls"))
{
  output <- match.arg(output)
  if (!nzchar(file) || !file.exists(file)) stop("Demo file not found: ", file)
  op <- options(width = 200); on.exit(options(op))

  cat("Running the onls gauge demo ...\n")
  env <- new.env()
  run <- function() utils::capture.output(source(file, echo = TRUE, print.eval = TRUE, max.deparse.length = Inf,
                                                 local = env, spaced = FALSE))
  grDevices::pdf(NULL)                                                ## plots of the demo are discarded
  TXT <- tryCatch(if (silent) suppressWarnings(suppressMessages(run())) else run(),
                  finally = grDevices::dev.off())
  RES <- env$.GAUGE
  SUM <- table(factor(RES$result, levels = c("PASS", "NOTE", "DIFF")))
  isCode <- grepl("^[>+] ", TXT)
  head1 <- paste0("onls ", utils::packageVersion("onls"), ": gauge demo")
  head2 <- paste0(format(Sys.time()), " - ", R.version.string)
  out <- NULL

  if (output == "console") {
    cat(head1, "\n", head2, "\n\n", sep = "")
    cat(TXT, sep = "\n")
  }

  if (output == "html") {
    esc <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)
    H <- esc(TXT)
    H[isCode] <- paste0("<span class='code'>", H[isCode], "</span>")
    H[!isCode] <- gsub("  PASS$", "  <b style='color:#00AA00'>PASS</b>", H[!isCode])
    H[!isCode] <- gsub("  NOTE$", "  <b style='color:#AA0000'>NOTE</b>", H[!isCode])
    H[!isCode] <- gsub("  DIFF$", "  <b style='color:#0000AA'>DIFF</b>", H[!isCode])
    H[!isCode] <- gsub(" (published|closed form|ODRPACK|NIST) ", " <b>\\1</b> ", H[!isCode])
    out <- file.path(outdir, "onls-gauge.html")
    writeLines(c("<!DOCTYPE html><html><head><meta charset='utf-8'>", paste0("<title>", head1, "</title><style>"),
                 "body { font-family: sans-serif; margin: 1.5em; }",
                 paste0("pre { font-size: ", fontsize, "px; line-height: 1.3; white-space: pre; overflow-x: auto; }"),
                 ".code { color: #666666; }", "</style></head><body>",
                 paste0("<h2>", head1, "</h2><p>", head2, "</p>"),
                 paste0("<p><b>PASS ", SUM[["PASS"]], " - NOTE ", SUM[["NOTE"]], " - DIFF ", SUM[["DIFF"]], "</b></p><pre>"),
                 H, "</pre></body></html>"), out)
  }

  if (output == "pdf") {
    out <- file.path(outdir, "onls-gauge.pdf")
    grDevices::pdf(out, width = 11.69, height = 8.27, family = "Courier")      ## A4 landscape, monospaced
    op2 <- graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
    cex <- min(fontsize/14 * 0.75, 0.94/graphics::strwidth(TXT[which.max(nchar(TXT))], cex = 1))
    lh <- 2 * graphics::strheight("M", cex = cex); nPage <- floor(0.88/lh)
    ## one line, drawn segment by segment: bold reference kinds and the coloured result
    drawLine <- function(s, y, code) {
      if (code) { graphics::text(0.03, y, s, adj = c(0, 0.5), cex = cex, col = "#666666"); return(invisible()) }
      seg <- regmatches(s, gregexpr("(PASS|NOTE|DIFF)$|.+?(?=(PASS|NOTE|DIFF)$)|.+$", s, perl = TRUE))[[1]]
      x <- 0.03
      for (g in seg) {
        col <- switch(g, PASS = "#00AA00", NOTE = "#AA0000", DIFF = "#0000AA", "black")
        font <- if (g %in% c("PASS", "NOTE", "DIFF")) 2 else 1
        graphics::text(x, y, g, adj = c(0, 0.5), cex = cex, col = col, font = font)
        x <- x + graphics::strwidth(g, cex = cex, font = font)
      }
    }
    pages <- split(seq_along(TXT), ceiling(seq_along(TXT)/nPage))
    for (p in seq_along(pages)) {
      if (p > 1) graphics::plot.new()
      graphics::text(0.03, 0.97, head1, adj = c(0, 0.5), font = 2, cex = 0.9)
      graphics::text(0.97, 0.97, paste0(head2, "  -  page ", p, "/", length(pages)), adj = c(1, 0.5), cex = 0.6, col = "grey40")
      y <- 0.93
      for (i in pages[[p]]) { drawLine(TXT[i], y, isCode[i]); y <- y - lh }
    }
    graphics::par(op2); grDevices::dev.off()
  }

  cat("\nonls gauge: PASS ", SUM[["PASS"]], ", NOTE ", SUM[["NOTE"]], ", DIFF ", SUM[["DIFF"]],
      if (!is.null(out)) paste0("\nWritten to ", out), "\n", sep = "")
  invisible(list(file = out, results = RES, summary = SUM))
}
