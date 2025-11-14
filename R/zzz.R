.onAttach <- function(lib, pkg)  {
  if (interactive()) {
    packageStartupMessage(
      "This is coenocliner ",
      utils::packageDescription(
        "coenocliner",
        field = "Version"
      ),
      appendLF = TRUE
    )
  }
}
