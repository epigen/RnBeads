% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Report-methods.R
\name{rnb.add.table}
\alias{rnb.add.table}
\title{rnb.add.table}
\usage{
rnb.add.table(
  report,
  tdata,
  row.names = TRUE,
  first.col.header = FALSE,
  indent = 0,
  tag.attrs = c(class = "tabdata"),
  thead = NULL,
  tcaption = NULL,
  na = "<span class=\\"disabled\\">n/a</span>"
)
}
\arguments{
\item{report}{Report to write the text to.}

\item{tdata}{Matrix or data frame to be presented in HTML form. Column names, if present, are used to
define table columns. If this table contains 0 (zero) rows or 0 columns, calling this
function has no effect.}

\item{row.names}{Flag indicating if row names should also be printed. If this parameter is \code{TRUE} and
\code{tdata} defines row names, these are printed in the left-most column and are displayed
as header cells. Keep in mind that \code{data.frame}s always define row names.}

\item{first.col.header}{Flag indicating if all cells in the first column must be displayed as header cells. Note
that, if both this parameter and \code{row.names} are \code{TRUE}, and \code{tdata} contains
row names, the constructed HTML table will have 2 columns of header cells.}

\item{indent}{Default indentation, in number of tabulation characters, to apply to HTML tags. This
indentation is also applied to \code{thead}.}

\item{tag.attrs}{Named \code{character} vector specifying the list of attributes to be set to the
\code{<table>} element. Setting this to \code{NULL} or an empty \code{character} vector
disables attributes.}

\item{thead}{\code{character} vector storing a table header to include. This can, for example, be a
\code{character} that defines column widths. Every element in this vector is written on a
separate line, applying the indentation given by \code{indent}.}

\item{tcaption}{Text to include as a caption below the table, or \code{NULL} if the table does not contain
caption.}

\item{na}{\code{character} to be used for printing \code{NA} values in the table. This parameter is
not considered when printing \code{thead} or the table's column names.}
}
\value{
The modified report, invisibly.
}
\description{
Generates HTML code for a table in the specified report.
}
\seealso{
\code{\link{rnb.add.tables}} for adding a listing of tables; \code{\linkS4class{Report}} for other functions
  adding contents to an HTML report
}
\author{
Yassen Assenov
}
