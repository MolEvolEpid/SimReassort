% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{internals}
\alias{internals}
\alias{current_popsize}
\alias{get_roottime}
\alias{get_tiptimes}
\alias{get_discrepancies}
\alias{drop_node}
\alias{nwk_recur}
\alias{write_nwk}
\title{internals}
\usage{
current_popsize(x)

get_roottime(x, prune = FALSE, init.include = FALSE)

get_tiptimes(x, prune = FALSE)

get_discrepancies(tr1, tr2)

drop_node(phy, node)

nwk_recur(p, df)

write_nwk(df)
}
\arguments{
\item{x}{\code{gpsim} object}

\item{prune}{logical}

\item{init.include}{logical}

\item{tr1}{character, newick}

\item{tr2}{character, newick}

\item{phy}{\code{phylo} object}

\item{p}{node label}

\item{df}{data frame}
}
\description{
internals
}
\keyword{internals}
