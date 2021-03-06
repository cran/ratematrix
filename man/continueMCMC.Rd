% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/continueMCMC.R
\name{continueMCMC}
\alias{continueMCMC}
\title{Continue unfinished MCMC chain or add more generations}
\usage{
continueMCMC(handle, add.gen = NULL, save.handle = TRUE, dir = NULL)
}
\arguments{
\item{handle}{the output of 'ratematrixMCMC'.}

\item{add.gen}{number of generations to be added to a finished chain. If 'NULL' (default), the function will only continue unfinished chains.}

\item{save.handle}{whether to save the updated 'handle' object to the directory. This can overwrite the previous handle file.}

\item{dir}{an optional path to the output files. See 'Details'.}
}
\value{
Function will write the parameter values for each generation and the log to files. The new generations will be appended to the same files created by 'ratematrixMCMC'.
}
\description{
Function to continue an unfinished MCMC chain or to append more generations to a previously finished MCMC. It works by reading the last state of the chain and the tunning parameters of the previous chain, then restarting it from this step.
}
\details{
The function will append the new generations to the same files created by the prior run of the 'ratematrixMCMC' function. The function will, by default, search for files in the same directory of the previous run (see 'handle$dir'). However, you can provide a new path (relative or absolute path) to the argument 'dir'. The path provided to 'dir' will override the path pointed by 'handle$dir'. The new 'handle' output from 'continueMCMC' will have an updated total number of generations and will also update the directory path, if required.
}
\examples{
\donttest{
## Continue unfinished run.
data(centrarchidae)
handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=10000
                         , dir=tempdir())
## Now add generations to the same MCMC chain.
handle.add <- continueMCMC(handle=handle, add.gen=10000)
}
}
\author{
Daniel S. Caetano and Luke J. Harmon
}
