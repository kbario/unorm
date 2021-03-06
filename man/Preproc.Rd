% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{Preproc}
\alias{Preproc}
\title{Standard 1D NMR Preprocessing}
\usage{
Preproc(
  X,
  ppm,
  meta,
  baseline = T,
  flip = F,
  cali = F,
  calibrant = "tsp",
  lineWid = 1,
  lowCut = 0.25,
  watCut = c(4.5, 5),
  ureCut = c(5.6, 6),
  uppCut = 9.5,
  noi_sh = c(9.5, 11)
)
}
\arguments{
\item{X}{A matrix containing the non-preprocessed NMR spectral data. The rows should containing all values of a single experiment, and the columns, the values of the chemical shift variables.}

\item{ppm}{An array of the chemical shift variables, column matched to X.}

\item{meta}{The matrix of metadata pertaining to the X matrix. This is crucial for the TSP calibration and line width calculation.}

\item{baseline}{Default is set to \code{TRUE}. This argument calls on \code{\link[metabom8:bline]{metabom8::bcor()}}, a baseline correcting function to smooth the spectral baselines and remove the influence of broad peaks.}

\item{flip}{Default is set to \code{FALSE}. \strong{NOT CRUCIAL FOR STANDARDLY ACQUIRED SPECTRA} This function checks the orientation of spectra. This is a particularly important function for spectra derived from a \strong{small number of scans}; the NMR orientates spectra based on what side of the y-axis has the most signal and water has a very large impact on this and without a large amount of metabolite signal on the correct side of the spectra (as seen in spectra derived from small scans) the spectra will be orientated the wrong way.}

\item{cali}{Default is set to \code{FALSE}. This calls on the \code{metabom8} function \code{\link[metabom8:calibrate]{metabom8::calibrate()}} which ensures peaks are aligned  based on the TSP signal. Calibration is standardly performed by the NMR but this can be performed again for certainty.}

\item{calibrant}{This is the signal you wish to use as the calibrant. This signal must be compatible with \code{\link[metabom8:calibrate]{metabom8::calibrate()}}. For urine spectra, \emph{tsp} is the main calibrant.}

\item{lineWid}{This argument provides import for \code{\link[metabom8:lw]{metabom8::lw()}}, a function that calculates the line width of peaks. Spectra with peaks that have small line widths have sharper and more precise results which is more desirable. A maximum cutoff of 1 ensures spectra contain robust results. Consider omitting spectra with line widths over 1.}

\item{lowCut}{A single floating point number defining the ppm value that the lower limit of the spectra are trimmed to.}

\item{watCut}{The lower and upper ppm values concatenated, from which the water region will be trimmed and omitted. Water regions provide no important information and should be removed prior to data analysis. Default is set to \code{c(4.5,5)}}

\item{ureCut}{The lower and upper ppm values concatenated, from which the urea region will be trimmed and omitted. Urea regions also provide no important information and should be removed prior to data analysis. Default is set to \code{c(5.6,6)}}

\item{uppCut}{A single floating point number defining the ppm value that the upper limit of the spectra are trimmed to.}

\item{noi_sh}{The shift of the noise region used to calculate estimation of noise}
}
\value{
This function returns a list with:
\enumerate{
\item The processed X matrix in the first element,
\item The processed ppm array in the second element, and
\item The line width results in a data frame in the third element.
}
\itemize{
\item Following the example below will extract the results quickly and easily.
}
}
\description{
\code{\link[=Preproc]{Preproc()}} is a function aimed at streamlining the preprocessing stage of analysing NMR spectra. It harnesses the power of the package \code{metabom8} to quickly and easily process spectra.
}
\details{
\code{\link[=Preproc]{Preproc()}} is powered by the \code{metabom8} package and is simply a tool to harness and streamline functions in \code{metabom8} when preprocessing NMR spectra, not replace them.
}
\section{The Pipeline}{
This function streamlines the preprocessing of NMR urine spectra by combining a range of functions. It:
\enumerate{
\item Orientates the spectra correctly,
\item Calibrates the spectra by a specific peak,
\item Calculates the line widths of the peaks and returns a warning with the spectra that exceed the specified threshold,
\item Removes the lower, upper, water and urea regions of the spectra,
\item Corrects the baseline of the spectra using asymmetric least squares
\item Verifies that the resulting X, ppm and meta objects match appropriately.
}
}

\examples{
path = system.file('extdata', package = 'unorm')
}
\author{
\email{kylebario1@gmail.com}
}
\concept{{preproc}}
