\name{aami}
\alias{aami3a}
\alias{aami3b}

\docType{data}
\title{ ANSI/AAMI EC13 Test Waveforms, 3a and 3b }
\description{
ANSI/AAMI EC13 Test Waveforms 3a and 3b, as obtained from the PhysioBank
database.
}
\usage{data(aami3a);
data(aami3b);}
\format{ 
 > str(aami3a)
 Time-Series [1:43081] from 0 to 59.8: 0.185 0.185 0.169 0.185 0.185
 0.185 0.185 0.185 0.185 0.2 ... \cr \cr
> str(aami3b)
 Time-Series [1:43142] from 0 to 59.9: 0.192 0.192 0.192 0.192 0.192 0.2 0.2 0.192 0.2 0.2 ...
} 

\details{
  

The following text is reproduced (abridged) from PhysioBank, page
\url{http://www.physionet.org/physiobank/database/aami-ec13/}.  Other
recordings belong to the same dataset and can be obtained from the same
page.



The files in this set can be used for testing a variety of devices that
monitor the electrocardiogram.  The recordings include both synthetic
and real waveforms. For details on these test waveforms and how to use
them, please refer to section 5.1.2.1, paragraphs (e) and (g) in the
reference below.  Each recording contains one ECG signal sampled at 720
Hz with 12-bit resolution.

}

\note{
   Timestamps in the datasets have been re-created at the indicated
  frequency of 720 Hz, whereas the original timestamps 
  only had three decimal digits' precision, and were therefore
  unevenly spaced. 
  }

\source{
  Cardiac monitors, heart rate meters, and alarms [American National Standard (ANSI/AAMI EC13:2002)]. Arlington, VA: Association for the Advancement of Medical Instrumentation, 2002
}
\references{

  Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG,
  Mietus JE, Moody GB, Peng CK, Stanley HE. \emph{PhysioBank, PhysioToolkit,
  and PhysioNet: Components of a New Research Resource for Complex
  Physiologic Signals.} Circulation 101(23):e215-e220 [Circulation
  Electronic Pages;
  \url{http://circ.ahajournals.org/cgi/content/full/101/23/e215}]; 2000
  (June 13).

}
\examples{

data(aami3a);
data(aami3b);

## Plot both as a multivariate TS object
##  only extract the fist 10 seconds

plot( main="ECG (mV)",
 window(
  cbind(aami3a,aami3b)   ,end=10)
)



}
\keyword{data}
