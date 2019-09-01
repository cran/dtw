
##
## Copyright (c) 2006-2019 of Toni Giorgino
##
## This file is part of the DTW package.
##
## DTW is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## DTW is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
## or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
## License for more details.
##
## You should have received a copy of the GNU General Public License
## along with DTW.  If not, see <http://www.gnu.org/licenses/>.
##


#' ANSI/AAMI EC13 Test Waveforms, 3a and 3b
#' 
#' ANSI/AAMI EC13 Test Waveforms 3a and 3b, as obtained from the PhysioBank
#' database.
#' 
#' 
#' The following text is reproduced (abridged) from PhysioBank, page
#' <https://www.physionet.org/content/aami-ec13/1.0.0/>.  Other
#' recordings belong to the dataset and can be obtained from the same page.
#' 
#' The files in this set can be used for testing a variety of devices that
#' monitor the electrocardiogram.  The recordings include both synthetic and
#' real waveforms. For details on these test waveforms and how to use them,
#' please refer to section 5.1.2.1, paragraphs (e) and (g) in the reference
#' below.  Each recording contains one ECG signal sampled at 720 Hz with 12-bit
#' resolution.
#' 
#' @note Timestamps in the datasets have been re-created at the indicated
#' frequency of 720 Hz, whereas the original timestamps in ms (at least in text
#' format) only had three decimal digits' precision, and were therefore
#' affected by substantial jittering.
#' @name aami
#' @aliases aami3a aami3b
#' @docType data
#' @format Time-series objects (class `ts`). 
#' @references 
#' * Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh,
#' Mark RG, Mietus JE, Moody GB, Peng CK, Stanley HE. *PhysioBank,
#' PhysioToolkit, and PhysioNet: Components of a New Research Resource for
#' Complex Physiologic Signals.* Circulation 101(23):e215-e220; 2000 (June
#' 13).
#' * Cardiac monitors, heart rate meters, and alarms; American National
#' Standard (ANSI/AAMI EC13:2002). Arlington, VA: Association for the
#' Advancement of Medical Instrumentation, 2002.
#' @source   <https://www.physionet.org/content/aami-ec13/1.0.0/>
#' @keywords datasets
#' @examples
#' data(aami3a);
#' data(aami3b);
#' 
#' ## Plot both as a multivariate TS object
#' ##  only extract the first 10 seconds
#' 
#' plot( main="ECG (mV)",
#'  window(
#'   cbind(aami3a,aami3b)   ,end=10)
#' )
#' 
#' 
#' 
#' 
NULL

