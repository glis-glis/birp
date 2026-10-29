## Resubmission 8

This is a resubmission. All issues of Resubmission 7 have been adressed:
 1) The necessary CMAKE-flag not to write to ~/.cmake/packages/armadillo/ has been set
 2) Only source-files are retained from any repository, the arma/NOTICE.txt file has been removed
 3) We added all authors of third party libraries as authors in DESCRIPTION with "cph" and added a inst/COPYRIGHT file.

## Resubmission 7

This is a resubmission. All issues of Resubmission 6 have been adressed:
 * More fallbacks for cmake path
 * Use R-defined CC, CXX and CXXFLAGS
 * Remove unnecessary test-files (which caused warnings)

 We changed the version number to be consistent within all our birp-repositories. We also changed the license to MPL-2.0, which is the license we use in all our programs.

## Resubmission 6

This is a resubmission. Local builds work now on Linux, Mac and Windows.


## Resubmission 5

This is a resubmission. The package was initially accepted by CRAN, 
but a compilation error occurred on one Fedora server because the stattools library 
build was skipped. I suspect this was caused by a temporary network issue or configuration setting
on that server, which prevented CMake’s FetchContent from populating the library.
Prof. Ripley suggested that the problem could be related to parallel builds. 
I have tested the package using the official fedora:40 Docker image with multiple threads 
(Sys.setenv(MAKEFLAGS = "-j8")), and it compiles successfully.
Although I could not reproduce the error locally, I have modified the package so that, 
on CRAN servers, no downloads via FetchContent are performed during compilation. 
This ensures that network configuration or connectivity issues will not affect the build.

## Resubmission 4

This is a resubmission. Although the package was initially accepted by CRAN, 
there was a compilation error on its Fedora server, which was unfortunately not solved
by the previous resubmission 3. I was still not able to reproduce the error locally, 
but I have now removed the step that locates the R library on the system, 
which failed on the Fedora server of CRAN. 
This likely solves the issue - but as I cannot reproduce the error, it's hard to be sure.

## Resubmission 3

This is a resubmission. Although the package was initially accepted by CRAN, 
there was a compilation error on its Fedora server, and the package was thus removed.
Despite not being able to reproduce the error on two local Fedora systems as well as 
Fedora on R-hub v2, I think I figured out the issue and fixed it accordingly. 
I apologize for contacting a CRAN team member directly. I was simply replying to the email and was unaware that it was a policy violation.
In this version, I now have:

* Likely fixed the compilation error on the Fedora CRAN server (as I cannot reproduce the error, it's hard to confirm this, but I'm confident this will work).
* Added a vignette.
* Improved the documentation.

## Resubmission 2
This is a resubmission. Thank you for your comments. In this version I have:

* Replaced instances of cat() by message() to make them suppressible.
* Removed gzstream.*, as this was an artifact that was never used.
* Added more contributors in DESCRIPTION.
* Added inst/AUTHORS file that explains in detail the contributions of all authors.
* inst/AUTHORS also describes the use of the fast-float library.
* inst/AUTHORS is referenced in DESCRIPTION.
* Added URL to DESCRIPTION.

## Resubmission 1
This is a resubmission. Thank you for your comments. In this version I have:

* Omitted the redundant "Provides functions to" from the description.
* Added more details about package functionality and implemented methods in the Description text.
* Added a reference to the description field.
* Replaced all T and F with TRUE and FALSE, respectively.
* Made sure all functions documentations define \value.
* Removed \dontrun{} and added small files needed for the examples to inst/extdata.
* Added a verbose argument that allows to turn off messages to the console completely.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
