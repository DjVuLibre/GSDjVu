@echo off
:: batch version by tumagonx <opensourcepack.blogspot.com>
:: assigned to the same license of original djvudigital

set infile=
set outfile=
set zip=
set ziparam=
set quiet=
set verbose=
set check=
set dryrun=
set help=
set lines=
set words=
set pdf=
set exact-color=
set threshold=
set bg-subsample=
set bg-slices=
set fg-colors=
set fg-image-colors=
set sepfile=
set epsf=
set psrotate=
set gsarg=
set cseparg=
set poppler=

::store whole arguments first
set argv= %* 

if "%~1" EQU "" (
echo Usage:  djvudigital [options] inputfile [outputfile]
echo More information is displayed by typing
echo.    djvudigital --help
exit /B 1
)

::mangling
set argv=%argv: --lines = --lines=1 %
set argv=%argv: --words = --words=1 %
set argv=%argv: --sepfile = --sepfile=1 %
set argv=%argv: --verbose = --verbose=1 %
set argv=%argv: --dryrun = --dryrun=1 %
set argv=%argv: --check = --check=1 %
set argv=%argv: --quiet = --quiet=1 %
set argv=%argv: --help = --help=1 %
set argv=%argv: -q = --quiet=1 %
set argv=%argv: -v = --verbose=1 %
set argv=%argv: --exact-color = --exact-color=1 %

::trims blankspace of options
set argv=%argv: --=---%>nul 2>&1
::trap remaining as dummy parameter
set argv=%argv: =---dummy=%>nul 2>&1
::save each argument
set argv=%argv:---=&set %>nul 2>&1

::7zsfx or plain mode
::if available switch to 8.3 path for executables location
if exist "%ROOT%" (
	for %%F in ("%ROOT%") do if exist "%%~sdpF" set ROOT=%%~sdpF
) else if exist "%~sdp0" (
	set ROOT=%~sdp0

) else set ROOT=%~dp0
::bundled files
set XML2DSED=%ROOT%xml2dsed.awk
set ZIP=%ROOT%gzip.exe
set AWK=%ROOT%gawk.exe

echo.

::help trap
if "%help%" EQU "1" goto :help

::find ghostscript
set GSDJVU=
for %%G in (gs.exe gsc.exe gswin32.exe gswin32c.exe) do if exist "%%~$PATH:G" set GSDJVU=%%~$PATH:G
if not  "%GS%" EQU "" if exist "%GS%" set GSDJVU=%GS%
if not  "%GSC%" EQU "" if exist "%GSC%" set GSDJVU=%GSC%
if not  "%GS_PROG%" EQU "" if exist "%GS_PROG%" set GSDJVU=%GS_PROG%
if not exist "%ROOT%gswin32c.exe" (
	if "%GSDJVU%" EQU "" (
		echo ERROR: Can't find Ghostscript
		exit /B 1
	)
) else set GSDJVU=%ROOT%gswin32c.exe

::find csepdjvu
set CSEPDJVU=
for %%G in (csepdjvu.exe) do if exist "%%~$PATH:G" set CSEPDJVU=%%~$PATH:G
if not exist "%ROOT%csepdjvu.exe" (
	if "%CSEPDJVU%" EQU "" (
		echo ERROR: Can't find csepdjuu
		exit /B 1
	)
) else set CSEPDJVU=%ROOT%csepdjvu.exe

::find pdftotext
set PDFTOTEXT=
for %%G in (pdftotext.exe) do if exist "%%~$PATH:G" set PDFTOTEXT=%%~$PATH:G
if exist "%ROOT%pdftotext.exe" set PDFTOTEXT=%ROOT%pdftotext.exe

::find djvused
set DJVUSED=
for %%G in (djvused.exe) do if exist "%%~$PATH:G" set DJVUSED=%%~$PATH:G
if exist "%ROOT%djvused.exe" set DJVUSED=%ROOT%djvused.exe

::check trap
if not "%check%" EQU "1" goto :findfile
set TEMPFILE=%TEMP%\%RANDOM%.txt
%GSDJVU% --help > %TEMPFILE%
if %ERRORLEVEL% EQU 0 (
	findstr "djvusep" %TEMPFILE%
) else (
	echo ERROR: Failed to execute Ghostscript
	del %TEMPFILE%>nul
	exit /B 1
)
if %ERRORLEVEL% EQU 1 (
	echo ERROR: %GSDJVU% have no djvu support
	del %TEMPFILE%>nul
	exit /B 1
) else del %TEMPFILE%>nul


:findfile
if "%~1" EQU "" (
	if "%infile%" EQU "" (
		echo ERROR: source file is missing
		exit /B 1
	)
) 
set buf=%~1
if "%~x1" EQU "" (
shift /1
goto :findfile
) else if "%buf:~0,1%" EQU "-" (
	shift /1
	goto :findfile
) else (
	if "%infile%" EQU "" (
		if not "%~3" EQU "" (
			shift /1
			goto :findfile
		)
		if not exist "%buf%" (
			echo ERROR: source file "%buf%" is not exist
			exit /B 1
		) else (
                        set inext=%~x1
			if /i not "%~x1" EQU ".pdf" if /i not "%~x1" EQU ".ps" if /i not "%~x1" EQU ".eps" if /i not "%~x1" EQU ".ai" if not "%zip%" EQU "" if /i not "%buf:~-7,7%" EQU ".pdf.gz" if /i not "%buf:~-6,6%" EQU ".ps.gz" if /i not "%buf:~-7,7%" EQU ".eps.gz" (
				echo ERROR: source filetype of "%~nx1" is not supported by Ghostscript
				exit /B 1
			)
			set infile=%~dpnx1
			if "%~2" EQU "" (
				set outfile=%~dpn1.djvu
				if /i "%~x1" EQU ".gz" for %%f in ("%~dpn1") do set outfile=%%~dpnf.djvu
				goto :checkargs
			)
			shift /1
			goto :findfile
		)
	) else if not exist "%~dp1" (
		echo ERROR: destination path "%~dp1" is not exist& exit /B 1
	) else (
		if /i not "%~x1" EQU ".djvu" echo ERROR: desination filetype of "%~nx1" is not supported by dejavu& exit /B 1
		set outfile=%buf%
		goto :checkargs
	)
)

:checkargs

::dpi trap
::TODO: numeric check and range
if "%dpi%" EQU "" set dpi=300

::verbose/quiet trap (quiet override verbose)
set csepverbosity=-v
set gsverbosity=
if "%verbose%" EQU "1" set csepverbosity=-vv
if "%quiet%" EQU "1" (
	set csepverbosity=
	set gsverbosity=-q
)

::poppler trap
set popplertext=0
set popplermeta=0
set dopoppler=0
if not "%poppler%" EQU "" for /f "delims=, tokens=1,2" %%F in ("%poppler%") do (
	if /i not "%%~F" EQU "" set poppler%%F=1
	if /i not "%%~G" EQU "" set poppler%%G=1
)
if %popplertext% EQU 1 set dopoppler=1
if %popplermeta% EQU 1 set dopoppler=1
if not "%poppler%" EQU "" if %dopoppler% EQU 0 echo ERROR: unrecognized option for --poppler& goto :help
if %dopoppler% EQU 1 (
        if /i not"%inext%" EQU ".pdf" (
		echo ERROR: input is not pdf while --poppler specified
		exit /B 1
	)
	if "%nopdftotext%" EQU "1" (
		echo ERROR: No pdftotext detected while --poppler specified
		exit /B 1
	)
)

::epsf trap
set gsepsf=
if "%epsf%" EQU "" (
set gsepsf=-dEPSCrop
) else (
if /i "%epsf%" EQU "no" set gsepsf=
if /i "%epsf%" EQU "ignore" set gsepsf=-dNOEPS
if /i "%epsf%" EQU "fit" set gsepsf=-dEPSFitPage
if /i "%epsf%" EQU "crop" set gsepsf=-dEPSCrop
)

::gsargs/csepargs trap
set csepargs=
set gsarg1=
set gsarg2=
set gsarg0=-sDEVICE=djvusep -dNOPAUSE -dBATCH -dSAFER
if "%words%" EQU "1" set gsarg0=%gsarg0% -dProvideUnicode -dExtractText
if "%lines%" EQU "1" (
	if "%words%" EQU "" set gsarg0=%gsarg0% -dProvideUnicode -dExtractText
	set csepargs=-t
)
if "%exact-color%" EQU "1" set gsarg0=-dUseCIEColor %gsarg0%
if not "%threshold%" EQU "" set gsarg0=-dThreshold=%threshold% %gsarg0%
if not "%bg-subsample%" EQU "" set gsarg0=-dBgSubsample=%bg-subsample% %gsarg0%
if not "%fg-colors%" EQU "" set gsarg0=-dFgColors=%fg-colors% %gsarg0%
if not "%fg-image-colors%" EQU "" set gsarg0=-dFgImgColors=%fg-image-colors% %gsarg0%
if not "%psrotate%" EQU "" if "%psrotate%" EQU "0" set gsarg2=-c "<< /Orientation 0 >> setpagedevice"
if not "%psrotate%" EQU "" if "%psrotate%" EQU "90" set gsarg2=-c "<< /Orientation 3 >> setpagedevice"
if not "%psrotate%" EQU "" if "%psrotate%" EQU "180" set gsarg2=-c "<< /Orientation 2 >> setpagedevice"
if not "%psrotate%" EQU "" if "%psrotate%" EQU "270" set gsarg2=-c "<< /Orientation 1 >> setpagedevice"
if not "%gsarg%" EQU "" set gsarg0=%gsarg:,= % %gsarg0%
if not "%bg-slices%" EQU "" set csepargs=-q %bg-slices% %csepargs%
if not "%cseparg%" EQU "" set csepargs=%csepargs% %cseparg:,= %

::gsprinted trap
if "%pdf%" EQU "" set gsprinted=-dPrinted
if not "%pdf%" EQU "" if "%pdf%" EQU "screen" set gsprinted=-dPrinted=false
if not "%pdf%" EQU "" if "%pdf%" EQU "printed" set gsprinted=-dPrinted=true

::sepfile trap
set backend=
if "%sepfile%" EQU "1" set outfile=%outfile:.djvu=.sep%
if "%sepfile%" EQU "1" (
set backend=-sOutputFile=%outfile%
) else set backend="-sOutputFile=|"%CSEPDJVU%" -d \"%dpi%\" -v - \"%outfile%\""
set ROOT=

if "%quiet%" EQU "" echo DJVUDIGITAL --- DjVuLibre-3.5
if /i "%infile:~-3,3%" EQU ".gz" (
	if "%dryrun%" EQU "1" (
		echo "%ZIP%" -c -d -q "%infile%" ^| "%GSDJVU%" -r%dpi% %gsverbosity% %gsprinted% %gsepsf% %backend% %gsarg0% %gsarg1% %gsarg2% -_ -c quit
	) else "%ZIP%"  -c -d -q "%infile%" | "%GSDJVU%" -r%dpi% %gsverbosity% %gsprinted% %gsepsf% %backend% %gsarg0% %gsarg1% %gsarg2% -_ -c quit
) else (
	if "%dryrun%" EQU "1" (
		echo "%GSDJVU%" -r%dpi% %gsverbosity% %gsprinted% %gsepsf% %backend% %gsarg0% %gsarg1% %gsarg2% -f "%infile%" -c quit
	) else "%GSDJVU%" -r%dpi% %gsverbosity% %gsprinted% %gsepsf% %backend% %gsarg0% %gsarg1% %gsarg2% -f "%infile%" -c quit
)
if %dopoppler% EQU 1 (
	if "%dryrun%" EQU "1" (
		echo "%PDFTOTEXT%" -bbox "%infile%" - ^| "%AWK%" -f "%XML2DSED%" dpi=%dpi% dometa=%popplermeta% dotext=%popplertext% ^| %DJVUSED% "%outfile%" -s
	) else "%PDFTOTEXT%" -bbox "%infile%" - | "%AWK%" -f "%XML2DSED%" dpi=%dpi% dometa=%popplermeta% dotext=%popplertext% | %DJVUSED% "%outfile%" -s
)
goto :EOF


:help
echo. Usage:  djvudigital [options] inputfile [outputfile]
echo. 
echo. OPTIONS
echo. --verbose, -v
echo.         Displays more informational messages while converting the file
echo. 
echo. --quiet, -q
echo.         Do not display informational messages while converting the file.
echo. 
echo. --dpi=resolution
echo.         Specify the desired resolution to resolution dots per inch.
echo.         The default is 300 dpi.
echo. 
echo. --psrotate=angle
echo.         Rotate the PostScript file by angle degrees clockwise.  
echo.         Only the values 0, 90, 180, and 270 are supported. 
echo.         This option only applies to PostScript files. PDF files are 
echo.         always converted according to their native orientation.
echo. 
echo. --epsf=disposition
echo.         Specify how to handle Encapsulated PostScript files.  
echo.         Argument disposition can take the values crop, fit, and ignore.
echo.         The default disposition crop creates a DjVu file whose size
echo.         matches the bounding box of the Encapsulated PostScript file.
echo.         Value fit rescales the graphics to the default page size.
echo.         Value ignore disables all Encapsulated PostScript specific code.
echo.         This option requires Ghostscript 7.07 or better.
echo. 
echo. --exact-color
echo.         Enables a more accurate rendering of the colors. 
echo.         This option requires Ghostscript 6.52 or better.
echo. 
echo. --threshold=thres
echo.         Specify a threshold for the foreground/background separation code.
echo.         Acceptable values of thres range from 0 to 100. Larger values 
echo.         place more information into the foreground layer.
echo.         The default threshold value is 80.
echo. 
echo. --bg-subsample=sub
echo.         Specify the background subsampling ratio. Argument sub must be
echo.         an integer between 1 and 6. The default value is 3.
echo. 
echo. --bg-slices=n+...+n
echo.         Specify the encoding quality of the background layer. The syntax
echo.         for the argument is similar to that described for the -slice 
echo.         option of command c44. The default is 72+11+10+10.
echo. 
echo. --fg-colors=ncolors
echo.         Specify the maximum number of distinct colors in the foreground
echo.         layer. Argument ncolors can take integer values between 
echo.         1 and 4000. The default value is 256.
echo. 
echo. --fg-image-colors=ncolors
echo.         Specify the  maximum number of distinct colors in an image 
echo.         for considering encoding it into the foreground layer.
echo.         Argument ncolors can take integer values between 1 and 4000.
echo.         The default value is 256.
echo. 
echo. --words
echo.         Extract the text from the PostScript code and incorporates 
echo.         this information into the DjVu file. This option records
echo.         the location of every word.
echo. 
echo. --lines
echo.         Extract the text from the PostScript code and incorporates
echo.         this  information into the DjVu file. This option saves
echo.         a few bytes by only recording the location of each line.
echo. 
echo. --poppler=arg1[,arg2]
echo.         Extract text overlay or metadata from existing pdf input using
echo.         Poppler's pdftotext and adding it into output using djvused
echo.         value can be text or meta.
echo. 
echo. --gsarg=arg1[,arg2,...,argN]
echo.         Insert extra arguments on the Ghostscript command line.
echo.
echo. --cseparg=arg1[,arg2,...,argN]
echo.         Insert extra arguments on the command line of program csepdjvu.
echo.
echo. --sepfile
echo.         Produces a separated data file instead of a DjVu file. Program
echo.         csepdjvu can then convert the separated data file into a DjVu file.
echo. 
echo. --check
echo.         Display the names of the two auxiliary programs found by 
echo.         djvudigital, namely a suitable Ghostscript interpreter and 
echo.         a suitable backend encoder.
echo. 
echo. --dryrun
echo.         Simply display the  Ghostscript command line generated by 
echo.         djvudigital without running it. No output file is produced.
echo. 
echo. --help  Display the manual page for djvudigital.
echo. 
