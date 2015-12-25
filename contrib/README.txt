
This is a Windows .BAT version of the djvudigital shell script.  It was
contributed by <tumagonx@users.sourceforge.net> and is offered here without
guarantees.

It searches for csepdjvu.exe and ghostcript along the command PATH or in the
same directory as the bat file.  The ghostscript program is assumed to contain
the gsdjvu driver (djvusep). If this is not the case, the script will fail to
work.

If you want to process compressed file, you need to install the executable
gzip.exe in the same directory as the bat file.  If you want to use the option
--poppler={text,meta}, you also need to have the file xml2dsed.awk and the
executable gawk.exe (from gnuwin32) in the same directory as the batch file,
and the executables djvused.exe and pdftodjvu.exe myst be found along the PATH
or again in the same directory as the batch file.

Reference:
http://opensourcepack.blogspot.com/2015/06/windows-version-of-djvudigital.html



