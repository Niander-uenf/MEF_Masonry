@echo off
setlocal
cd /d "%~dp0"
title MEF_Masonry - build the executable
echo Building the MEF_Masonry executable (dist\MEF_Masonry). It takes a few minutes.
echo.
set "PY="
where py >nul 2>&1
if not errorlevel 1 (
  py --version
  if not errorlevel 1 set "PY=py"
)
if not defined PY (
  where python >nul 2>&1
  if not errorlevel 1 (
    python --version
    if not errorlevel 1 set "PY=python"
  )
)
if not defined PY (
  echo Python was not found. Install Python 3.13 and run this file again.
  pause
  goto end
)
echo Checking the libraries...
%PY% -m pip install -r requirements.txt
%PY% build_exe.py
if errorlevel 1 (
  echo.
  echo The build failed. Copy the messages above and send them to the support.
) else (
  echo.
  echo Done. Copy the whole dist\MEF_Masonry folder to distribute the program.
)
pause
:end
endlocal
