@echo off
setlocal
cd /d "%~dp0"
title MEF_Masonry
echo ==================================================
echo   MEF_Masonry - block masonry at high temperature
echo ==================================================
echo.
echo [1/3] Looking for Python...
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
if not defined PY goto no_python
echo.
echo [2/3] Checking the numpy, scipy and matplotlib libraries...
%PY% -c "import numpy, scipy, matplotlib" >nul 2>&1
if not errorlevel 1 goto start
echo Installing the libraries. This only happens the first time and takes a few minutes...
%PY% -m pip install -r requirements.txt
if errorlevel 1 (
  echo Trying again in the user folder...
  %PY% -m ensurepip --upgrade >nul 2>&1
  %PY% -m pip install --user -r requirements.txt
)
%PY% -c "import numpy, scipy, matplotlib"
if errorlevel 1 goto pip_error
:start
echo.
echo [3/3] Opening the program in the browser.
echo Keep this window open while you use the program. To quit, close this window.
echo.
%PY% app.py %*
if errorlevel 1 goto app_error
goto end

:no_python
echo.
echo Python was not found on this computer.
echo  1. If you downloaded the python-manager .msix, double-click it and click Install.
echo     Then open the Command Prompt and run:  py install 3.13
echo  2. Or install it from python.org, ticking "Add python.exe to PATH".
echo Then close this window and double-click this file again.
echo.
pause
goto end

:pip_error
echo.
echo The numpy, scipy and matplotlib libraries could not be installed.
echo Check the internet connection and try again. If it still fails,
echo copy the messages above and send them to the support.
echo.
pause
goto end

:app_error
echo.
echo The program stopped with an error. Copy the messages above and send them to the support.
echo.
pause

:end
endlocal
