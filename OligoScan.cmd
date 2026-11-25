@echo off
REM  OligoScan desktop application. Plain Python, no packaging, no bundling,
REM  nothing for endpoint protection to remove.
setlocal
cd /d C:\WebUI\oligoscan
where python >nul 2>&1 || ( echo Python not found on PATH. & pause & exit /b 1 )
python -c "import tkinter" 2>nul || (
  echo tkinter is missing. Re-run the Python installer and tick "tcl/tk and IDLE".
  pause & exit /b 1
)
python -c "import numpy" 2>nul || (
  echo Installing numpy ^(needed only for the ASO GTF annotation^)...
  python -m pip install --quiet numpy
)
python oligoscan_app.py
if errorlevel 1 pause
endlocal
