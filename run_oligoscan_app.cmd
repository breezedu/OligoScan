@echo off
REM  Launch the OligoScan desktop tool. No build, no packaging, nothing for
REM  endpoint protection to object to -- it is just Python running locally.
setlocal
cd /d C:\WebUI\oligoscan
where python >nul 2>&1 || ( echo Python not found on PATH. & pause & exit /b 1 )
python -c "import flask" 2>nul || (
  echo Installing Flask ^(one-time^)...
  python -m pip install --quiet flask || ( echo pip failed & pause & exit /b 1 )
)
python web\desktop.py
pause
endlocal
