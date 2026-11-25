@echo off
REM  Build OligoScan.exe. Logs everything to build_log.txt.
REM
REM  Uses ONE interpreter throughout, resolved from sys.executable rather than
REM  from `where python`: on this machine `where python` returns a shim
REM  (AppData\Local\Python\bin\python.exe) while code actually runs under
REM  pythoncore-3.14-64\python.exe. Installing with one and running with the
REM  other loses the package silently.
setlocal enabledelayedexpansion
cd /d "%~dp0"
set LOG=build_log.txt
echo OligoScan build %DATE% %TIME% > %LOG%
echo Working directory: %CD% >> %LOG%

for /f "delims=" %%P in ('python -c "import sys;print(sys.executable)"') do set PY=%%P
if not defined PY ( echo Python not found. & echo python not found >> %LOG% & goto :done )
echo Interpreter: !PY! & echo Interpreter: !PY! >> %LOG%
"!PY!" -c "import sys;print(sys.version)" >> %LOG% 2>&1

echo Installing build dependencies...
"!PY!" -m pip install --upgrade pyinstaller flask >> %LOG% 2>&1
if errorlevel 1 ( echo   pip install FAILED - see %LOG% & goto :done )

echo Confirming PyInstaller is importable by THIS interpreter...
"!PY!" -m PyInstaller --version >> %LOG% 2>&1
if errorlevel 1 (
  echo   PyInstaller still not importable - see %LOG%
  echo PyInstaller not importable after install >> %LOG%
  goto :done
)
for /f "delims=" %%V in ('"!PY!" -m PyInstaller --version 2^>nul') do echo   PyInstaller %%V

echo Verifying the scan engines...
build_win\oligoscan-v2.exe --selftest >> %LOG% 2>&1
if errorlevel 1 ( echo   v2 FAILED SELFTEST - not shipping this & goto :done )
echo   v2 OK
build_win\oligoscan-v3.exe --selftest >> %LOG% 2>&1
if errorlevel 1 ( echo   v3 failed - expected without AVX2, app falls back to v2 ) else ( echo   v3 OK )

echo Building...
"!PY!" -m PyInstaller --clean --noconfirm --workpath build_pyi --distpath dist web\OligoScan.spec >> %LOG% 2>&1
echo PyInstaller exit code: !ERRORLEVEL! >> %LOG%

echo. >> %LOG%
echo --- dist --- >> %LOG%
dir dist >> %LOG% 2>&1

if exist dist\OligoScan.exe (
  if not exist dist\references mkdir dist\references
  echo.
  for %%A in (dist\OligoScan.exe) do echo   BUILT: dist\OligoScan.exe  (%%~zA bytes^)
  echo   Put a .osr reference in dist\references\ then double-click it.
) else (
  echo.
  echo   BUILD DID NOT PRODUCE dist\OligoScan.exe - see %LOG%
)

:done
echo.
echo  Log: %~dp0build_log.txt
endlocal
