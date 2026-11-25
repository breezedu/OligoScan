@echo off
REM  Prove the Windows engine produces the same answers as the verified Linux
REM  build -- not just correct arithmetic (--selftest) but correct file I/O,
REM  .osr loading and output formatting, which selftest never touches.
setlocal
cd /d "%~dp0"

echo Running the demo scan from a FASTA...
build_win\oligoscan-v2.exe examples\demo_ref.fa examples\demo_oligos.tsv 1 3 4 anti 1 1 > build_win\demo_win_fa.tsv 2>nul

echo Packing the demo reference...
build_win\oligoscan-v2.exe pack examples\demo_ref.fa build_win\demo.osr >nul 2>&1

echo Running the same scan from the packed .osr...
build_win\oligoscan-v2.exe build_win\demo.osr examples\demo_oligos.tsv 1 3 4 anti 1 1 > build_win\demo_win_osr.tsv 2>nul

echo Same scan on the AVX2 build...
build_win\oligoscan-v3.exe build_win\demo.osr examples\demo_oligos.tsv 1 3 4 anti 1 1 > build_win\demo_win_v3.tsv 2>nul

echo.
echo  Row counts (all three must match, and match Linux):
for %%F in (build_win\demo_win_fa.tsv build_win\demo_win_osr.tsv build_win\demo_win_v3.tsv) do (
  for /f %%C in ('find /c /v "" ^< %%F') do echo    %%F : %%C lines
)
echo.
echo  FASTA vs .osr identical?
fc /b build_win\demo_win_fa.tsv build_win\demo_win_osr.tsv >nul && echo    YES || echo    NO - DIFFERENT
echo  v2 vs v3 identical?
fc /b build_win\demo_win_osr.tsv build_win\demo_win_v3.tsv >nul && echo    YES || echo    NO - DIFFERENT
echo.
echo  Paste the output above back, plus the first few lines of:
echo    build_win\demo_win_osr.tsv
endlocal
