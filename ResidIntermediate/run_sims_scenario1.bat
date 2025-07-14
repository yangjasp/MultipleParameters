:: run all scenarios for Multiple Outcome Two-Phase sampling
@echo on
setlocal enabledelayedexpansion

:: set R version for your machine
:: set r_version=4.2.2
set r_version=4.3.1
For /F "Delims=" %%0 In ('where /r "C:\Program Files\R" Rscript.exe') do set scriptpath="%%~0"
echo Using R executable !scriptpath!
::set "scriptpath = C:\Users\T551275\Downloads\R\R-4.2.2\bin\Rscript.exe"

:: set up arguments to pass
set phase2size=1000
set phase1size=10000

:: send output to directory
set outdir=%cd%\rout

if not exist %outdir% mkdir %outdir%

:: Add
if not exist %outdir%\files mkdir %outdir%\files

:: run the simulation
for %%S in (1 2 3 4) do (
   set this_outfile=!outdir!\output_n%%S.out
   echo Running sim = %%S
   "C:\Program Files\R\R-4.2.2\bin\x64\Rscript.exe" run_simulations_scenario1.R --scenario %%S --pathname %cd%\rout\files --n !phase2size! --N !phase1size! 1>!this_outfile! 2>&1
  )
)

pause
