@echo off
rem Build the cfmsDrape CPU DLL out of tree for the engine-link gates L0 and L2 (docs/GATE_LINK_L0.md, GATE_LINK_L2.md).
rem Source: %CFMS_DRAPE_DIR% (default ..\cfmsDrape_dev_base378d20c if present, else ..\cfmsDrape_dev), its drape\ CMake project. The engine repo is only read.
rem Output: %TEMP%\cfmsDrape_build_<short commit>\Release\cfmsDrape.dll, or %TEMP%\cfmsDrape_build_<tag> when
rem CFMS_DRAPE_BUILD_TAG is set (the sxjacobian-link worktree uses tag sxjlink so that an uncommitted branch at
rem 378d20c never shares a build folder with main). Only the cfmsDrape target is built;
rem the engine's cfms_deploy target (copy to the shared predist folder) is never invoked.
setlocal
set "ROOT=%~dp0.."
rem The baseline engine (hand kernels) is commit 378d20c. Engine main moved on from it on 2026-11-13, so a worktree
rem pinned to that commit (..\cfmsDrape_dev_base378d20c) is used when it exists.
if not defined CFMS_DRAPE_DIR if exist "%ROOT%\..\cfmsDrape_dev_base378d20c\drape\CMakeLists.txt" set "CFMS_DRAPE_DIR=%ROOT%\..\cfmsDrape_dev_base378d20c"
if not defined CFMS_DRAPE_DIR set "CFMS_DRAPE_DIR=%ROOT%\..\cfmsDrape_dev"
if not exist "%CFMS_DRAPE_DIR%\drape\CMakeLists.txt" (
  echo engine repo not found: %CFMS_DRAPE_DIR%
  exit /b 1
)
for /f %%c in ('git -C "%CFMS_DRAPE_DIR%" rev-parse --short HEAD') do set "COMMIT=%%c"
set "B=%TEMP%\cfmsDrape_build_%COMMIT%"
if defined CFMS_DRAPE_BUILD_TAG set "B=%TEMP%\cfmsDrape_build_%CFMS_DRAPE_BUILD_TAG%"
if not defined VCINSTALLDIR (
  call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
)
cmake -S "%CFMS_DRAPE_DIR%\drape" -B "%B%" -G "Visual Studio 17 2022" -A x64 || exit /b 1
cmake --build "%B%" --config Release --target cfmsDrape -- /m || exit /b 1
echo %B%\Release\cfmsDrape.dll
exit /b 0
