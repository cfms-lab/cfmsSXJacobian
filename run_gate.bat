@echo off
rem cfmsSXJacobian_dev: gates (docs/MODEL_PROBLEM.md = start, docs/GATE_SPEED.md = speed, docs/GATE_ASSEMBLY.md = assembly,
rem   docs/GATE_ENGINE.md = engine hand-code comparison, docs/GATE_STRETCH_SHEAR.md = stretch and shear).
rem   1) generate C++ from the specs  2) samples  3) Mathematica reference (skipped if absent)
rem   4) speed gate S3 (flat functions unchanged)  4b) engine gate E1: extract engine code byte for byte
rem      (tests\engine_ref\extract_engine.py; skipped if the engine repo is absent)
rem   5) MSVC build (Release, C++14 switch, code is C++11; engine code compiled without /W4 /WX)
rem   6) start gate G1-G4  7) speed gate S1  8) assembly gate R1'-R2
rem   9) timing gates S2 and R3' over 8 fixed load addresses (tests\timing_bases.py, about 2 min)
rem  10) engine gate E2 (test_engine.exe), timing gate E4 (tests\timing_bases.py engine)
rem  11) stretch and shear gate T1-T5 (test_stretch_shear.exe, runs right after the assembly gate)
rem  12) Bridson bending gate (docs/GATE_BRIDSON.md): spec B1-B4 always, engine B5-B7 and timing gate B8
rem      (tests\timing_bases.py bridson) only when the engine repo is present
rem  13) Gauss-Newton gate (docs/GATE_GAUSS_NEWTON.md): N1 always (test_gauss_newton.exe), timing gate N2
rem      (tests\timing_bases.py gn) only when the engine repo is present
rem  14) integration comparison (docs/GATE_INTEGRATION.md): I0-I2 (test_integration.exe, about 40 s), then timing (c) and the
rem      recommendation per force (tests\timing_bases.py integ, about 1 min). Engine Bridson is included when the repo is present
rem  15) stretch tension-only gate (docs/GATE_TENSION.md): P0 (tests\check_tension_p0.py), P1-P2 (test_tension.exe),
rem      integration (test_integration.exe tension, about 80 s), timing (c') and the recommendation (tests\timing_bases.py tension)
rem  16) shear Gauss-Newton re-comparison with stretch T (docs/GATE_SHEAR_GN.md): integration (test_integration.exe shear_gn),
rem      timing (c'') and the recommendation (tests\timing_bases.py shear_gn)
rem  17) engine link L1 (docs/GATE_LINK_L1.md): L1-0 (tests\check_link_l1_p0.py), L1-1 - L1-3 (test_link_l1.exe)
rem  18) engine link L0 (docs/GATE_LINK_L0.md), only when the engine repo is there: builds the cfmsDrape CPU DLL
rem      out of tree (tests\build_engine_dll.bat -> %TEMP%\cfmsDrape_build_<commit>) and probes it (tests\link_l0.py)
rem  19) engine link L2 (docs/GATE_LINK_L2.md), only when the engine repo and the sxjacobian-link worktree are there
rem      (%CFMS_DRAPE_LINK_DIR%, default ..\cfmsDrape_dev_sxjacobian-link): builds the branch DLL out of tree
rem      (-> %TEMP%\cfmsDrape_build_sxjlink) and compares it with the step 18 main DLL (tests\link_l2.py)
rem  20) engine link L4 bending (docs/GATE_LINK_L4.md), same condition as 19: integrates two cloth scenes with the
rem      step 19 branch DLL (engine Oh bending vs generated Bridson GN) and times the engine default setting on
rem      8 rebased DLL copies (tests\link_l4.py, uses editbin); since the 2nd (c) attempt (docs/GATE_LINK_L4.md)
rem      it writes results\gate_link_l4_2.json and compares each run with the 1st attempt results\gate_link_l4.json
rem  21) engine link L4R bending (docs/GATE_LINK_L4R.md), same condition as 19, right after 20: the same runs measured
rem      against small-h reference solutions (tests\link_l4r.py, 30 reference runs, about 12 min); time (c) is the
rem      fixed 2nd-attempt value, it is not measured again
rem  22) Provot shear gate (docs/GATE_PROVOT_SHEAR.md): spec provot_shear.py, V0 = the unchanged-variant checks above,
rem      V1-V7 (test_provot_shear.exe, runs right after the stretch and shear gate). First step of engine link L3
rem  23) engine link L3 tension and shear (docs/GATE_LINK_L3.md), same condition as 19, right after 19: the branch's
rem      generated tension / Provot shear kernels (CFMS_CPU_SXJ_MEMBRANE=on) against the engine hand kernels of the
rem      same DLL, and switch-off against the step 18 main DLL (tests\link_l3.py); reads the step 17 and 19 results
rem  24) Poisson force gate (docs/GATE_POISSON.md): spec sul_poisson.py (Sul, Kim, Park 2010), Q0 = the unchanged-variant
rem      checks above and tests\check_poisson_p0.py, Q1-Q7 (test_poisson.exe, runs right after the Provot shear gate)
rem  25) Poisson uniaxial-tension integration (docs/GATE_POISSON_INTEGRATION.md): static balance U0a, global assembly U0b,
rem      96 runs each of the roller-end scene (judged) and the clamped scene (reported), Poisson stiffness exact vs GN
rem      (test_poisson_integration.exe, a few seconds, runs right after the Poisson force gate)
rem Build output goes to %TEMP%\cfmsSXJacobian_build. Numbers go to results\*.json.
rem   run_gate.bat           full run
rem   run_gate.bat nomma     reuse results\mma_ref.txt and results\mma_ref_tri.txt (Mathematica takes about 30 s)
setlocal
set "ROOT=%~dp0"
set "OUT=%TEMP%\cfmsSXJacobian_build"
if not exist "%OUT%" mkdir "%OUT%"
pushd "%ROOT%"

python codegen\sxgen.py codegen\specs\baraff_bend.py --out generated --results results || goto :fail
python codegen\sxgen.py codegen\specs\baraff_bend_compact.py --out generated --results results || goto :fail
python codegen\sxgen.py codegen\specs\baraff_stretch.py --out generated --results results || goto :fail
python codegen\sxgen.py codegen\specs\baraff_shear.py --out generated --results results || goto :fail
python codegen\sxgen.py codegen\specs\bridson_bend.py --out generated --results results || goto :fail
python codegen\sxgen.py codegen\specs\provot_shear.py --out generated --results results || goto :fail
python codegen\sxgen.py codegen\specs\sul_poisson.py --out generated --results results || goto :fail
python model\make_hinges.py results\hinges.txt || goto :fail
python model\make_triangles.py results\triangles.txt || goto :fail
python model\make_provot_diag.py results\provot_diag.txt || goto :fail
python model\make_poisson_beta.py results\poisson_beta.txt || goto :fail

set "WS=C:\Program Files\Wolfram Research\WolframScript\wolframscript.exe"
if /I "%~1"=="nomma" goto :s3
if exist "%WS%" (
  "%WS%" -file model\baraff_bend_ref.wl 200 || goto :fail
  "%WS%" -file model\baraff_stretch_shear_ref.wl 200 || goto :fail
  "%WS%" -file model\bridson_bend_ref.wl 200 || goto :fail
  "%WS%" -file model\provot_shear_ref.wl 200 || goto :fail
  "%WS%" -file model\sul_poisson_ref.wl 200 || goto :fail
  "%WS%" -file model\sul_poisson_uniaxial_ref.wl || goto :fail
) else (
  echo wolframscript not found - G1, T1, B1, V1 and Q1 use the committed results\mma_ref*.txt
)

:s3
python tests\check_flat_unchanged.py || goto :fail
python tests\engine_ref\extract_engine.py
set "EXT=%ERRORLEVEL%"

if not defined VCINSTALLDIR (
  call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
)
set FLAGS=/nologo /EHsc /W4 /WX /permissive- /utf-8 /Zc:__cplusplus /std:c++14 /O2 /MD /DNDEBUG /bigobj
cl %FLAGS% /c third_party\sxmesh\SXMesh.cpp /Fo"%OUT%\SXMesh.obj" || goto :fail
cl %FLAGS% tests\test_baraff_bend.cpp "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_baraff_bend.exe" || goto :fail
cl %FLAGS% tests\test_speed.cpp "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_speed.exe" || goto :fail
cl %FLAGS% tests\test_assembly.cpp "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_assembly.exe" || goto :fail
cl %FLAGS% tests\test_stretch_shear.cpp "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_stretch_shear.exe" || goto :fail
cl %FLAGS% tests\test_tension.cpp "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_tension.exe" || goto :fail
cl %FLAGS% tests\test_provot_shear.cpp "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_provot_shear.exe" || goto :fail
cl %FLAGS% tests\test_poisson.cpp "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_poisson.exe" || goto :fail
cl %FLAGS% tests\test_poisson_integration.cpp "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_poisson_integration.exe" || goto :fail
cl %FLAGS% tests\test_link_l1.cpp "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_link_l1.exe" || goto :fail
cl %FLAGS% tests\link_l0_ref.cpp /Fo"%OUT%\\" /Fe"%OUT%\link_l0_ref.exe" || goto :fail
if not "%EXT%"=="0" goto :specbuild
set EFLAGS=/nologo /EHsc /W3 /utf-8 /O2 /MD /DNDEBUG /bigobj
cl %EFLAGS% /c tests\engine_ref\engine_macri.cpp /Fo"%OUT%\engine_macri.obj" || goto :fail
cl %EFLAGS% /Itests\engine_ref /c generated\engine_ref\SDEmysse.cpp /Fo"%OUT%\SDEmysse.obj" || goto :fail
cl %FLAGS% tests\test_engine.cpp "%OUT%\engine_macri.obj" "%OUT%\SDEmysse.obj" "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_engine.exe" || goto :fail
cl %FLAGS% /DSX_WITH_ENGINE tests\test_bridson.cpp "%OUT%\engine_macri.obj" "%OUT%\SDEmysse.obj" "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_bridson.exe" || goto :fail
cl %FLAGS% /DSX_WITH_ENGINE tests\test_gauss_newton.cpp "%OUT%\engine_macri.obj" "%OUT%\SDEmysse.obj" "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_gauss_newton.exe" || goto :fail
cl %FLAGS% /DSX_WITH_ENGINE tests\test_integration.cpp "%OUT%\engine_macri.obj" "%OUT%\SDEmysse.obj" "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_integration.exe" || goto :fail
goto :skipbuild
:specbuild
cl %FLAGS% tests\test_bridson.cpp "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_bridson.exe" || goto :fail
cl %FLAGS% tests\test_gauss_newton.cpp "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_gauss_newton.exe" || goto :fail
cl %FLAGS% tests\test_integration.cpp "%OUT%\SXMesh.obj" /Fo"%OUT%\\" /Fe"%OUT%\test_integration.exe" || goto :fail
:skipbuild
rem Every gate runs even if an earlier one fails; the summary line lists the failed ones.
set "FAILS="
"%OUT%\test_baraff_bend.exe" || set "FAILS=%FAILS% start"
"%OUT%\test_speed.exe" || set "FAILS=%FAILS% speed-S1"
"%OUT%\test_assembly.exe" || set "FAILS=%FAILS% assembly-R1'-R2"
"%OUT%\test_stretch_shear.exe" || set "FAILS=%FAILS% stretch-shear-T1-T5"
"%OUT%\test_provot_shear.exe" || set "FAILS=%FAILS% provot-shear-V1-V7"
python tests\check_poisson_p0.py || set "FAILS=%FAILS% poisson-Q0"
"%OUT%\test_poisson.exe" || set "FAILS=%FAILS% poisson-Q1-Q7"
"%OUT%\test_poisson_integration.exe" || set "FAILS=%FAILS% poisson-integration"
"%OUT%\test_bridson.exe" || set "FAILS=%FAILS% bridson-B1-B5"
"%OUT%\test_gauss_newton.exe" || set "FAILS=%FAILS% gauss-newton-N1"
"%OUT%\test_integration.exe" || set "FAILS=%FAILS% integration-I0"
python tests\check_tension_p0.py || set "FAILS=%FAILS% tension-P0"
"%OUT%\test_tension.exe" || set "FAILS=%FAILS% tension-P1-P2"
"%OUT%\test_integration.exe" tension || set "FAILS=%FAILS% tension-integration"
"%OUT%\test_integration.exe" shear_gn || set "FAILS=%FAILS% shear-gn-integration"
python tests\timing_bases.py "%OUT%" || set "FAILS=%FAILS% timing-S2-R3'"
python tests\timing_bases.py "%OUT%" integ || set "FAILS=%FAILS% timing-integ"
python tests\timing_bases.py "%OUT%" tension || set "FAILS=%FAILS% timing-tension"
python tests\timing_bases.py "%OUT%" shear_gn || set "FAILS=%FAILS% timing-shear-gn"
python tests\check_link_l1_p0.py || set "FAILS=%FAILS% link-L1-0"
"%OUT%\test_link_l1.exe" || set "FAILS=%FAILS% link-L1"
if "%EXT%"=="3" goto :skipengine
if "%EXT%"=="0" goto :runengine
set "FAILS=%FAILS% engine-E1"
goto :summary
:runengine
"%OUT%\test_engine.exe" || set "FAILS=%FAILS% engine-E2"
python tests\timing_bases.py "%OUT%" engine || set "FAILS=%FAILS% timing-E4"
python tests\timing_bases.py "%OUT%" bridson || set "FAILS=%FAILS% timing-B8"
python tests\timing_bases.py "%OUT%" gn || set "FAILS=%FAILS% timing-N2"
set "EDIR=%CFMS_DRAPE_DIR%"
if not defined EDIR if exist "%ROOT%..\cfmsDrape_dev_base378d20c\drape\CMakeLists.txt" set "EDIR=%ROOT%..\cfmsDrape_dev_base378d20c"
if not defined EDIR set "EDIR=%ROOT%..\cfmsDrape_dev"
for /f %%c in ('git -C "%EDIR%" rev-parse --short HEAD') do set "ECOMMIT=%%c"
call tests\build_engine_dll.bat > "%OUT%\engine_dll_build.log" 2>&1 || set "FAILS=%FAILS% link-L0-build"
python tests\link_l0.py "%TEMP%\cfmsDrape_build_%ECOMMIT%\Release\cfmsDrape.dll" "%OUT%\link_l0_ref.exe" || set "FAILS=%FAILS% link-L0"
set "LDIR=%CFMS_DRAPE_LINK_DIR%"
if not defined LDIR set "LDIR=%ROOT%..\cfmsDrape_dev_sxjacobian-link"
if not exist "%LDIR%\drape\CMakeLists.txt" goto :skiplink
call :buildlink || set "FAILS=%FAILS% link-L2-build"
python tests\link_l2.py "%TEMP%\cfmsDrape_build_%ECOMMIT%\Release\cfmsDrape.dll" "%TEMP%\cfmsDrape_build_sxjlink\Release\cfmsDrape.dll" "%OUT%\link_l0_ref.exe" || set "FAILS=%FAILS% link-L2"
python tests\link_l3.py "%TEMP%\cfmsDrape_build_%ECOMMIT%\Release\cfmsDrape.dll" "%TEMP%\cfmsDrape_build_sxjlink\Release\cfmsDrape.dll" || set "FAILS=%FAILS% link-L3"
python tests\link_l4.py "%TEMP%\cfmsDrape_build_sxjlink\Release\cfmsDrape.dll" --out results\gate_link_l4_2.json || set "FAILS=%FAILS% link-L4"
python tests\link_l4r.py "%TEMP%\cfmsDrape_build_sxjlink\Release\cfmsDrape.dll" || set "FAILS=%FAILS% link-L4R"
goto :summary
:skiplink
echo engine link worktree not found (%LDIR%) - engine link L2, L3, L4 and L4R skipped
goto :summary
:skipengine
echo engine repo not found - engine gate E1-E4, Bridson engine part B5-B8, Gauss-Newton timing N2 and engine link L0 - L4 skipped
:summary
if defined FAILS goto :gatefail

popd
echo GATE OK
exit /b 0

:gatefail
popd
echo GATE FAILED:%FAILS%
exit /b 1

:fail
popd
echo GATE FAILED
exit /b 1

rem Build the sxjacobian-link worktree DLL without leaking CFMS_DRAPE_DIR into the later steps.
:buildlink
setlocal
set "CFMS_DRAPE_DIR=%LDIR%"
set "CFMS_DRAPE_BUILD_TAG=sxjlink"
call tests\build_engine_dll.bat > "%OUT%\engine_link_dll_build.log" 2>&1
endlocal & exit /b %ERRORLEVEL%
