@echo off

REM set compiler_flags=-Od -Z7 -FC -W4 -WX -wd4214 -wd4152 -wd4221 -wd4996 -wd4201 -wd4100 -wd4204 -nologo -D_DEBUG
set compiler_flags=-O2 -Z7 -FC -W4 -WX -wd4214 -wd4152 -wd4221 -wd4996 -wd4201 -wd4100 -wd4204 -nologo
set linker_flags=-incremental:no

if not exist build mkdir build
pushd build

cl %compiler_flags% -TC -LD -DEXPORT ..\no_gl.c /link %linker_flags% kernel32.lib user32.lib
cl %compiler_flags% -TC ..\viewer.c /link /out:viewer.exe %linker_flags% kernel32.lib user32.lib gdi32.lib comdlg32.lib no_gl.lib /subsystem:windows /ENTRY:mainCRTStartup

popd
