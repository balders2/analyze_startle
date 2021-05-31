

@echo off

echo Installing python...
pause
.\python-2.7.15.amd64.msi

echo Installing Matlab Runtime Compiler...
pause
.\MyAppInstaller_web.exe


echo Installing bioread...
pause
set PATH=%PATH%;C:\Python27
set PATH=%PATH%;C:\Python27\Scripts
cd \
cd c:\python27\Scripts\
pip install numpy
pip install scipy
pip install bioread

pause

