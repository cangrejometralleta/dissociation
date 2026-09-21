@echo off
rem Builds the Package, then Runs the Demo. The logic lives in test.sh.
call "%~dp0run-bash.cmd" test.sh %*
exit /b %errorlevel%
