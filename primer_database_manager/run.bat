
@echo off
set "SCRIPT_DIR=%~dp0"
cd /d "%SCRIPT_DIR%"
title Primer Database Manager
color 0A
echo ============================================================
echo   🧬 Primer Database Manager
echo ============================================================
echo.
echo Starting server...
echo Access the app at: http://127.0.0.1:5001
echo Press Ctrl+C to stop the server
echo.
python app.py
pause
