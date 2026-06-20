
@echo off
set "SCRIPT_DIR=%~dp0"
cd /d "%SCRIPT_DIR%"
title Primer Database Manager - Installer
color 0A

echo ============================================================
echo   🧬 Primer Database Manager - Installation
echo ============================================================
echo.

echo [1/3] Checking Python installation...
python --version >nul 2>&1
if errorlevel 1 (
    color 0C
    echo.
    echo [ERROR] Python is not installed or not in PATH.
    echo Please install Python 3.8 or higher from:
    echo https://www.python.org/downloads/
    echo.
    echo Make sure to check "Add Python to PATH" during installation.
    echo.
    pause
    exit /b 1
)
echo [OK] Python found:
python --version
echo.

echo [2/3] Installing required packages...
echo This may take a few minutes...
echo.

pip install flask werkzeug -i https://pypi.tuna.tsinghua.edu.cn/simple --trusted-host pypi.tuna.tsinghua.edu.cn >nul 2>&1

if errorlevel 1 (
    echo [WARNING] Mirror failed. Trying without mirror...
    pip install flask werkzeug
)

if errorlevel 1 (
    color 0C
    echo [ERROR] Failed to install packages.
    echo Please check your internet connection and try again.
    pause
    exit /b 1
)

echo [OK] All packages installed successfully.
echo.

echo [3/3] Creating run script...
(
echo @echo off
echo set "SCRIPT_DIR=%%~dp0"
echo cd /d "%%SCRIPT_DIR%%"
echo title Primer Database Manager
echo color 0A
echo echo ============================================================
echo echo   🧬 Primer Database Manager
echo echo ============================================================
echo echo.
echo echo Starting server...
echo echo Access the app at: http://127.0.0.1:5001
echo echo Press Ctrl+C to stop the server
echo echo.
echo python app.py
echo pause
) > run.bat

echo [OK] run.bat created.
echo.

echo ============================================================
echo   ✅ Installation completed successfully!
echo ============================================================
echo.
echo To run the program:
echo   1. Double-click on "run.bat"
echo   2. Open browser and go to http://127.0.0.1:5001
echo.
echo Default login:
echo   Username: admin
echo   Password: admin123
echo.
echo ============================================================
pause
