@echo off
echo ========================================
echo    QuickNGS - Installation Script
echo    Powered by Pourdad Panahi
echo ========================================
echo.

echo [1/2] Installing Python dependencies...
pip install flask>=2.3.0 biopython==1.79 pyyaml==6.0 plotly==5.9.0

if %errorlevel% neq 0 (
    echo.
    echo ERROR: Failed to install dependencies.
    echo Please check your internet connection and try again.
    pause
    exit /b 1
)

echo.
echo [2/2] Creating required folders...
if not exist "uploads" mkdir uploads
if not exist "results" mkdir results

echo.
echo ========================================
echo    Installation Complete!
echo.
echo    Run: python app.py
echo    Open: http://127.0.0.1:5002
echo ========================================
pause
