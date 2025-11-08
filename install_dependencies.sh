#!/bin/bash
# Helper script to install GDAL Python bindings on macOS with Homebrew

# Check if we're in a virtual environment
if [ -z "$VIRTUAL_ENV" ]; then
    echo "Warning: No virtual environment detected."
    echo "Please activate your virtual environment first:"
    echo "  source venv/bin/activate"
    exit 1
fi

echo "Installing GDAL Python bindings for Homebrew GDAL..."

# Set environment variables for Homebrew GDAL
export GDAL_CONFIG=/opt/homebrew/bin/gdal-config
export CPLUS_INCLUDE_PATH=/opt/homebrew/include
export C_INCLUDE_PATH=/opt/homebrew/include
export LIBRARY_PATH=/opt/homebrew/lib

# Install dependencies
pip install -r requirements.txt

echo "Installation complete!"

