#!/bin/bash
# Install GDAL Python bindings in virtual environment

if [ -z "$VIRTUAL_ENV" ]; then
    echo "Error: Please activate your virtual environment first"
    echo "Example: source venv/bin/activate"
    exit 1
fi

echo "Installing GDAL Python bindings..."
echo "Virtual environment: $VIRTUAL_ENV"

# Set environment variables for Homebrew GDAL
export GDAL_CONFIG=/opt/homebrew/bin/gdal-config
export CPLUS_INCLUDE_PATH=/opt/homebrew/include
export C_INCLUDE_PATH=/opt/homebrew/include
export LIBRARY_PATH=/opt/homebrew/lib

# Install gdal and pyproj
pip install gdal pyproj

echo "Done! Try running your script again."
