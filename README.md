# Garmin Custom Map Converter

A Python tool to convert GeoTIFF files into Garmin-compatible KMZ files for use with Garmin devices (such as Rino series).

## Features

- Converts GeoTIFF (.tif) files to JPEG images
- Automatically resizes images to meet Garmin requirements (max 1024x1024 pixels, max 3MB)
- Automatically detects and converts coordinate systems to WGS84 (required for KML)
- Extracts georeferencing information and creates KML files
- Transforms KML format to Garmin-compatible structure
- Packages everything into KMZ files ready for Garmin import

## Requirements

- Python 3.6 or higher
- GDAL library (includes coordinate transformation capabilities - see installation instructions below)
- Pillow library (for image resizing - installs automatically with requirements.txt)

## Installation

### 1. Install GDAL

GDAL installation varies by operating system:

#### macOS
```bash
# Using Homebrew (recommended)
brew install gdal

# Then install Python bindings in your virtual environment
# Make sure your virtual environment is activated first
source venv/bin/activate  # or your venv path

# Set environment variables for Homebrew GDAL
export GDAL_CONFIG=/opt/homebrew/bin/gdal-config
export CPLUS_INCLUDE_PATH=/opt/homebrew/include
export C_INCLUDE_PATH=/opt/homebrew/include
export LIBRARY_PATH=/opt/homebrew/lib

# Install dependencies
pip install -r requirements.txt

# Or use the helper script:
# ./install_dependencies.sh
```

#### Linux (Ubuntu/Debian)
```bash
sudo apt-get update
sudo apt-get install gdal-bin libgdal-dev
pip install gdal
```

#### Linux (Fedora/CentOS/RHEL)
```bash
sudo dnf install gdal gdal-devel
pip install gdal
```

#### Windows
Download GDAL binaries from [OSGeo4W](https://trac.osgeo.org/osgeo4w/) or use conda:
```bash
conda install -c conda-forge gdal
```

### 2. Install Python Dependencies

```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

Convert a GeoTIFF file to KMZ:
```bash
python garmin-map-converter.py input.tif
```

This will create `input.kmz` in the same directory as the input file.

### Specify Output File

```bash
python garmin-map-converter.py input.tif -o output.kmz
```

### Command-Line Options

```
positional arguments:
  input                 Path to input GeoTIFF file

optional arguments:
  -h, --help            Show help message and exit
  -o OUTPUT, --output OUTPUT
                        Path to output KMZ file (default: input filename with .kmz extension)
  --temp-dir TEMP_DIR   Temporary directory for intermediate files (default: system temp)
```

## How It Works

1. **TIF to JPG Conversion**: Uses GDAL to convert the GeoTIFF to a JPEG image
2. **Image Resizing**: Automatically resizes images to meet Garmin requirements:
   - Maximum dimensions: 1024x1024 pixels (maintains aspect ratio)
   - Maximum file size: 3MB (adjusts quality and size as needed)
3. **Coordinate System Detection**: Reads the coordinate reference system (CRS) from the GeoTIFF metadata using GDAL
4. **Coordinate Transformation**: Converts coordinates from the source CRS (e.g., UTM) to WGS84 (EPSG:4326) using GDAL's built-in `osr.CoordinateTransformation` - no additional packages needed!
5. **Georeferencing Extraction**: Extracts geographic bounds (north, south, east, west) in WGS84 decimal degrees
6. **KML Generation**: Creates an initial KML file with the georeferencing information
7. **KML Format Fix**: Transforms the KML to Garmin-compatible format:
   - Changes encoding to iso-8859-1
   - Updates namespace to kml 2.1
   - Wraps in Document element
   - Ensures proper structure for Garmin devices
8. **KMZ Packaging**: Packages the KML (named `doc.kml`) and JPG image into a KMZ (ZIP) file

## Output Format

The generated KMZ file contains:
- A KML file named `doc.kml` with Garmin-compatible format (coordinates in WGS84)
- A JPEG image file referenced by the KML

The KMZ file can be directly imported into Garmin BaseCamp or transferred to compatible Garmin devices.

## Coordinate System Support

The script automatically handles coordinate system conversion using GDAL's built-in `osr` module:
- **Detects CRS**: Reads the coordinate reference system from GeoTIFF metadata (WKT format)
- **Converts to WGS84**: Transforms coordinates to WGS84 (EPSG:4326) required for KML using `osr.CoordinateTransformation`
- **Handles projections**: Works with UTM, State Plane, and other projected coordinate systems
- **No extra dependencies**: Uses only GDAL - no need for pyproj or other packages!

If your GeoTIFF uses a projected coordinate system (like UTM), the script will automatically convert it to geographic coordinates (latitude/longitude) for the KML file.

## Troubleshooting

### GDAL Import Error
If you get `ImportError: No module named 'osgeo'`, make sure GDAL is properly installed:
- Verify GDAL installation: `gdalinfo --version`
- Reinstall Python bindings: `pip install --upgrade gdal`

### Coordinate Transformation Errors
If you encounter errors related to coordinate transformation:
- Ensure your GeoTIFF has proper CRS metadata (check with `gdalinfo yourfile.tif`)
- The script uses GDAL's built-in `osr` module for coordinate transformation - no additional packages needed
- If transformation fails, the script will warn you and return coordinates as-is (they may not be in WGS84)

### File Not Found Error
Ensure the input TIF file path is correct and the file exists.

### Permission Errors
Make sure you have write permissions in the output directory.

## License

This project is open source and available for use.

