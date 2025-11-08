# Garmin Custom Map Converter

Made by Maiko for the love of the game

A Python tool to convert GeoTIFF files into Garmin-compatible KMZ or IMG files for use with Garmin devices (such as Rino series).

## Features

### Option 1: KMZ Conversion
- Converts GeoTIFF (.tif) files to JPEG images
- Automatically resizes images to meet Garmin requirements (max 1024x1024 pixels, max 3MB)
- Automatically detects and converts coordinate systems to WGS84 (required for KML)
- Extracts georeferencing information and creates KML files
- Transforms KML format to Garmin-compatible structure
- Packages everything into KMZ files ready for Garmin import

### Option 2: IMG Conversion
- Converts GeoTIFF (.tif) files to Garmin IMG format
- Uses gdal_polygonize to convert raster to vector (OSM format)
- Uses mkgmap to create native Garmin IMG files
- Produces IMG files that can be directly loaded onto Garmin devices

## Requirements

### For Both Conversion Types
- Python 3.6 or higher
- GDAL library (includes coordinate transformation and image resizing capabilities - see installation instructions below)

### Additional Requirements for IMG Conversion
- Java Runtime Environment (JRE) - required to run mkgmap.jar
- mkgmap.jar - download from [mkgmap.org.uk](https://www.mkgmap.org.uk/)

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

### 3. Install mkgmap (Required for IMG Conversion Only)

**Easiest method:**
1. Download the mkgmap zip file from [mkgmap.org.uk](https://www.mkgmap.org.uk/download/mkgmap.html)
2. Unzip it in the current directory (it will create a folder like `mkgmap-r4923`)
3. That's it! The script will automatically find it and all dependencies.

**Alternative methods:**
- Place `mkgmap.jar` directly in the current directory
- Use the `--mkgmap-path` argument to specify a custom location

## Usage

### Interactive Mode (Recommended)

When you run the script without specifying a conversion type, you'll be prompted to choose:

```bash
python garmin-map-converter.py input.tif
```

You'll see:
```
Select conversion type: 1) KMZ  2) IMG [1]:
```

- Press Enter or type `1` for KMZ conversion
- Type `2` for IMG conversion

### KMZ Conversion

Convert a GeoTIFF file to KMZ:
```bash
python garmin-map-converter.py input.tif --type kmz
```

Or use the interactive prompt (default):
```bash
python garmin-map-converter.py input.tif
```

This will create `input.kmz` in the same directory as the input file.

### IMG Conversion

Convert a GeoTIFF file to IMG:
```bash
python garmin-map-converter.py input.tif --type img
```

Or use the interactive prompt:
```bash
python garmin-map-converter.py input.tif
# Then select option 2
```

This will create `input.img` in the same directory as the input file.

### Specify Output File

For KMZ conversion:
```bash
python garmin-map-converter.py input.tif --type kmz -o output.kmz
```

For IMG conversion:
```bash
python garmin-map-converter.py input.tif --type img -o output.img
```

### Specify mkgmap.jar Location

If mkgmap.jar is not in the current directory or PATH:
```bash
python garmin-map-converter.py input.tif --type img --mkgmap-path /path/to/mkgmap.jar
```

### Command-Line Options

```
positional arguments:
  input                 Path to input GeoTIFF file

optional arguments:
  -h, --help            Show help message and exit
  -o OUTPUT, --output OUTPUT
                        Path to output file (default: input filename with .kmz or .img extension)
  --type {kmz,img}      Conversion type: kmz or img (if not provided, will prompt interactively)
  --mkgmap-path MKGMAP_PATH
                        Path to mkgmap.jar (required for IMG conversion, will auto-detect if not provided)
  --temp-dir TEMP_DIR   Temporary directory for intermediate files (default: system temp)
```

## How It Works

### KMZ Conversion Process

1. **TIF to JPG Conversion**: Uses GDAL to convert the GeoTIFF to a JPEG image
2. **Image Resizing**: Automatically resizes images to meet Garmin requirements using GDAL:
   - Maximum dimensions: 1024x1024 pixels (maintains aspect ratio)
   - Maximum file size: 3MB (adjusts JPEG quality and size as needed)
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

### IMG Conversion Process

1. **TIF to Shapefile Conversion**: Uses `gdal_polygonize.py` to convert the GeoTIFF raster to ESRI Shapefile format
   - This converts the raster image into polygons/vectors
2. **Shapefile to OSM Conversion**: Uses GDAL Python bindings to read the Shapefile and generate OSM XML format
   - Automatically transforms coordinates to WGS84 (lat/lon) if needed
   - Creates proper OSM XML structure with nodes and ways
   - OSM format is required as input for mkgmap
3. **OSM to IMG Conversion**: Uses `mkgmap.jar` (Java tool) to convert the OSM file to Garmin IMG format
   - Generates a unique mapname (8-digit numeric ID) based on the input filename
   - Creates a native Garmin IMG file that can be directly loaded onto Garmin devices

## Output Format

### KMZ Output

The generated KMZ file contains:
- A KML file named `doc.kml` with Garmin-compatible format (coordinates in WGS84)
- A JPEG image file referenced by the KML

The KMZ file can be directly imported into Garmin BaseCamp or transferred to compatible Garmin devices.

### IMG Output

The generated IMG file is a native Garmin map file that can be:
- Directly copied to a Garmin device's memory card
- Loaded into Garmin BaseCamp
- Used with Garmin MapSource or other Garmin mapping software

The IMG file uses a unique mapname (8-digit numeric ID) derived from the input filename to ensure compatibility with Garmin devices.

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

### IMG Conversion Errors

#### mkgmap.jar Not Found
If you get an error that mkgmap.jar is not found:
- Download mkgmap.jar from [mkgmap.org.uk](https://www.mkgmap.org.uk/)
- Place it in the current directory, or
- Use `--mkgmap-path` to specify its location

#### Java Not Found
If you get an error that Java is not found:
- Install Java Runtime Environment (JRE)
- Verify installation: `java -version`
- Make sure Java is in your PATH

#### gdal_polygonize.py Not Found
If you get an error that gdal_polygonize.py is not found:
- Make sure GDAL is properly installed (gdal_polygonize.py comes with GDAL)
- Verify GDAL installation: `gdalinfo --version`
- Verify tool is available: `gdal_polygonize.py --version`
- On some systems, gdal_polygonize.py may be in a different location - check your GDAL installation
- The Shapefile to OSM conversion uses GDAL Python bindings (ogr module), which should be available if GDAL Python package is installed

#### mkgmap Conversion Fails
If mkgmap fails during conversion:
- Check that your GeoTIFF has valid georeferencing information
- Verify the OSM file was created successfully (check temp directory if using `--temp-dir`)
- Some complex GeoTIFF files may need preprocessing before conversion

## License

This project is open source and available for use.

