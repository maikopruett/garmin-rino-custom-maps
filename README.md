# Garmin Custom Map Converter

A Python tool to convert GeoTIFF files into Garmin-compatible KMZ files for use with Garmin devices (such as Rino series).

## Features

- Converts GeoTIFF (.tif) files to JPEG images
- Extracts georeferencing information and creates KML files
- Transforms KML format to Garmin-compatible structure
- Packages everything into KMZ files ready for Garmin import

## Requirements

- Python 3.6 or higher
- GDAL library (see installation instructions below)

## Installation

### 1. Install GDAL

GDAL installation varies by operating system:

#### macOS
```bash
# Using Homebrew (recommended)
brew install gdal

# Then install Python bindings
pip install gdal
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
2. **Georeferencing Extraction**: Extracts geographic bounds (north, south, east, west) from the GeoTIFF metadata
3. **KML Generation**: Creates an initial KML file with the georeferencing information
4. **KML Format Fix**: Transforms the KML to Garmin-compatible format:
   - Changes encoding to iso-8859-1
   - Updates namespace to kml 2.1
   - Wraps in Document element
   - Ensures proper structure for Garmin devices
5. **KMZ Packaging**: Packages the KML and JPG image into a KMZ (ZIP) file

## Output Format

The generated KMZ file contains:
- A KML file with Garmin-compatible format
- A JPEG image file referenced by the KML

The KMZ file can be directly imported into Garmin BaseCamp or transferred to compatible Garmin devices.

## Troubleshooting

### GDAL Import Error
If you get `ImportError: No module named 'osgeo'`, make sure GDAL is properly installed:
- Verify GDAL installation: `gdalinfo --version`
- Reinstall Python bindings: `pip install --upgrade gdal`

### File Not Found Error
Ensure the input TIF file path is correct and the file exists.

### Permission Errors
Make sure you have write permissions in the output directory.

## License

This project is open source and available for use.

