#!/usr/bin/env python3
"""
Garmin Custom Map Converter
Converts GeoTIFF files to Garmin-compatible KMZ files.
"""

import argparse
import os
import shutil
import tempfile
import zipfile
from osgeo import gdal
from osgeo import osr
from xml.etree import ElementTree as ET
try:
    from PIL import Image
    PIL_AVAILABLE = True
except ImportError:
    PIL_AVAILABLE = False


def convert_coordinates_to_wgs84(minx, miny, maxx, maxy, projection_wkt=None):
    """
    Convert coordinates from source CRS to WGS84 (EPSG:4326) using GDAL's osr.
    
    Args:
        minx, miny, maxx, maxy: Bounding box coordinates in source CRS
        projection_wkt: WKT projection string from GDAL
        
    Returns:
        tuple: (west, south, east, north) in WGS84 decimal degrees
    """
    # Check if coordinates look like they're already in WGS84
    if -180 <= minx <= 180 and -180 <= maxx <= 180 and -90 <= miny <= 90 and -90 <= maxy <= 90:
        # Already in WGS84, return as-is
        return minx, miny, maxx, maxy
    
    if not projection_wkt or not projection_wkt.strip():
        print("Warning: No projection information found. Assuming coordinates are already in WGS84.")
        return minx, miny, maxx, maxy
    
    try:
        # Create source spatial reference from WKT
        source_srs = osr.SpatialReference()
        source_srs.ImportFromWkt(projection_wkt)
        
        # Create target spatial reference (WGS84)
        target_srs = osr.SpatialReference()
        target_srs.ImportFromEPSG(4326)  # WGS84
        
        # Check if already WGS84
        if source_srs.IsSame(target_srs):
            return minx, miny, maxx, maxy
        
        # Create coordinate transformation
        transform = osr.CoordinateTransformation(source_srs, target_srs)
        
        # Transform all four corners to ensure we get the correct bounding box
        # This handles cases where the bounding box doesn't align with coordinate axes
        corners = [
            (minx, miny),  # Southwest
            (maxx, miny),  # Southeast
            (maxx, maxy),  # Northeast
            (minx, maxy)   # Northwest
        ]
        
        transformed_corners = []
        for x, y in corners:
            # TransformPoint returns (latitude, longitude, z) for WGS84, so we need to swap
            lat, lon, _ = transform.TransformPoint(x, y)
            transformed_corners.append((lon, lat))
        
        # Extract all longitudes and latitudes
        lons = [corner[0] for corner in transformed_corners]
        lats = [corner[1] for corner in transformed_corners]
        
        # Get bounding box
        west = min(lons)
        east = max(lons)
        south = min(lats)
        north = max(lats)
        
        return west, south, east, north
        
    except Exception as e:
        print(f"Warning: Could not transform coordinates to WGS84: {str(e)}")
        print("Returning coordinates as-is. They may not be in WGS84.")
        # Try to get EPSG code for better error message
        try:
            source_srs = osr.SpatialReference()
            source_srs.ImportFromWkt(projection_wkt)
            epsg_code = source_srs.GetAuthorityCode(None)
            if epsg_code:
                print(f"Source CRS appears to be EPSG:{epsg_code}")
        except:
            pass
        return minx, miny, maxx, maxy


def resize_image_for_garmin(image_path, max_size=1024, max_file_size_mb=3):
    """
    Resize and optimize image to meet Garmin requirements:
    - Max dimensions: 1024x1024 pixels
    - Max file size: 3MB
    
    Args:
        image_path: Path to input JPG image
        max_size: Maximum width or height in pixels (default: 1024)
        max_file_size_mb: Maximum file size in MB (default: 3)
    
    Returns:
        Path to the resized image (may be the same file if already compliant)
    """
    if not PIL_AVAILABLE:
        print("Warning: PIL/Pillow not available. Cannot resize image.")
        print("Install with: pip install Pillow")
        return image_path
    
    max_file_size_bytes = max_file_size_mb * 1024 * 1024
    temp_path = None
    
    try:
        # Open the image
        img = Image.open(image_path)
        original_size = img.size
        original_file_size = os.path.getsize(image_path)
        
        # Check if resizing is needed
        needs_resize = False
        if original_size[0] > max_size or original_size[1] > max_size:
            needs_resize = True
            print(f"Image size {original_size[0]}x{original_size[1]} exceeds {max_size}x{max_size}, resizing...")
        
        # Resize if needed (maintain aspect ratio)
        if needs_resize:
            # Calculate new size maintaining aspect ratio
            ratio = min(max_size / original_size[0], max_size / original_size[1])
            new_size = (int(original_size[0] * ratio), int(original_size[1] * ratio))
            img = img.resize(new_size, Image.Resampling.LANCZOS)
            print(f"Resized to {new_size[0]}x{new_size[1]}")
        
        # Save with quality adjustment to meet file size requirement
        quality = 95
        temp_path = image_path + ".tmp"
        
        # Try different quality levels until file size is acceptable
        for q in range(95, 40, -5):
            img.save(temp_path, "JPEG", quality=q, optimize=True)
            file_size = os.path.getsize(temp_path)
            
            if file_size <= max_file_size_bytes:
                break
            quality = q
        
        # If still too large, resize further
        if os.path.getsize(temp_path) > max_file_size_bytes:
            print(f"File size still too large ({os.path.getsize(temp_path) / 1024 / 1024:.2f}MB), resizing further...")
            while os.path.getsize(temp_path) > max_file_size_bytes and img.size[0] > 256 and img.size[1] > 256:
                # Reduce size by 10%
                new_size = (int(img.size[0] * 0.9), int(img.size[1] * 0.9))
                img = img.resize(new_size, Image.Resampling.LANCZOS)
                img.save(temp_path, "JPEG", quality=85, optimize=True)
                print(f"Resized to {new_size[0]}x{new_size[1]}, size: {os.path.getsize(temp_path) / 1024 / 1024:.2f}MB")
        
        # Replace original with resized version
        if needs_resize or os.path.getsize(temp_path) != original_file_size:
            shutil.move(temp_path, image_path)
            final_size = os.path.getsize(image_path)
            print(f"Final image: {img.size[0]}x{img.size[1]}, {final_size / 1024 / 1024:.2f}MB")
        else:
            if temp_path and os.path.exists(temp_path):
                os.remove(temp_path)
        
        return image_path
        
    except Exception as e:
        print(f"Warning: Could not resize image: {str(e)}")
        if temp_path and os.path.exists(temp_path):
            os.remove(temp_path)
        return image_path


def tif_to_kml_img(input_tif, output_dir):
    """
    Convert GeoTIFF to JPG image and generate initial KML file.
    
    Args:
        input_tif: Path to input GeoTIFF file
        output_dir: Directory to save output files
        
    Returns:
        tuple: (output_jpg_path, output_kml_path, bounds_dict)
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    base_name = os.path.splitext(os.path.basename(input_tif))[0]
    output_img = os.path.join(output_dir, base_name + ".jpg")
    output_kml = os.path.join(output_dir, base_name + ".kml")
    
    # Open the GeoTIFF
    dataset = gdal.Open(input_tif)
    if not dataset:
        raise ValueError("Could not open the TIF file")
    
    # Export to JPEG (convert raster to image)
    gdal.Translate(output_img, dataset, format="JPEG")
    
    # Resize and optimize image to meet Garmin requirements (1024x1024 max, 3MB max)
    output_img = resize_image_for_garmin(output_img)
    
    # Get georeferencing info
    gt = dataset.GetGeoTransform()
    width = dataset.RasterXSize
    height = dataset.RasterYSize
    
    # Calculate bounds in source coordinate system
    minx = gt[0]
    maxy = gt[3]
    maxx = minx + gt[1] * width
    miny = maxy + gt[5] * height
    
    # Get projection WKT string
    projection_wkt = dataset.GetProjection()
    
    # Convert coordinates to WGS84 if needed
    west, south, east, north = convert_coordinates_to_wgs84(minx, miny, maxx, maxy, projection_wkt)
    
    # Store bounds for later use
    bounds = {
        'north': north,
        'south': south,
        'east': east,
        'west': west
    }
    
    # Create the initial KML file referencing the image (using WGS84 coordinates)
    kml_content = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <GroundOverlay>
    <name>{base_name}</name>
    <Icon>
      <href>{os.path.basename(output_img)}</href>
    </Icon>
    <LatLonBox>
      <north>{north}</north>
      <south>{south}</south>
      <east>{east}</east>
      <west>{west}</west>
    </LatLonBox>
  </GroundOverlay>
</kml>
"""
    
    with open(output_kml, "w") as f:
        f.write(kml_content)
    
    dataset = None  # Close dataset
    
    return output_img, output_kml, bounds


def fix_kml_format(input_kml, output_kml, image_filename, source_tif_path):
    """
    Transform KML from generated format to Garmin-compatible format.
    
    Args:
        input_kml: Path to input KML file
        output_kml: Path to output fixed KML file
        image_filename: Name of the JPG image file
        source_tif_path: Path to source TIF file for description
    """
    # Parse the input KML
    tree = ET.parse(input_kml)
    root = tree.getroot()
    
    # Extract GroundOverlay data
    ground_overlay = root.find('.//{http://www.opengis.net/kml/2.2}GroundOverlay')
    if ground_overlay is None:
        raise ValueError("Could not find GroundOverlay in KML file")
    
    name_elem = ground_overlay.find('{http://www.opengis.net/kml/2.2}name')
    name = name_elem.text if name_elem is not None else os.path.splitext(os.path.basename(source_tif_path))[0]
    
    latlonbox = ground_overlay.find('{http://www.opengis.net/kml/2.2}LatLonBox')
    north = latlonbox.find('{http://www.opengis.net/kml/2.2}north').text
    south = latlonbox.find('{http://www.opengis.net/kml/2.2}south').text
    east = latlonbox.find('{http://www.opengis.net/kml/2.2}east').text
    west = latlonbox.find('{http://www.opengis.net/kml/2.2}west').text
    
    # Create the Garmin-compatible KML structure
    kml_content = f"""<?xml version="1.0" encoding="iso-8859-1"?>
<kml xmlns="http://earth.google.com/kml/2.1">
<Document>
  <name>{name}</name>
  <GroundOverlay>
    <name>{os.path.basename(source_tif_path)}</name>
    <description>Source Image: {os.path.abspath(source_tif_path)}</description>
    <Icon>
      <href>{image_filename}</href>
    </Icon>
    <LatLonBox>
      <north>{north}</north>
      <south>{south}</south>
      <east>{east}</east>
      <west>{west}</west>
    </LatLonBox>
  </GroundOverlay>
</Document>
</kml>
"""
    
    with open(output_kml, "w", encoding="iso-8859-1") as f:
        f.write(kml_content)


def create_kmz(kml_path, image_path, output_kmz):
    """
    Package KML and JPG image into a KMZ file.
    
    Args:
        kml_path: Path to KML file
        image_path: Path to JPG image file
        output_kmz: Path to output KMZ file
    """
    with zipfile.ZipFile(output_kmz, 'w', zipfile.ZIP_DEFLATED) as kmz:
        # Add KML file (must be named doc.kml for Garmin compatibility)
        kmz.write(kml_path, "doc.kml")
        
        # Add image file
        image_basename = os.path.basename(image_path)
        kmz.write(image_path, image_basename)


def main():
    parser = argparse.ArgumentParser(
        description='Convert GeoTIFF files to Garmin-compatible KMZ files'
    )
    parser.add_argument(
        'input',
        help='Path to input GeoTIFF file'
    )
    parser.add_argument(
        '-o', '--output',
        help='Path to output KMZ file (default: input filename with .kmz extension)',
        default=None
    )
    parser.add_argument(
        '--temp-dir',
        help='Temporary directory for intermediate files (default: system temp)',
        default=None
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' does not exist.")
        return 1
    
    # Determine output path
    if args.output is None:
        base_name = os.path.splitext(os.path.basename(args.input))[0]
        output_dir = os.path.dirname(args.input) or '.'
        args.output = os.path.join(output_dir, base_name + ".kmz")
    
    # Create temporary directory
    temp_dir = args.temp_dir
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp(prefix="garmin_converter_")
        cleanup_temp = True
    else:
        os.makedirs(temp_dir, exist_ok=True)
        cleanup_temp = False
    
    try:
        # Step 1: Convert TIF to JPG + initial KML
        print(f"Converting {args.input} to JPG and generating KML...")
        jpg_path, initial_kml_path, bounds = tif_to_kml_img(args.input, temp_dir)
        
        # Step 2: Fix KML format
        print("Fixing KML format for Garmin compatibility...")
        doc_kml_path = os.path.join(temp_dir, "doc.kml")
        image_filename = os.path.basename(jpg_path)
        fix_kml_format(initial_kml_path, doc_kml_path, image_filename, args.input)
        
        # Step 3: Package into KMZ
        print(f"Creating KMZ file: {args.output}")
        create_kmz(doc_kml_path, jpg_path, args.output)
        
        print(f"Successfully created {args.output}")
        return 0
        
    except Exception as e:
        print(f"Error: {str(e)}")
        return 1
        
    finally:
        # Clean up temporary directory if we created it
        if cleanup_temp and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)


if __name__ == "__main__":
    exit(main())

