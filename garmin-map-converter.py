#!/usr/bin/env python3
"""
Garmin Custom Map Converter
Converts GeoTIFF files to Garmin-compatible KMZ or IMG files.
"""

import argparse
import os
import shutil
import subprocess
import tempfile
import zipfile
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
from xml.etree import ElementTree as ET


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
    Resize and optimize image to meet Garmin requirements using GDAL:
    - Max dimensions: 1024x1024 pixels
    - Max file size: 3MB
    
    Args:
        image_path: Path to input JPG image
        max_size: Maximum width or height in pixels (default: 1024)
        max_file_size_mb: Maximum file size in MB (default: 3)
    
    Returns:
        Path to the resized image (may be the same file if already compliant)
    """
    max_file_size_bytes = max_file_size_mb * 1024 * 1024
    temp_path = None
    
    try:
        # Open the image to get dimensions
        dataset = gdal.Open(image_path)
        if not dataset:
            print("Warning: Could not open image for resizing")
            return image_path
        
        width = dataset.RasterXSize
        height = dataset.RasterYSize
        original_file_size = os.path.getsize(image_path)
        
        # Check if resizing is needed
        needs_resize = False
        if width > max_size or height > max_size:
            needs_resize = True
            print(f"Image size {width}x{height} exceeds {max_size}x{max_size}, resizing...")
        
        # Calculate new size maintaining aspect ratio
        if needs_resize:
            ratio = min(max_size / width, max_size / height)
            new_width = int(width * ratio)
            new_height = int(height * ratio)
        else:
            new_width = width
            new_height = height
        
        # If file size is acceptable and no resize needed, return early
        if not needs_resize and original_file_size <= max_file_size_bytes:
            dataset = None
            return image_path
        
        # Create temp path for resized image
        temp_path = image_path + ".tmp"
        
        # Try different JPEG quality levels and sizes to meet file size requirement
        quality_levels = [95, 85, 75, 65, 55, 45]
        current_width = new_width
        current_height = new_height
        
        for quality in quality_levels:
            # Set JPEG quality through creation options
            translate_options = gdal.TranslateOptions(
                format='JPEG',
                width=current_width,
                height=current_height,
                creationOptions=[f'JPEG_QUALITY={quality}']
            )
            
            # Resize the image
            gdal.Translate(temp_path, dataset, options=translate_options)
            file_size = os.path.getsize(temp_path)
            
            if file_size <= max_file_size_bytes:
                break
            
            # If still too large, reduce size further
            if file_size > max_file_size_bytes and current_width > 256 and current_height > 256:
                current_width = int(current_width * 0.9)
                current_height = int(current_height * 0.9)
                print(f"File size {file_size / 1024 / 1024:.2f}MB too large, reducing to {current_width}x{current_height}...")
        
        # Final check - if still too large, keep reducing size
        while os.path.getsize(temp_path) > max_file_size_bytes and current_width > 256 and current_height > 256:
            current_width = int(current_width * 0.9)
            current_height = int(current_height * 0.9)
            translate_options = gdal.TranslateOptions(
                format='JPEG',
                width=current_width,
                height=current_height,
                creationOptions=['JPEG_QUALITY=75']
            )
            gdal.Translate(temp_path, dataset, options=translate_options)
            print(f"Resized to {current_width}x{current_height}, size: {os.path.getsize(temp_path) / 1024 / 1024:.2f}MB")
        
        # Replace original with resized version
        if needs_resize or os.path.getsize(temp_path) != original_file_size:
            shutil.move(temp_path, image_path)
            final_size = os.path.getsize(image_path)
            print(f"Final image: {current_width}x{current_height}, {final_size / 1024 / 1024:.2f}MB")
        else:
            if os.path.exists(temp_path):
                os.remove(temp_path)
        
        dataset = None
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


def find_mkgmap_jar(mkgmap_path=None):
    """
    Find mkgmap.jar file.
    
    Args:
        mkgmap_path: Optional user-specified path to mkgmap.jar
        
    Returns:
        Path to mkgmap.jar if found
        
    Raises:
        FileNotFoundError: If mkgmap.jar cannot be found
    """
    # If user provided a path, use it
    if mkgmap_path:
        if os.path.exists(mkgmap_path):
            return mkgmap_path
        else:
            raise FileNotFoundError(f"mkgmap.jar not found at specified path: {mkgmap_path}")
    
    current_dir = os.getcwd()
    
    # Check current directory for mkgmap.jar
    current_dir_jar = os.path.join(current_dir, "mkgmap.jar")
    if os.path.exists(current_dir_jar):
        return current_dir_jar
    
    # Check for mkgmap folders in current directory (common when unzipped)
    # Look for folders starting with "mkgmap" (e.g., mkgmap-r4923, mkgmap-0.0.0, etc.)
    if os.path.exists(current_dir):
        for item in os.listdir(current_dir):
            item_path = os.path.join(current_dir, item)
            # Check if it's a directory starting with "mkgmap"
            if os.path.isdir(item_path) and item.lower().startswith("mkgmap"):
                # Look for mkgmap.jar inside this directory
                jar_in_folder = os.path.join(item_path, "mkgmap.jar")
                if os.path.exists(jar_in_folder):
                    return jar_in_folder
                # Also check for mkgmap.jar in subdirectories (some distributions have it nested)
                for subitem in os.listdir(item_path):
                    subitem_path = os.path.join(item_path, subitem)
                    if os.path.isdir(subitem_path):
                        nested_jar = os.path.join(subitem_path, "mkgmap.jar")
                        if os.path.exists(nested_jar):
                            return nested_jar
    
    # Check PATH
    path_jar = shutil.which("mkgmap.jar")
    if path_jar and os.path.exists(path_jar):
        return path_jar
    
    # Check if mkgmap.jar is in PATH as a command (unlikely but possible)
    mkgmap_cmd = shutil.which("mkgmap")
    if mkgmap_cmd:
        # If mkgmap command exists, check if it's a jar wrapper
        # For now, we'll just raise an error and let user specify
        pass
    
    raise FileNotFoundError(
        "mkgmap.jar not found. Please:\n"
        "  1. Download mkgmap zip from https://www.mkgmap.org.uk/\n"
        "  2. Unzip it in the current directory (it will create a folder like 'mkgmap-rXXXX')\n"
        "  3. Or place mkgmap.jar directly in the current directory\n"
        "  4. Or use --mkgmap-path to specify the location"
    )


def generate_mapname(input_filename):
    """
    Generate a unique mapname (numeric ID) from input filename.
    
    Args:
        input_filename: Path to input file
        
    Returns:
        String of 8 digits derived from filename hash
    """
    import hashlib
    base_name = os.path.splitext(os.path.basename(input_filename))[0]
    # Create a hash from the filename and take first 8 digits
    hash_obj = hashlib.md5(base_name.encode())
    hash_hex = hash_obj.hexdigest()
    # Convert hex to numeric string (take first 8 characters, convert to int, then back to string)
    numeric_id = str(int(hash_hex[:8], 16))[:8]
    # Pad with zeros if needed to ensure 8 digits
    return numeric_id.zfill(8)


def shapefile_to_osm(shapefile_path, osm_path):
    """
    Convert Shapefile to OSM XML format using GDAL Python bindings.
    
    Args:
        shapefile_path: Path to input Shapefile (.shp)
        osm_path: Path to output OSM XML file
    """
    # Open the shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    datasource = driver.Open(shapefile_path, 0)
    if datasource is None:
        raise RuntimeError(f"Could not open shapefile: {shapefile_path}")
    
    layer = datasource.GetLayer()
    
    # Get spatial reference and set up coordinate transformation if needed
    source_srs = layer.GetSpatialRef()
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(4326)  # WGS84
    
    transform = None
    if source_srs and not source_srs.IsSame(target_srs):
        transform = osr.CoordinateTransformation(source_srs, target_srs)
    
    # Create OSM XML structure
    osm_root = ET.Element("osm", version="0.6", generator="garmin-map-converter")
    
    node_id = 1
    way_id = 1
    
    def transform_point(x, y):
        """Transform a point to WGS84 if needed."""
        if transform:
            try:
                # TransformPoint returns (lat, lon, z) but we need (lon, lat)
                lat, lon, _ = transform.TransformPoint(x, y)
                return lon, lat
            except:
                # If transformation fails, return as-is
                return x, y
        return x, y
    
    # Process each feature (polygon) in the shapefile
    for feature in layer:
        geometry = feature.GetGeometryRef()
        if geometry is None:
            continue
        
        # Get geometry type
        geom_type = geometry.GetGeometryType()
        
        if geom_type == ogr.wkbPolygon or geom_type == ogr.wkbMultiPolygon:
            # Handle polygons
            if geom_type == ogr.wkbPolygon:
                polygons = [geometry]
            else:  # MultiPolygon
                polygons = []
                for i in range(geometry.GetGeometryCount()):
                    polygons.append(geometry.GetGeometryRef(i))
            
            for polygon in polygons:
                # Get exterior ring
                exterior_ring = polygon.GetGeometryRef(0)
                if exterior_ring is None:
                    continue
                
                # Create nodes for each point in the ring
                node_ids = []
                point_count = exterior_ring.GetPointCount()
                
                for i in range(point_count):
                    point = exterior_ring.GetPoint(i)
                    lon, lat = transform_point(point[0], point[1])
                    
                    # Create node element
                    node_elem = ET.SubElement(osm_root, "node", 
                                             id=str(node_id),
                                             lat=str(lat),
                                             lon=str(lon))
                    node_ids.append(node_id)
                    node_id += 1
                
                # Create way element connecting the nodes
                if len(node_ids) > 0:
                    way_elem = ET.SubElement(osm_root, "way", id=str(way_id))
                    
                    # Add node references
                    for nid in node_ids:
                        ET.SubElement(way_elem, "nd", ref=str(nid))
                    
                    # Close the way by adding first node again if it's a polygon
                    if len(node_ids) > 2:
                        ET.SubElement(way_elem, "nd", ref=str(node_ids[0]))
                    
                    # Add tag for area
                    ET.SubElement(way_elem, "tag", k="area", v="yes")
                    
                    way_id += 1
        
        elif geom_type == ogr.wkbLineString or geom_type == ogr.wkbMultiLineString:
            # Handle linestrings
            if geom_type == ogr.wkbLineString:
                linestrings = [geometry]
            else:  # MultiLineString
                linestrings = []
                for i in range(geometry.GetGeometryCount()):
                    linestrings.append(geometry.GetGeometryRef(i))
            
            for linestring in linestrings:
                node_ids = []
                point_count = linestring.GetPointCount()
                
                for i in range(point_count):
                    point = linestring.GetPoint(i)
                    lon, lat = transform_point(point[0], point[1])
                    
                    node_elem = ET.SubElement(osm_root, "node",
                                             id=str(node_id),
                                             lat=str(lat),
                                             lon=str(lon))
                    node_ids.append(node_id)
                    node_id += 1
                
                if len(node_ids) > 1:
                    way_elem = ET.SubElement(osm_root, "way", id=str(way_id))
                    for nid in node_ids:
                        ET.SubElement(way_elem, "nd", ref=str(nid))
                    way_id += 1
    
    datasource = None
    
    # Write OSM XML to file
    tree = ET.ElementTree(osm_root)
    # ET.indent is available in Python 3.9+, use it if available
    if hasattr(ET, 'indent'):
        ET.indent(tree, space="  ")
    tree.write(osm_path, encoding="utf-8", xml_declaration=True)


def tif_to_img(input_tif, output_path, mkgmap_path=None):
    """
    Convert GeoTIFF to Garmin IMG format using gdal_polygonize and mkgmap.
    
    Args:
        input_tif: Path to input GeoTIFF file
        output_path: Path to output IMG file
        mkgmap_path: Optional path to mkgmap.jar (will auto-detect if not provided)
        
    Returns:
        Path to the created IMG file
    """
    # Find mkgmap.jar
    mkgmap_jar = find_mkgmap_jar(mkgmap_path)
    
    # Create temporary directory for intermediate files
    temp_dir = tempfile.mkdtemp(prefix="garmin_img_")
    
    try:
        base_name = os.path.splitext(os.path.basename(input_tif))[0]
        shapefile_path = os.path.join(temp_dir, base_name + ".shp")
        osm_path = os.path.join(temp_dir, base_name + ".osm")
        mkgmap_output_dir = os.path.join(temp_dir, "mkgmap_output")
        os.makedirs(mkgmap_output_dir, exist_ok=True)
        
        # Step 1: Convert GeoTIFF to Shapefile using gdal_polygonize
        print(f"Converting {input_tif} to Shapefile format...")
        try:
            result = subprocess.run(
                ["gdal_polygonize.py", input_tif, "-f", "ESRI Shapefile", shapefile_path],
                check=True,
                capture_output=True,
                text=True
            )
        except FileNotFoundError:
            raise FileNotFoundError(
                "gdal_polygonize.py not found. Make sure GDAL is properly installed and "
                "gdal_polygonize.py is in your PATH."
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"gdal_polygonize.py failed: {e.stderr}")
        
        # Check if shapefile was created (it creates multiple files, so check for .shp)
        if not os.path.exists(shapefile_path):
            raise RuntimeError("gdal_polygonize.py did not create Shapefile")
        
        # Step 2: Convert Shapefile to OSM format using Python/GDAL
        print(f"Converting Shapefile to OSM format...")
        try:
            shapefile_to_osm(shapefile_path, osm_path)
        except Exception as e:
            raise RuntimeError(f"Failed to convert Shapefile to OSM: {str(e)}")
        
        if not os.path.exists(osm_path):
            raise RuntimeError("Failed to create OSM file")
        
        # Step 3: Convert OSM to IMG using mkgmap
        print(f"Converting OSM to IMG format using mkgmap...")
        mapname = generate_mapname(input_tif)
        
        # Check if Java is available
        java_cmd = shutil.which("java")
        if not java_cmd:
            raise FileNotFoundError(
                "Java runtime not found. Please install Java (JRE) to use IMG conversion."
            )
        
        # Build classpath for mkgmap - check for additional JAR dependencies
        mkgmap_dir = os.path.dirname(os.path.abspath(mkgmap_jar))
        mkgmap_jar_name = os.path.basename(mkgmap_jar)
        
        # Collect all JAR files in the mkgmap directory
        jar_files = [os.path.abspath(mkgmap_jar)]
        jar_paths_set = {os.path.abspath(mkgmap_jar)}
        
        if os.path.exists(mkgmap_dir):
            # Check mkgmap directory for JAR files
            for file in os.listdir(mkgmap_dir):
                if file.endswith('.jar') and file != mkgmap_jar_name:
                    jar_path = os.path.abspath(os.path.join(mkgmap_dir, file))
                    if os.path.isfile(jar_path) and jar_path not in jar_paths_set:
                        jar_files.append(jar_path)
                        jar_paths_set.add(jar_path)
            
            # Check for lib subdirectory (common pattern)
            lib_dir = os.path.join(mkgmap_dir, "lib")
            if os.path.exists(lib_dir) and os.path.isdir(lib_dir):
                for file in os.listdir(lib_dir):
                    if file.endswith('.jar'):
                        jar_path = os.path.abspath(os.path.join(lib_dir, file))
                        if jar_path not in jar_paths_set:
                            jar_files.append(jar_path)
                            jar_paths_set.add(jar_path)
        
        # Also check current directory for common dependency names
        current_dir = os.getcwd()
        common_deps = ['osmosis-core.jar', 'osmosis-pbf.jar', 'osmosis-xml.jar', 
                       'osmosis-osm-binary.jar', 'osmosis-pbfmigrate.jar']
        for dep in common_deps:
            dep_path = os.path.join(current_dir, dep)
            if os.path.exists(dep_path):
                dep_abs_path = os.path.abspath(dep_path)
                if dep_abs_path not in jar_paths_set:
                    jar_files.append(dep_abs_path)
                    jar_paths_set.add(dep_abs_path)
        
        # Check for lib directory in current directory
        current_lib_dir = os.path.join(current_dir, "lib")
        if os.path.exists(current_lib_dir) and os.path.isdir(current_lib_dir):
            for file in os.listdir(current_lib_dir):
                if file.endswith('.jar'):
                    jar_path = os.path.abspath(os.path.join(current_lib_dir, file))
                    if jar_path not in jar_paths_set:
                        jar_files.append(jar_path)
                        jar_paths_set.add(jar_path)
        
        # Build Java command
        # Use classpath method if we have multiple JARs, otherwise try -jar
        if len(jar_files) > 1:
            # Multiple JARs - use classpath
            classpath = os.pathsep.join(jar_files)
            print(f"Found {len(jar_files)} JAR files, using classpath method")
            java_args = [
                "java", "-cp", classpath,
                "uk.me.parabola.mkgmap.Main",
                f"--output-dir={mkgmap_output_dir}",
                f"--mapname={mapname}",
                osm_path
            ]
        else:
            # Single JAR - use -jar method
            print(f"Using single JAR: {mkgmap_jar}")
            java_args = [
                "java", "-jar", mkgmap_jar,
                f"--output-dir={mkgmap_output_dir}",
                f"--mapname={mapname}",
                osm_path
            ]
        
        try:
            result = subprocess.run(
                java_args,
                check=True,
                capture_output=True,
                text=True
            )
        except subprocess.CalledProcessError as e:
            error_msg = f"mkgmap failed: {e.stderr}"
            if e.stdout:
                error_msg += f"\nstdout: {e.stdout}"
            
            # Check if it's a missing dependency error
            if "NoClassDefFoundError" in error_msg or "ClassNotFoundException" in error_msg:
                error_msg += (
                    "\n\nThis error suggests mkgmap is missing required dependencies (likely Osmosis JARs).\n"
                    "To fix this:\n"
                    "  1. Download Osmosis from https://github.com/openstreetmap/osmosis/releases\n"
                    "  2. Extract the Osmosis JAR files (osmosis-core.jar, osmosis-pbf.jar, etc.)\n"
                    "  3. Place them in the same directory as mkgmap.jar, or in a 'lib' subdirectory\n"
                    "  4. Alternatively, place them in the current directory or a 'lib' subdirectory"
                )
            
            raise RuntimeError(error_msg)
        
        # Find the generated IMG file (mkgmap creates it with the mapname)
        img_filename = f"{mapname}.img"
        generated_img = os.path.join(mkgmap_output_dir, img_filename)
        
        if not os.path.exists(generated_img):
            # Sometimes mkgmap creates files with different names, search for .img files
            img_files = [f for f in os.listdir(mkgmap_output_dir) if f.endswith('.img')]
            if img_files:
                generated_img = os.path.join(mkgmap_output_dir, img_files[0])
            else:
                raise RuntimeError(f"mkgmap did not create IMG file in {mkgmap_output_dir}")
        
        # Move IMG file to final output location
        output_dir = os.path.dirname(output_path) or '.'
        os.makedirs(output_dir, exist_ok=True)
        shutil.move(generated_img, output_path)
        
        print(f"Successfully created {output_path}")
        return output_path
        
    finally:
        # Clean up temporary directory
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)


def main():
    parser = argparse.ArgumentParser(
        description='Convert GeoTIFF files to Garmin-compatible KMZ or IMG files'
    )
    parser.add_argument(
        'input',
        help='Path to input GeoTIFF file'
    )
    parser.add_argument(
        '-o', '--output',
        help='Path to output file (default: input filename with .kmz or .img extension)',
        default=None
    )
    parser.add_argument(
        '--type',
        choices=['kmz', 'img'],
        default=None,
        help='Conversion type: kmz or img (if not provided, will prompt interactively)'
    )
    parser.add_argument(
        '--mkgmap-path',
        help='Path to mkgmap.jar (required for IMG conversion, will auto-detect if not provided)',
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
    
    # Determine conversion type
    conversion_type = args.type
    if conversion_type is None:
        # Interactive prompt
        while True:
            try:
                user_input = input("Select conversion type: 1) KMZ  2) IMG [1]: ").strip()
                if user_input == '' or user_input == '1':
                    conversion_type = 'kmz'
                    break
                elif user_input == '2':
                    conversion_type = 'img'
                    break
                else:
                    print("Invalid choice. Please enter 1 or 2.")
            except (EOFError, KeyboardInterrupt):
                print("\nCancelled.")
                return 1
    
    # Determine output path based on conversion type
    if args.output is None:
        base_name = os.path.splitext(os.path.basename(args.input))[0]
        output_dir = os.path.dirname(args.input) or '.'
        extension = '.img' if conversion_type == 'img' else '.kmz'
        args.output = os.path.join(output_dir, base_name + extension)
    
    # Branch based on conversion type
    if conversion_type == 'kmz':
        # KMZ conversion (existing flow)
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
    
    elif conversion_type == 'img':
        # IMG conversion (new flow)
        try:
            tif_to_img(args.input, args.output, args.mkgmap_path)
            return 0
        except Exception as e:
            print(f"Error: {str(e)}")
            return 1


if __name__ == "__main__":
    exit(main())

