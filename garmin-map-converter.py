#!/usr/bin/env python3
"""
Garmin Custom Map Converter
Converts GeoTIFF files to Garmin-compatible KMZ, IMG, or GMAPI files.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
from xml.etree import ElementTree as ET

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_TEMP_ROOT = SCRIPT_DIR / ".garmin_temp"


def _ensure_project_temp_root():
    PROJECT_TEMP_ROOT.mkdir(parents=True, exist_ok=True)
    return PROJECT_TEMP_ROOT


def create_project_temp_dir(prefix):
    """
    Create a temporary directory inside the project folder instead of /var.
    """
    temp_root = _ensure_project_temp_root()
    return tempfile.mkdtemp(prefix=prefix, dir=str(temp_root))


def create_project_temp_file(prefix, suffix=".tmp", mode="w+", encoding="utf-8"):
    """
    Create a temporary file inside the project folder.
    """
    temp_root = _ensure_project_temp_root()
    return tempfile.NamedTemporaryFile(
        mode=mode,
        encoding=encoding,
        delete=False,
        dir=str(temp_root),
        prefix=prefix,
        suffix=suffix
    )

# Suppress GDAL warning about UseExceptions
gdal.UseExceptions()


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
    Find mkgmap.jar file in the project directory.
    
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
    
    # Check current directory for mkgmap.jar
    current_dir_jar = os.path.join(os.getcwd(), "mkgmap.jar")
    if os.path.exists(current_dir_jar):
        return current_dir_jar
    
    raise FileNotFoundError(
        "mkgmap.jar not found in current directory. Please:\n"
        "  1. Place mkgmap.jar in the project directory\n"
        "  2. Or use --mkgmap-path to specify the location"
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


def _process_feature_batch(feature_data_batch, transform_wkt):
    """
    Process a batch of features and return simplified way coordinate data.
    Designed for multiprocessing without managing global node/way IDs.
    
    Args:
        feature_data_batch: List of feature data dicts with geometry coordinates
        transform_wkt: WKT string for coordinate transformation (or None)
    
    Returns:
        List of dicts with keys:
            - 'coords': list of (lon, lat) tuples
            - 'is_area': bool indicating if geometry should be closed and tagged as area
    """
    # Set up coordinate transformation if needed
    transform = None
    if transform_wkt:
        try:
            source_srs = osr.SpatialReference()
            source_srs.ImportFromWkt(transform_wkt)
            target_srs = osr.SpatialReference()
            target_srs.ImportFromEPSG(4326)  # WGS84
            if not source_srs.IsSame(target_srs):
                transform = osr.CoordinateTransformation(source_srs, target_srs)
        except Exception:
            transform = None
    
    def transform_point(x, y):
        """Transform a point to WGS84 if needed."""
        if transform:
            try:
                lat, lon, _ = transform.TransformPoint(x, y)
                return lon, lat
            except Exception:
                return x, y
        return x, y
    
    ways_list = []
    
    for feature_data in feature_data_batch:
        geom_type = feature_data['type']
        geometries = feature_data['geometries']
        is_area_geom = geom_type in ['polygon', 'multipolygon']
        
        for geometry in geometries:
            min_points = 3 if is_area_geom else 2
            if len(geometry) < min_points:
                continue
            
            coords = []
            for x, y in geometry:
                lon, lat = transform_point(x, y)
                coords.append((lon, lat))
            
            if len(coords) >= min_points:
                ways_list.append({
                    'coords': coords,
                    'is_area': is_area_geom
                })
    
    return ways_list


def _extract_points_fast(geometry_obj):
    """
    Efficiently extract points from a geometry object.
    Uses ExportToWkb and manual parsing for better performance.
    
    Args:
        geometry_obj: OGR Geometry object
    
    Returns:
        List of (x, y) tuples
    """
    point_count = geometry_obj.GetPointCount()
    if point_count == 0:
        return []
    
    # Use GetPoint() but optimize by pre-allocating list
    points = [None] * point_count
    for i in range(point_count):
        point = geometry_obj.GetPoint(i)
        points[i] = (point[0], point[1])
    return points


def _process_single_geometry(geometry):
    """
    Process a single cloned geometry object and extract feature data.
    This function is designed to be called by parallel executors.
    
    Args:
        geometry: Cloned OGR Geometry object (independent of feature)
    
    Returns:
        Feature data dict or None if geometry is invalid
    """
    if geometry is None:
        return None
    
    geom_type = geometry.GetGeometryType()
    geometries = []
    
    try:
        if geom_type == ogr.wkbPolygon:
            # Single polygon
            exterior_ring = geometry.GetGeometryRef(0)
            if exterior_ring:
                points = _extract_points_fast(exterior_ring)
                if len(points) > 0:
                    geometries.append(points)
            if len(geometries) > 0:
                return {
                    'type': 'polygon',
                    'geometries': geometries
                }
        
        elif geom_type == ogr.wkbMultiPolygon:
            # MultiPolygon
            geom_count = geometry.GetGeometryCount()
            for i in range(geom_count):
                polygon = geometry.GetGeometryRef(i)
                if polygon:
                    exterior_ring = polygon.GetGeometryRef(0)
                    if exterior_ring:
                        points = _extract_points_fast(exterior_ring)
                        if len(points) > 0:
                            geometries.append(points)
            if len(geometries) > 0:
                return {
                    'type': 'multipolygon',
                    'geometries': geometries
                }
        
        elif geom_type == ogr.wkbLineString:
            # Single linestring
            points = _extract_points_fast(geometry)
            if len(points) > 0:
                geometries.append(points)
                return {
                    'type': 'linestring',
                    'geometries': geometries
                }
        
        elif geom_type == ogr.wkbMultiLineString:
            # MultiLineString
            geom_count = geometry.GetGeometryCount()
            for i in range(geom_count):
                linestring = geometry.GetGeometryRef(i)
                if linestring:
                    points = _extract_points_fast(linestring)
                    if len(points) > 0:
                        geometries.append(points)
            if len(geometries) > 0:
                return {
                    'type': 'multilinestring',
                    'geometries': geometries
                }
    except Exception:
        # If processing fails, return None (will be skipped)
        return None
    
    return None


def _process_geometry_batch(geometry_batch):
    """
    Process a batch of geometries in parallel.
    
    Args:
        geometry_batch: List of OGR Geometry objects or WKB byte strings
    
    Returns:
        List of feature data dicts (None values filtered out)
    """
    results = []
    for geometry in geometry_batch:
        if geometry is None:
            continue
        if isinstance(geometry, (bytes, bytearray, memoryview)):
            try:
                geometry = ogr.CreateGeometryFromWkb(bytes(geometry))
            except Exception:
                continue
        result = _process_single_geometry(geometry)
        if result is not None:
            results.append(result)
    return results


def _extract_feature_data(shapefile_path, max_workers=None):
    """
    Extract feature geometries from a shapefile as a streaming generator.
    Uses multiprocessing to parallelize geometry extraction without holding the
    entire dataset in memory at once.
    
    Args:
        shapefile_path: Path to shapefile
        max_workers: Maximum number of worker processes (None = auto-detect CPU count)
    
    Returns:
        Tuple of (feature_batch_iterator, source_srs_wkt, feature_count)
            feature_batch_iterator: generator yielding lists of feature data dicts
    """
    import multiprocessing
    
    if max_workers is None:
        max_workers = max(1, multiprocessing.cpu_count())
    
    driver = ogr.GetDriverByName("ESRI Shapefile")
    datasource = driver.Open(shapefile_path, 0)
    if datasource is None:
        raise RuntimeError(f"Could not open shapefile: {shapefile_path}")
    
    layer = datasource.GetLayer()
    feature_count = layer.GetFeatureCount()
    layer.ResetReading()
    
    source_srs = layer.GetSpatialRef()
    source_srs_wkt = source_srs.ExportToWkt() if source_srs else None
    
    batch_size = max(100, feature_count // (max_workers * 10)) if feature_count > 1000 else 50
    progress_interval = max(1, feature_count // 100) if feature_count > 1000 else 10
    
    def feature_batch_iterator():
        nonlocal datasource
        
        geometry_batch = []
        processed_features = 0
        processed_geometries = 0
        total_geometries = 0
        processed_batches = 0
        batch_index = 0
        next_batch_to_yield = 0
        
        pending_results = {}
        
        outstanding_limit = max(2, max_workers * 2)
        
        print(f"Reading features from shapefile (using {max_workers} processes for extraction)...")
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            
            def submit_current_batch():
                nonlocal geometry_batch, batch_index, total_geometries
                if not geometry_batch:
                    return
                batch = geometry_batch
                geometry_batch = []
                batch_len = len(batch)
                total_geometries += batch_len
                future = executor.submit(_process_geometry_batch, batch)
                futures[future] = (batch_index, batch_len)
                batch_index += 1
            
            def drain_one():
                nonlocal processed_geometries, processed_batches, next_batch_to_yield
                if not futures:
                    return
                for future in as_completed(list(futures.keys())):
                    batch_idx, batch_len = futures.pop(future)
                    try:
                        batch_results = future.result()
                    except Exception as e:
                        print(f"\nError processing batch: {e}")
                        raise
                    pending_results[batch_idx] = batch_results
                    processed_batches += 1
                    processed_geometries += batch_len
                    
                    if total_geometries > 0:
                        geometry_progress = (processed_geometries / total_geometries) * 100
                        batch_progress = (processed_batches / max(1, batch_index)) * 100
                        sys.stdout.write(
                            f"\rExtracting features: {processed_geometries}/{total_geometries} geometries "
                            f"({geometry_progress:.1f}%) | {processed_batches}/{max(1, batch_index)} batches "
                            f"({batch_progress:.1f}%)"
                        )
                        sys.stdout.flush()
                    
                    while next_batch_to_yield in pending_results:
                        result = pending_results.pop(next_batch_to_yield)
                        next_batch_to_yield += 1
                        yield result
                    break
            
            feature = layer.GetNextFeature()
            while feature is not None:
                geometry_ref = feature.GetGeometryRef()
                
                if geometry_ref is None:
                    feature = None
                    processed_features += 1
                    if feature_count > 0 and (processed_features % progress_interval == 0 or processed_features == feature_count):
                        progress = (processed_features / feature_count) * 100
                        sys.stdout.write(f"\rReading features: {processed_features}/{feature_count} ({progress:.1f}%)")
                        sys.stdout.flush()
                    feature = layer.GetNextFeature()
                    continue
                
                geometry_wkb = geometry_ref.ExportToWkb()
                feature = None  # Release feature immediately
                
                geometry_batch.append(bytes(geometry_wkb))
                processed_features += 1
                
                if feature_count > 0 and (processed_features % progress_interval == 0 or processed_features == feature_count):
                    progress = (processed_features / feature_count) * 100
                    sys.stdout.write(f"\rReading features: {processed_features}/{feature_count} ({progress:.1f}%)")
                    sys.stdout.flush()
                
                if len(geometry_batch) >= batch_size:
                    submit_current_batch()
                    if len(futures) >= outstanding_limit:
                        yield from drain_one()
                
                feature = layer.GetNextFeature()
            
            print()  # New line after reading progress
            
            # Submit any remaining geometries
            submit_current_batch()
            
            if batch_index > 0:
                print(f"Processing geometries using {max_workers} processes...")
            
            # Drain all remaining futures
            while futures:
                yield from drain_one()
            
            print()  # New line after extraction progress
        
        datasource = None
    
    return feature_batch_iterator(), source_srs_wkt, feature_count


def shapefile_to_osm(shapefile_path, osm_path, max_workers=None):
    """
    Convert Shapefile to OSM XML format using GDAL Python bindings.
    Streams data through multiprocessing so that memory use stays bounded even
    for extremely large datasets.
    
    Args:
        shapefile_path: Path to input Shapefile (.shp)
        osm_path: Path to output OSM XML file
        max_workers: Maximum number of worker processes (None = auto-detect CPU count)
    """
    import multiprocessing
    import tempfile
    
    # Auto-detect CPU count if not specified
    if max_workers is None:
        max_workers = max(1, multiprocessing.cpu_count() - 1)
    
    print(f"Extracting features from shapefile...")
    feature_batch_iter, source_srs_wkt, feature_count = _extract_feature_data(shapefile_path)
    
    if feature_count == 0:
        # Create empty OSM file
        osm_root = ET.Element("osm", version="0.6", generator="garmin-map-converter")
        tree = ET.ElementTree(osm_root)
        if hasattr(ET, 'indent'):
            ET.indent(tree, space="  ")
        tree.write(osm_path, encoding="utf-8", xml_declaration=True)
        return
    
    def format_coord(value):
        text = f"{value:.10f}".rstrip('0').rstrip('.')
        return text if text else "0"
    
    ways_temp = create_project_temp_file(prefix="garmin_osm_ways_", suffix=".tmp")
    ways_temp_path = ways_temp.name
    
    try:
        with open(osm_path, "w", encoding="utf-8") as osm_file:
            osm_file.write('<?xml version="1.0" encoding="utf-8"?>\n')
            osm_file.write('<osm version="0.6" generator="garmin-map-converter">\n')
            
            node_coord_map = {}
            current_node_id = 1
            current_way_id = 1
            
            def ensure_node_id(coord):
                nonlocal current_node_id
                if coord in node_coord_map:
                    return node_coord_map[coord]
                lon, lat = coord
                node_id = current_node_id
                node_coord_map[coord] = node_id
                node_line = (
                    f'  <node id="{node_id}" lat="{format_coord(lat)}" '
                    f'lon="{format_coord(lon)}" />\n'
                )
                osm_file.write(node_line)
                current_node_id += 1
                return node_id
            
            processed_features = 0
            outstanding_limit = max(2, max_workers * 2)
            pending_batches = {}
            next_batch_to_write = 0
            
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = {}
                
                def submit_batch(batch_idx, feature_batch):
                    future = executor.submit(_process_feature_batch, feature_batch, source_srs_wkt)
                    futures[future] = (batch_idx, len(feature_batch))
                
                def drain_one():
                    nonlocal processed_features, next_batch_to_write, current_way_id
                    if not futures:
                        return
                    for future in as_completed(list(futures.keys())):
                        batch_idx, batch_len = futures.pop(future)
                        try:
                            batch_ways = future.result()
                        except Exception as e:
                            print(f"\nError processing batch {batch_idx}: {e}")
                            raise
                        pending_batches[batch_idx] = batch_ways
                        processed_features += batch_len
                        progress = (processed_features / feature_count) * 100 if feature_count else 100.0
                        sys.stdout.write(
                            f"\rConverting features to OSM: {processed_features}/{feature_count} ({progress:.1f}%)"
                        )
                        sys.stdout.flush()
                        
                        while next_batch_to_write in pending_batches:
                            ways_list = pending_batches.pop(next_batch_to_write)
                            next_batch_to_write += 1
                            for way_data in ways_list:
                                coords = way_data['coords']
                                if way_data['is_area'] and len(coords) < 3:
                                    continue
                                if not way_data['is_area'] and len(coords) < 2:
                                    continue
                                
                                node_ids = []
                                for lon, lat in coords:
                                    coord_key = (lon, lat)
                                    node_ids.append(ensure_node_id(coord_key))
                                
                                if way_data['is_area'] and len(node_ids) > 2 and node_ids[0] != node_ids[-1]:
                                    node_ids.append(node_ids[0])
                                
                                if len(node_ids) < 2:
                                    continue
                                
                                ways_temp.write(f'  <way id="{current_way_id}">\n')
                                for node_id in node_ids:
                                    ways_temp.write(f'    <nd ref="{node_id}" />\n')
                                if way_data['is_area'] and len(node_ids) > 2:
                                    ways_temp.write('    <tag k="area" v="yes" />\n')
                                ways_temp.write('  </way>\n')
                                current_way_id += 1
                        break
                
                for batch_idx, feature_batch in enumerate(feature_batch_iter):
                    submit_batch(batch_idx, feature_batch)
                    if len(futures) >= outstanding_limit:
                        drain_one()
                
                while futures:
                    drain_one()
            
            print()  # New line after progress meter
            
            ways_temp.flush()
            ways_temp.seek(0)
            shutil.copyfileobj(ways_temp, osm_file)
            osm_file.write('</osm>\n')
    
    finally:
        ways_temp.close()
        if os.path.exists(ways_temp_path):
            os.remove(ways_temp_path)
    
    print(f"Successfully created OSM file at {osm_path}")


def run_mkgmap_for_gmapi(input_tif, mkgmap_output_dir, mkgmap_path=None, java_max_memory=None):
    """
    Run mkgmap to generate all map files (.img, .tdb, .typ, .mdx) needed for GMAPI.
    
    Args:
        input_tif: Path to input GeoTIFF file
        mkgmap_output_dir: Directory where mkgmap should write output files
        mkgmap_path: Optional path to mkgmap.jar (will auto-detect if not provided)
        java_max_memory: Optional string to pass to Java -Xmx (examples: '6G', '8192m')
        
    Returns:
        dict: Dictionary with paths to generated files and mapname
    """
    # Find mkgmap.jar
    mkgmap_jar = find_mkgmap_jar(mkgmap_path)
    
    # Create temporary directory for intermediate files
    temp_dir = create_project_temp_dir("garmin_gmapi_")
    
    try:
        base_name = os.path.splitext(os.path.basename(input_tif))[0]
        shapefile_path = os.path.join(temp_dir, base_name + ".shp")
        osm_path = os.path.join(temp_dir, base_name + ".osm")
        
        os.makedirs(mkgmap_output_dir, exist_ok=True)
        
        # Step 1: Convert GeoTIFF to Shapefile using gdal_polygonize
        print(f"Step 1/3: Converting GeoTIFF to Shapefile format...")
        print(f"  Running gdal_polygonize.py (this may take a while for large files)...")
        try:
            result = subprocess.run(
                ["gdal_polygonize.py", input_tif, "-f", "ESRI Shapefile", shapefile_path],
                check=True,
                capture_output=True,
                text=True
            )
            print(f"  ✓ Shapefile created successfully")
        except FileNotFoundError:
            raise FileNotFoundError(
                "gdal_polygonize.py not found. Make sure GDAL is properly installed and "
                "gdal_polygonize.py is in your PATH."
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"gdal_polygonize.py failed: {e.stderr}")
        
        # Check if shapefile was created
        if not os.path.exists(shapefile_path):
            raise RuntimeError("gdal_polygonize.py did not create Shapefile")
        
        # Step 2: Convert Shapefile to OSM format using Python/GDAL
        print(f"Step 2/3: Converting Shapefile to OSM format...")
        try:
            shapefile_to_osm(shapefile_path, osm_path)
        except Exception as e:
            raise RuntimeError(f"Failed to convert Shapefile to OSM: {str(e)}")
        
        if not os.path.exists(osm_path):
            raise RuntimeError("Failed to create OSM file")
        
        # Step 3: Convert OSM to IMG using mkgmap (with all necessary outputs)
        print(f"Step 3/3: Running mkgmap to generate map files...")
        print(f"  Processing OSM data with mkgmap (this may take a while)...")
        mapname = generate_mapname(input_tif)
        
        # Check if Java is available
        java_cmd = shutil.which("java")
        if not java_cmd:
            raise FileNotFoundError(
                "Java runtime not found. Please install Java (JRE) to use GMAPI conversion."
            )
        
        # Build classpath for mkgmap
        current_dir = os.getcwd()
        mkgmap_jar_abs = os.path.abspath(mkgmap_jar)
        
        jar_files = [mkgmap_jar_abs]
        jar_paths_set = {mkgmap_jar_abs}
        
        lib_dir = os.path.join(current_dir, "lib")
        if os.path.exists(lib_dir) and os.path.isdir(lib_dir):
            for file in os.listdir(lib_dir):
                if file.endswith('.jar'):
                    jar_path = os.path.abspath(os.path.join(lib_dir, file))
                    if jar_path not in jar_paths_set:
                        jar_files.append(jar_path)
                        jar_paths_set.add(jar_path)
        
        # Build Java command - use --tdbfile and other options to ensure all files are generated
        java_base = ["java"]
        if java_max_memory:
            java_base.append(f"-Xmx{java_max_memory}")
        java_args_jar = java_base + [
            "-jar", mkgmap_jar,
            f"--output-dir={mkgmap_output_dir}",
            f"--mapname={mapname}",
            f"--tdbfile={base_name}.tdb",
            osm_path
        ]
        
        try:
            result = subprocess.run(
                java_args_jar,
                check=True,
                capture_output=True,
                text=True
            )
        except subprocess.CalledProcessError as e:
            # Try classpath method if -jar fails
            if "NoClassDefFoundError" in e.stderr or "ClassNotFoundException" in e.stderr:
                print(f"  -jar method failed, trying with classpath (found {len(jar_files)} JAR files)")
                jar_files_clean = [j for j in jar_files if j != mkgmap_jar_abs]
                jar_files_clean.insert(0, mkgmap_jar_abs)
                classpath = os.pathsep.join(jar_files_clean)
                java_args_cp = java_base + [
                    "-cp", classpath,
                    "uk.me.parabola.mkgmap.Main",
                    f"--output-dir={mkgmap_output_dir}",
                    f"--mapname={mapname}",
                    f"--tdbfile={base_name}.tdb",
                    osm_path
                ]
                try:
                    result = subprocess.run(
                        java_args_cp,
                        check=True,
                        capture_output=True,
                        text=True
                    )
                    print(f"  ✓ mkgmap completed successfully")
                except subprocess.CalledProcessError as e2:
                    error_msg = f"mkgmap failed: {e2.stderr}"
                    if e2.stdout:
                        error_msg += f"\nstdout: {e2.stdout}"
                    raise RuntimeError(error_msg)
            else:
                error_msg = f"mkgmap failed: {e.stderr}"
                if e.stdout:
                    error_msg += f"\nstdout: {e.stdout}"
                raise RuntimeError(error_msg)
        else:
            # Success with -jar method
            print(f"  ✓ mkgmap completed successfully")
        
        # Collect all generated files
        generated_files = {}
        
        # Find .img file
        img_filename = f"{mapname}.img"
        img_path = os.path.join(mkgmap_output_dir, img_filename)
        if not os.path.exists(img_path):
            img_files = [f for f in os.listdir(mkgmap_output_dir) if f.endswith('.img')]
            if img_files:
                img_path = os.path.join(mkgmap_output_dir, img_files[0])
                img_filename = img_files[0]
            else:
                raise RuntimeError(f"mkgmap did not create IMG file in {mkgmap_output_dir}")
        generated_files['img'] = img_path
        generated_files['img_name'] = img_filename
        
        # Find .tdb file
        tdb_filename = f"{base_name}.tdb"
        tdb_path = os.path.join(mkgmap_output_dir, tdb_filename)
        if not os.path.exists(tdb_path):
            tdb_files = [f for f in os.listdir(mkgmap_output_dir) if f.endswith('.tdb')]
            if tdb_files:
                tdb_path = os.path.join(mkgmap_output_dir, tdb_files[0])
                tdb_filename = tdb_files[0]
        if os.path.exists(tdb_path):
            generated_files['tdb'] = tdb_path
            generated_files['tdb_name'] = tdb_filename
        
        # Find .typ file
        typ_files = [f for f in os.listdir(mkgmap_output_dir) if f.endswith('.typ')]
        if typ_files:
            typ_path = os.path.join(mkgmap_output_dir, typ_files[0])
            generated_files['typ'] = typ_path
            generated_files['typ_name'] = typ_files[0]
        
        # Find .mdx file
        mdx_files = [f for f in os.listdir(mkgmap_output_dir) if f.endswith('.mdx')]
        if mdx_files:
            mdx_path = os.path.join(mkgmap_output_dir, mdx_files[0])
            generated_files['mdx'] = mdx_path
            generated_files['mdx_name'] = mdx_files[0]
        
        generated_files['mapname'] = mapname
        generated_files['base_name'] = base_name
        
        return generated_files
        
    finally:
        # Clean up temporary directory
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)


def tif_to_img(input_tif, output_path, mkgmap_path=None, java_max_memory=None):
    """
    Convert GeoTIFF to Garmin IMG format using gdal_polygonize and mkgmap.
    
    Args:
        input_tif: Path to input GeoTIFF file
        output_path: Path to output IMG file
        mkgmap_path: Optional path to mkgmap.jar (will auto-detect if not provided)
        java_max_memory: Optional string for Java -Xmx (e.g., '4G')
        
    Returns:
        Path to the created IMG file
    """
    # Find mkgmap.jar
    mkgmap_jar = find_mkgmap_jar(mkgmap_path)
    
    # Create temporary directory for intermediate files
    temp_dir = create_project_temp_dir("garmin_img_")
    
    try:
        base_name = os.path.splitext(os.path.basename(input_tif))[0]
        shapefile_path = os.path.join(temp_dir, base_name + ".shp")
        osm_path = os.path.join(temp_dir, base_name + ".osm")
        mkgmap_output_dir = os.path.join(temp_dir, "mkgmap_output")
        os.makedirs(mkgmap_output_dir, exist_ok=True)
        
        # Step 1: Convert GeoTIFF to Shapefile using gdal_polygonize
        print(f"Step 1/3: Converting GeoTIFF to Shapefile format...")
        print(f"  Running gdal_polygonize.py (this may take a while for large files)...")
        try:
            result = subprocess.run(
                ["gdal_polygonize.py", input_tif, "-f", "ESRI Shapefile", shapefile_path],
                check=True,
                capture_output=True,
                text=True
            )
            print(f"  ✓ Shapefile created successfully")
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
        print(f"Step 2/3: Converting Shapefile to OSM format...")
        try:
            shapefile_to_osm(shapefile_path, osm_path)
        except Exception as e:
            raise RuntimeError(f"Failed to convert Shapefile to OSM: {str(e)}")
        
        if not os.path.exists(osm_path):
            raise RuntimeError("Failed to create OSM file")
        
        # Step 3: Convert OSM to IMG using mkgmap
        print(f"Step 3/3: Converting OSM to IMG format using mkgmap...")
        print(f"  Processing OSM data with mkgmap (this may take a while)...")
        mapname = generate_mapname(input_tif)
        
        # Check if Java is available
        java_cmd = shutil.which("java")
        if not java_cmd:
            raise FileNotFoundError(
                "Java runtime not found. Please install Java (JRE) to use IMG conversion."
            )
        
        # Build classpath for mkgmap - collect JAR files from project directory
        current_dir = os.getcwd()
        mkgmap_jar_abs = os.path.abspath(mkgmap_jar)
        
        # Start with mkgmap.jar
        jar_files = [mkgmap_jar_abs]
        jar_paths_set = {mkgmap_jar_abs}
        
        # Check for lib directory in current directory
        lib_dir = os.path.join(current_dir, "lib")
        if os.path.exists(lib_dir) and os.path.isdir(lib_dir):
            for file in os.listdir(lib_dir):
                if file.endswith('.jar'):
                    jar_path = os.path.abspath(os.path.join(lib_dir, file))
                    if jar_path not in jar_paths_set:
                        jar_files.append(jar_path)
                        jar_paths_set.add(jar_path)
        
        # Build Java command
        # According to mkgmap docs: java -jar mkgmap.jar [options] [input files]
        # Try -jar first (standard method), fall back to -cp if needed
        
        # First try the standard -jar method (works even with dependencies in lib/)
        print(f"Using mkgmap.jar (standard method)")
        java_base = ["java"]
        if java_max_memory:
            java_base.append(f"-Xmx{java_max_memory}")
        java_args_jar = java_base + [
            "-jar", mkgmap_jar,
            f"--output-dir={mkgmap_output_dir}",
            f"--mapname={mapname}",
            osm_path
        ]
        
        try:
            # Try -jar method first
            result = subprocess.run(
                java_args_jar,
                check=True,
                capture_output=True,
                text=True
            )
        except subprocess.CalledProcessError as e:
            # If -jar fails with missing dependencies, try classpath method
            if "NoClassDefFoundError" in e.stderr or "ClassNotFoundException" in e.stderr:
                print(f"  -jar method failed, trying with classpath (found {len(jar_files)} JAR files)")
                # Ensure mkgmap.jar is first in the classpath
                jar_files_clean = [j for j in jar_files if j != mkgmap_jar_abs]
                jar_files_clean.insert(0, mkgmap_jar_abs)
                classpath = os.pathsep.join(jar_files_clean)
                java_args_cp = java_base + [
                    "-cp", classpath,
                    "uk.me.parabola.mkgmap.Main",
                    f"--output-dir={mkgmap_output_dir}",
                    f"--mapname={mapname}",
                    osm_path
                ]
                try:
                    result = subprocess.run(
                        java_args_cp,
                        check=True,
                        capture_output=True,
                        text=True
                    )
                    print(f"  ✓ mkgmap completed successfully")
                except subprocess.CalledProcessError as e2:
                    # Both methods failed, format error message
                    error_msg = f"mkgmap failed: {e2.stderr}"
                    if e2.stdout:
                        error_msg += f"\nstdout: {e2.stdout}"
                    
                    # Check for main class not found error
                    if "Could not find or load main class" in error_msg:
                        error_msg += (
                            f"\n\nMain class 'uk.me.parabola.mkgmap.Main' not found.\n"
                            f"Troubleshooting:\n"
                            f"  1. Verify mkgmap.jar is valid: try 'java -jar {mkgmap_jar_abs} --version'\n"
                            f"  2. If that works, the issue is with the classpath setup\n"
                            f"  3. Check that all JAR files in lib/ are valid"
                        )
                    raise RuntimeError(error_msg)
            else:
                # Re-raise if it's a different error
                error_msg = f"mkgmap failed: {e.stderr}"
                if e.stdout:
                    error_msg += f"\nstdout: {e.stdout}"
                raise RuntimeError(error_msg)
        else:
            # Success with -jar method
            print(f"  ✓ mkgmap completed successfully")
        
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


def tif_to_gmapi(input_tif, output_path, mkgmap_path=None, map_name=None, 
                  family_id=None, product_id=1, data_version="27", create_zip=False,
                  java_max_memory=None):
    """
    Convert GeoTIFF to Garmin GMAPI format (BaseCamp-compatible).
    
    Args:
        input_tif: Path to input GeoTIFF file
        output_path: Path to output GMAPI directory or .gmapi file
        mkgmap_path: Optional path to mkgmap.jar (will auto-detect if not provided)
        map_name: Name for the map (default: derived from input filename)
        family_id: Family ID for the map (default: derived from filename hash)
        product_id: Product ID (default: 1)
        data_version: Data version string (default: "27")
        create_zip: If True, create a .gmapi.zip file instead of directory
        java_max_memory: Optional string for Java -Xmx (e.g., '6G')
        
    Returns:
        Path to the created GMAPI directory or file
    """
    import hashlib
    
    # Determine output directory
    if output_path.endswith('.gmapi.zip'):
        # Remove both .gmapi.zip extensions to get directory name
        output_dir = Path(output_path[:-10])  # Remove '.gmapi.zip'
    elif output_path.endswith('.gmapi'):
        # Remove .gmapi extension to get directory name
        output_dir = Path(output_path[:-6])  # Remove '.gmapi'
    else:
        output_dir = Path(output_path)
    
    # Clean and create gmapi directory
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)
    
    # Create subproduct folder (Garmin expects "Product1")
    subproduct_dir = output_dir / "Product1"
    subproduct_dir.mkdir()
    
    # Generate map name if not provided
    if map_name is None:
        base_name = os.path.splitext(os.path.basename(input_tif))[0]
        map_name = base_name.replace("_", " ").title()
    
    # Generate family_id if not provided (8-digit number from filename hash)
    if family_id is None:
        base_name = os.path.splitext(os.path.basename(input_tif))[0]
        hash_obj = hashlib.md5(base_name.encode())
        hash_hex = hash_obj.hexdigest()
        family_id = int(hash_hex[:8], 16) % 100000000  # Ensure it's 8 digits max
    
    # Create temporary directory for mkgmap output
    temp_mkgmap_dir = create_project_temp_dir("mkgmap_gmapi_")
    
    try:
        # Run mkgmap to generate all map files
        print("Generating map files with mkgmap...")
        generated_files = run_mkgmap_for_gmapi(input_tif, temp_mkgmap_dir, mkgmap_path, java_max_memory)
        
        # Copy map tiles and support files to appropriate locations
        print("Packaging GMAPI files...")
        
        # Determine file names for XML
        base_name = generated_files.get('base_name', os.path.splitext(os.path.basename(input_tif))[0])
        tdb_name = generated_files.get('tdb_name', f"{base_name}.tdb")
        mdx_name = generated_files.get('mdx_name', f"{base_name}.mdx")
        typ_name = generated_files.get('typ_name', "100D8.TYP")
        
        # BaseMap name should match the base part of TDB/MDX filenames (without extension)
        # Extract base name from TDB filename (remove .tdb extension)
        if tdb_name and tdb_name.endswith('.tdb'):
            base_map_name = tdb_name[:-4]  # Remove .tdb extension
        elif mdx_name and mdx_name.endswith('.mdx'):
            base_map_name = mdx_name[:-4]  # Remove .mdx extension
        else:
            # Fallback: use base_name with spaces replaced by underscores
            base_map_name = base_name.replace(" ", "_")
        
        # Copy .mdx and .typ files to root level (not Product1)
        if 'mdx' in generated_files:
            shutil.copy(generated_files['mdx'], output_dir / mdx_name)
            print(f"  ✓ Copied MDX file to root: {mdx_name}")
        
        if 'typ' in generated_files:
            shutil.copy(generated_files['typ'], output_dir / typ_name)
            print(f"  ✓ Copied TYP file to root: {typ_name}")
        
        # Copy .tdb file to Product1 directory
        if 'tdb' in generated_files:
            shutil.copy(generated_files['tdb'], subproduct_dir / tdb_name)
            print(f"  ✓ Copied TDB file to Product1: {tdb_name}")
        
        # Copy .img files to Product1 directory (if any)
        if 'img' in generated_files:
            shutil.copy(generated_files['img'], subproduct_dir / generated_files['img_name'])
            print(f"  ✓ Copied IMG file to Product1: {generated_files['img_name']}")
        
        # Copy all other files and subdirectories from mkgmap output to Product1
        # (excluding .mdx, .typ which go to root, and .tdb which we already copied)
        mkgmap_output = Path(temp_mkgmap_dir)
        excluded_files = set()
        if 'mdx' in generated_files:
            excluded_files.add(Path(generated_files['mdx']).name)
        if 'typ' in generated_files:
            excluded_files.add(Path(generated_files['typ']).name)
        if 'tdb' in generated_files:
            excluded_files.add(Path(generated_files['tdb']).name)
        if 'img' in generated_files:
            excluded_files.add(Path(generated_files['img']).name)
        
        # Copy all files and directories from mkgmap output
        for item in mkgmap_output.iterdir():
            item_name = item.name
            if item_name not in excluded_files:
                dest = subproduct_dir / item_name
                if item.is_dir():
                    if dest.exists():
                        shutil.rmtree(dest)
                    shutil.copytree(item, dest)
                    print(f"  ✓ Copied directory to Product1: {item_name}/")
                elif item.is_file():
                    shutil.copy(item, dest)
                    print(f"  ✓ Copied file to Product1: {item_name}")
        
        # Create Info.xml (not MapProduct.xml) at root level
        print("Creating Info.xml...")
        root = ET.Element("MapProduct", xmlns="http://www.garmin.com/xmlschemas/MapProduct/v1")
        
        name_elem = ET.SubElement(root, "Name")
        name_elem.text = map_name
        
        dataversion_elem = ET.SubElement(root, "DataVersion")
        dataversion_elem.text = data_version
        
        dataformat_elem = ET.SubElement(root, "DataFormat")
        dataformat_elem.text = "Original"
        
        fid_elem = ET.SubElement(root, "ID")
        fid_elem.text = str(family_id)
        
        idx_elem = ET.SubElement(root, "IDX")
        idx_elem.text = mdx_name
        
        typ_elem = ET.SubElement(root, "TYP")
        typ_elem.text = typ_name
        
        sub_elem = ET.SubElement(root, "SubProduct")
        ET.SubElement(sub_elem, "Name").text = map_name
        ET.SubElement(sub_elem, "ID").text = str(product_id)
        ET.SubElement(sub_elem, "BaseMap").text = base_map_name
        ET.SubElement(sub_elem, "TDB").text = tdb_name
        ET.SubElement(sub_elem, "Directory").text = "Product1"
        
        # Save XML with formatting matching the example
        xml_path = output_dir / "Info.xml"
        tree = ET.ElementTree(root)
        # ET.indent is available in Python 3.9+
        if hasattr(ET, 'indent'):
            ET.indent(tree, space="  ")
        # Write with XML declaration matching example format
        tree.write(xml_path, encoding="UTF-8", xml_declaration=True)
        
        # Post-process to match example formatting (add newlines between elements)
        with open(xml_path, 'r', encoding='UTF-8') as f:
            content = f.read()
        
        # Add newlines between major elements to match example format
        # Replace closing/opening tag pairs with newlines
        content = content.replace('</Name>\n    <DataVersion>', '</Name>\n\n  <DataVersion>')
        content = content.replace('</DataVersion>\n    <DataFormat>', '</DataVersion>\n\n  <DataFormat>')
        content = content.replace('</DataFormat>\n    <ID>', '</DataFormat>\n\n  <ID>')
        content = content.replace('</ID>\n    <IDX>', '</ID>\n\n  <IDX>')
        content = content.replace('</IDX>\n    <TYP>', '</IDX>\n\n  <TYP>')
        content = content.replace('</TYP>\n    <SubProduct>', '</TYP>\n\n  <SubProduct>')
        content = content.replace('</SubProduct>\n</MapProduct>', '</SubProduct>\n\n</MapProduct>')
        
        # Fix the opening tag if needed
        content = content.replace('<MapProduct xmlns=', '\n<MapProduct xmlns=')
        
        # Add standalone="no" to XML declaration to match example
        if 'standalone=' not in content:
            content = content.replace('<?xml version="1.0" encoding="UTF-8"?>', 
                                     '<?xml version="1.0" encoding="UTF-8" standalone="no" ?>')
        
        # Clean up any double newlines at the start
        content = content.lstrip('\n')
        
        with open(xml_path, 'w', encoding='UTF-8') as f:
            f.write(content)
        
        print(f"  ✓ Created Info.xml")
        
        # Optionally create zip file
        final_output = output_dir
        if create_zip or output_path.endswith('.gmapi.zip') or output_path.endswith('.gmapi'):
            print("Creating GMAPI archive...")
            zip_path = Path(str(output_dir) + ".zip")
            shutil.make_archive(str(output_dir), "zip", output_dir)
            if zip_path.exists():
                # Rename to .gmapi or .gmapi.zip if requested
                if output_path.endswith('.gmapi'):
                    gmapi_path = Path(output_path)
                    zip_path.rename(gmapi_path)
                    final_output = gmapi_path
                    print(f"  ✓ Created GMAPI package: {gmapi_path.resolve()}")
                elif output_path.endswith('.gmapi.zip'):
                    gmapi_zip_path = Path(output_path)
                    zip_path.rename(gmapi_zip_path)
                    final_output = gmapi_zip_path
                    print(f"  ✓ Created GMAPI zip package: {gmapi_zip_path.resolve()}")
                else:
                    final_output = zip_path
                    print(f"  ✓ Created GMAPI zip package: {zip_path.resolve()}")
        else:
            print(f"  ✓ GMAPI package created at: {output_dir.resolve()}")
        
        return str(final_output)
        
    finally:
        # Clean up temporary mkgmap directory
        if os.path.exists(temp_mkgmap_dir):
            shutil.rmtree(temp_mkgmap_dir)


def find_tif_files(directory):
    """
    Find all TIF/TIFF files in the specified directory.
    
    Args:
        directory: Directory to search for TIF files
        
    Returns:
        List of paths to TIF files
    """
    tif_files = []
    if not os.path.exists(directory) or not os.path.isdir(directory):
        return tif_files
    
    for file in os.listdir(directory):
        file_lower = file.lower()
        if file_lower.endswith('.tif') or file_lower.endswith('.tiff'):
            full_path = os.path.join(directory, file)
            if os.path.isfile(full_path):
                tif_files.append(full_path)
    
    return sorted(tif_files)


def select_tif_file_interactive(main_folder):
    """
    Display TIF files in the main folder and allow user to select one.
    
    Args:
        main_folder: Directory to search for TIF files
        
    Returns:
        Path to selected TIF file, or None if cancelled
    """
    tif_files = find_tif_files(main_folder)
    
    if not tif_files:
        print(f"No TIF files found in {main_folder}")
        return None
    
    print(f"\nFound {len(tif_files)} TIF file(s) in {main_folder}:\n")
    for i, tif_file in enumerate(tif_files, 1):
        filename = os.path.basename(tif_file)
        file_size = os.path.getsize(tif_file)
        size_mb = file_size / (1024 * 1024)
        print(f"  {i}) {filename} ({size_mb:.2f} MB)")
    
    print()
    while True:
        try:
            user_input = input(f"Select a file (1-{len(tif_files)}) or 'q' to quit: ").strip()
            if user_input.lower() == 'q':
                print("Cancelled.")
                return None
            
            try:
                selection = int(user_input)
                if 1 <= selection <= len(tif_files):
                    selected_file = tif_files[selection - 1]
                    print(f"Selected: {os.path.basename(selected_file)}\n")
                    return selected_file
                else:
                    print(f"Invalid selection. Please enter a number between 1 and {len(tif_files)}.")
            except ValueError:
                print("Invalid input. Please enter a number or 'q' to quit.")
        except (EOFError, KeyboardInterrupt):
            print("\nCancelled.")
            return None


def main():
    parser = argparse.ArgumentParser(
        description='Convert GeoTIFF files to Garmin-compatible KMZ, IMG, or GMAPI files'
    )
    parser.add_argument(
        'input',
        nargs='?',
        help='Path to input GeoTIFF file (optional - if not provided, will search for TIF files in the main folder)'
    )
    parser.add_argument(
        '-o', '--output',
        help='Path to output file (default: input filename with .kmz, .img, or .gmapi extension)',
        default=None
    )
    parser.add_argument(
        '--type',
        choices=['kmz', 'img', 'gmapi'],
        default=None,
        help='Conversion type: kmz, img, or gmapi (if not provided, will prompt interactively)'
    )
    parser.add_argument(
        '--mkgmap-path',
        help='Path to mkgmap.jar (required for IMG and GMAPI conversion, will auto-detect if not provided)',
        default=None
    )
    parser.add_argument(
        '--temp-dir',
        help='Temporary directory for intermediate files (default: system temp)',
        default=None
    )
    parser.add_argument(
        '--java-max-memory',
        default='4G',
        help="Maximum Java heap size passed to mkgmap via -Xmx (examples: '4G', '8192m')."
    )
    
    args = parser.parse_args()
    
    # If no input file provided, search for TIF files in the main folder
    if args.input is None:
        # Get the directory where the script is located
        script_dir = os.path.dirname(os.path.abspath(__file__))
        if script_dir == '':
            script_dir = os.getcwd()
        
        selected_file = select_tif_file_interactive(script_dir)
        if selected_file is None:
            return 1
        args.input = selected_file
    
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
                user_input = input("Select conversion type: 1) KMZ  2) IMG  3) GMAPI [1]: ").strip()
                if user_input == '' or user_input == '1':
                    conversion_type = 'kmz'
                    break
                elif user_input == '2':
                    conversion_type = 'img'
                    break
                elif user_input == '3':
                    conversion_type = 'gmapi'
                    break
                else:
                    print("Invalid choice. Please enter 1, 2, or 3.")
            except (EOFError, KeyboardInterrupt):
                print("\nCancelled.")
                return 1
    
    # Determine output path based on conversion type
    if args.output is None:
        base_name = os.path.splitext(os.path.basename(args.input))[0]
        output_dir = os.path.dirname(args.input) or '.'
        if conversion_type == 'img':
            extension = '.img'
        elif conversion_type == 'gmapi':
            extension = '.gmapi'
        else:
            extension = '.kmz'
        args.output = os.path.join(output_dir, base_name + extension)
    
    # Branch based on conversion type
    if conversion_type == 'kmz':
        # KMZ conversion (existing flow)
        # Create temporary directory
        temp_dir = args.temp_dir
        if temp_dir is None:
            temp_dir = create_project_temp_dir("garmin_converter_")
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
        # IMG conversion
        try:
            tif_to_img(args.input, args.output, args.mkgmap_path, args.java_max_memory)
            return 0
        except Exception as e:
            print(f"Error: {str(e)}")
            return 1
    
    elif conversion_type == 'gmapi':
        # GMAPI conversion (BaseCamp-compatible)
        try:
            tif_to_gmapi(args.input, args.output, args.mkgmap_path, java_max_memory=args.java_max_memory)
            return 0
        except Exception as e:
            print(f"Error: {str(e)}")
            return 1


if __name__ == "__main__":
    exit(main())

