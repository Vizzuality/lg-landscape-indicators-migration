#!/usr/bin/env python
# coding: utf-8

# # 1. Landscape indicator computation notebook
#
# Scripts for preprocessing and generating landscape indicators using the 50km sLUC methodology
#  - Deforestation,
#  - Deforestation carbon, and
#  - Cropland expansion in natural lands
#
# General approach:
#  - First, we create or resample source layers to 100m.
#  - Then we do the overlay analysis to identify deforestation, cropland expansion, etc.
#  - Then we use these outputs to compute the final kerneled layers at 1000m.
#
# Assets will be saved to your `GEE_PROJECT` in the `WORKING_FOLDER` defined below.
# Generates ~150GB of data.

# In[1]:


# Imports

import os
import time

import ee
import eeUtil
import geemap.foliumap as gmap
from ee import EEException

# Initialize Earth Engine
PROJECT = os.getenv("GEE_PROJECT")
GEE_JSON = os.getenv("GEE_JSON")

assert PROJECT is not None, "Please set GEE_PROJECT environment variable"
assert (
  GEE_JSON is not None
), "Please set GEE_JSON environment variable with service account credentials"

eeUtil.init()

# In[ ]:


# Constants / options for the landscape indicators analysis
PROJECTION = "EPSG:4326"
WORKING_FOLDER = f"projects/{PROJECT}/assets/landscape_indicators_20230821"
SRC_FOLDER = "projects/ee-fgassert/assets/landscape_indicators_source"
WORLD_GEOM = ee.Geometry.Polygon(
  coords=[[[-180, -85], [-180, 85], [180, 85], [180, -85], [-180, -85]]],
  proj=PROJECTION,
  geodesic=False,
)
EXPORT_BUCKET = os.getenv("GEE_BUCKET")
EXPORT_PREFIX = "landscape_indicators_20231023"

TARGET_YEAR = 2022  # the target year for analysis
START_YEAR = TARGET_YEAR - 20  # the start year for deforestation analysis

# scales in meters
ANALYSIS_SCALE = 100  # the working scale for overlays in the main analysis
KERNEL_SCALE = 1000  # the scale for the kernel analysis
KERNEL_RADIUS = 50000

# additional radii to compute for validation
RADII = [500, 5000, 10000, 25000, 50000]

ALL_INPUT_LAYERS = [
  f"tree_loss_to{TARGET_YEAR}",
  "nonnat_forest",
  "intact_forests",
  "primary_tropical_forests",
  "burned_area_mask",
  "disturbance",
  "carbon_filled_2000",
  "esri_cropland_2020",
  f"esri_cropland_{TARGET_YEAR}",
  "natural_lands",
  "forest_loss_fire",
]

# layers that must be computed at a higher resolution before downsampling for the kernel
PRE_KERNEL_LAYERS = [
  f"tree_loss_to{TARGET_YEAR}",
  "deforest",
  "deforest_carbon",
  "natcrop_expansion",
  "natcrop_reduction",
  "natcrop_flii_loss",
  "natcrop_bii_loss",
  "nonnatural_excl_builtwater",
]

#### ASSETS ########

####################
# derived from RESOLVE ecoregions
NON_SUBTROPIC_ASSET = f"{SRC_FOLDER}/not_subtropic"
# Noon et al. 2022 vulnerable carbon total for 2010
CARBON_ASSET = f"{SRC_FOLDER}/Vulnerable_C_Total_2010"
# ESA CCI land cover for 2010
ESA_CCI_ASSET = f"{SRC_FOLDER}/ESACCI_LC_300m_P1Y_2010_v207"
# Forest Landscape Integrity Index for 2020
FLII_ASSET = f"{SRC_FOLDER}/flii_earth"

####################
# 3rd party hosted assets
HANSEN_ASSET = "UMD/hansen/global_forest_change_2022_v1_10"
FOREST_DYNAMICS_ASSET = "projects/glad/GLCLU2020/Forest_type"
INTACT_FORESTS_ASSET = "users/potapovpeter/IFL_2000"
PRIMARY_TROPICAL_ASSET = "UMD/GLAD/PRIMARY_HUMID_TROPICAL_FORESTS/v1/2001"

ESRI_LULC_IC = "projects/sat-io/open-datasets/landcover/ESRI_Global-LULC_10m_TS"
CROPLAND_2019_IC = "users/potapovpeter/Global_cropland_2019"
SBTN_ASSET = "projects/wri-datalab/SBTN/natLands_beta/naturalLands_allClasses_20230516"

FOREST_LOSS_FIRE_IC = "users/sashatyu/2001-2022_fire_forest_loss"
FIRECCI_IC = "ESA/CCI/FireCCI/5_1"
MODIS_BA_IC = "MODIS/006/MCD64A1"

BII_IC = "projects/ebx-data/assets/earthblox/IO/BIOINTACT"

# In[3]:


# create working folder if not exists
eeUtil.createFolder(WORKING_FOLDER)


# In[4]:


def get_sbtn_layers():
  """Get relevant SBTN layers.

  All Classes values:
  Natural:     2 = natural forests, 3 = natural short vegetation, 4 = natural water, 5 = mangroves
               6 = bare, 7 = snow, 8 = wet natural forests, 9 = natural peat forests,
               10 = wet natural short vegetation, 11 = natural peat short vegetation
  Non-Natural: 12 = cropland, 13 = built, 14 = non-natural tree cover, 15 = non-natural short
               vegetation, 16 = non-natural water, 17 = wet non-natural forests,
               18 = non-natural peat forests, 19 = wet non-natural short vegetation,
               20 = non-natural peat short vegetation.
  """
  sbtn_all_classes = ee.Image(SBTN_ASSET)
  return {
    "natural_lands": sbtn_all_classes.lte(11),
    "nonnat_forest": sbtn_all_classes.eq(ee.Image([14, 17, 18])).reduce("anyNonZero"),
    "nonnatural_excl_builtwater": sbtn_all_classes.gte(12)
    .And(sbtn_all_classes.neq(13))
    .And(sbtn_all_classes.neq(16)),
    "available_area_mask": sbtn_all_classes.eq([4, 13, 16]).reduce("anyNonZero").Not(),
  }


def get_mask():
  """Get mask of available area at target resolution."""
  return ee.Image(SBTN_ASSET).eq([4, 13, 16]).reduce("anyNonZero").Not()


def get_tree_loss():
  """Gets Hansen tree loss ensuring proportional reduction in reduced resoultions."""
  hanson_data = ee.Image(HANSEN_ASSET)
  loss_year = hanson_data.select(3)
  loss_portion = hanson_data.select(1)

  return {
    "tree_loss_allyrs": loss_portion,
    f"tree_loss_to{TARGET_YEAR}": loss_portion.multiply(loss_year.gte(START_YEAR - 2000)),
  }


def get_forest_loss_fire():
  return {
    "forest_loss_fire": ee.ImageCollection(FOREST_LOSS_FIRE_IC).mosaic().gt(1).unmask(),
  }


def get_esri_cropland(start_year=2020, end_year=TARGET_YEAR, lag_years=3):
  """Gets Esri cropland at reduced resolution using mean pyramiding.

  Each year represents the maximum cropland area in the previous `lag_years`.
  """

  esri_lulc = ee.ImageCollection(ESRI_LULC_IC)
  esri_crop_annual = ee.List.sequence(start_year - lag_years, end_year).map(
    lambda y: (
      esri_lulc.filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y, 12, 31))
      .mosaic()
      .eq(5)
    )
  )
  esri_crop_lagged = ee.List.sequence(0, esri_crop_annual.length().subtract(lag_years + 1)).map(
    lambda i: (
      ee.ImageCollection(esri_crop_annual.slice(i, ee.Number(i).add(lag_years + 1))).reduce(
        "max"
      )
    )
  )

  layers = {}
  for i in range(end_year - start_year + 1):
    y = start_year + i
    layers[f"esri_cropland_{y}"] = ee.Image(esri_crop_lagged.get(i))

  return layers


def get_forest_dynamics():
  """Gets Potapov forest dynamics at reduced resolution using mean pyramiding."""
  forest_dynamics = ee.Image(FOREST_DYNAMICS_ASSET)
  return {"gain": forest_dynamics.eq(3), "disturbance": forest_dynamics.eq(4)}


def get_burned_area(kernel_radius=250):
  """Gets FireCII and Modis Burned Area products

  FireCCI is a bit more accurate and reported at 250m but only runs through 2020
  MODIS Burned Area is reported at 500m but is operationally near real-time

  We run both though a kernel to smooth out the edges and expand the filter region
  """
  firecci = ee.ImageCollection(FIRECCI_IC)
  modis_ba = ee.ImageCollection(MODIS_BA_IC)

  firecci = (
    firecci.filterDate(ee.Date.fromYMD(START_YEAR, 1, 1), ee.Date.fromYMD(TARGET_YEAR, 12, 31))
    .select(0)
    .mosaic()
    .setDefaultProjection(firecci.first().select(0).projection())
  )
  modis_ba = (
    modis_ba.filterDate(ee.Date.fromYMD(START_YEAR, 1, 1), ee.Date.fromYMD(TARGET_YEAR, 12, 31))
    .select(0)
    .mosaic()
    .setDefaultProjection(modis_ba.first().select(0).projection())
  )

  layers = {"firecci": firecci, "modis_ba": modis_ba, "burned_area_mask": firecci.Or(modis_ba)}

  return {
    k: image.reduceNeighborhood(
      "sum", ee.Kernel.square(kernel_radius, "meters"), None, False, "boxcar"
    ).gt(0)
    for k, image in layers.items()
  }


def get_vulnerable_carbon_filled():
  """Computes a filled version of the Noon et al. vulnerable carbon layer for 2000.

  The original layer is at 300m for the year 2010. We fill in forest areas for the year 2000
  by looking at the difference in forest cover, and filling with the mean carbon value for forests
  in a 10km radius around the forest pixel.
  """
  carbon = ee.Image(CARBON_ASSET).unmask()
  tree_cover_2000 = ee.Image(HANSEN_ASSET).select(1).unmask()
  esa_cci = ee.Image(ESA_CCI_ASSET).unmask()

  cci_tree_cover_2010 = (
    esa_cci.gte(40).And(esa_cci.lt(100)).Or(esa_cci.eq(160)).Or(esa_cci.eq(170))
  )
  tree2000_not_2010 = tree_cover_2000.And(cci_tree_cover_2010.Not())
  carbon_kernel_mean = (
    carbon.multiply(cci_tree_cover_2010)
    .selfMask()
    .reduceNeighborhood("mean", ee.Kernel.square(10000, "meters"), "mask", False, "boxcar")
  )
  carbon_filled_2000 = carbon.unmask().where(tree2000_not_2010, carbon_kernel_mean).selfMask()

  return {"carbon": carbon, "carbon_filled_2000": carbon_filled_2000}


def get_intact_forests():
  """Get intact and primary humid tropical forest layers."""
  return {
    "intact_forests": ee.Image(INTACT_FORESTS_ASSET),
    "primary_tropical_forests": ee.Image(PRIMARY_TROPICAL_ASSET),
  }


def get_cropland_2019():
  """Get potapov cropland layer for 2019."""
  return {"cropland_2019": ee.ImageCollection(CROPLAND_2019_IC).mosaic()}


def get_flii():
  """Get forest landsacpe integrity index."""
  return {"flii": ee.Image(FLII_ASSET).divide(10000)}


def get_bii():
  """Get forest biodiversity intactness for 2020."""
  return {
    "bii": ee.ImageCollection(BII_IC)
    .filterDate("2020-01-01", "2020-12-31")
    .mosaic()
    .setDefaultProjection(PROJECTION, None, 100)
  }


def get_all_input_layers(cache=False):
  """Get all layers used in the analysis."""
  return {
    **get_sbtn_layers(),
    **get_esri_cropland(),
    **get_tree_loss(),
    **get_intact_forests(),
    **get_vulnerable_carbon_filled(),
    **get_forest_loss_fire(),
    **get_burned_area(),
    **get_forest_dynamics(),
    **get_cropland_2019(),
    **get_flii(),
    **get_bii(),
  }


# In[5]:


# functions for logical operations on 0-1 floating point rasters


def and_assume_overlap(*args):
  return ee.Image([*args]).reduce("min")


def and_assume_uniform(*args):
  return ee.Image([*args]).reduce("product")


def and_assume_disjoint(*args):
  im = ee.Image([*args]).reduce("sum")
  return im.subtract(im.mask()).max(0)


def or_assume_overlap(*args):
  return ee.Image([*args]).reduce("max")


def or_assume_uniform(*args):
  return inverse(inverse(args).reduce("product"))


def or_assume_disjoint(*args):
  return ee.Image([*args]).reduce("sum").min(1)


def inverse(a):
  return ee.Image(a).subtract(1).multiply(-1)


# In[6]:


def compute_deforest(layers=None):
  """Calculate deforestation from forest loss excluding likely types of non-deforestation loss.

  Returns the input layers plus the new 'deforest' layer
  """
  layers = layers or get_all_input_layers()

  tree_loss = layers[f"tree_loss_to{TARGET_YEAR}"]
  non_nat_forest = layers["nonnat_forest"]
  not_primary_forest = inverse(
    or_assume_overlap(layers["intact_forests"], layers["primary_tropical_forests"])
  )
  disturbance = layers["disturbance"]
  layers["burned_area_mask"]
  forest_loss_fire = layers["forest_loss_fire"]
  not_subtropic = ee.Image(NON_SUBTROPIC_ASSET).unmask()

  forest_ok = and_assume_overlap(non_nat_forest, not_primary_forest)
  disturbance_ok = and_assume_overlap(disturbance, not_primary_forest)
  # burned_area_ok = and_assume_uniform(burned_area, not_subtropic)
  burned_area_ok = and_assume_uniform(forest_loss_fire, not_subtropic)
  loss_ok = or_assume_overlap(forest_ok, disturbance_ok, burned_area_ok)

  deforest = tree_loss.subtract(loss_ok).max(0)

  return {**layers, "deforest": deforest}


def compute_deforest_carbon(layers=None):
  """Calculate the carbon loss from deforestation.

  Returns the input layers plus the new 'deforest_carbon' layer
  """
  layers = layers or compute_deforest()

  deforest = layers["deforest"]
  deforest_carbon = deforest.multiply(layers["carbon_filled_2000"])

  return {**layers, "deforest_carbon": deforest_carbon}


def compute_cropland_expansion(layers=None):
  """Compute cropland expansion and reduction layers from Esri cropland data

  Mask to natural lands to get natural cropland expansion and reduction layers.
  """
  layers = layers or get_all_input_layers()

  cropland_end = layers[f"esri_cropland_{TARGET_YEAR}"]
  cropland_2020 = layers["esri_cropland_2020"]
  natural = layers["natural_lands"]
  cropland_expansion = cropland_end.subtract(cropland_2020).max(0)
  cropland_reduction = cropland_2020.subtract(cropland_end).max(0)
  natcrop_expansion = and_assume_disjoint(cropland_expansion, natural)
  natcrop_reduction = and_assume_disjoint(cropland_reduction, natural)

  return {
    **layers,
    "cropland_expansion": cropland_expansion,
    "cropland_reduction": cropland_reduction,
    "natcrop_expansion": natcrop_expansion,
    "natcrop_reduction": natcrop_reduction,
  }


def compute_biodiversity_loss(layers=None):
  """Calulate the forest landscape integrity index and biodiversity intactness index
  in areas of cropland expansion
  """
  layers = layers or compute_cropland_expansion()

  cropland_expansion = layers["natcrop_expansion"]
  flii = layers["flii"]
  bii = layers["bii"]

  return {
    **layers,
    "natcrop_flii_loss": flii.multiply(cropland_expansion),
    "natcrop_bii_loss": bii.multiply(cropland_expansion),
  }


def compute_derived_layers(input_layers=None):
  """Compute the derived land use change layers before kernels"""
  layers = input_layers or get_all_input_layers()
  layers = compute_deforest(layers)
  layers = compute_deforest_carbon(layers)
  layers = compute_cropland_expansion(layers)
  layers = compute_biodiversity_loss(layers)

  return {k: layers[k] for k in PRE_KERNEL_LAYERS}


# In[7]:


def compute_kernels(derived_layers=None, kernel_radius=KERNEL_RADIUS):
  """computed kerneled layers"""
  layers = derived_layers or compute_derived_layers()

  # mask all layers to only non-urban non-water areas
  mask = get_mask()
  kernel_layers = {
    "tree_loss_kernel": layers[f"tree_loss_to{TARGET_YEAR}"],
    "deforest_kernel": layers["deforest"],
    "deforest_carbon_kernel": layers["deforest_carbon"],
    "natcrop_expansion_kernel": layers["natcrop_expansion"],
    "natcrop_reduction_kernel": layers["natcrop_reduction"],
    "natcrop_flii_loss_kernel": layers["natcrop_flii_loss"],
    "natcrop_bii_loss_kernel": layers["natcrop_bii_loss"],
    "nonnatural_kernel": layers["nonnatural_excl_builtwater"],
  }

  kernel = ee.Kernel.circle(radius=kernel_radius, units="meters")
  for k, image in kernel_layers.items():
    image = image.updateMask(mask)
    kernel_layers[k] = image.reduceNeighborhood("mean", kernel, "mask", False)

  return kernel_layers


def compute_indicators(kernel_layers=None):
  """Compute final landscape indicators"""
  layers = kernel_layers or compute_kernels()

  natural_crop_reduction_kernel = layers["natcrop_reduction_kernel"]
  natural_crop_conversion_kernel = layers["natcrop_expansion_kernel"]
  natural_crop_flii_loss_kernel = layers["natcrop_flii_loss_kernel"]
  natural_crop_bii_loss_kernel = layers["natcrop_bii_loss_kernel"]
  deforest_kernel = layers["deforest_kernel"]
  deforest_carbon_kernel = layers["deforest_carbon_kernel"]
  non_nat_kernel = layers["nonnatural_kernel"]
  forest_loss_kernel = layers["tree_loss_kernel"]

  # avoid divide by zero errors and cap at 1ha/ha
  human_lu = non_nat_kernel.add(0.000001)

  tree_loss_by_human_lu = forest_loss_kernel.divide(human_lu).min(1).max(0)
  deforest_by_human_lu = deforest_kernel.divide(human_lu)

  # scale carbon by excess deforest to cap at 1ha/ha
  excess_deforest_per_ha_human_lu = deforest_by_human_lu.where(deforest_by_human_lu.lt(1), 1)
  deforest_by_human_lu = deforest_by_human_lu.min(1)
  deforest_carbon_by_human_lu = deforest_carbon_kernel.divide(human_lu).divide(
    excess_deforest_per_ha_human_lu
  )

  # allow shifting ag
  natcrop_net_conversion = natural_crop_conversion_kernel.subtract(
    natural_crop_reduction_kernel
  ).max(0)
  natcrop_conversion_by_human_lu = natural_crop_conversion_kernel.divide(human_lu).min(1).max(0)
  natcrop_net_conversion_by_human_lu = natcrop_net_conversion.divide(human_lu).min(1).max(0)

  # flii loss in cropland expansion areas
  natcrop_flii_loss_by_human_lu = natural_crop_flii_loss_kernel.divide(human_lu).min(1).max(0)
  # bii loss in cropland expansion areas
  natcrop_bii_loss_by_human_lu = natural_crop_bii_loss_kernel.divide(human_lu).min(1).max(0)

  # normalize to annual values
  tree_loss_by_human_lu = tree_loss_by_human_lu.divide(20)
  deforest_by_human_lu = deforest_by_human_lu.divide(20)
  # convert to CO2
  deforest_carbon_by_human_lu = deforest_carbon_by_human_lu.divide(20).multiply(3.66)

  conversion_years = TARGET_YEAR - 2020  # 2020 is baseline year for natural land conversion
  natcrop_conversion_by_human_lu = natcrop_conversion_by_human_lu.divide(conversion_years)
  natcrop_net_conversion_by_human_lu = natcrop_net_conversion_by_human_lu.divide(conversion_years)
  natcrop_flii_loss_by_human_lu = natcrop_flii_loss_by_human_lu.divide(conversion_years)

  return {
    "deforest_by_human_lu": deforest_by_human_lu,
    "deforest_carbon_by_human_lu": deforest_carbon_by_human_lu,
    "tree_loss_by_human_lu": tree_loss_by_human_lu,
    "natural_crop_conversion_by_human_lu": natcrop_conversion_by_human_lu,
    "natural_crop_net_conversion_by_human_lu": natcrop_net_conversion_by_human_lu,
    "natural_crop_flii_loss_by_human_lu": natcrop_flii_loss_by_human_lu,
    "natural_crop_bii_loss_by_human_lu": natcrop_bii_loss_by_human_lu,
  }


def compute_viz_layers(indicator_layers=None):
  indicator_layers = indicator_layers or compute_indicators()
  cropland = get_cropland_2019()["cropland_2019"]
  nonnat = get_sbtn_layers()["nonnatural_excl_builtwater"]
  get_mask()

  layers = {}
  for k, image in indicator_layers.items():
    layers[f"cropland_{k}"] = cropland.multiply(image)
    layers[f"nonnat_{k}"] = nonnat.multiply(image)

  return layers


# In[8]:


def _cache_layers(layers, scale, reduce_resolution=False, default_scale=30):
  """Save or load cached layers, optionally downsampling with reduceResolution."""
  for k, image in layers.items():
    asset_id = f"{WORKING_FOLDER}/{k}_{scale}m"
    image = image.unmask()

    if reduce_resolution:
      _orig_scale = image.projection().nominalScale().getInfo()
      if _orig_scale > 100000:  # GEE lost track of the scale, set it manually
        image = image.setDefaultProjection(
          PROJECTION, None, default_scale
        ).reduceResolution("mean")
      elif _orig_scale < scale:  #
        image = image.reduceResolution("mean")
      # else: original scale is coarser than target scale, don't reduce resolution
      _wait_for_ee_internal_cache(image)

    layers[k] = eeUtil.findOrSaveImage(
      image, asset_id, region=WORLD_GEOM, scale=scale, crs=PROJECTION, dtype="float"
    )
  return layers


def _wait_for_ee_internal_cache(image, retries=3):
  """Sometimes reduce resolution will fail due to there not being a default projection set,
  even though we set one manually when it needs to be. This seems to be related to an
  internal GEE caching issue.
  """
  try:
    image.getInfo()
  except EEException as e:
    if e.args[0].startswith("Image.reduceResolution"):
      time.sleep(15)
      if retries > 0:
        _wait_for_ee_internal_cache(image, retries - 1)
      else:
        raise e


def _export_layers(layers, scale=KERNEL_SCALE):
  """Export layers to GCS"""
  for k, image in layers.items():
    blob = f"{EXPORT_PREFIX}/{k}_{scale}m.tif"
    eeUtil.exportImage(
      image,
      blob,
      bucket=EXPORT_BUCKET,
      region=WORLD_GEOM,
      crs=PROJECTION,
      scale=scale,
      dtype="float",
    )


def _check_layers_precomputed(layers):
  """Check if derived layers have been precomputed at ANALYSIS_SCALE"""
  cached_assets = list(eeUtil.ls(f"{WORKING_FOLDER}"))
  precomputed = all([f"{k}_{ANALYSIS_SCALE}m" in cached_assets for k in layers])
  if not precomputed:
    print("WARNING: Cannot compute next step until the following layers are computed.")
    print([k for k in layers if f"{k}_{ANALYSIS_SCALE}m" not in cached_assets])
    print("Wait for current tasks to complete and then re-run this cell.")
    print("Current tasks:", [t["description"] for t in eeUtil.getTasks(True)])
  return precomputed


def _wait_for_precomputed(layers):
  """Wait for layers to be computed"""
  cached_assets = list(eeUtil.ls(f"{WORKING_FOLDER}"))
  precomputed = [f"{k}_{ANALYSIS_SCALE}m" in cached_assets for k in layers]
  elapsed = 0
  while not all(precomputed):
    print(
      f"Waiting for layers. ({sum(precomputed)}/{len(layers)} complete), {elapsed}m", end=""
    )
    time.sleep(60)
    cached_assets = list(eeUtil.ls(f"{WORKING_FOLDER}"))
    precomputed = [f"{k}_{ANALYSIS_SCALE}m" in cached_assets for k in layers]
    elapsed += 1
  return True


def compute_all_landscape_indicators(kernel_radius=KERNEL_RADIUS):
  """Main function to compute all landscape indicators."""

  input_layers = get_all_input_layers()
  _cache_layers(input_layers, ANALYSIS_SCALE, reduce_resolution=True)
  _wait_for_precomputed(ALL_INPUT_LAYERS)
  input_layers = _cache_layers(input_layers, ANALYSIS_SCALE, reduce_resolution=True)

  derived_layers = compute_derived_layers(input_layers)
  _cache_layers(derived_layers, ANALYSIS_SCALE)
  _wait_for_precomputed(PRE_KERNEL_LAYERS)
  derived_layers = _cache_layers(derived_layers, ANALYSIS_SCALE)

  kernel_layers = compute_kernels(derived_layers, kernel_radius)
  indicators = compute_indicators(kernel_layers)
  viz_layers = compute_viz_layers(indicators)

  rad_km = int(kernel_radius / 1000)
  kernel_layers = _cache_layers(
    {f"{k}_{rad_km}km": i for k, i in kernel_layers.items()}, KERNEL_SCALE
  )
  indicators = _cache_layers({f"{k}_{rad_km}km": i for k, i in indicators.items()}, KERNEL_SCALE)
  viz_layers = _cache_layers({f"{k}_{rad_km}km": i for k, i in viz_layers.items()}, KERNEL_SCALE)

  return {
    **input_layers,
    **derived_layers,
    **kernel_layers,
    **indicators,
    **viz_layers,
  }


def preview_all_landscape_indicators(kernel_radius=KERNEL_RADIUS):
  """"""
  input_layers = get_all_input_layers()
  derived_layers = compute_derived_layers(input_layers)
  kernel_layers = compute_kernels(derived_layers, kernel_radius)
  indicators = compute_indicators(kernel_layers)
  viz_layers = compute_viz_layers(indicators)

  return {
    **input_layers,
    **derived_layers,
    **kernel_layers,
    **indicators,
    **viz_layers,
  }


# In[9]:


### Preview all layers
# You can preview all the calculations and intermediate layers below

# The map will render quickly(ish) but is inaccurate because many input datasets don't
# have nearest-neighbor resampling or other masking issues in their pyramids/overviews.


def visualize(layers):
  map = gmap.Map()
  map.add_basemap("SATELLITE")
  for k, layer in layers.items():
    vmax = 1
    if "kernel" in k:
      vmax *= 0.2
    if "carbon" in k:
      vmax *= 50
    if "human_lu" in k:
      vmax *= 0.05
    map.add_layer(
      layer.updateMask(layer.divide(vmax)),
      {"min": 0, "max": vmax * 2, "palette": ["black", "red", "orange", "yellow"]},
      k,
      False,
    )
  return map


visualize(preview_all_landscape_indicators())

# In[ ]:


# Compute and visualize final layers

# When caching the layers, we'll call reduceResolution to accurately
# resample the input layers to lower resolutions.

# NOTE: Computation needs to be done in three steps:
#       First, we resample to 100m.
#       Then we do the overlay analysis to identify deforestation, cropland expansion, etc.
#       Then we'll use these outputs to comput the final kerneled layers at 1000m.

# NOTE: Sometimes this will raise an exception due to there not being a default
#       projection set, even though we set one manually when it needs to be.
#       This seems to be related to an internal GEE caching issue. If this happens,
#       just run the cell again

layers = compute_all_landscape_indicators()

visualize(layers)


# In[ ]:


# Compute layers with varying radii if necessary


def cache_radii_layers():
  for r in RADII:
    compute_all_landscape_indicators(kernel_radius=r)


# cache_radii_layers()


# In[ ]:


# make assets public
eeUtil.setAcl(WORKING_FOLDER, "public", recursive=True)
