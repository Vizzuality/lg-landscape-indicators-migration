#!/usr/bin/env python
# coding: utf-8

# # Extract country level statistics from global layers
#
# Computes statistics for each indicator and exports the results to Google Cloud Storage.
#

# In[1]:


# Imports
import os

import ee
import eeUtil
from _01_compute_landscape_indicators import visualize

# Initialize Earth Engine
PROJECT = os.getenv("GEE_PROJECT")
GEE_JSON = os.getenv("GEE_JSON")

assert PROJECT is not None, "Please set GEE_PROJECT environment variable"
assert (
  GEE_JSON is not None
), "Please set GEE_JSON environment variable with service account credentials"

eeUtil.init()

# In[4]:


PROJECTION = "EPSG:4326"
SRC_FOLDER = "projects/ee-fgassert/assets/landscape_indicators_source"
WORKING_FOLDER = f"projects/{PROJECT}/assets/landscape_indicators_20231023"
WORLD_GEOM = ee.Geometry.Polygon(
  coords=[[[-180, -85], [-180, 85], [180, 85], [180, -85], [-180, -85]]],
  proj=PROJECTION,
  geodesic=False,
)
EXPORT_BUCKET = os.getenv("GEE_BUCKET")
EXPORT_PREFIX = "landscape_indicators_20231023"

TARGET_YEAR = 2022

ANALYSIS_SCALE = 100
KERNEL_SCALE = 1000

RADII = [500, 5000, 10000, 25000, 50000]

GAUL_FC = "FAO/GAUL/2015/level0"

COUNTRY_ANALYSIS_LAYERS = [
  "cropland_deforest_by_human_lu",
  "cropland_tree_loss_by_human_lu",
  "cropland_deforest_carbon_by_human_lu",
  "cropland_natural_crop_net_conversion_by_human_lu",
  "nonnat_deforest_by_human_lu",
  "nonnat_tree_loss_by_human_lu",
  "nonnat_deforest_carbon_by_human_lu",
  "nonnat_natural_crop_net_conversion_by_human_lu",
]


# In[5]:


def stack_layers_to_bands(layers):
  """Stacks a list of layers into a single image"""
  im = ee.Image(1).rename("const")
  for k, v in layers.items():
    im = im.addBands(v.rename(k))
  return im


def compute_country_stats_by_radii():
  gaul = ee.FeatureCollection(GAUL_FC)
  layers = {}
  for r in RADII:
    for k in COUNTRY_ANALYSIS_LAYERS:
      rad_km = int(r / 1000)
      key = f"{k}_{rad_km}km_{KERNEL_SCALE}m"
      layers[key] = ee.Image(f"{WORKING_FOLDER}/{key}")

  im = stack_layers_to_bands(layers)
  im = im.multiply(ee.Image.pixelArea()).divide(1e6)  # convert to km2
  table = im.reduceRegions(gaul, "sum", KERNEL_SCALE).select([".*"], None, False)
  asset_id = f"{EXPORT_PREFIX}/sluc_kernel_country_stats_{KERNEL_SCALE}m"
  desc = eeUtil.eeutil._getExportDescription(asset_id)
  if desc in [t["description"] for t in eeUtil.getTasks(True)]:
    print(f"{asset_id} in progress")
  else:
    ee.batch.Export.table.toCloudStorage(table, desc, EXPORT_BUCKET, asset_id).start()


def and_assume_disjoint(*args):
  im = ee.Image([*args]).reduce("sum")
  return im.subtract(im.mask()).max(0)


def get_dluc_layers():
  mask = ee.Image(f"{WORKING_FOLDER}/available_area_mask_{ANALYSIS_SCALE}m")
  cropland = ee.Image(f"{WORKING_FOLDER}/cropland_2019_{ANALYSIS_SCALE}m").updateMask(mask)
  nonnat = ee.Image(f"{WORKING_FOLDER}/nonnatural_excl_builtwater_{ANALYSIS_SCALE}m").updateMask(
    mask
  )
  treeloss = (
    ee.Image(f"{WORKING_FOLDER}/tree_loss_to{TARGET_YEAR}_{ANALYSIS_SCALE}m")
    .divide(20)
    .updateMask(mask)
  )
  deforest = ee.Image(f"{WORKING_FOLDER}/deforest_{ANALYSIS_SCALE}m").divide(20).updateMask(mask)
  deforest_carbon = (
    ee.Image(f"{WORKING_FOLDER}/deforest_carbon_{ANALYSIS_SCALE}m").divide(20).updateMask(mask)
  )

  nat_treeloss = treeloss.multiply(nonnat.subtract(1).multiply(-1))
  return {
    "cropland": cropland,
    "nonnat": nonnat,
    "deforest": deforest,
    "deforest_carbon": deforest_carbon,
    "treeloss": treeloss,
    "nat_treeloss": nat_treeloss,
    "dluc_cropland_treeloss": cropland.multiply(treeloss),
    "dluc_nonnat_treeloss": nonnat.multiply(treeloss),
    "dluc_cropland_nattreeloss": cropland.multiply(nat_treeloss),
    "dluc_nonnat_nattreeloss": nonnat.multiply(nat_treeloss),
    "dluc_cropland_deforest": cropland.multiply(deforest),
    "dluc_nonnat_deforest": nonnat.multiply(deforest),
    "dluc_cropland_deforest_carbon": cropland.multiply(deforest_carbon),
    "dluc_nonnat_deforest_carbon": nonnat.multiply(deforest_carbon),
  }


def compute_country_stats_dluc():
  """Compute the DLUC indicator for each country"""
  gaul = ee.FeatureCollection(GAUL_FC)
  dluc_layers = get_dluc_layers()

  im = stack_layers_to_bands(dluc_layers)
  im = im.multiply(ee.Image.pixelArea()).divide(1e6)  # convert to km2
  table = im.reduceRegions(gaul, "sum", ANALYSIS_SCALE, tileScale=4).select([".*"], None, False)
  asset_id = f"{EXPORT_PREFIX}/country_stats_{ANALYSIS_SCALE}m"
  desc = eeUtil.eeutil._getExportDescription(asset_id)
  if desc in [t["description"] for t in eeUtil.getTasks(True)]:
    print(f"{asset_id} in progress")
  else:
    ee.batch.Export.table.toCloudStorage(table, desc, EXPORT_BUCKET, asset_id).start()


# In[6]:


visualize(get_dluc_layers())

# In[7]:


# run country-level analysis

compute_country_stats_by_radii()
compute_country_stats_dluc()

# In[ ]:


[f"{t['state']} {t['description']}" for t in eeUtil.getTasks()]
# [ee.data.cancelTask(t['id']) for t in eeUtil.getTasks(True)]
