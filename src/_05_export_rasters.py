# Imports
import os

import eeUtil

# Initialize Earth Engine
PROJECT = os.getenv("GEE_PROJECT")
GEE_JSON = os.getenv("GEE_JSON")
eeUtil.init()

WORKING_FOLDER = f"projects/{PROJECT}/assets/landscape_indicators_20230821"
EXPORT_BUCKET = os.getenv("GEE_BUCKET")
EXPORT_PREFIX = "landscape_indicators_20231023"

export_layers = [
  "deforest_100m",
  "deforest_carbon_100m",
  "natcrop_expansion_100m",
  "natcrop_reduction_100m",
  "natcrop_flii_loss_100m",
  "natcrop_bii_loss_100m",
]


def export_final_layers():
  return eeUtil.export(
    [
      a
      for a in eeUtil.ls(WORKING_FOLDER, True)
      if a.endswith("by_human_lu_50km_1000m") or os.path.basename(a) in export_layers
    ],
    bucket=EXPORT_BUCKET,
    prefix=EXPORT_PREFIX,
    overwrite=True,
  )


export_final_layers()
