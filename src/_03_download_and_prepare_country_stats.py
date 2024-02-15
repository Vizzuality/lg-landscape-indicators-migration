#!/usr/bin/env python
# coding: utf-8

# # 3. Download and prepare country stats
#
# Downloads statistics from google cloud storage and combines them with FAO STAT data
# on emissions land use and forests. Creates a single csv file with all data.
#

# In[5]:


import os

import pandas as pd

# Set up credentials
GEE_JSON = os.getenv("GEE_JSON")
PROJECT = os.getenv("PROJECT")
EXPORT_BUCKET = os.getenv("GEE_BUCKET")
EXPORT_PREFIX = "landscape_indicators_20231023"

DATA_FOLDER = "../data/raw"

# In[6]:


drop_columns = [
  "system:index",
  "ADM0_NAME",
  "DISP_AREA",
  "EXP0_YEAR",
  "STATUS",
  "STR0_YEAR",
  "Shape_Area",
  "Shape_Leng",
  ".geo",
]

country_stats = pd.read_csv(f"gs://{EXPORT_BUCKET}/{EXPORT_PREFIX}/country_stats_100m.csv")
country_kernel_stats = pd.read_csv(
  f"gs://{EXPORT_BUCKET}/{EXPORT_PREFIX}/sluc_kernel_country_stats_1000m.csv"
)
country_stats.set_index("ADM0_CODE", inplace=True)
country_stats.drop(columns=drop_columns, inplace=True)
country_kernel_stats.set_index("ADM0_CODE", inplace=True)
country_kernel_stats.drop(columns=drop_columns, inplace=True)
country_stats = country_stats.join(country_kernel_stats, rsuffix="2")

print(country_stats.columns)
country_stats.head()

# In[7]:


# Load FAO data

country_df = pd.read_csv(f"{DATA_FOLDER}/country-codes.csv")
fao_forest_data = pd.read_csv(
  f"{DATA_FOLDER}/Emissions_Land_Use_Forests_E_All_Data_(Normalized).csv", encoding="latin-1"
)
fao_landuse_data = pd.read_csv(
  f"{DATA_FOLDER}/Inputs_LandUse_E_All_Data_(Normalized).csv", encoding="latin-1"
)

fao_forest_data = fao_forest_data[
  ((fao_forest_data["Year"] <= 2022) & (fao_forest_data["Year"] > 2002))
]
fao_forest_data["M49"] = fao_forest_data["Area Code (M49)"].apply(lambda x: x[1:]).astype(float)
fao_landuse_data = fao_landuse_data[
  ((fao_landuse_data["Year"] <= 2022) & (fao_landuse_data["Year"] > 2002))
]
fao_landuse_data["M49"] = fao_landuse_data["Area Code (M49)"].apply(lambda x: x[1:]).astype(float)

fao_area = (
  fao_forest_data[
    (
      (fao_forest_data["Source"] == "FAO TIER 1")
      & (fao_forest_data["Item"] == "Net Forest conversion")
      & (fao_forest_data["Element"] == "Area")
    )
  ]
  .groupby("M49")["Value"]
  .mean()
  * 10
)  # convert to km2
fao_area.name = "fao_forestloss"

fao_emissions = (
  fao_forest_data[
    (
      (fao_forest_data["Source"] == "FAO TIER 1")
      & (fao_forest_data["Item"] == "Net Forest conversion")
      & (fao_forest_data["Element"] == "Net emissions/removals (CO2) (Forest land)")
    )
  ]
  .groupby("M49")["Value"]
  .mean()
  * 1000
  # kT to tons
)
fao_emissions.name = "fao_emissions"

fao_cropland = (
  fao_landuse_data[((fao_landuse_data["Item"] == "Cropland"))].groupby("M49")["Value"].mean() * 10
)  # convert to km2
fao_cropland.name = "fao_cropland"

fao_agland = (
  fao_landuse_data[((fao_landuse_data["Item"] == "Agricultural land"))]
  .groupby("M49")["Value"]
  .mean()
  * 10
)  # convert to km2
fao_agland.name = "fao_agland"

all_fao = pd.concat([fao_area, fao_emissions, fao_cropland, fao_agland], axis=1)
all_fao.head()

# In[8]:


c = country_df.drop(211)
c.GAUL = c.GAUL.astype(float)
c = c.filter(["GAUL", "CLDR display name", "Sub-region Name", "M49", "ISO3166-1-Alpha-3"])
c.rename(columns={"CLDR display name": "country_name", "Sub-region Name": "region"}, inplace=True)

all_stats = (
  country_stats.groupby("ADM0_CODE").sum().join(c.set_index("GAUL")).join(all_fao, on="M49")
)

all_stats["deforest_carbon"] *= 100 * 3.66  # ha->km2 and tC -> tCO2e
all_stats["nonnat_deforest_carbon_by_human_lu_0km_1000m"] *= 100  # ha->km2
all_stats["nonnat_deforest_carbon_by_human_lu_5km_1000m"] *= 100  # ha->km2
all_stats["nonnat_deforest_carbon_by_human_lu_10km_1000m"] *= 100  # ha->km2
all_stats["nonnat_deforest_carbon_by_human_lu_25km_1000m"] *= 100  # ha->km2
all_stats["nonnat_deforest_carbon_by_human_lu_50km_1000m"] *= 100  # ha->km2
all_stats["cropland_deforest_carbon_by_human_lu_0km_1000m"] *= 100  # ha->km2
all_stats["cropland_deforest_carbon_by_human_lu_5km_1000m"] *= 100  # ha->km2
all_stats["cropland_deforest_carbon_by_human_lu_10km_1000m"] *= 100  # ha->km2
all_stats["cropland_deforest_carbon_by_human_lu_25km_1000m"] *= 100  # ha->km2
all_stats["cropland_deforest_carbon_by_human_lu_50km_1000m"] *= 100  # ha->km2

# annualize
all_stats["sluc_treeloss_nonnat_ha"] = all_stats["treeloss"] / all_stats["nonnat"]
all_stats["sluc_treeloss_cropland_ha"] = all_stats["treeloss"] / all_stats["cropland"]
all_stats["sluc_treeloss_cropland"] = all_stats["cropland"] * all_stats["sluc_treeloss_nonnat_ha"]
all_stats["sluc_nat_treeloss_nonnat_ha"] = all_stats["nat_treeloss"] / all_stats["nonnat"]
all_stats["sluc_nat_treeloss_cropland_ha"] = all_stats["nat_treeloss"] / all_stats["cropland"]
all_stats["sluc_nat_treeloss_cropland"] = (
  all_stats["cropland"] * all_stats["sluc_nat_treeloss_nonnat_ha"]
)
all_stats["sluc_deforest_nonnat_ha"] = all_stats["deforest"] / all_stats["nonnat"]
all_stats["sluc_deforest_cropland_ha"] = all_stats["deforest"] / all_stats["cropland"]
all_stats["sluc_deforest_cropland"] = all_stats["cropland"] * all_stats["sluc_deforest_nonnat_ha"]
all_stats["sluc_deforest_carbon_nonnat_ha"] = all_stats["deforest_carbon"] / all_stats["nonnat"]
all_stats["sluc_deforest_carbon_cropland_ha"] = all_stats["deforest_carbon"] / all_stats["cropland"]
all_stats["sluc_deforest_carbon_cropland"] = (
  all_stats["cropland"] * all_stats["sluc_deforest_carbon_nonnat_ha"]
)
all_stats["sluc_fao_agland_ha"] = all_stats["fao_forestloss"] / all_stats["fao_agland"]
all_stats["sluc_fao_cropland_ha"] = all_stats["fao_forestloss"] / all_stats["fao_cropland"]
all_stats["sluc_fao_cropland"] = all_stats["fao_cropland"] * all_stats["sluc_fao_agland_ha"]
all_stats["sluc_emissions_fao_agland_ha"] = all_stats["fao_emissions"] / all_stats["fao_agland"]
all_stats["sluc_emissions_fao_cropland_ha"] = all_stats["fao_emissions"] / all_stats["fao_cropland"]
all_stats["sluc_emissions_fao_cropland"] = (
  all_stats["fao_cropland"] * all_stats["sluc_emissions_fao_agland_ha"]
)

print(all_stats.columns)
all_stats.to_csv("../data/processed/all_country_stats.csv")
all_stats.head()
