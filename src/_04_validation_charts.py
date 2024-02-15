#!/usr/bin/env python
# coding: utf-8

# # Validation charts and statistics
#
# Generates charts and statistics to validate the data in the country comparison dataset.

# In[16]:


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# In[17]:


all_stats = pd.read_csv("../data/processed/all_country_stats.csv")
all_stats.head()

# In[18]:


# Comparison of deforestation definitions
N = 20
SORT_COL = "fao_forestloss"
title = "Forest loss area (km2) / yr"
cols = {
  "country_name": "Country",
  "treeloss": "Tree loss (Hansen)",
  "nat_treeloss": "Natural forest loss (Hansen, Mazur)",
  "deforest": 'Deforest (treeloss excl. nonnat, fire, distrub) "LG Method"',
  "fao_forestloss": "FAO forest conversion",
}
all_stats.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.bar(x="Country", title=title)

all_stats["fao_ratio"] = np.log(all_stats["deforest"] / all_stats["fao_forestloss"])
all_stats[all_stats["fao_forestloss"] > 0].sort_values("fao_ratio", ascending=False).head(N).filter(
  cols.keys()
).rename(columns=cols).plot.bar(x="Country", title=title)
all_stats["fao_ratio"] = np.log(all_stats["deforest"] / all_stats["fao_forestloss"])
all_stats[all_stats["fao_forestloss"] > 0].sort_values("fao_ratio").head(N).filter(
  cols.keys()
).rename(columns=cols).plot.bar(x="Country", title=title)

REGION_LIST = [
  "Southern Asia",  # 'Southern Europe', 'Northern Africa',
  "Polynesia",
  "Sub-Saharan Africa",
  "Latin America and the Caribbean",  # 'Western Asia',
  # 'Australia and New Zealand', 'Western Europe', 'Eastern Europe',
  # 'Northern America',
  "South-eastern Asia",
  "Eastern Asia",
  # 'Northern Europe',
  "Melanesia",
  "Micronesia",  # 'Central Asia'
]
all_stats["tropical"] = all_stats.region.isin(REGION_LIST)
all_stats_trop = all_stats[all_stats["tropical"]]
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.bar(x="Country", title=title)

# In[19]:


# Comparison of sLUC
c = ["#0c1063", "#0d6cac"]

fig, axs = plt.subplots(2, 1, sharex=True)
title = "Deforestation total (km2/yr)"
cols = {
  "country_name": "Country",
  "deforest": "Deforest (LG method)",
  "fao_forestloss": "FAO forest conversion",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.bar(x="Country", title=title, color=c, ax=axs[0])

title = "sLUC deforestation attributable to cropland"
cols = {
  "country_name": "Country",
  "sluc_deforest_cropland": "Deforest attr. to cropland (LG, Potapov)",
  "sluc_fao_cropland": "FAO forest conversion attr. to cropland",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.bar(x="Country", title=title, color=c, ax=axs[1])
fig.show()

# In[20]:


# Comparison of sLUC
c = ["#0c1063", "#0d6cac"]

fig, axs = plt.subplots(2, 1, sharex=True)
title = "Deforestation GHG total (tCO2e/yr)"
cols = {
  "country_name": "Country",
  "deforest_carbon": "Deforest GHG (LG method)",
  "fao_emissions": "FAO forest conversion GHG",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.bar(x="Country", title=title, color=c, ax=axs[0])

title = "sLUC deforestation GHG attributable to cropland"
cols = {
  "country_name": "Country",
  "sluc_deforest_carbon_cropland": "Deforest attr. to cropland (LG, Potapov)",
  "sluc_emissions_fao_cropland": "FAO emissions attr. to cropland",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.bar(x="Country", title=title, color=c, ax=axs[1])
fig.show()

# In[21]:


# dLUC
c = ["#fe6598", "#333333"]
N = 10

fig, axs = plt.subplots(1, 2, sharey=True)
fig.suptitle("sLUC vs dLUC treeloss (km2/yr)")
title = "attr. to all human land (km2)"
cols = {
  "country_name": "Country",
  "treeloss": "sLUC",
  "dluc_nonnat_treeloss": "dLUC",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.barh(x="Country", title=title, color=c, ax=axs[1], width=0.8).invert_yaxis()
title = "attr. to cropland (km2)"
cols = {
  "country_name": "Country",
  "sluc_treeloss_cropland": "sLUC",
  "dluc_cropland_treeloss": "dLUC",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.barh(x="Country", title=title, color=c, ax=axs[0], width=0.8).invert_yaxis()
fig.show()

fig, axs = plt.subplots(1, 2, sharey=True)
fig.suptitle("sLUC vs dLUC natural forest loss (km2/yr)")
title = "attr. to all human land (km2)"
cols = {
  "country_name": "Country",
  "nat_treeloss": "sLUC",
  "dluc_nonnat_nattreeloss": "dLUC",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.barh(x="Country", title=title, color=c, ax=axs[1], width=0.8).invert_yaxis()
title = "attr. to cropland (km2)"
cols = {
  "country_name": "Country",
  "sluc_nat_treeloss_cropland": "sLUC",
  "dluc_cropland_nattreeloss": "dLUC",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.barh(x="Country", title=title, color=c, ax=axs[0], width=0.8).invert_yaxis()
fig.show()

fig, axs = plt.subplots(1, 2, sharey=True)
fig.suptitle("sLUC vs dLUC deforestation (km2/yr)")
title = "Attr. to all human land"
cols = {
  "country_name": "Country",
  "deforest": "sLUC",
  "dluc_nonnat_deforest": "dLUC",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.barh(x="Country", title=title, color=c, ax=axs[1], width=0.8).invert_yaxis()
title = "Attr. to cropland"
cols = {
  "country_name": "Country",
  "sluc_deforest_cropland": "sLUC",
  "dluc_cropland_deforest": "dLUC",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.barh(x="Country", title=title, color=c, ax=axs[0], width=0.8).invert_yaxis()
fig.show()

# In[22]:


# spatial sLUC

c = ["#fe6598", "#0C1063", "#1C218A", "#2C33B1", "#3C44D8", "#4C55FF", "#333333"]

title = "sLUC vs dLUC for cropland (km2/yr)"
cols = {
  "country_name": "Country",
  "sluc_deforest_cropland": "sLUC countrywide",
  "cropland_deforest_by_human_lu_50km_1000m": "sLUC 50km kernel",
  "cropland_deforest_by_human_lu_25km_1000m": "sLUC 25km kernel",
  "cropland_deforest_by_human_lu_10km_1000m": "sLUC 10km kernel",
  "cropland_deforest_by_human_lu_5km_1000m": "sLUC 5km kernel",
  "cropland_deforest_by_human_lu_0km_1000m": "sLUC 0.5km kernel",
  "dluc_cropland_deforest": "dLUC",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.barh(x="Country", color=c, title=title, width=0.8).invert_yaxis()

title = "sLUC vs dLUC for all non-natural areas (km2/yr)"
cols = {
  "country_name": "Country",
  "deforest": "sLUC countrywide",
  "nonnat_deforest_by_human_lu_50km_1000m": "sLUC 50km kernel",
  "nonnat_deforest_by_human_lu_25km_1000m": "sLUC 25km kernel",
  "nonnat_deforest_by_human_lu_10km_1000m": "sLUC 10km kernel",
  "nonnat_deforest_by_human_lu_5km_1000m": "sLUC 5km kernel",
  "nonnat_deforest_by_human_lu_0km_1000m": "sLUC 0.5km kernel",
  "dluc_nonnat_deforest": "dLUC",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.barh(x="Country", color=c, title=title, width=0.8).invert_yaxis()

# In[23]:


# spatial sLUC

c = ["#fe6598", "#0C1063", "#1C218A", "#2C33B1", "#3C44D8", "#4C55FF", "#333333"]

title = "sLUC vs dLUC for cropland (tCO2/yr)"
cols = {
  "country_name": "Country",
  "sluc_deforest_carbon_cropland": "sLUC countrywide",
  "cropland_deforest_carbon_by_human_lu_50km_1000m": "sLUC 50km kernel",
  "cropland_deforest_carbon_by_human_lu_25km_1000m": "sLUC 25km kernel",
  "cropland_deforest_carbon_by_human_lu_10km_1000m": "sLUC 10km kernel",
  "cropland_deforest_carbon_by_human_lu_5km_1000m": "sLUC 5km kernel",
  "cropland_deforest_carbon_by_human_lu_0km_1000m": "sLUC 0.5km kernel",
  "dluc_cropland_deforest_carbon": "dLUC",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.barh(x="Country", color=c, title=title, width=0.8).invert_yaxis()

title = "sLUC vs dLUC for all non-natural areas (tCO2/yr)"
cols = {
  "country_name": "Country",
  "deforest_carbon": "sLUC countrywide",
  "nonnat_deforest_carbon_by_human_lu_50km_1000m": "sLUC 50km kernel",
  "nonnat_deforest_carbon_by_human_lu_25km_1000m": "sLUC 25km kernel",
  "nonnat_deforest_carbon_by_human_lu_10km_1000m": "sLUC 10km kernel",
  "nonnat_deforest_carbon_by_human_lu_5km_1000m": "sLUC 5km kernel",
  "nonnat_deforest_carbon_by_human_lu_0km_1000m": "sLUC 0.5km kernel",
  "dluc_nonnat_deforest_carbon": "dLUC",
}
all_stats_trop.sort_values(SORT_COL, ascending=False).head(N).filter(cols.keys()).rename(
  columns=cols
).plot.barh(x="Country", color=c, title=title, width=0.8).invert_yaxis()

# In[24]:


print("Cropland")
print("Tropical countries")
denom = "sluc_fao_cropland"
denom = "sluc_deforest_cropland"

all_stats_f = all_stats[all_stats[denom] > 0]
all_stats_trop_f = all_stats_trop[all_stats_trop[denom] > 0]

cols = {
  "country_name": "Country",
  "sluc_deforest_cropland": "sLUC countrywide",
  "cropland_deforest_by_human_lu_50km_1000m": "sLUC 50km kernel",
  "cropland_deforest_by_human_lu_25km_1000m": "sLUC 25km kernel",
  "cropland_deforest_by_human_lu_10km_1000m": "sLUC 10km kernel",
  "cropland_deforest_by_human_lu_5km_1000m": "sLUC 5km kernel",
  "cropland_deforest_by_human_lu_0km_1000m": "sLUC 0.5km kernel",
  "dluc_cropland_deforest": "dLUC",
}

for col, lyr in cols.items():
  if col != "country_name":
    print(lyr, round(all_stats_trop_f[col].sum() / all_stats_trop_f[denom].sum(), 3))
print("All countries")
for col, lyr in cols.items():
  if col != "country_name":
    print(lyr, round(all_stats_f[col].sum() / all_stats_f[denom].sum(), 3))

print("All HLU")
print("Tropical countries")
denom = "fao_forestloss"
denom = "deforest"

all_stats_f = all_stats[all_stats[denom] > 0]
all_stats_trop_f = all_stats_trop[all_stats_trop[denom] > 0]

cols = {
  "country_name": "Country",
  "deforest": "sLUC countrywide",
  "nonnat_deforest_by_human_lu_50km_1000m": "sLUC 50km kernel",
  "nonnat_deforest_by_human_lu_25km_1000m": "sLUC 25km kernel",
  "nonnat_deforest_by_human_lu_10km_1000m": "sLUC 10km kernel",
  "nonnat_deforest_by_human_lu_5km_1000m": "sLUC 5km kernel",
  "nonnat_deforest_by_human_lu_0km_1000m": "sLUC 0.5km kernel",
  "dluc_nonnat_deforest": "dLUC",
}

for col, lyr in cols.items():
  if col != "country_name":
    print(lyr, round(all_stats_trop_f[col].sum() / all_stats_trop_f[denom].sum(), 3))
print("All countries")
for col, lyr in cols.items():
  if col != "country_name":
    print(lyr, round(all_stats_f[col].sum() / all_stats_f[denom].sum(), 3))

# In[25]:


# ALL
# all_stats_f = all_stats_f[~all_stats_f["country_name"].isin(('Canada', 'US'))]
print("ALL COUNTRIES", all_stats_f["fao_forestloss"].count())
print("DEFORESTATION")
print("r2", all_stats_f["fao_forestloss"].corr(all_stats_f["deforest"]) ** 2)
print(
  "r2 nonnat",
  all_stats_f["fao_forestloss"].corr(all_stats_f["nonnat_deforest_by_human_lu_50km_1000m"]) ** 2,
)
print(
  "r2 cropland",
  all_stats_f["sluc_fao_cropland"].corr(all_stats_f["cropland_deforest_by_human_lu_50km_1000m"])
  ** 2,
)
print("this/fao", all_stats_f["deforest"].sum() / all_stats_f["fao_forestloss"].sum())
print(
  "this/fao nonnat",
  all_stats_f["nonnat_deforest_by_human_lu_50km_1000m"].sum()
  / all_stats_f["fao_forestloss"].sum(),
)
print(
  "this/fao cropland",
  all_stats_f["cropland_deforest_by_human_lu_50km_1000m"].sum()
  / all_stats_f["sluc_fao_cropland"].sum(),
)
print("DEFORESTATION TOTALS km2/yr")
print("this", all_stats_f["deforest"].sum())
print("fao", all_stats_f["fao_forestloss"].sum())
print("this cropland", all_stats_f["cropland_deforest_by_human_lu_50km_1000m"].sum())
print("fao cropland", all_stats_f["sluc_fao_cropland"].sum())

print("\n\nCO2")
print("r2", all_stats_f["fao_emissions"].corr(all_stats_f["deforest_carbon"]) ** 2)
print(
  "r2 nonnat",
  all_stats_f["fao_emissions"].corr(all_stats_f["nonnat_deforest_carbon_by_human_lu_50km_1000m"])
  ** 2,
)
print(
  "r2 cropland",
  all_stats_f["sluc_emissions_fao_cropland"].corr(
    all_stats_f["cropland_deforest_carbon_by_human_lu_50km_1000m"]
  )
  ** 2,
)
print("this/fao", all_stats_f["deforest_carbon"].sum() / all_stats_f["fao_emissions"].sum())
print(
  "this/fao nonnat",
  all_stats_f["nonnat_deforest_carbon_by_human_lu_50km_1000m"].sum()
  / all_stats_f["fao_emissions"].sum(),
)
print(
  "this/fao cropland",
  all_stats_f["cropland_deforest_carbon_by_human_lu_50km_1000m"].sum()
  / all_stats_f["sluc_emissions_fao_cropland"].sum(),
)

print("CO2 TOTALS tCO2/yr")
print("this", all_stats_f["deforest_carbon"].sum())
print("fao", all_stats_f["fao_emissions"].sum())
print("this cropland", all_stats_f["cropland_deforest_carbon_by_human_lu_50km_1000m"].sum())
print("fao cropland", all_stats_f["sluc_emissions_fao_cropland"].sum())

# In[26]:


lbl_countries = ["Canada", "US", "Botswana", "Pakistan"]

all_vars = [
  {
    "deforest": "This study",
    "fao_forestloss": "FAO",
  },
  {
    "cropland_deforest_by_human_lu_50km_1000m": "This study",
    "sluc_fao_cropland": "FAO",
  },
  {
    "deforest_carbon": "This study",
    "fao_emissions": "FAO",
  },
  {
    "cropland_deforest_carbon_by_human_lu_50km_1000m": "This study",
    "sluc_emissions_fao_cropland": "FAO",
  },
]

lims = [[1, 1e5], [1, 1e5], [1e3, 1e10], [1e3, 1e9]]

fig = plt.figure(figsize=(6, 6.5), dpi=150, layout="constrained")
subfigs = fig.subfigures(2, 1)
_axs0 = subfigs[0].subplots(1, 2)
subfigs[0].suptitle("Area (km2/yr)")
_axs1 = subfigs[1].subplots(1, 2)
subfigs[1].suptitle("GHG Emissions (tCO2/yr)")
axs = np.concatenate([_axs0, _axs1])
for i in range(4):
  ax = axs[i]
  vars = all_vars[i]
  lim = lims[i]
  ax.plot(lim, lim, color="k", linestyle="--", alpha=0.5)

  x = list(vars.keys())[0]
  y = list(vars.keys())[1]
  color = "tab:blue"
  if i > 1:
    color = "tab:orange"

  # all_stats_f.filter(vars.keys()).rename(columns=vars).
  ax.scatter(x=all_stats_f[x], y=all_stats_f[y], color=color, s=15)
  ax.set_xlim(lim)
  ax.set_ylim(lim)
  ax.set_yscale("log")
  ax.set_xscale("log")
  if i == 0:
    ax.set_title("Total")
    ax.set_ylabel("FAO")
  elif i == 1:
    ax.set_title("Cropland")
  elif i == 2:
    ax.set_ylabel("FAO")
    ax.set_xlabel("This study")
  elif i == 3:
    ax.set_xlabel("This study")

  for c in lbl_countries:
    i = all_stats_f[all_stats_f["country_name"] == c].index[0]
    txt = all_stats_f["ISO3166-1-Alpha-3"][i]
    _x = all_stats_f[x][i]
    _y = all_stats_f[y][i]
    ax.annotate(txt, (_x, _y), fontsize=10)

plt.show()

# In[27]:


# TROPICAL
print("TROPICAL COUNTRIES", all_stats_trop_f["fao_forestloss"].count())
print("DEFORESTATION")
print("r2", all_stats_trop_f["fao_forestloss"].corr(all_stats_trop_f["deforest"]) ** 2)
print(
  "r2 nonnat",
  all_stats_trop_f["fao_forestloss"].corr(
    all_stats_trop_f["nonnat_deforest_by_human_lu_50km_1000m"]
  )
  ** 2,
)
print(
  "r2 cropland",
  all_stats_trop_f["sluc_fao_cropland"].corr(
    all_stats_trop_f["cropland_deforest_by_human_lu_50km_1000m"]
  )
  ** 2,
)
print("this/fao", all_stats_trop_f["deforest"].sum() / all_stats_trop_f["fao_forestloss"].sum())
print(
  "this/fao nonnat",
  all_stats_trop_f["nonnat_deforest_by_human_lu_50km_1000m"].sum()
  / all_stats_trop_f["fao_forestloss"].sum(),
)
print(
  "this/fao cropland",
  all_stats_trop_f["cropland_deforest_by_human_lu_50km_1000m"].sum()
  / all_stats_trop_f["sluc_fao_cropland"].sum(),
)
print("DEFORESTATION TOTALS km2/yr")
print("this", all_stats_trop_f["deforest"].sum())
print("fao", all_stats_trop_f["fao_forestloss"].sum())
print("this cropland", all_stats_trop_f["cropland_deforest_by_human_lu_50km_1000m"].sum())
print("fao cropland", all_stats_trop_f["sluc_fao_cropland"].sum())

print("\n\nCO2")
print("r2", all_stats_trop_f["fao_emissions"].corr(all_stats_trop_f["deforest_carbon"]) ** 2)
print(
  "r2 nonnat",
  all_stats_trop_f["fao_emissions"].corr(
    all_stats_trop_f["nonnat_deforest_carbon_by_human_lu_50km_1000m"]
  )
  ** 2,
)
print(
  "r2 cropland",
  all_stats_trop_f["sluc_emissions_fao_cropland"].corr(
    all_stats_trop_f["cropland_deforest_carbon_by_human_lu_50km_1000m"]
  )
  ** 2,
)
print(
  "this/fao", all_stats_trop_f["deforest_carbon"].sum() / all_stats_trop_f["fao_emissions"].sum()
)
print(
  "this/fao nonnat",
  all_stats_trop_f["nonnat_deforest_carbon_by_human_lu_50km_1000m"].sum()
  / all_stats_trop_f["fao_emissions"].sum(),
)
print(
  "this/fao cropland",
  all_stats_trop_f["cropland_deforest_carbon_by_human_lu_50km_1000m"].sum()
  / all_stats_trop_f["sluc_emissions_fao_cropland"].sum(),
)

print("CO2 TOTALS tCO2/yr")
print("this", all_stats_trop_f["deforest_carbon"].sum())
print("fao", all_stats_trop_f["fao_emissions"].sum())
print("this cropland", all_stats_trop_f["cropland_deforest_carbon_by_human_lu_50km_1000m"].sum())
print("fao cropland", all_stats_trop_f["sluc_emissions_fao_cropland"].sum())
