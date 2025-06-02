# replicate_figures.py

import os
import pandas as pd
import plotly.express as px
import math
import json
from urllib.request import urlopen

### Set directories ###
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')
IMG_DIR = os.path.join(BASE_DIR, 'images')

os.makedirs(IMG_DIR, exist_ok=True)

### Load county data (geojson) ###
with urlopen('https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json') as response:
    counties = json.load(response)

### Load county-level results ###
df = pd.read_csv(os.path.join(DATA_DIR, "county_results.csv"))
df.fips = df.fips.apply(lambda x: "{:05d}".format(x))
df['LogFedCattle'] = df['FedCattle'].apply(lambda x: math.log10(x))

### Function to generate choropleth maps ###
def generate_maps(col, label, range_, ticks, png_size):
    scale = px.colors.sequential.OrRd
    fig = px.choropleth(df, geojson=counties, locations="fips", color=df[col],
                        color_continuous_scale=scale, range_color=range_)

    fig.update_traces(marker_line_width=0)
    fig.update_geos(scope='usa', fitbounds="locations", landcolor="LightGray")
    fig.update_layout(
        coloraxis_colorbar=dict(
            title=label,
            orientation="h",
            yanchor="top",
            y=0.00,
            title_side="top",
            len=0.7,
            thickness=10,
            tickvals=ticks[0],
            ticktext=ticks[1],
            tickfont=dict(size=16)
        ),
        margin={"l":10,"r":10,"t":30,"b":10,"autoexpand":True},
        width=png_size[0],
        height=png_size[1],
        title_font_size=18
    )
    fig.write_image(os.path.join(IMG_DIR, f"{col}.png"), scale=4.5)

### Generate Figure 1 ###
generate_maps('LogFedCattle', 'Fed Cattle (head)', [2,6],
              [[2,3,4,5,6], ['100','1,000','10,000','100,000','1,000,000']],
              [900, 600])

### Generate Figure 4 ###
col_list = ['Price','Markdown','MarkdownSpatial','MarkdownMulti','MarkdownContract','CapEffect']
col_dict = {
    'Price':           {'range':[137,145], 'label':'Price ($/cwt)',         'ticks':[[137,139,141,143,145],[137,139,141,143,145]]},
    'Markdown':        {'range':[2.5,5],   'label':'Markdown ($/cwt)',      'ticks':[[2.5,3,3.5,4,4.5,5],['2.50','3.00','3.50','4.00','4.50','5.00']]},
    'MarkdownSpatial': {'range':[1.5,3],   'label':'Markdown ($/cwt)',      'ticks':[[1.5,2,2.5,3],['1.50','2.00','2.50','3.00']]},
    'MarkdownMulti':   {'range':[0,0.5],   'label':'Markdown ($/cwt)',      'ticks':[[0,0.1,0.2,0.3,0.4,0.5],['0','0.10','0.20','0.30','0.40','0.50']]},
    'MarkdownContract':{'range':[0,2],     'label':'Markdown ($/cwt)',      'ticks':[[0,0.5,1,1.5,2],['0','0.50','1.00','1.50','2.00']]},
    'CapEffect':       {'range':[0,2],     'label':'Capacity effect ($/cwt)','ticks':[[0,0.5,1,1.5,2],['0','0.50','1.00','1.50','2.00']]}
}
for col in col_list:
    generate_maps(col, col_dict[col]['label'], col_dict[col]['range'], col_dict[col]['ticks'], [600, 400])

### Generate Figure 2: Plant locations ###
plants = pd.read_csv(os.path.join(DATA_DIR, "plant_data.csv"))
plants['firm_name'] = plants['firm_id'].map({1: "JBS", 2: "Tyson", 3: "Cargill", 4: "National Beef", 5: "Other"}).fillna("Other")

fig = px.scatter_geo(plants, lat="lat", lon="lon",
                     size="capacity", color="firm_name",
                     labels={'firm_name': 'Firm Ownership'},
                     color_discrete_sequence=px.colors.qualitative.G10)

fig.update_geos(scope='usa', landcolor="LightGray")
fig.update_layout(
    margin={"l":0,"r":0,"t":0,"b":0,"pad":4,"autoexpand":True},
    legend=dict(xanchor="right", x=0.85, yanchor="bottom", y=0.1, font_size=14),
    width=1200,
    height=500,
    title_font_size=18
)
fig.write_image(os.path.join(IMG_DIR, "plants.png"), scale=1.5)
