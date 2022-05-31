import scipy as scp
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
import csv
from scipy.cluster.hierarchy import dendrogram, linkage
import base64
from enum import Enum
from io import BytesIO, StringIO
from typing import Union
import streamlit as st
from clustergrammer_widget import *
net = Network(clustergrammer_widget)
from clustergrammer2 import net
st.title('Babylon: A Stream-Lit On-The-Fly Hierarchal Clustering Application for User Gene Lists')

@st.cache
def load_data():
    df = pd.read_csv("../data/data.csv")
    return df
@st.cache
def load_covariates():
    covariates=pd.read_csv("../data/covariates.csv")
    return covariates
# Create a text element and let the reader know the data is loading.
data_load_state = st.text('Loading data...')
#load_data()
#df = pd.read_csv("../data/Hugo_data.csv")
df = load_data()
covariates=load_covariates()
# Notify the reader that the data was successfully loaded.
data_load_state.text("Done!")

# USER INPUT GENES
gene_set_option = st.selectbox('Gene Set', ("KEGG AMYOTROPHIC LATERAL SCLEROSIS ALS", "KEGG UBIQUITIN MEDIATED PROTEOLYSIS","KEGG RNA DEGRADATION","KEGG CALCIUM SIGNALING PATHWAY","SYNAPTIC VESICLE","REGULATION OF NEURON APOPTOSIS","PROTEASOME"),key = "gene_choice")
if gene_set_option=="KEGG AMYOTROPHIC LATERAL SCLEROSIS ALS":
    user_input = st.text_area("Gene Set", "ALS2 APAF1 BAD BAX BCL2 BCL2L1 BID CASP1 CASP3 CASP9 CAT CCS CHP1 CHP2 CYCS DAXX DERL1 GPX1 GRIA1 GRIA2 GRIN1 GRIN2A GRIN2B GRIN2C GRIN2D MAP2K3 MAP2K6 MAP3K5 MAPK11 MAPK12 MAPK13 MAPK14 NEFH NEFL NEFM NOS1 PPP3CA PPP3CB PPP3CC PPP3R1 PPP3R2 PRPH PRPH2 RAB5A RAC1 SLC1A2 SOD1 TNF TNFRSF1A TNFRSF1B TOMM40 TOMM40L TP53", key="KEGG_ALS")
if gene_set_option=="KEGG UBIQUITIN MEDIATED PROTEOLYSIS":
    user_input= st.text_area("Gene Set", "AIRE ANAPC1 ANAPC10 ANAPC11 ANAPC13 ANAPC2 ANAPC4 ANAPC5 ANAPC7 BIRC2 BIRC3 BIRC6 BRCA1 BTRC CBL CBLB CBLC CDC16 CDC20 CDC23 CDC26 CDC27 CDC34 COP1 CUL1 CUL2 CUL3 CUL4A CUL4B CUL5 CUL7 DDB1 DDB2 DET1 ELOB ELOC ERCC8 FANCL FBXO2 FBXO4 FBXW11 FBXW7 FBXW8 FZR1 HERC1 HERC2 HERC3 HERC4 HUWE1 ITCH KEAP1 KLHL13 KLHL9 MAP3K1 MDM2 MGRN1 MID1 NEDD4 NEDD4L NHLRC1 PIAS1 PIAS2 PIAS3 PIAS4 PML PPIL2 PRKN PRPF19 RBX1 RCHY1 RHOBTB2 RNF7 SAE1 SIAH1 SKP1 SKP1P2 SKP2 SMURF1 SMURF2 SOCS1 SOCS3 STUB1 SYVN1 TRAF6 TRIM32 TRIM37 TRIP12 UBA1 UBA2 UBA3 UBA6 UBA7 UBE2A UBE2B UBE2C UBE2D1 UBE2D2 UBE2D3 UBE2D4 UBE2E1 UBE2E2 UBE2E3 UBE2F UBE2G1 UBE2G2 UBE2H UBE2I UBE2J1 UBE2J2 UBE2K UBE2L3 UBE2L6 UBE2M UBE2N UBE2NL UBE2O UBE2Q1 UBE2Q2 UBE2QL1 UBE2R2 UBE2S UBE2U UBE2W UBE2Z UBE3A UBE3B UBE3C UBE4A UBE4B UBOX5 UBR5 VHL WWP1 WWP2 XIAP",key="KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS")
if gene_set_option == "KEGG RNA DEGRADATION":
    user_input= st.text_area("Gene Set", "C1D C1DP2 C1DP3 CNOT1 CNOT10 CNOT2 CNOT3 CNOT4 CNOT6 CNOT6L CNOT7 CNOT8 CNOT9 DCP1A DCP1B DCP2 DCPS DDX6 DIS3 EDC3 EDC4 ENO1 ENO2 ENO3 EXOSC1 EXOSC10 EXOSC2 EXOSC3 EXOSC4 EXOSC5 EXOSC6 EXOSC7 EXOSC8 EXOSC9 HSPA9 HSPD1 LSM1 LSM2 LSM3 LSM4 LSM5 LSM6 LSM7 LSM8 MPHOSPH6 MTREX PAPOLA PAPOLB PAPOLG PARN PATL1 PNPT1 SKIV2L TENT4A TTC37 WDR61 XRN1 XRN2 ZCCHC7",key="KEGG_RNA_DEGRADATION")
if gene_set_option == "KEGG CALCIUM SIGNALING PATHWAY":
    user_input= st.text_area("Gene Set","ADCY1 ADCY2 ADCY3 ADCY4 ADCY7 ADCY8 ADCY9 ADORA2A ADORA2B ADRA1A ADRA1B ADRA1D ADRB1 ADRB2 ADRB3 AGTR1 ATP2A1 ATP2A2 ATP2A3 ATP2B1 ATP2B2 ATP2B3 ATP2B4 AVPR1A AVPR1B BDKRB1 BDKRB2 BST1 CACNA1A CACNA1B CACNA1C CACNA1D CACNA1E CACNA1F CACNA1G CACNA1H CACNA1I CACNA1S CALM1 CALM2 CALM3 CALML3 CALML5 CALML6 CAMK2A CAMK2B CAMK2D CAMK2G CAMK4 CCKAR CCKBR CD38 CHP1 CHP2 CHRM1 CHRM2 CHRM3 CHRM5 CHRNA7 CYSLTR1 CYSLTR2 DRD1 DRD5 EDNRA EDNRB EGFR ERBB2 ERBB3 ERBB4 F2R GNA11 GNA14 GNA15 GNAL GNAQ GNAS GRIN1 GRIN2A GRIN2C GRIN2D GRM1 GRM5 GRPR HRH1 HRH2 HTR2A HTR2B HTR2C HTR4 HTR5A HTR6 HTR7 ITPKA ITPKB ITPR1 ITPR2 ITPR3 LHCGR LTB4R2 MYLK MYLK2 MYLK3 NOS1 NOS2 NOS3 NTSR1 OXTR P2RX1 P2RX2 P2RX3 P2RX4 P2RX5 P2RX6 P2RX7 PDE1A PDE1B PDE1C PDGFRA PDGFRB PHKA1 PHKA2 PHKB PHKG1 PHKG2 PLCB1 PLCB2 PLCB3 PLCB4 PLCD1 PLCD3 PLCD4 PLCE1 PLCG1 PLCG2 PLCZ1 PLN PPID PPP3CA PPP3CB PPP3CC PPP3R1 PPP3R2 PRKACA PRKACB PRKACG PRKCA PRKCB PRKCG PRKX PTAFR PTGER1 PTGER3 PTGFR PTK2B RYR1 RYR2 RYR3 SLC25A31 SLC25A4 SLC25A5 SLC25A6 SLC8A1 SLC8A2 SLC8A3 SPHK1 SPHK2 TACR1 TACR2 TACR3 TBXA2R TNNC1 TNNC2 TRHR TRPC1 VDAC1 VDAC2 VDAC2P5 VDAC3",key="KEGG_CALCIUM_SIGNALING_PATHWAY")
if gene_set_option == "SYNAPTIC VESICLE":
    user_input= st.text_area("Gene Set","AMPH APBA1 CLN3 DMXL2 HCRT ICA1 MT3 PPT1 RAB3A SEPTIN5 SLC30A3 SNAPIN SVOP SYN3 TRAPPC4",key="SYNAPTIC_VESICLE")
if gene_set_option == "REGULATION OF NEURON APOPTOSIS":
    user_input= st.text_area("Gene Set","AKT1S1 BAX CDK5 CDK5R1 GDNF GRM4 NF1 PCSK9 PPT1 RASA1 SOD1 TGFB2",key="REGULATION_OF_NEURON_APOPTOSIS")
if gene_set_option == "PROTEASOME":
    user_input= st.text_area("Gene Set","IFNG POMP PSMA1 PSMA2 PSMA3 PSMA4 PSMA5 PSMA6 PSMA6P4 PSMA7 PSMA8 PSMB1 PSMB10 PSMB11 PSMB2 PSMB3 PSMB4 PSMB5 PSMB6 PSMB7 PSMB8 PSMB9 PSMC1 PSMC1P4 PSMC2 PSMC3 PSMC4 PSMC5 PSMC6 PSMD1 PSMD11 PSMD12 PSMD13 PSMD14 PSMD2 PSMD3 PSMD4 PSMD6 PSMD7 PSMD8 PSME1 PSME2 PSME3 PSME4 PSMF1 SEM1",key="KEGG_PROTEASOME")
#else:
    #user_input = st.text_area("Gene Set", "ALS2 APAF1 BAD BAX BCL2 BCL2L1 BID CASP1 CASP3 CASP9 CAT CCS CHP1 CHP2 CYCS DAXX DERL1 GPX1 GRIA1 GRIA2 GRIN1 GRIN2A GRIN2B GRIN2C GRIN2D MAP2K3 MAP2K6 MAP3K5 MAPK11 MAPK12 MAPK13 MAPK14 NEFH NEFL NEFM NOS1 PPP3CA PPP3CB PPP3CC PPP3R1 PPP3R2 PRPH PRPH2 RAB5A RAC1 SLC1A2 SOD1 TNF TNFRSF1A TNFRSF1B TOMM40 TOMM40L TP53",key="default_KEGG")
#user_input = st.text_input("Gene List", 0)
gene_list = user_input.split()

#Grab selected genes
#st.write(gene_list)
df_selected = df.query('gene_id in @gene_list')
#df_selected=df[df['gene_id'].isin(gene_list)]
#st.write(df_selected)
df_selected = df_selected.set_index('gene_id')
#df_selected=df_selected.astype('float16')

# Log Scale selected data Option
st.write("Normalizing Options")
if st.checkbox('Log scale'):
    #Log Scale the selected data
    df_selected = (1+df_selected)/2 # (-1,1] -> (0,1]
    df_selected=np.log(df_selected)

#Covariate Selection
st.write("Covariate Selection")
sex=st.checkbox('Sex')
race=st.checkbox('Race')
ethnicity=st.checkbox('Ethnicity')
subject_group=st.checkbox('Subject Group')
site_of_onset=st.checkbox('Site of Onset')
#Color Options
#color_option = st.selectbox('Color', ("mako", "RdBu","vlag"),key = "color_choice")
color_option = st.selectbox('Color', ("aggrnyl", "agsunset", "algae", "amp", "armyrose", "balance","blackbody", "bluered", "blues", "blugrn", "bluyl", "brbg","brwnyl", "bugn", "bupu", "burg", "burgyl", "cividis", "curl", "darkmint", "deep", "delta", "dense", "earth", "edge", "electric","emrld", "fall", "geyser", "gnbu", "gray", "greens", "greys","haline", "hot", "hsv", "ice", "icefire", "inferno", "jet","magenta", "magma", "matter", "mint", "mrybm", "mygbm", "oranges","orrd", "oryel", "oxy", "peach", "phase", "picnic", "pinkyl","piyg", "plasma", "plotly3", "portland", "prgn", "pubu", "pubugn","puor", "purd", "purp", "purples", "purpor", "rainbow", "rdbu","rdgy", "rdpu", "rdylbu", "rdylgn", "redor", "reds", "solar","spectral", "speed", "sunset", "sunsetdark", "teal", "tealgrn","tealrose", "tempo", "temps", "thermal", "tropic", "turbid","turbo", "twilight", "viridis", "ylgn", "ylgnbu", "ylorbr","ylorrd"), key="color_choice")
# PLOT FIGURE
def plot():
    net.load_df(df_selected)
    net.cluster()
    fig=net.widget()
    #fig
plot()

# if sex:
#     #Select Sex
#     sex_dict=pd.Series(covariates.Sex.values,index=covariates["Participant_ID"]).to_dict()
#     df_selected.rename(columns=sex_dict,inplace=True)
#     df_selected = df_selected.loc[:, df_selected.columns.notnull()]
#     lut = dict(zip(df_selected.columns[0:len(df_selected.columns)].unique(), "rb"))
#     col_colors = df_selected.columns[0:len(df_selected.columns)].map(lut)
#     def plot():
#         fig=dash_bio.Clustergram(data=df_selected,column_labels=list(df_selected.columns.values),row_labels=list(df_selected.index),height=1000,width=5000,color_map=color_option,line_width=1)
#         #fig=sns.clustermap(df_selected, metric="euclidean", standard_scale=1, method="ward",cmap=color_option,col_colors=col_colors)
#         fig
#     plot()
# if race:
#     #Select Race data
#     race_dict=pd.Series(covariates.Race.values,index=covariates["Participant_ID"]).to_dict()
#     df_selected.rename(columns=race_dict, inplace=True)
#     df_selected = df_selected.loc[:, df_selected.columns.notnull()]
#     lut = dict(zip(df_selected.columns[0:len(df_selected.columns)].unique(), "rbgcmykw"))
#     col_colors = df_selected.columns[0:len(df_selected.columns)].map(lut)
#     def plot():
#         fig=dash_bio.Clustergram(data=df_selected,column_labels=list(df_selected.columns.values),row_labels=list(df_selected.index),height=1000,width=5000,color_map=color_option,line_width=1)
#         #fig=sns.clustermap(df_selected, metric="euclidean", standard_scale=1, method="ward",cmap=color_option,col_colors=col_colors)
#         fig
#     plot()
# if ethnicity:
#     #Select Ethnicity data
#     ethnicity_dict=pd.Series(covariates.Ethnicity.values,index=covariates["Participant_ID"]).to_dict()
#     df_selected.rename(columns=ethnicity_dict, inplace=True)
#     df_selected = df_selected.loc[:, df_selected.columns.notnull()]
#     lut = dict(zip(df_selected.columns[0:len(df_selected.columns)].unique(), "rb"))
#     col_colors = df_selected.columns[0:len(df_selected.columns)].map(lut)
#     def plot():
#         fig=dash_bio.Clustergram(data=df_selected,column_labels=list(df_selected.columns.values),row_labels=list(df_selected.index),height=1000,width=5000,color_map=color_option,line_width=1)
#         #fig=sns.clustermap(df_selected, metric="euclidean", standard_scale=1, method="ward",cmap=color_option,col_colors=col_colors)
#         fig
#     plot()
# if subject_group:
#     subject_group_dict=pd.Series(covariates["Subject Group"].values,index=covariates["Participant_ID"]).to_dict()
#     df_selected.rename(columns=subject_group_dict, inplace=True)
#     df_selected = df_selected.loc[:, df_selected.columns.notnull()]
#     lut = dict(zip(df_selected.columns[0:len(df_selected.columns)].unique(), "rbgc"))
#     col_colors = df_selected.columns[0:len(df_selected.columns)].map(lut)
#     def plot():
#         fig=dash_bio.Clustergram(data=df_selected,column_labels=list(df_selected.columns.values),row_labels=list(df_selected.index),height=1000,width=5000,color_map=color_option,line_width=1)
#         #fig=sns.clustermap(df_selected, metric="euclidean", standard_scale=1, method="ward",cmap=color_option,col_colors=col_colors)
#         fig
#     plot()
# if site_of_onset:
#     site_of_onset_dict=pd.Series(covariates["Site of Onset"].values,index=covariates["Participant_ID"]).to_dict()
#     df_selected.rename(columns=site_of_onset_dict, inplace=True)
#     df_selected = df_selected.loc[:, df_selected.columns.notnull()]
#     lut = dict(zip(df_selected.columns[0:len(df_selected.columns)].unique(), "rbgcm"))
#     col_colors = df_selected.columns[0:len(df_selected.columns)].map(lut)
#     def plot():
#         fig=dash_bio.Clustergram(data=df_selected,column_labels=list(df_selected.columns.values),row_labels=list(df_selected.index),height=1000,width=5000,color_map=color_option,line_width=1)
#         #fig=sns.clustermap(df_selected, metric="euclidean", standard_scale=1, method="ward",cmap=color_option,col_colors=col_colors)
#         fig
#     plot()
# else:
#     def plot():
#         fig=dash_bio.Clustergram(data=df_selected,column_labels=list(df_selected.columns.values),row_labels=list(df_selected.index),height=1000,width=1000,color_map=color_option,line_width=1)
#         #fig=sns.clustermap(df_selected, metric="euclidean", standard_scale=1, method="ward",cmap=color_option)
#         fig
#     plot()

#fig=sns.clustermap(df_selected, metric="euclidean", standard_scale=1, method="ward",cmap=color_option)
#st.pyplot(plot())

# Download clustermap
    # look into kwargs arg on pyplot for plt.savefig option
# Download selected data_table Option
@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')
csv = convert_df(df_selected)
st.download_button(
    label="Download selected data as CSV",
    data=csv,
    file_name='selected_data.csv',
    mime='text/csv',)


#csv = df_selected.to_csv(index=True)
#b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
#href = f'<a href="data:file/csv;base64,{b64}">Download Selected Data</a> (right-click and save as &lt;some_name&gt;.csv)'
#st.markdown(href, unsafe_allow_html=True)
