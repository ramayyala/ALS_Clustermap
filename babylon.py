import scipy as scp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import csv
from clustergrammer_widget import *
net = Network(clustergrammer_widget)
from clustergrammer2 import net
import ipywidgets as widgets
from ipywidgets import interact
import panel as pn
pn.extension('ipywidgets')
# set a title for your dashboard
title = '### ALS Clustermap'
subtitle = 'test'

#gene_list = open("../data/geneset.txt").read().splitlines()
#gene_list
user_input = pn.widgets.input.TextAreaInput(value='ALS2 APAF1 BAD BAX BCL2 BCL2L1 BID CASP1 CASP3 CASP9 CAT CCS CHP1 CHP2 CYCS DAXX DERL1 GPX1 GRIA1 GRIA2 GRIN1 GRIN2A GRIN2B GRIN2C GRIN2D MAP2K3 MAP2K6 MAP3K5 MAPK11 MAPK12 MAPK13 MAPK14 NEFH NEFL NEFM NOS1 PPP3CA PPP3CB PPP3CC PPP3R1 PPP3R2 PRPH PRPH2 RAB5A RAC1 SLC1A2 SOD1 TNF TNFRSF1A TNFRSF1B TOMM40 TOMM40L TP53',name='Gene List Input')
#user_input.value='ALS2 APAF1 BAD BAX BCL2 BCL2L1 BID CASP1 CASP3 CASP9 CAT CCS CHP1 CHP2 CYCS DAXX DERL1 GPX1 GRIA1 GRIA2 GRIN1 GRIN2A GRIN2B GRIN2C GRIN2D MAP2K3 MAP2K6 MAP3K5 MAPK11 MAPK12 MAPK13 MAPK14 NEFH NEFL NEFM NOS1 PPP3CA PPP3CB PPP3CC PPP3R1 PPP3R2 PRPH PRPH2 RAB5A RAC1 SLC1A2 SOD1 TNF TNFRSF1A TNFRSF1B TOMM40 TOMM40L TP53'
#gene_list = open("../data/geneset.txt").read().splitlines()
#user_dropdown = pn.widgets.Select(name='Select Gene Set', options=['KEGG AMYOTROPHIC LATERAL SCLEROSIS ALS', 'KEGG UBIQUITIN MEDIATED PROTEOLYSIS', 'KEGG RNA DEGRADATION'],value='KEGG AMYOTROPHIC LATERAL SCLEROSIS ALS')
#user_dropdown.value='KEGG AMYOTROPHIC LATERAL SCLEROSIS ALS'
#if user_dropdown.value == 'KEGG AMYOTROPHIC LATERAL SCLEROSIS ALS':
#    user_input = 'ALS2 APAF1 BAD BAX BCL2 BCL2L1 BID CASP1 CASP3 CASP9 CAT CCS CHP1 CHP2 CYCS DAXX DERL1 GPX1 GRIA1 GRIA2 GRIN1 GRIN2A GRIN2B GRIN2C GRIN2D MAP2K3 MAP2K6 MAP3K5 MAPK11 MAPK12 MAPK13 MAPK14 NEFH NEFL NEFM NOS1 PPP3CA PPP3CB PPP3CC PPP3R1 PPP3R2 PRPH PRPH2 RAB5A RAC1 SLC1A2 SOD1 TNF TNFRSF1A TNFRSF1B TOMM40 TOMM40L TP53'
#if user_dropdown.value == 'KEGG UBIQUITIN MEDIATED PROTEOLYSIS':
#    user_input = 'AIRE ANAPC1 ANAPC10 ANAPC11 ANAPC13 ANAPC2 ANAPC4 ANAPC5 ANAPC7 BIRC2 BIRC3 BIRC6 BRCA1 BTRC CBL CBLB CBLC CDC16 CDC20 CDC23 CDC26 CDC27 CDC34 COP1 CUL1 CUL2 CUL3 CUL4A CUL4B CUL5 CUL7 DDB1 DDB2 DET1 ELOB ELOC ERCC8 FANCL FBXO2 FBXO4 FBXW11 FBXW7 FBXW8 FZR1 HERC1 HERC2 HERC3 HERC4 HUWE1 ITCH KEAP1 KLHL13 KLHL9 MAP3K1 MDM2 MGRN1 MID1 NEDD4 NEDD4L NHLRC1 PIAS1 PIAS2 PIAS3 PIAS4 PML PPIL2 PRKN PRPF19 RBX1 RCHY1 RHOBTB2 RNF7 SAE1 SIAH1 SKP1 SKP1P2 SKP2 SMURF1 SMURF2 SOCS1 SOCS3 STUB1 SYVN1 TRAF6 TRIM32 TRIM37 TRIP12 UBA1 UBA2 UBA3 UBA6 UBA7 UBE2A UBE2B UBE2C UBE2D1 UBE2D2 UBE2D3 UBE2D4 UBE2E1 UBE2E2 UBE2E3 UBE2F UBE2G1 UBE2G2 UBE2H UBE2I UBE2J1 UBE2J2 UBE2K UBE2L3 UBE2L6 UBE2M UBE2N UBE2NL UBE2O UBE2Q1 UBE2Q2 UBE2QL1 UBE2R2 UBE2S UBE2U UBE2W UBE2Z UBE3A UBE3B UBE3C UBE4A UBE4B UBOX5 UBR5 VHL WWP1 WWP2 XIAP'
#if user_dropdown.value == 'KEGG RNA DEGRADATION':
#    user_input ='C1D C1DP2 C1DP3 CNOT1 CNOT10 CNOT2 CNOT3 CNOT4 CNOT6 CNOT6L CNOT7 CNOT8 CNOT9 DCP1A DCP1B DCP2 DCPS DDX6 DIS3 EDC3 EDC4 ENO1 ENO2 ENO3 EXOSC1 EXOSC10 EXOSC2 EXOSC3 EXOSC4 EXOSC5 EXOSC6 EXOSC7 EXOSC8 EXOSC9 HSPA9 HSPD1 LSM1 LSM2 LSM3 LSM4 LSM5 LSM6 LSM7 LSM8 MPHOSPH6 MTREX PAPOLA PAPOLB PAPOLG PARN PATL1 PNPT1 SKIV2L TENT4A TTC37 WDR61 XRN1 XRN2 ZCCHC7'
gene_list=list(user_input.value.split())
df=pd.read_csv("https://media.githubusercontent.com/media/ramayyala/ALS_Clustermap/master/data/data.csv")
# Load and format the data
# tell Panel what your plot "depends" on.
# This defines what should trigger a change in the chart.
# both values in depends() will be used in our below Altair chart as filters
@pn.depends(user_input)
def get_plot(gene_set): # start function
    df_selected=df[df['gene_id'].isin(list(gene_set.split()))]
    df_selected = df_selected.set_index('gene_id')
    net.load_df(df_selected)
    net.cluster()
    fig=net.widget()
    return fig

# create the Panel object, passing in all smaller objects
#text_area_input.controls(jslink=True), text_area_input
pn.Row(
   pn.Column(title, subtitle,user_input),get_plot).servable()
