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
import requests, json
pn.extension('ipywidgets')
material = pn.template.BootstrapTemplate(site_url="https://dataportal.answerals.org/search",logo="https://raw.githubusercontent.com/ramayyala/ALS_Clustermap/master/misc/logo.png",title='ANSWER ALS CLUSTERMAP',header_background="#204cac",sidebar_width=410)

#User Input
_user_input = {
    'KEGG AMYOTROPHIC LATERAL SCLEROSIS ALS':['ALS2 APAF1 BAD BAX BCL2 BCL2L1 BID CASP1 CASP3 CASP9 CAT CCS CHP1 CHP2 CYCS DAXX DERL1 GPX1 GRIA1 GRIA2 GRIN1 GRIN2A GRIN2B GRIN2C GRIN2D MAP2K3 MAP2K6 MAP3K5 MAPK11 MAPK12 MAPK13 MAPK14 NEFH NEFL NEFM NOS1 PPP3CA PPP3CB PPP3CC PPP3R1 PPP3R2 PRPH PRPH2 RAB5A RAC1 SLC1A2 SOD1 TNF TNFRSF1A TNFRSF1B TOMM40 TOMM40L TP53'],
    'KEGG UBIQUITIN MEDIATED PROTEOLYSIS':['AIRE ANAPC1 ANAPC10 ANAPC11 ANAPC13 ANAPC2 ANAPC4 ANAPC5 ANAPC7 BIRC2 BIRC3 BIRC6 BRCA1 BTRC CBL CBLB CBLC CDC16 CDC20 CDC23 CDC26 CDC27 CDC34 COP1 CUL1 CUL2 CUL3 CUL4A CUL4B CUL5 CUL7 DDB1 DDB2 DET1 ELOB ELOC ERCC8 FANCL FBXO2 FBXO4 FBXW11 FBXW7 FBXW8 FZR1 HERC1 HERC2 HERC3 HERC4 HUWE1 ITCH KEAP1 KLHL13 KLHL9 MAP3K1 MDM2 MGRN1 MID1 NEDD4 NEDD4L NHLRC1 PIAS1 PIAS2 PIAS3 PIAS4 PML PPIL2 PRKN PRPF19 RBX1 RCHY1 RHOBTB2 RNF7 SAE1 SIAH1 SKP1 SKP1P2 SKP2 SMURF1 SMURF2 SOCS1 SOCS3 STUB1 SYVN1 TRAF6 TRIM32 TRIM37 TRIP12 UBA1 UBA2 UBA3 UBA6 UBA7 UBE2A UBE2B UBE2C UBE2D1 UBE2D2 UBE2D3 UBE2D4 UBE2E1 UBE2E2 UBE2E3 UBE2F UBE2G1 UBE2G2 UBE2H UBE2I UBE2J1 UBE2J2 UBE2K UBE2L3 UBE2L6 UBE2M UBE2N UBE2NL UBE2O UBE2Q1 UBE2Q2 UBE2QL1 UBE2R2 UBE2S UBE2U UBE2W UBE2Z UBE3A UBE3B UBE3C UBE4A UBE4B UBOX5 UBR5 VHL WWP1 WWP2 XIAP'],
    'KEGG RNA DEGRADATION':['C1D C1DP2 C1DP3 CNOT1 CNOT10 CNOT2 CNOT3 CNOT4 CNOT6 CNOT6L CNOT7 CNOT8 CNOT9 DCP1A DCP1B DCP2 DCPS DDX6 DIS3 EDC3 EDC4 ENO1 ENO2 ENO3 EXOSC1 EXOSC10 EXOSC2 EXOSC3 EXOSC4 EXOSC5 EXOSC6 EXOSC7 EXOSC8 EXOSC9 HSPA9 HSPD1 LSM1 LSM2 LSM3 LSM4 LSM5 LSM6 LSM7 LSM8 MPHOSPH6 MTREX PAPOLA PAPOLB PAPOLG PARN PATL1 PNPT1 SKIV2L TENT4A TTC37 WDR61 XRN1 XRN2 ZCCHC7'],
    'KEGG CALCIUM SIGNALING PATHWAY':['ADCY1 ADCY2 ADCY3 ADCY4 ADCY7 ADCY8 ADCY9 ADORA2A ADORA2B ADRA1A ADRA1B ADRA1D ADRB1 ADRB2 ADRB3 AGTR1 ATP2A1 ATP2A2 ATP2A3 ATP2B1 ATP2B2 ATP2B3 ATP2B4 AVPR1A AVPR1B BDKRB1 BDKRB2 BST1 CACNA1A CACNA1B CACNA1C CACNA1D CACNA1E CACNA1F CACNA1G CACNA1H CACNA1I CACNA1S CALM1 CALM2 CALM3 CALML3 CALML5 CALML6 CAMK2A CAMK2B CAMK2D CAMK2G CAMK4 CCKAR CCKBR CD38 CHP1 CHP2 CHRM1 CHRM2 CHRM3 CHRM5 CHRNA7 CYSLTR1 CYSLTR2 DRD1 DRD5 EDNRA EDNRB EGFR ERBB2 ERBB3 ERBB4 F2R GNA11 GNA14 GNA15 GNAL GNAQ GNAS GRIN1 GRIN2A GRIN2C GRIN2D GRM1 GRM5 GRPR HRH1 HRH2 HTR2A HTR2B HTR2C HTR4 HTR5A HTR6 HTR7 ITPKA ITPKB ITPR1 ITPR2 ITPR3 LHCGR LTB4R2 MYLK MYLK2 MYLK3 NOS1 NOS2 NOS3 NTSR1 OXTR P2RX1 P2RX2 P2RX3 P2RX4 P2RX5 P2RX6 P2RX7 PDE1A PDE1B PDE1C PDGFRA PDGFRB PHKA1 PHKA2 PHKB PHKG1 PHKG2 PLCB1 PLCB2 PLCB3 PLCB4 PLCD1 PLCD3 PLCD4 PLCE1 PLCG1 PLCG2 PLCZ1 PLN PPID PPP3CA PPP3CB PPP3CC PPP3R1 PPP3R2 PRKACA PRKACB PRKACG PRKCA PRKCB PRKCG PRKX PTAFR PTGER1 PTGER3 PTGFR PTK2B RYR1 RYR2 RYR3 SLC25A31 SLC25A4 SLC25A5 SLC25A6 SLC8A1 SLC8A2 SLC8A3 SPHK1 SPHK2 TACR1 TACR2 TACR3 TBXA2R TNNC1 TNNC2 TRHR TRPC1 VDAC1 VDAC2 VDAC2P5 VDAC3'],
    'SYNAPTIC VESICLE':['AMPH APBA1 CLN3 DMXL2 HCRT ICA1 MT3 PPT1 RAB3A SEPTIN5 SLC30A3 SNAPIN SVOP SYN3 TRAPPC4'],
    'REGULATION OF NEURON APOPTOSIS':['AKT1S1 BAX CDK5 CDK5R1 GDNF GRM4 NF1 PCSK9 PPT1 RASA1 SOD1 TGFB2'],
    'PROTEASOME':['IFNG POMP PSMA1 PSMA2 PSMA3 PSMA4 PSMA5 PSMA6 PSMA6P4 PSMA7 PSMA8 PSMB1 PSMB10 PSMB11 PSMB2 PSMB3 PSMB4 PSMB5 PSMB6 PSMB7 PSMB8 PSMB9 PSMC1 PSMC1P4 PSMC2 PSMC3 PSMC4 PSMC5 PSMC6 PSMD1 PSMD11 PSMD12 PSMD13 PSMD14 PSMD2 PSMD3 PSMD4 PSMD6 PSMD7 PSMD8 PSME1 PSME2 PSME3 PSME4 PSMF1 SEM1'],
    'NFKAPPAB65_01':['AAMDC ABI3 ACAN ACTN3 AKT1S1 ALG6 AMOTL1 APPL1 ARHGEF2 ASCL3 ASH1L ATP1B1 BAZ2B BCL3 BCL6B BDNF BIRC3 BLCAP BMF BMP2K C1QL1 CACNG3 CALCOCO1 CASKIN2 CCDC107 CCL5 CCM2L CD40 CD69 CD70 CD86 CDC42SE1 CDK6 CHD4 CHD6 CLCN1 CLCN2 CLDN5 CLOCK COL11A2 COL16A1 COQ8B CREB1 CSF1R CSF2RB CTDSP1 CTDSPL2 CUEDC1 CXCL10 CXCL11 CXCL16 CXCR5 CYLD CYP2D6 DAP3 DCLK1 DDR1 DOCK4 DSC2 E2F3 EBF1 EHF EIF4A2 EIF4G1 EIF5A ENO3 ERN1 FAM117A FAM43B FGF1 FGF12 FGF17 FLOT1 FOXS1 FTHL17 FUT7 G3BP1 GADD45B GATA4 GDPD5 GNG4 GNGT2 GPBP1 GPM6A GREM1 GRK5 HIVEP1 HNRNPR HSD3B7 HSP90B1 HTR3B ICAM1 IER3 IER5 IFNB1 IL13 IL17C IL1RAPL1 IL27 IL6ST ILK ITPKC JAK3 KANSL1L KAT7 KCNN2 KCNT2 KLK9 KRT23 KRT36 KY LAMA1 LINC01138 LIX1L LRCH1 LTB MADCAM1 MAML2 MAP3K11 MAP3K8 MAPK6 MIA MIDEAS MIR17HG MITF MLLT11 MLLT6 MMP9 MOB3C MSC MSX1 NDUFB9 NFAT5 NFKB2 NFKBIA NFKBIB NFKBID NLK NR2F2 NXPH4 ORAI1 PAN2 PARP8 PCBP4 PCDH10 PCDH12 PCSK2 PFN1 PLXNB1 PNKD POU2F3 PPP1R13B PRDM12 PRRT2 PTGES PTHLH PURG RANBP10 RAP2C RASGRP4 RASSF2 RBMS1 REL RELB RFX5 RIN2 RND1 RNF43 RPS19 RPS6KA4 RRAS RRP8 RSF1 S1PR2 SDC4 SEC63 SH2B3 SIN3A SIRT2 SIX4 SIX5 SLAMF8 SLC12A2 SLC16A6 SLC44A1 SLC6A12 SMOC1 SMPD3 SOX10 SOX3 SOX5 SP6 STAT6 STX19 STX4 SUCO TATDN1 TBC1D17 TCEA2 TFE3 TIAL1 TJAP1 TLX1 TLX3 TNFRSF1B TNFRSF9 TNFSF15 TNFSF18 TNIP1 TP53 TP63 TRAF4 TRIB2 TRIM47 TRPC4 TSEN54 TSLP TSNAXIP1 TSPEAR TUT1 UACA UBD UBE2D3 UBE2H UBE2I UPF2 VEZF1 WNT10A WNT10B WNT4 WRAP53 WRN YWHAQ YWHAZ YY1AP1 ZDHHC24 ZDHHC8 ZEB1 ZFHX3 ZIC4 ZMYND15'],
    'NFKAPPAB_01':['AAMDC ABHD8 ACAN ACTN1 ACTN3 ADGRB2 AGPAT1 ALG6 ANKFN1 AP1S2 APPL1 ARHGAP44 ARHGAP5 ARHGAP8 ARPC2 ASCL3 ASH1L ATOH1 ATP1B1 AZIN1 BCKDK BCL3 BCL6B BDNF BFSP1 BIRC3 BMF BMP2K BNC2 C1QL1 CALCOCO1 CBX2 CCDC107 CCL5 CCN2 CD40 CD70 CD83 CDC14A CDK6 CFAP69 CHD4 CHD6 CLCN1 COL11A2 COL16A1 CPD CREB1 CTAGE4 CXCL10 CXCL11 CXCL16 CXCL2 CXCL5 CXCL6 CXCL9 CXCR5 CYB5A CYLD DAP3 DCLK1 DDR1 EBI3 EDN2 EGF EHF EIF5A ENO3 ERN1 ETV6 FAM117A FBXL12 FGF1 FGF17 FLRT1 FOXJ2 FOXS1 FTHL17 FUT7 G3BP1 GABRB1 GADD45B GATA4 GDPD5 GEN1 GGNBP2 GNAO1 GNB1 GNG4 GPBP1 GPHN GREM1 GRIN2D GRK5 HCFC1 HCST HIVEP1 HOXA11 HOXB9 HSD11B2 HSD3B7 HSP90B1 ICAM1 IER5 IFNB1 IGDCC3 IL13 IL1A IL1RAPL1 IL1RN IL23A IL27 IL27RA IL2RA IL4I1 INO80D IRF1 ITGB4 JAK3 JARID2 KANSL1L KCNH3 KCNN3 KCNS3 KCNT2 KLK9 KRTAP13-1 KY LAMA1 LIG1 LINC01138 LIX1L LRCH1 LTA LTB MADCAM1 MAG MAML2 MAP3K8 MAP4K2 MED1 MIDEAS MIR17HG MMP9 MOB3C MSC MSN MSX1 NFAT5 NFKB2 NFKBIB NFKBID NIBAN1 NIPBL NLK NOD2 NOL4 NR2F2 NTN1 NUFIP2 OPCML ORAI1 P2RY10 PAN2 PARP8 PCBP4 PCSK2 PFN1 PLAU PNKD PRDM12 PRX PTMS RALGDS RAP2C RBM14 RBMS1 RELB RIMS2 RIN2 RND1 RPS6KA4 RRAS RSF1 S1PR2 SALL1 SAMSN1 SDC4 SEC14L2 SESN2 SH2B3 SHOX2 SIN3A SLAMF8 SLC11A2 SLC12A2 SLC44A1 SLC6A12 SMC6 SMOC1 SMPD3 SOBP SOX10 SOX5 SP6 SPTB STAT6 STC2 STON2 SUN2 TAFAZZIN TAL1 TASL TCEA2 TIFA TJAP1 TLX1 TNF TNFRSF1B TNFRSF9 TNFSF15 TNFSF18 TNIP1 TP53 TP63 TRAF4 TRIB2 TRIM47 TSLP TSPEAR TUBGCP4 UACA UBD UBE2D3 UBE2H UBE2I UBE4B VCAM1 VEZF1 VSX2 WNT10A WNT10B WRAP53 XKR8 YY1AP1 ZBTB11 ZBTB5 ZBTB9 ZFP36L2 ZMYND15 ZNF232 ZNF384 ZNF821 ZSCAN29 ZSWIM9'],
    'TNF Signaling Pathway (GO)':['ACTN4 ADAM17 ADIPOQ AIM2 APOA1 BIRC7 CARD14 CASP1 CCDC3 CD70 CDIP1 CDIP1 CHUK CLDN18 COMMD7 COMMD7 CPNE1 CPNE1 CYLD CYLD EDA2R EXT1 F2RL1 FOXO3 GAS6 GPS2 HIPK1 IKBKB IKBKB IKBKB ILK JAK2 KRT8 KRT8 KRT18 LAPTM5 LIMS1 NAIP1 NAIP2 NAIP5 NAIP6 NAIP7 NFKBIA NKIRAS1 NKIRAS2 NOL3 NR1H4 OTULIN OTULIN OTULIN PELI3 PELI3 PIAS3 PIAS4 PLVAP PPP2CB PRKN PTK2B PTPN2 PYCARD PYCARD RELA RFFL RIPK1 RRAGA SHARPIN SHARPIN SPATA2 SPATA2 SPATA2 SPHK1 SPHK1 ST18 STAT1 SYK TJP2 TNF TNF TNF TNFRSF1A TNFRSF1A TNFRSF1A TNFRSF11A TNFRSF11A TNFRSF11A TNFRSF11A TNFRSF13C TNFRSF17 TNFSF11 TNFSF11 TNFSF11 TNFSF11 TNFSF13B TNFSF18 TNFSF18 TRAF1 TRAF2 TRAF2 TRAF2 TRAF2 TRAF3 TRAF3 TRAF3 TRAF3IP2 TRAF5 TRAF6 TRAIP TRIM32 TRP53 TXNDC17 UBE2K UMOD'],
    'Sex-Linked':['XIST SRY UTY UTX ZFY USP9X USP9Y PRKY EIF1AY DDX3Y DDX3X VCX3B TSIX JPX ZFX PUDP KDM5D']
}
user_dropdown = pn.widgets.Select(value='KEGG AMYOTROPHIC LATERAL SCLEROSIS ALS',name='Select Gene Set', options=['KEGG AMYOTROPHIC LATERAL SCLEROSIS ALS', 'KEGG UBIQUITIN MEDIATED PROTEOLYSIS', 'KEGG RNA DEGRADATION','KEGG CALCIUM SIGNALING PATHWAY','SYNAPTIC VESICLE','REGULATION OF NEURON APOPTOSIS','PROTEASOME','NFKAPPAB65_01','NFKAPPAB_01','TNF Signaling Pathway (GO)','Sex-Linked'])
user_input = pn.widgets.input.TextAreaInput(value=_user_input[user_dropdown.value][0],name='Gene List Input',width=300, height=400)
@pn.depends(user_dropdown.param.value, watch=True)
def _update_input(user_dropdown):
    gene_list=_user_input[user_dropdown]
    user_input.value = gene_list[0]

#Participant List
participant_input = pn.widgets.input.TextAreaInput(value='',name='Participant List',width=300, height=300,max_length=100000)
try:
    url=requests.get(str(pn.state.session_args.get('endpoint')[0].decode('utf-8'))+'?'+'searchId='+str(pn.state.session_args.get('searchId')[0].decode('utf-8'))+'&'+'releaseName='+str(pn.state.session_args.get('releaseName')[0].decode('utf-8')))
    participants=list(url.json())
    participant_input.value=' '.join(participants)
except Exception:
    participant_input.value = 'CASE-NEUAA599TMX CASE-NEUAB000NKC CASE-NEUAF553MJ3 CASE-NEUAG603XLK CASE-NEUAM655HF7 CASE-NEUAT234RK6 CASE-NEUAW157NMJ CASE-NEUAX665ZHY CASE-NEUBA169GXD CASE-NEUBC998WWB CASE-NEUBD218YR3 CASE-NEUBD288RXQ CASE-NEUBK117YXL CASE-NEUBY734PFR CASE-NEUCB613CA9 CASE-NEUCE965ZGK CASE-NEUCF538BRM CASE-NEUCN596RR3 CASE-NEUCU245GBQ CASE-NEUDD018KAH CASE-NEUDE902GCT CASE-NEUDG000ZG5 CASE-NEUDG272XWC CASE-NEUDP155HFH CASE-NEUDT709YHN CASE-NEUEB422WW0 CASE-NEUEK191WYC CASE-NEUEM720BUU CASE-NEUEN017PCJ CASE-NEUEN476CLW CASE-NEUEU558MNK CASE-NEUFB989DT1 CASE-NEUFH122WN7 CASE-NEUFV237VCZ CASE-NEUFY342UNG CASE-NEUGD965XVD CASE-NEUGE540TC4 CASE-NEUGJ081HKR CASE-NEUGP781PDU CASE-NEUGR121TFD CASE-NEUGW326BRV CASE-NEUGW340YEB CASE-NEUHB491NGF CASE-NEUHG644RYB CASE-NEUHG791RV5 CASE-NEUHL999WWT CASE-NEUHV216EHA CASE-NEUHW530WD8 CASE-NEUHW590KBH CASE-NEUHZ302VU4 CASE-NEUJA217MTJ CASE-NEUJA666PYD CASE-NEUJC191RLE CASE-NEUJH197AK2 CASE-NEUJX780XXK CASE-NEUJY426MBU CASE-NEUJY536DKF CASE-NEUKD025JPF CASE-NEULA777DBT CASE-NEULD354RZB CASE-NEULL588ALG CASE-NEULP998KDJ CASE-NEULY177TTN CASE-NEULZ548ZXV CASE-NEUME498PCJ CASE-NEUMN012EVP CASE-NEUMN922PMZ CASE-NEUMY871DGF CASE-NEUNA248WXL CASE-NEUNG326MFP CASE-NEUNL415AFW CASE-NEUPK546ZLD CASE-NEUPN525XEW CASE-NEUPR357AUF CASE-NEUPR600MBU CASE-NEURA639MK8 CASE-NEURJ362MXH CASE-NEURR881FKY CASE-NEURX315DYH CASE-NEURX909UL6 CASE-NEUTA057AF6 CASE-NEUTA689LN5 CASE-NEUTB230DA3 CASE-NEUTD713DE3 CASE-NEUTJ613AH9 CASE-NEUTN952DDG CASE-NEUTU360YJY CASE-NEUUE852BB1 CASE-NEUVF888UHM CASE-NEUVJ560JGZ CASE-NEUVL876PUV CASE-NEUVM674HUA CASE-NEUVN746WKV CASE-NEUVR636VTP CASE-NEUVW680LPK CASE-NEUVW999FP9 CASE-NEUVX902YNL CASE-NEUWD538KT3 CASE-NEUWH380YHV CASE-NEUWJ389FC7 CASE-NEUXP289KRC CASE-NEUXP595FFP CASE-NEUXV122FN3 CASE-NEUXX361RUG CASE-NEUYG298CA9 CASE-NEUYH924UCE CASE-NEUYK029WUU CASE-NEUYL008GM4 CASE-NEUYL149PRF CASE-NEUYP235ZLD CASE-NEUYY225MNZ CASE-NEUYY878JGP CASE-NEUZN836GME CASE-NEUZN936HJ9 CASE-NEUZP278MR4 CASE-NEUZT557DHF CASE-NEUZV656DD1 CASE-NEUZX521TKK CASE-NEUZX847VWV CTRL-NEUAJ025JC3 CTRL-NEUAJ928PAA CTRL-NEUCA748GF2 CTRL-NEUCV136DHM CTRL-NEUCV809LL4 CTRL-NEUDE949BP3 CTRL-NEUDM126GNG CTRL-NEUEU392AE8 CTRL-NEUEY565NWT CTRL-NEUFE306EFY CTRL-NEUFZ500KDB CTRL-NEUFZ508VBV CTRL-NEUHE723FGT CTRL-NEUJH290RH7 CTRL-NEUJX341NDP CTRL-NEUKW131XJ2 CTRL-NEULL933JXY CTRL-NEUMN061ATZ CTRL-NEUNW343RXP CTRL-NEUPH301NNX CTRL-NEURJ861MMD CTRL-NEUUV825HYF CTRL-NEUVZ050YX7 CTRL-NEUWT164JRQ CTRL-NEUYM205MRL CASE-NEUAE228FF6 CASE-NEUAE993EPR CASE-NEUAG241NUD CASE-NEUAG766ULB CASE-NEUAL076FCE CASE-NEUAP285GGU CASE-NEUAW717TN6 CASE-NEUAY067UTB CASE-NEUAZ394JEZ CASE-NEUBC901KL3 CASE-NEUBD141MGD CASE-NEUBM949LA5 CASE-NEUBN979ZZ5 CASE-NEUBT273FH3 CASE-NEUBW008RJ5 CASE-NEUCD502BFU CASE-NEUCG511MJ6 CASE-NEUCP218XVV CASE-NEUCT022VTK CASE-NEUCT842RJV CASE-NEUCU076ADN CASE-NEUCV578DFJ CASE-NEUDH063DEA CASE-NEUDM949BAR CASE-NEUDY534KGP CASE-NEUDZ473EGM CASE-NEUEC400NYR CASE-NEUEM029XXZ CASE-NEUET072VDG CASE-NEUEU318NY2 CASE-NEUEY478NZP CASE-NEUFH461AAN CASE-NEUFU395MJD CASE-NEUFX386KMB CASE-NEUGH995TFK CASE-NEUGL543NJ1 CASE-NEUGR539MVT CASE-NEUHK991AVP CASE-NEUHL096UF0 CASE-NEUHL814WMU CASE-NEUHM532NDD CASE-NEUHT569HT1 CASE-NEUHY206ZEQ CASE-NEUJG311WGV CASE-NEUJG885PY7 CASE-NEUJK720JCL CASE-NEUJL547WVQ CASE-NEUJL595JAZ CASE-NEUJP935AVF CASE-NEUJX990GR5 CASE-NEUKA860NUG CASE-NEUKP160HX8 CASE-NEUKR376CW3 CASE-NEUKV547CFA CASE-NEUKY704JVA CASE-NEULA694DDC CASE-NEULH729AU2 CASE-NEULL648LJ1 CASE-NEULM691PBK CASE-NEULT851ENP CASE-NEULY328CRJ CASE-NEUMB242CLN CASE-NEUMH634MKT CASE-NEUMW598YMT CASE-NEUNJ155KYL CASE-NEUNJ938DUL CASE-NEUNL303HLF CASE-NEUNR020KV6 CASE-NEUNZ171UYF CASE-NEUPG593UTV CASE-NEUPM937TMY CASE-NEUPY050ANK CASE-NEURF720KHA CASE-NEUTC791PVE CASE-NEUTD314VFT CASE-NEUTD866EAH CASE-NEUTL257PNR CASE-NEUTL699WJ0 CASE-NEUTM934BPY CASE-NEUUC458RFT CASE-NEUUD158CJF CASE-NEUUF613PHL CASE-NEUUH658EYH CASE-NEUUJ507EAB CASE-NEUUL256UC9 CASE-NEUUL292XRC CASE-NEUUL311NRQ CASE-NEUUP280KXL CASE-NEUUU506WZF CASE-NEUUZ216XY8 CASE-NEUVD687XD8 CASE-NEUVP060WFG CASE-NEUVR814YCY CASE-NEUVU735HU6 CASE-NEUVV225XKB CASE-NEUWD946RDA CASE-NEUWH955FJF CASE-NEUWM344ZLM CASE-NEUWP426NBR CASE-NEUWX167DAH CASE-NEUWY079WUH CASE-NEUWZ614ARQ CASE-NEUXD985ZU1 CASE-NEUXE491RMM CASE-NEUXG265ME9 CASE-NEUXH833LL4 CASE-NEUXP495XE2 CASE-NEUYC055XPZ CASE-NEUYC198DCR CASE-NEUYC303AJY CASE-NEUYD306HAB CASE-NEUYG208KVV CASE-NEUYJ705VN9 CASE-NEUYP226CPW CASE-NEUYT193FFG CASE-NEUYY614DN8 CASE-NEUZD231YZ6 CASE-NEUZE432DZM CASE-NEUZF321EW4 CASE-NEUZF473UEP CASE-NEUZJ053JGZ CASE-NEUZK054DP5 CASE-NEUZT902WVB CASE-NEUZW701NNF CASE-NEUZY128BJ2 CASE-NEUZY975XKL CTRL-NEUAA485DZL CTRL-NEUDA782GW3 CTRL-NEUDT762KUL CTRL-NEUEB210XRC CTRL-NEUFL733GX5 CTRL-NEUHZ716BZ2 CTRL-NEUMA002VLD CTRL-NEUMF089KLV CTRL-NEUML507PFJ CTRL-NEUMT184NWC CTRL-NEUNC876ZB2 CTRL-NEUNN472ACB CTRL-NEUPL878MTL CTRL-NEUPW536ZKZ CTRL-NEURV546WMW CTRL-NEUWN092BVG CTRL-NEUXC258VTR CTRL-NEUXP955XW7 CTRL-NEUXW311EFC CTRL-NEUZL045YD3' ##All 290 participants selected by default


#Covariate Selector
covariate_selection = pn.widgets.RadioButtonGroup(
    name='Covariates', options=['Sex', 'Ethnicity', 'Race','Subject Group','Case_Control'], button_type='success')

# Color Pickers
positive_col = pn.widgets.ColorPicker(name='Positive Value Color', value='#ff1900')
negative_col = pn.widgets.ColorPicker(name='Negative Value Color', value='#0055ff')

#Load Covariates and Data
covariates=pd.read_csv("covariates.csv")
df=pd.read_csv("data.csv")

#covariates=pd.read_csv("https://media.githubusercontent.com/media/ramayyala/ALS_Clustermap/master/data/covariates.csv")
#df=pd.read_csv("https://media.githubusercontent.com/media/ramayyala/ALS_Clustermap/master/data/data.csv")


@pn.depends(user_input,covariate_selection.param.value,positive_col,negative_col,participant_input)
def get_plot(gene_set,covariate,positive,negative,participant): # start function
    gene_list=list(user_input.value.split())
    df_selected=df[df['gene_id'].isin(list(gene_set.split()))]
    df_selected = df_selected.set_index('gene_id')
    df_selected = (1+df_selected)/2 # (-1,1] -> (0,1]
    df_selected=np.log(df_selected)
    df_selected=df_selected[df_selected.columns.intersection(list(participant.split()))]
    covariate_dict=pd.Series(covariates[covariate].values,index=covariates["Participant_ID"]).to_dict()
    columns=list(df_selected.columns)
    for i in range( len(columns)):
        columns[i]='Sample: '+columns[i],covariate+": "+covariate_dict.get(columns[i])
    df_selected.columns = tuple(columns)
    net.load_df(df_selected)
    net.set_matrix_colors(pos=positive, neg=negative)
    net.cluster()
    fig=net.widget()
    return fig

material.sidebar.append(user_dropdown)
material.sidebar.append(user_input)
material.sidebar.append(participant_input)
material.sidebar.append(covariate_selection)
material.sidebar.append(positive_col)
material.sidebar.append(negative_col)

material.main.append(
    pn.Row(get_plot,height=1000)
    )
material.servable();
