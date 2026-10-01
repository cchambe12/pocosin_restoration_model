
1) Dataset Title: Data from: Assessment of soil carbon sequestration or losses from drained short pocosins located in Hyde County, NC, during the years 2020 and 2021

2) Investigators:
	

	Neal Flanagan, Ph.D.
	Research Professor
	Duke University Wetland Center
	Box 90328
	Durham, NC 27708
	email: nflanaga@duke.edu

	Curt Richardson, Ph.D.
	Professor
	Duke University Wetland Center
	Box 90333
	Durham, NC 27708
	email: curtr@duke.edu

3) Date of Data Collection July 2019 to March 2021

4) Locations: Hyde County, NC, USA

6) Language: English

7) Funding: Duke University, Grantham Foundation


DATA AND FILE OVERVIEW


These files contain the data used to generate the figures contained in the following journal article:

Article ID: GCB16366
Article DOI: 10.1111/gcb.16366
Article: Annual Carbon Sequestration and Loss Rates Under Altered Hydrology and Fire Regimes in Southeastern USA Pocosin Peatlands 
Journal: Global Change Biology 



FILES A-D:


A.  	File: C_FLUX_TOVI_OUTPUT.csv
	Description:30 minute carbon fluxes calculated Li-Cor TOVI software from Eddy Covariance tower at site G11N2 as umols (not C or CO2 mass).  Data is gapfilled after QC, 	U* threshold processing, and footprint analysis, quality flags.
	Column descriptions:
	1) Date and Time
	2) GPP - Gross Primary Productivity (umol.s^-1.m^-2)			
	3) NEE - Net Ecosystem Exchange (umol.s^-1.m^-2)	
	4) RECO - Ecosystem Respiration (umol.s^-1.m^-2)			

B.	File: STAT_DATA_NOT_GAPFILLED.csv
	Description: SUbset of actual observations (not gapfilled), and filtered to remove data with quality flags or low instrument signal.  NEE processed using R package 	REddyProc 	https://cran.rstudio.com/web/packages/REddyProc/ 
	Column descriptions:
	1) RAD_GLOB - Global Radiation (W.m^-2)
	2) WTD - Depth of watertable below soil surface (cm)
	3) NEE - Net Ecosystem Exchange (umol.s^-1.m^-2)


C.	File: WELL_WATERTABLE_DEPTHS.csv
	Description: Depth of watertable below soil surface (cm).  Negative values indicate position below soil surace, postive value indicates depth of inundation above the 	soil. Symbol "----" indicated USGS well was not yet operational at site PLNWR_DNL.
	Column description:
	1) TimeStamp- Date and Time
	2) PD_G11 - Well located in block G11 of Private Drained Pococsin Site (35.690526°N, 76.410167°W)
	3) PD_F11 - Well located in block F11 of Private Drained Pococsin Site (35.696953°N, 76.426731°W)
	4) PD_E11 - Well located in block E11 of Private Drained Pococsin Site (35.697869°N, 76.443342°W)
	5) PLNWR_D5 - Well located in restored block D5 of the Pocosin Lakes National Wildlife Refuge (35.702722°N, 76.454694°W)
	6) PLNWR_DNL - Well located in reference block of the Pocosin Lakes National Wildlife Refuge  (35.690344°N, 76.528300°W)

D.	File: CHAMBER_CH4_TEMP.csv
	Description: The temperature coefficients (Q10) for soil incubations.
	Column descriptions:
      1) Soil moisture (v/v)
      2) Soil temperature (°C)                              
      3) CH4 (nmol m-2 s-1)         




SHARING AND ACCESS INFORMATION

The data that support the findings of this study are openly available under a CC0 waiver.


Citation:

Richardson, C.J., Flanagan, N.E., Wang, H., M. Ho. 2022.  Dataset for assessment of soil carbon sequestration or losses from drained short pocosins located in Hyde County, NC, during the years 2020 and 2021.  Duke Research Data Repository. https://doi.org/10.7924/r46w9hp7v. 