# Determining-Critical-Periods-for-Thermal-Accliamtion-Using-a-DLNM-approach

This repository contains data and code used in the paper "Determining Critical Periods for Thermal Acclimation Using a Distributed Lag Non-linear Model approach". The paper test the use of Distributed Lag Non-linear Models (DLNM) to model the effect of water temperaure on variation of mayfly Critical Thermal maxima (CTmax) accounting for delayed variation in time and potenital linearity of these realtions. 
The scripts and data avaible here drive throug the full pipeline presetned in the paper, from water temperature data managing, data preparation and model implementation.
The piepline is organized as following:

-STEP 1: water temperaure exposures
In this step water temperature data recorded in filed are used to estiamate water temperature non recorded, to have for each CTmax sample the preceding 28 days of water temperature exposures on an hourly basis.
  
Datasets: "temp_data_16"(original field measured water temperaures with HOBO temperaure loggers for 2016).
          "temp_data_17"(original field measured water temperaures with HOBO temperaure loggers for 2017).
          "tulloch_bridge_2016" (MetOffice hourly air temperaure data used to estimate  missing water temperaure in Loch Morie)
           "strath_air_2016"(MetOffice hourly air temperaure data used to estimate missing water temperaure in Glenquey reservoir 2016)
           "strath_air_2017"(MetOffice hourly air temperaure data used to estimate missing water temperaure in Glenquey reservoir 2017)
           "cairnwell_2017"(MetOffice hourly air temperaure data used to estimate missing water temperaure in Loch Rescobie and Loch of Lintrathen 2017)
            
codes to run: "2016_T_pred", "2017_T_pred"
  
  
  -STEP 2: data managing. Temperaure data (i.e.exposures) are standardized with the three methods presetnted and for each mean CTmax is produced the raltive of exposures profile.
  
  datasets:"allCT" (original dataset of all CTmax measurement for each individual sampled)
  codes to run: "comp_variables"
  
  -STEP 3: DLNM implementation. Three distinct DLNM are implemented (one ofr each exposures: Z-score, Distance from trend and temperaure variation rate).
  
  codes to run:"dlnm", "heatplot_function"(support function to produce estiamted DLNM effect with the filtering of non-significative areas)
