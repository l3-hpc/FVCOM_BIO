cgem_vars = c('A1','Qn1','Qp1','Zoo1','Zoo2','NO3','NH4','PO4','DIC','O2',
              'OM1CA','OM1NA','OM1PA','OM2CA','OM2NA','OM2PA','OM1CZ','OM1NZ',
              'OM1PZ','OM2CZ','OM2NZ','OM2PZ','OM1R','OM2R','CDOM','Si','OM1BC',
              'OM2BC','Alk','Tr')

cgem_units = c('cells/m3','mmol-N/cell','mmol-P/cell','individuals/m3','individuals/m3','mmol-N/m3','mmol-N/m3','mmol-P/m3',
              'mmol-C/m3','mmol-O2/m3',
              'mmol-C/m3','mmol-N/m3','mmol-P/m3','mmol-C/m3','mmol-N/m3','mmol-P/m3','mmol-C/m3','mmol-N/m3',
              'mmol-P/m3','mmol-C/m3','mmol-N/m3','mmol-P/m3','mmol-C/m3','mmol-C/m3','CDOM','Si','mmol-C/m3',
              'mmol-C/m3','mmol-HCO3/m3','none')

cgem_init_values = c(6.e7,0.30649887E-8,0.19438481E-9,150.,1505.,5.,1.,2.,2134.,172.,
                     1.E-5,1.5E-6,1.E-7,1.E-5,1.5E-6,1.E-7,1.E-5,1.5E-6,
                     1.E-7,1.E-5,1.5E-6,1.E-7,0.00001,0.00001,2.,15.,0.00001,
                     0.00001,2134.,1.)


for(i in 1:length(cgem_vars)){
  filename = paste0("NUTRIENT_INI_",i,".dat")
  writeLines(c("CONSTANT",cgem_init_values[i]),filename)
}
