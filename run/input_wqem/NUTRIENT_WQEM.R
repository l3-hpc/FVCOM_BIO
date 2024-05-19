wqem_vars = c('DOC','DIA','GRE','ZOO','LOC','ROC','SRP','DOP','LOP','ROP',
              'NH4','NO3','DON','LON','RON','SA','SU','DO2','TR')
#,'JDIAN','JDIAP','JGREN','JGREP')

wqem_units = c('C kg/m3','C kg/m3','C kg/m3','C kg/m3','C kg/m3','C kg/m3','P kg/m3','P kg/m3',
              'P kg/m3','P kg/m3','N kg/m3','N kg/m3','N kg/m3','N kg/m3','N kg/m3','Si kg/m3',
              'Si kg/m3','O2 kg/m3','none')
#,'N kg/m3','P kg/m3','N kg/m3','P kg/m3')

wqem_init_values = c(2.5e-3,8e-5,8e-5,1.6e-5,2e-2,2e-4,3e-6,4e-6,7e-6,6e-6,
                     2e-5,3e-4,5.5e-5,1.1e-4,1.1e-4,9e-4,1e-4,8e-3,1.e0)
#,4e-3,4e-3,4e-3,4e-3)


for(i in 1:length(wqem_vars)){
  filename = paste0("NUTRIENT_INI_",i,".dat")
  writeLines(c("CONSTANT",wqem_init_values[i]),filename)
}
