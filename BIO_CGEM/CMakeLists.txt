set(schismmodlibs cgem ${schismmodlibs} PARENT_SCOPE)
add_library(cgem grid.F90 grid_setup.F90 cgem.F90 cgem_setup.F90 cgem_run.F90 cgem_step.F90
cgem_sink.F90
surface_flux.F90
calc_solar_zenith.F90
get_solar.F90
call_iop_par.F90
IOP_PARattenuation.F90
calcAgrow.F90
cgem_utils.F90
date_time.F90
moc_src/constants.F90
moc_src/depth2press.F90
moc_src/f2pCO2.F90
moc_src/gasx.F90
moc_src/p2fCO2.F90
moc_src/p80.F90
moc_src/phsolvers.F90
moc_src/rho.F90
moc_src/rhoinsitu.F90
moc_src/singledouble.F90
moc_src/sw_adtg.F90
moc_src/sw_ptmp.F90
moc_src/sw_temp.F90
moc_src/vars.F90
moc_src/varsolver.F90
 )
add_dependencies(cgem core ${schismmodlibs})
target_link_libraries(cgem core ${schismmodlibs})
