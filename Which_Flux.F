Module Which_Flux
! =========================================================
! Define which fluxes shall be used
! =========================================================

! Which_Flux(iO2surf,iDICsurf,iSOC,iMPB,iNutEx,iCMAQ,iInRemin)

INTEGER, PARAMETER :: iO2surf  = 1 !O2 surface flux
INTEGER, PARAMETER :: iDICsurf = 2 !DIC surface flux
INTEGER, PARAMETER :: iSOC     = 3 !Sediment Oxygen Consumption
INTEGER, PARAMETER :: iMPB     = 4 !Microphytobethos
INTEGER, PARAMETER :: iNutEx   = 5 !Sediment Nutrient Fluxes
INTEGER, PARAMETER :: iCMAQ    = 6 !CMAQ surface deposition of NH4 and NO3
INTEGER, PARAMETER :: iInRemin = 7 !Instant Remineralization in bottom layer
INTEGER, PARAMETER :: iSDM     = 8 !Sediment Diagenesis Model
INTEGER, PARAMETER :: i_Si      = 9 !Silica (SA, SRP) Fluxes

END Module Which_Flux
