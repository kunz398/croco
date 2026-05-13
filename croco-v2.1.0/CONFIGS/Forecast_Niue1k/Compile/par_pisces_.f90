










!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!







                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      

                      
                      
                      
                      
                      


                      





                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      

                      
                     








                      
                      
                      
                      
                      


!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA,
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!















































































































!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!





















!









!-# define float dfloat
!-# define FLoaT dfloat
!-# define FLOAT dfloat
!-# define sqrt dsqrt
!-# define SQRT dsqrt
!-# define exp dexp
!-# define EXP dexp
!-# define dtanh dtanh
!-# define TANH dtanh









MODULE par_pisces
   !!======================================================================
   !!                        ***  par_pisces  ***
   !! TOP :   set the PISCES parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   IMPLICIT NONE
   PUBLIC

   !                                                                !!** Floating point **
   INTEGER, PUBLIC, PARAMETER ::   sp = SELECTED_REAL_KIND( 6, 37)   !: single precision (real 4)
   INTEGER, PUBLIC, PARAMETER ::   dp = SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
   INTEGER, PUBLIC, PARAMETER ::   wp = dp                              !: working precision
   !                                                                !!** Integer **
   INTEGER, PUBLIC, PARAMETER ::   i4 = SELECTED_INT_KIND( 9)        !: single precision (integer 4)
   INTEGER, PUBLIC, PARAMETER ::   i8 = SELECTED_INT_KIND(14)        !: double precision (integer 8)

   !                                                                !!** Integer **
   INTEGER, PUBLIC, PARAMETER ::   lc  = 256                          !: Lenght of Character strings
   INTEGER, PUBLIC, PARAMETER ::   lca = 400                          !: Lenght of Character arrays


   LOGICAL, PUBLIC, PARAMETER ::   lk_pisces     = .FALSE.  !: PISCES flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces     = 0      !: number of PISCES passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_2d  = 0      !: additional 2d output ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_pisces_3d  = 0      !: additional 3d output ('key_trc_diaadd')
   !
   INTEGER, PUBLIC ::   jpdic     !: dissolved inoganic carbon concentration 
   INTEGER, PUBLIC ::   jptal     !: total alkalinity 
   INTEGER, PUBLIC ::   jpoxy     !: oxygen carbon concentration 
   INTEGER, PUBLIC ::   jppoc     !: small particulate organic phosphate concentration
   INTEGER, PUBLIC ::   jpphy     !: phytoplancton concentration 
   INTEGER, PUBLIC ::   jpzoo     !: zooplancton concentration
   INTEGER, PUBLIC ::   jpdoc    !: dissolved organic carbon concentration 
   INTEGER, PUBLIC ::   jpno3    !: Nitrates Concentration
   INTEGER, PUBLIC ::   jpfer    !: Iron Concentration
   INTEGER, PUBLIC ::   jpcal     !: calcite  concentration 
   INTEGER, PUBLIC ::   jppo4     !: phosphate concentration 
   INTEGER, PUBLIC ::   jpsil     !: silicate concentration
   INTEGER, PUBLIC ::   jpdia     !: Diatoms Concentration
   INTEGER, PUBLIC ::   jpmes     !: Mesozooplankton Concentration
   INTEGER, PUBLIC ::   jpgsi     !: (big) Silicate Concentration
   INTEGER, PUBLIC ::   jpbfe     !: Big iron particles Concentration
   INTEGER, PUBLIC ::   jpgoc     !: big particulate organic phosphate concentration
   INTEGER, PUBLIC ::   jpsfe     !: Small iron particles Concentration
   INTEGER, PUBLIC ::   jpdfe     !: Diatoms iron Concentration
   INTEGER, PUBLIC ::   jpdsi     !: Diatoms Silicate Concentration
   INTEGER, PUBLIC ::   jpnfe     !: Nano iron Concentration
   INTEGER, PUBLIC ::   jpnch     !: Nano Chlorophyll Concentration
   INTEGER, PUBLIC ::   jpdch     !: Diatoms Chlorophyll Concentration
   INTEGER, PUBLIC ::   jpnh4     !: Ammonium Concentration
   INTEGER, PUBLIC ::   jplgw    !: Ammonium Concentration
   INTEGER, PUBLIC ::   jpdon    !: DON concentration 
   INTEGER, PUBLIC ::   jpdop    !: DOP concentration 
   INTEGER, PUBLIC ::   jppon    !: PON concentration
   INTEGER, PUBLIC ::   jppop    !: POP concentration
   INTEGER, PUBLIC ::   jpnph     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppph     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpndi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppdi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppic     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpnpi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpppi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppfe     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppch     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpgon    !: GON concentration
   INTEGER, PUBLIC ::   jpgop    !: GOP concentration

   INTEGER, PUBLIC ::   jp_flxco2  
   INTEGER, PUBLIC ::   jp_flxo2   
   INTEGER, PUBLIC ::   jp_kgco2   
   INTEGER, PUBLIC ::   jp_dpco2   
   INTEGER, PUBLIC ::   jp_sinkco2 
   INTEGER, PUBLIC ::   jp_sinkfer 
   INTEGER, PUBLIC ::   jp_sinksil 
   INTEGER, PUBLIC ::   jp_sinkcal 
   INTEGER, PUBLIC ::   jp_heup    
   INTEGER, PUBLIC ::   jp_sildep   
   INTEGER, PUBLIC ::   jp_po4dep   
   INTEGER, PUBLIC ::   jp_no3dep   
   INTEGER, PUBLIC ::   jp_nh4dep   
   INTEGER, PUBLIC ::   jp_nitrpot 

   INTEGER, PUBLIC ::   jp_hi      
   INTEGER, PUBLIC ::   jp_co3     
   INTEGER, PUBLIC ::   jp_co3sat  
   INTEGER, PUBLIC ::   jp_etot    
   INTEGER, PUBLIC ::   jp_pphy    
   INTEGER, PUBLIC ::   jp_pphy2   
   INTEGER, PUBLIC ::   jp_pnew    
   INTEGER, PUBLIC ::   jp_pnew2   
   INTEGER, PUBLIC ::   jp_pbsi    
   INTEGER, PUBLIC ::   jp_pfed    
   INTEGER, PUBLIC ::   jp_pfen    
   INTEGER, PUBLIC ::   jp_pnewo2  
   INTEGER, PUBLIC ::   jp_prego2  
   INTEGER, PUBLIC ::   jp_grapoc   
   INTEGER, PUBLIC ::   jp_grapoc2   
   INTEGER, PUBLIC ::   jp_mico2  
   INTEGER, PUBLIC ::   jp_meso2  
   INTEGER, PUBLIC ::   jp_nitrifo2 
   INTEGER, PUBLIC ::   jp_remino2 
   INTEGER, PUBLIC ::   jp_nfixo2  
   INTEGER, PUBLIC ::   jp_irondep  
   INTEGER, PUBLIC ::   jp_ironsed  

   ! Starting/ending PISCES do-loop indices (N.B. no PISCES : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jptra       = jp_pisces                  !: First index of PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs0     = 1                  !: First index of PISCES tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_pcs1     = jp_pisces          !: Last  index of PISCES tracers

   REAL(wp), PUBLIC ::  mMass_C  = 12.00      ! Molar mass of carbon
   REAL(wp), PUBLIC ::  mMass_N  = 14.00      ! Molar mass of nitrogen
   REAL(wp), PUBLIC ::  mMass_P  = 31.00      ! Molar mass of phosphorus
   REAL(wp), PUBLIC ::  mMass_Fe = 55.85      ! Molar mass of iron
   REAL(wp), PUBLIC ::  mMass_Si = 28.00      ! Molar mass of silver

   !!======================================================================
END MODULE par_pisces
