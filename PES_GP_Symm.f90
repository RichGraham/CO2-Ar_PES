

 ! Test program
implicit none

integer k,i, choice

call load_GP_CO2_Ar_Data
call fixedAngleSlice

end
!




subroutine fixedAngleSlice()
  use PES_CO2_Ar_details
  implicit none
  double precision rab(3)
  integer i, itot
  double precision  r, beta1, e, e_GP, asymp_CO2_Ar, PES_CO2_Ar
    
  itot=500
  beta1 =  0   /180.0*3.14159265359

  open (unit=15, file="PES_CO2_Ar_Out.dat ", status='replace')
  
  do i=0, itot

     ! specify centre-to-centre separation
     r = (  0.5 + 15.0*i/(1.0*itot) ) 

     call computeDistances(r,beta1,rab)
     
     
     e=PES_CO2_Ar( rab)
     !e_GP = PES_CO2_Ar_GP( xStar)
     write(15,*) r , e 
     
  enddo
  write(6,*)'Written to file: PES_CO2_Ar_Out.dat '
  close(15)

end subroutine fixedAngleSlice
  
  
  
subroutine computeDistances(r,beta1, rab)
  use PES_CO2_Ar_details
  implicit none
  double precision  rab(3), r, beta1
  integer ia, ib, ir,k

  rab(1)=r
  rab(2)= SQRT( (lCO*SIN(beta1))**2 +(lCO*COS(beta1)-r)**2  )
  rab(3)= SQRT( (lCO*SIN(beta1))**2 +(lCO*COS(beta1)+r)**2  )
   
end subroutine
