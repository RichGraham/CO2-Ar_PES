!!CO2-Ar PES
!!Distances in Angstrom and energies in Hartrees
module PES_details
  double precision :: gpRMax = 0.0  
  double precision :: gpRMin = -100.0  
  double precision :: lCO = 1.1632
  double precision :: AngToBohr =  1.8897259885789
  double precision :: gpEmax = 0.456 !! Needs recoding so that this is read from training data
  !! Maximum energy in Hartree seen within the geometric constraint
  interface PES_GP 
     function PES_GP(xStar) 
       implicit none 
       double precision:: PES_GP
       double precision, dimension(:) ::  xStar
     end function PES_GP
  end interface PES_GP
end module PES_details


module GP_variables
  double precision, allocatable :: alpha (:), lScale(:), xTraining(:,:), xTrainingPerm(:,:)
  double precision expVar,NuggVar
  integer :: nDim=6
  integer :: nTraining=135 
end module GP_variables


 ! Test program
implicit none

integer k,i, choice

!call load_GP_Data
call fixedAngleSlice

end
!

subroutine fixedAngleSlice()
  use PES_details
  implicit none
  double precision rab(3)
  integer i, itot
  double precision  r, beta1, e, e_GP, asymp, PES
    
  itot=500
  beta1 =  90   /180.0*3.14159265359

  open (unit=15, file="PES_Out.dat ", status='replace')
  
  do i=0, itot

     ! specify centre-to-centre separation
     r = (  2.4 + 15.0*i/(1.0*itot) ) 

     call computeDistances(r,beta1,rab)
     
     
     e=PES( rab)
     !e_GP = PES_GP( xStar)
     write(15,*) r , e 
     
  enddo
  write(6,*)'Written to file: PES_Out.dat '
  close(15)

end subroutine fixedAngleSlice
  
  
  
subroutine computeDistances(r,beta1, rab)
  use PES_details
  implicit none
  double precision  rab(3), r, beta1
  integer ia, ib, ir,k

  rab(1)=r
  rab(2)= SQRT( (lCO*SIN(beta1))**2 +(lCO*COS(beta1)-r)**2  )
  rab(3)= SQRT( (lCO*SIN(beta1))**2 +(lCO*COS(beta1)+r)**2  )
   
end subroutine

double precision function asymp(rab)
  use PES_details
  implicit none
  double precision rab(3)
  double precision c1,c2,c3, cAr, o1Ar, o2Ar

  c1= (  rab(1)**2 + lCO**2 - rab(2)**2)/2.0/rab(1)/lCO
  
  c2= (  rab(2)**2 + lCO**2 - rab(1)**2)/2.0/rab(2)/lCO
  c3= (  rab(3)**2 + lCO**2 - rab(1)**2)/2.0/rab(3)/lCO


  cAr  = - ( 4.64 * (1 + 3*c1**2)  +  3.30 * (5 - 3*c1**2) ) / (rab(1)*AngToBohr)**6
  o1Ar = - ( 8.69 * (1 + 3*c2**2)  +  4.76 * (5 - 3*c2**2) ) / (rab(2)*AngToBohr)**6
  o2Ar = - ( 8.69 * (1 + 3*c3**2)  +  4.76 * (5 - 3*c3**2) ) / (rab(3)*AngToBohr)**6
  
  asymp= cAr + o1Ar + o2Ar
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Gaussian Process Code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pom


  
subroutine load_GP_Data
  use GP_variables
  use PES_details
  implicit none
  
  double precision, allocatable::  xStar(:)
  integer i,j
  double precision :: dum
  character (len=90) :: filename

  allocate (alpha(nTraining), lScale(nDim), xTraining(nDim,nTraining),xTrainingPerm(nDim,nTraining), xStar(nDim))

  !====Load hyperparameters====
  write (filename, '( "N", I4.4, "/HyperParams_Symm.dat" )' )  nTraining
  open (unit = 7, file = filename)
  !Only need to read some as others are tied.
  read (7,*) lScale(4), lScale(3), lScale(2), lScale(1),expVar,NuggVar
  !Copy over the tied values
  lScale(5) = lScale(3)
  lScale(6) = lScale(4)
  !print *,"HyperParams",lScale(1), lScale(2), lScale(3), lScale(4), lScale(5), lScale(6),expVar,NuggVar
  close(7)
  
  !====Load alpha coefficients====
  write (filename, '( "N", I4.4, "/alpha_Symm.dat" )' )  nTraining
  open (unit = 7, file = filename)
  do i=1,nTraining
     read (7,*) alpha(i)
     !print *,"alpha ",i, alpha(i)
  end do
  close(7)


  !====Load training data x values ====
  write (filename, '( "N", I4.4, "/xTraining.dat" )' )  nTraining
  open (unit = 7, file = filename)
    
  do i=1,nTraining
     read (7,*) xTraining(1,i), xTraining(2,i), xTraining(3,i), xTraining(4,i), xTraining(5,i), xTraining(6,i)  
  end do
  close(7)
  
  !! Permute the training vectors
  xTrainingPerm = xTraining
  do i=1,nTraining
     xTrainingPerm(3,i)=xTraining(5,i)
     xTrainingPerm(5,i)=xTraining(3,i)
     xTrainingPerm(4,i)=xTraining(6,i)
     xTrainingPerm(6,i)=xTraining(4,i)
  end do

end subroutine load_GP_Data
  
function PES_GP(xStar)
  use GP_variables
  implicit none
  double precision, dimension(:) :: xStar
  double precision:: PES_GP
  integer beta,i
  double precision kSqExp, kSqExpPerm, kKern

  kKern=0

  do i=1,nTraining
     kSqExp=1;
     kSqExpPerm=1;
     do beta=1,nDim
        kSqExp  =  kSqExp * ( exp( - (xStar(beta)-xTraining(beta,i))**2 /2.0/lScale(beta)**2) )
        kSqExpPerm  =  kSqExpPerm * ( exp( - (xStar(beta)-xTrainingPerm(beta,i))**2 /2.0/lScale(beta)**2) ) 
     end do
     kKern = kKern + alpha(i) * (kSqExp + kSqExpPerm)
  end do
  
  PES_GP=kKern * expVar
end function PES_GP



function PES( rab )
  !! Takes in rab in Angstrom
  use PES_details
  implicit none
  double precision rab(3), xStar(3), asymp
  double precision  PES

  
  if( rab(1) > gpRMax  .AND.  rab(2) > gpRMax .AND.  rab(3) > gpRMax &
       ) then !!Use asymptotic function
     PES = asymp(rab)
     
  else if (rab(1) < gpRMin  .OR.  rab(2) < gpRMin  .OR.  rab(3) < gpRMin &
       ) then !! Use repulsive approximation function
     PES=gpEmax
     
  else !! Use the Guassian Process function
     xStar(:) = 1/rab(:)
     PES = PES_GP( xStar)
  end if

end function PES
